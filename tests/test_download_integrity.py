# SPDX-License-Identifier: GPL-3.0-or-later
"""Intégrité des téléchargements : réponse tronquée, page d'erreur reçue en HTTP 200, cache invalide.

Les tests utilisent un vrai serveur HTTP local (127.0.0.1), pas de simulation de la bibliothèque :
c'est le comportement réel de urllib et http.client face à une connexion coupée qui est vérifié."""
import http.server
import os
import socketserver
import sys
import threading
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import download  # noqa: E402
import sources  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402

XML_ERROR = b'<?xml version="1.0"?><ServiceExceptionReport>serveur surcharge</ServiceExceptionReport>'
PNG = b"\x89PNG\r\n\x1a\n" + b"\x00" * 64


@pytest.fixture
def server():
    """Serveur local dont le comportement est piloté par `state` : liste de réponses à servir dans
    l'ordre ('ok', 'truncated', 'xml200'), la dernière étant répétée. `state['hits']` compte les requêtes."""
    state = {"script": ["ok"], "hits": 0, "payload": b"II*\x00" + bytes(range(256)) * 40}

    class Handler(http.server.BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, *a):
            pass

        def do_GET(self):
            i = min(state["hits"], len(state["script"]) - 1)
            mode = state["script"][i]
            state["hits"] += 1
            payload = state["payload"]
            if mode == "truncated":
                self.send_response(200)
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload[: len(payload) // 2])
                self.wfile.flush()
                self.close_connection = True
            elif mode == "xml200":
                self.send_response(200)
                self.send_header("Content-Length", str(len(XML_ERROR)))
                self.end_headers()
                self.wfile.write(XML_ERROR)
            else:
                self.send_response(200)
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload)

    srv = socketserver.TCPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    state["url"] = f"http://127.0.0.1:{srv.server_address[1]}/tile"
    yield state
    srv.shutdown()
    srv.server_close()


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    monkeypatch.setattr(sources.time, "sleep", lambda s: None)


def test_http_get_retries_a_truncated_response_then_succeeds(server):
    server["script"] = ["truncated", "ok"]
    data = sources._http_get(server["url"])
    assert data == server["payload"] and server["hits"] == 2


def test_http_get_persistent_truncation_is_source_unavailable_not_a_raw_exception(server):
    server["script"] = ["truncated"]
    with pytest.raises(sources.SourceUnavailable):
        sources._http_get(server["url"])
    assert server["hits"] == sources.HTTP_RETRIES


def test_download_rejects_http_200_error_page_and_does_not_cache_it(server, tmp_path):
    server["script"] = ["xml200"]
    dest = str(tmp_path / "c" / "tile.tif")
    with pytest.raises(sources.SourceUnavailable) as e:
        sources._download(server["url"], dest, validate=sources._is_tiff)
    assert "ServiceExceptionReport" in str(e.value)       # le début de la réponse aide au diagnostic
    assert not os.path.exists(dest) and not os.path.exists(dest + ".part")


def test_download_discards_an_invalid_cached_file_and_downloads_again(server, tmp_path):
    dest = tmp_path / "c" / "tile.tif"
    dest.parent.mkdir()
    dest.write_bytes(XML_ERROR)                            # cache empoisonné par une exécution antérieure
    out = sources._download(server["url"], str(dest), validate=sources._is_tiff)
    assert out == str(dest) and dest.read_bytes() == server["payload"] and server["hits"] == 1


def test_download_keeps_a_valid_cached_file_without_network(server, tmp_path):
    dest = tmp_path / "tile.tif"
    dest.write_bytes(server["payload"])
    sources._download(server["url"], str(dest), validate=sources._is_tiff)
    assert server["hits"] == 0


def test_download_without_validation_keeps_previous_behaviour(server, tmp_path):
    server["script"] = ["xml200"]
    dest = tmp_path / "x.json"
    assert sources._download(server["url"], str(dest)) == str(dest) and dest.read_bytes() == XML_ERROR


def test_gsi_tile_rejects_non_png_and_does_not_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", lambda *a, **k: XML_ERROR)
    with pytest.raises(sources.SourceUnavailable):
        sources._gsi_tile(14, 100, 50, str(tmp_path))
    assert not list(tmp_path.rglob("*.png*"))


def test_gsi_tile_discards_invalid_cached_png_and_downloads_again(tmp_path, monkeypatch):
    good = None
    from test_sources import make_png
    good = make_png(np.full((256, 256), 100.0))
    cached = tmp_path / "gsi" / "dem1a_png" / "17" / "100_50.png"
    cached.parent.mkdir(parents=True)
    cached.write_bytes(XML_ERROR)
    monkeypatch.setattr(sources, "_http_get", lambda *a, **k: good)
    h, layer = sources._gsi_tile(17, 100, 50, str(tmp_path))
    assert layer == "dem1a_png" and h.shape == (256, 256) and cached.read_bytes() == good


def test_gsi_tile_unreadable_png_is_removed_from_cache(tmp_path, monkeypatch):
    bad = PNG                                              # signature PNG valide, contenu illisible
    monkeypatch.setattr(sources, "_http_get", lambda *a, **k: bad)
    with pytest.raises(sources.SourceUnavailable):
        sources._gsi_tile(17, 100, 50, str(tmp_path))
    assert not list(tmp_path.rglob("100_50.png"))


def test_copernicus_truncated_download_is_not_cached(server, tmp_path, monkeypatch):
    server["script"] = ["truncated"]
    monkeypatch.setattr(download, "tile_url", lambda *a, **k: server["url"])
    with pytest.raises(RuntimeError, match="incomplet"):
        download.download_copernicus(35.36, 138.73, 500, str(tmp_path / "cache"), log=lambda *_: None)
    assert not list((tmp_path / "cache").glob("*.tif"))   # avant correctif : un .tif tronqué y restait


def test_copernicus_error_page_is_not_cached(server, tmp_path, monkeypatch):
    server["script"] = ["xml200"]
    monkeypatch.setattr(download, "tile_url", lambda *a, **k: server["url"])
    with pytest.raises(RuntimeError):
        download.download_copernicus(35.36, 138.73, 500, str(tmp_path / "cache"), log=lambda *_: None)
    assert not list((tmp_path / "cache").glob("*.tif"))


def test_copernicus_invalid_cached_tile_is_downloaded_again(server, tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    cache.mkdir()
    name = download.tile_name(35, 138)
    (cache / f"{name}.tif").write_bytes(XML_ERROR)
    good = tmp_path / "good.tif"
    write_synthetic_tif(good, np.full((10, 10), 100.0), 138.0, 36.0, 0.1, 4326)
    server["payload"] = good.read_bytes()
    monkeypatch.setattr(download, "tile_url", lambda *a, **k: server["url"])
    paths = download.download_copernicus(35.36, 138.73, 500, str(cache), log=lambda *_: None)
    assert (cache / f"{name}.tif").read_bytes() == server["payload"] and server["hits"] == 1 and len(paths) == 1
