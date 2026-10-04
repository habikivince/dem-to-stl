# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests de geocode.py (--place). Aucun accès réseau : Nominatim est simulé."""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import circle_dem_to_stl  # noqa: E402
import geocode  # noqa: E402
import sources  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402
import numpy as np  # noqa: E402

quiet = lambda *a, **k: None  # noqa: E731

FUJI = [
    {"display_name": "Mont Fuji, Shizuoka, Japon", "lat": "35.36049", "lon": "138.72735", "category": "natural", "type": "volcano"},
    {"display_name": "Fuji, Shizuoka, Japon", "lat": "35.16190", "lon": "138.67620", "category": "boundary", "type": "administrative"},
]


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch):
    monkeypatch.setattr(geocode.time, "sleep", lambda s: None)


def fake_nominatim(payload, calls):
    def _get(url, timeout=None, retries=None):
        calls.append(url)
        return json.dumps(payload).encode("utf-8")
    return _get


def test_resolve_first_result_pick_and_cache(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(sources, "_http_get", fake_nominatim(FUJI, calls))
    lat, lon = geocode.resolve("Mont Fuji", 1, str(tmp_path), log=quiet)
    assert (lat, lon) == pytest.approx((35.36049, 138.72735))
    assert "q=Mont+Fuji" in calls[0] and "format=jsonv2" in calls[0] and "limit=5" in calls[0]
    # 2e appel (même requête, casse/espaces différents) : cache disque, aucun réseau
    monkeypatch.setattr(sources, "_http_get", lambda *a, **k: pytest.fail("réseau malgré le cache"))
    assert geocode.resolve("  mont   FUJI ", 2, str(tmp_path), log=quiet) == pytest.approx((35.16190, 138.67620))


def test_resolve_warns_on_administrative_centre(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", fake_nominatim(FUJI, []))
    msgs = []
    geocode.resolve("Fuji", 2, str(tmp_path), log=msgs.append)
    assert any("zone administrative" in m for m in msgs)
    assert any("autre candidat 1" in m for m in msgs)


def test_resolve_errors(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", fake_nominatim([], []))
    with pytest.raises(geocode.GeocodeError, match="aucun résultat"):
        geocode.resolve("zzzzzz", 1, str(tmp_path), log=quiet)
    monkeypatch.setattr(sources, "_http_get", fake_nominatim(FUJI, []))
    with pytest.raises(geocode.GeocodeError, match="hors limites"):
        geocode.resolve("Mont Fuji", 3, str(tmp_path), log=quiet)

    def down(url, **k):
        raise sources.SourceUnavailable("hors ligne")
    monkeypatch.setattr(sources, "_http_get", down)
    with pytest.raises(geocode.GeocodeError, match="indisponible"):
        geocode.resolve("Everest", 1, str(tmp_path / "other"), log=quiet)
    monkeypatch.setattr(sources, "_http_get", lambda *a, **k: b"<html>oops</html>")
    with pytest.raises(geocode.GeocodeError, match="invalide"):
        geocode.resolve("Everest2", 1, str(tmp_path / "other"), log=quiet)
    with pytest.raises(geocode.GeocodeError, match="vide"):
        geocode.resolve("   ", 1, str(tmp_path), log=quiet)


def test_rate_limit_waits_between_requests(tmp_path, monkeypatch):
    waits = []
    monkeypatch.setattr(geocode.time, "sleep", waits.append)
    monkeypatch.setattr(sources, "_http_get", fake_nominatim(FUJI, []))
    geocode._last_request[0] = geocode.time.monotonic()
    geocode.search("Mont Fuji", str(tmp_path), log=quiet)
    assert waits and 0 < waits[0] <= geocode.MIN_INTERVAL_S


def test_cli_place_implies_auto_source_and_is_exclusive(tmp_path, monkeypatch):
    lat, lon = 46.0, 7.6
    x, y = sources._to_epsg(2056, lat, lon)
    tif = tmp_path / "ch.tif"
    write_synthetic_tif(tif, np.full((1500, 1500), 1500.0), x - 1500.0, y + 1500.0, 2.0, 2056)
    seen = {}

    def fake_fetch(lat_, lon_, radius_m, cache_dir, pixel_size_m, log=print):
        seen["pos"] = (lat_, lon_)
        return [str(tif)]

    monkeypatch.setitem(sources.SOURCES, "swisstopo", sources.Source("swisstopo", "x", [(5, 45, 11, 48)], fake_fetch))
    monkeypatch.setattr(geocode, "resolve", lambda q, pick, cache_dir, log=print: (lat, lon))
    out = tmp_path / "o.stl"
    monkeypatch.setattr(sys, "argv", ["x", "--place", "Cervin", "--radius", "500", "--diameter-mm", "50",
                                      "--vexag", "1", "--base-mm", "2", "--print-spacing-mm", "0.5",
                                      "--output", str(out), "--cache-dir", str(tmp_path / "c")])
    circle_dem_to_stl.main()
    assert out.stat().st_size > 84 and seen["pos"] == (lat, lon)

    monkeypatch.setattr(sys, "argv", ["x", "--place", "Cervin", "--lat", "46", "--lon", "7", "--radius", "500",
                                      "--output", str(tmp_path / "p.stl")])
    with pytest.raises(SystemExit):
        circle_dem_to_stl.main()


def test_cli_place_not_found_exits_cleanly(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", fake_nominatim([], []))
    monkeypatch.setattr(sys, "argv", ["x", "--place", "zzzz", "--radius", "500", "--output", str(tmp_path / "o.stl"),
                                      "--cache-dir", str(tmp_path / "c")])
    with pytest.raises(SystemExit) as e:
        circle_dem_to_stl.main()
    assert "aucun résultat" in str(e.value)
