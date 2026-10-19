# SPDX-License-Identifier: GPL-3.0-or-later
"""Retry-After : sur 429 et 503, on attend au moins la durée demandée par le serveur (plafonnée)."""
import email.utils
import http.server
import socketserver
import sys
import threading
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import sources  # noqa: E402


@pytest.fixture
def server():
    """Serveur local : state['script'] liste des réponses (code, valeur de Retry-After ou None), la dernière est répétée."""
    state = {"script": [(200, None)], "hits": 0}

    class Handler(http.server.BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, *a):
            pass

        def do_GET(self):
            code, ra = state["script"][min(state["hits"], len(state["script"]) - 1)]
            state["hits"] += 1
            body = b"ok" if code == 200 else b"trop de requetes"
            self.send_response(code)
            if ra is not None:
                self.send_header("Retry-After", ra)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

    srv = socketserver.TCPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    state["url"] = f"http://127.0.0.1:{srv.server_address[1]}/x"
    yield state
    srv.shutdown()
    srv.server_close()


@pytest.fixture
def sleeps(monkeypatch):
    calls = []
    monkeypatch.setattr(sources.time, "sleep", lambda s: calls.append(s))
    return calls


def test_waits_at_least_the_retry_after_given_in_seconds(server, sleeps):
    server["script"] = [(429, "7"), (200, None)]
    assert sources._http_get(server["url"]) == b"ok" and sleeps == [7.0]


def test_503_is_treated_like_429(server, sleeps):
    server["script"] = [(503, "3"), (200, None)]
    assert sources._http_get(server["url"]) == b"ok" and sleeps == [3.0]


def test_without_retry_after_the_exponential_backoff_is_unchanged(server, sleeps):
    server["script"] = [(429, None), (429, None), (200, None)]
    assert sources._http_get(server["url"]) == b"ok" and sleeps == [1, 2]


def test_a_short_retry_after_does_not_shorten_the_backoff(server, sleeps):
    server["script"] = [(429, "0"), (200, None)]
    assert sources._http_get(server["url"]) == b"ok" and sleeps == [1]


def test_the_wait_is_capped(server, sleeps):
    server["script"] = [(429, "86400"), (200, None)]
    assert sources._http_get(server["url"]) == b"ok" and sleeps == [float(sources.HTTP_RETRY_AFTER_MAX)]


def test_http_date_format(server, sleeps):
    when = email.utils.format_datetime(datetime.now(timezone.utc) + timedelta(seconds=10), usegmt=True)
    server["script"] = [(429, when), (200, None)]
    assert sources._http_get(server["url"]) == b"ok" and 8.0 <= sleeps[0] <= 10.0


def test_unreadable_retry_after_falls_back_to_the_backoff(server, sleeps):
    server["script"] = [(429, "bientôt".encode("utf-8").decode("latin-1")), (200, None)]
    assert sources._http_get(server["url"]) == b"ok" and sleeps == [1]


def test_persistent_429_ends_in_source_unavailable_after_all_attempts(server, sleeps):
    server["script"] = [(429, "2")]
    with pytest.raises(sources.SourceUnavailable):
        sources._http_get(server["url"])
    assert server["hits"] == sources.HTTP_RETRIES and sleeps == [2.0] * (sources.HTTP_RETRIES - 1)


def test_404_is_still_none_without_waiting(server, sleeps):
    server["script"] = [(404, "9")]
    assert sources._http_get(server["url"]) is None and sleeps == []


@pytest.mark.parametrize("value,expected", [(None, 0.0), ("", 0.0), ("5", 5.0), (" 12 ", 12.0), ("-3", 0.0), ("abc", 0.0), ("999999", 60.0)])
def test_retry_after_parsing(value, expected):
    assert sources._retry_after_seconds(value) == expected


def test_past_http_date_means_no_extra_wait():
    now = datetime(2026, 10, 11, 12, 0, 0, tzinfo=timezone.utc)
    past = email.utils.format_datetime(now - timedelta(seconds=30), usegmt=True)
    future = email.utils.format_datetime(now + timedelta(seconds=15), usegmt=True)
    assert sources._retry_after_seconds(past, now=now) == 0.0
    assert sources._retry_after_seconds(future, now=now) == 15.0
