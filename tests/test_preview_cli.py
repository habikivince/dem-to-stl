# SPDX-License-Identifier: GPL-3.0-or-later
"""Option --preview de circle_dem_to_stl.py, de bout en bout (main()), sur le cône synthétique de test_report_cli
(la lecture GDAL est court-circuitée : seul le câblage est testé)."""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import circle_dem_to_stl  # noqa: E402
from test_preview import decode_png  # noqa: E402
from test_report_cli import cone, run  # noqa: E402,F401


def test_preview_is_written_before_meshing_and_matches_the_grid(cone, monkeypatch, capsys, tmp_path):
    out = run(monkeypatch, capsys, tmp_path, "--preview")
    png = next(tmp_path.glob("*.png"))
    img = decode_png(png.read_bytes())
    assert img.shape == (113, 113, 2)                                  # grille de travail : 113 x 113 points
    assert img[56, 56, 1] == 255 and img[0, 0, 1] == 0                 # centre de l'île opaque, coin transparent
    assert out.index("aperçu :") < out.index("triangles (dessus)")     # écrit avant l'étape longue
    assert "nord en haut" in out and str(png) in out.replace("\n", "")


def test_preview_is_listed_in_the_report(cone, monkeypatch, capsys, tmp_path):
    run(monkeypatch, capsys, tmp_path, "--preview", "--report")
    rep = json.loads(next(tmp_path.glob("*.json")).read_text(encoding="utf-8"))
    assert Path(rep["outputs"]["preview"]) == next(tmp_path.glob("*.png")).resolve()


def test_no_preview_without_the_option_and_dry_run_writes_nothing(cone, monkeypatch, capsys, tmp_path):
    run(monkeypatch, capsys, tmp_path)
    assert not list(tmp_path.glob("*.png"))
    other = tmp_path / "vide"
    out = run(monkeypatch, capsys, other, "--dry-run", "--preview", str(other / "a.png"))
    assert "a.png" in out
    assert not other.exists() or not list(other.iterdir())


def test_unwritable_preview_does_not_stop_the_model(cone, monkeypatch, capsys, tmp_path):
    def boom(*a, **k):
        raise OSError("disque plein")
    monkeypatch.setattr(circle_dem_to_stl.preview, "write_preview", boom)
    out = run(monkeypatch, capsys, tmp_path, "--preview", "--report")
    assert "ATTENTION : aperçu non écrit (disque plein)" in out and list(tmp_path.glob("*.stl"))
    rep = json.loads(next(tmp_path.glob("*.json")).read_text(encoding="utf-8"))
    assert "preview" not in rep["outputs"]
