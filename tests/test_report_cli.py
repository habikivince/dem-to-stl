# SPDX-License-Identifier: GPL-3.0-or-later
"""Ligne de relief du journal et option --report de circle_dem_to_stl.py, de bout en bout (main()).

Le raster est remplacé par un cône synthétique renvoyé directement à la place de raster.warp_window : aucune lecture
GDAL n'intervient, seul le câblage est testé (journal, calculs, écriture du JSON). Le relief annoncé est confronté à la
hauteur réelle du STL écrit."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import circle_dem_to_stl  # noqa: E402
import formats  # noqa: E402
import raster  # noqa: E402

CX, CY = 1000.0, 2000.0


@pytest.fixture
def cone(monkeypatch):
    """Île conique (80 m, rayon 100 m) sur nodata, et court-circuit de la lecture GDAL."""
    def warp(input_path, dst, cx, cy, radius, pixel, progress_cb=None, fill=True, inside_fn=None, resample="auto",
             extrapolate_edges=False, stats=None):
        n = int(round(2 * radius / pixel)) + 1
        xs = cx - radius + np.arange(n) * pixel
        ys = cy + radius - np.arange(n) * pixel
        X, Y = np.meshgrid(xs, ys)
        dist = np.hypot(X - cx, Y - cy)
        return np.where(dist <= 100.0, 80.0 - dist * 0.6, -9999.0), xs, ys, -9999.0

    monkeypatch.setattr(raster, "resolve_target_crs_and_center", lambda p, lat, lon, cx, cy: ("WKT", cx, cy))
    monkeypatch.setattr(raster, "warp_window", warp)
    monkeypatch.setattr(raster, "check_readable", lambda p: None)
    monkeypatch.setattr(formats, "prepare", lambda paths, opts=None, log=None: ("raster", ["synthetique.tif"]))


def run(monkeypatch, capsys, tmp_path, *extra):
    argv = ["circle_dem_to_stl.py", "--input", "synthetique.tif", "--cx", str(CX), "--cy", str(CY), "--radius", "140",
            "--diameter-mm", "56", "--base-mm", "3", "--vexag", "3", "--sea-level", "0",
            "--print-spacing-mm", "0.5", "--output-dir", str(tmp_path), *extra]
    monkeypatch.setattr(sys, "argv", argv)
    circle_dem_to_stl.main()
    out = capsys.readouterr()
    return out.out + out.err


def test_relief_line_is_logged_without_any_option(cone, monkeypatch, capsys, tmp_path):
    out = run(monkeypatch, capsys, tmp_path)
    assert "relief : altitudes de 0.0 à 80.0 m (dénivelé 80 m), soit 48.0 mm au-dessus du socle" in out
    assert "16.0 mm à l'échelle vraie" in out
    assert not list(tmp_path.glob("*.json"))


def test_report_matches_the_written_model(cone, monkeypatch, capsys, tmp_path):
    run(monkeypatch, capsys, tmp_path, "--format", "stl,3mf", "--report")
    stl = next(tmp_path.glob("*.stl"))
    rep = json.loads(next(tmp_path.glob("*.json")).read_text(encoding="utf-8"))
    assert rep["outputs"].keys() == {"stl", "3mf"} and Path(rep["outputs"]["stl"]) == stl.resolve()
    tris = check_stl.read_stl(str(stl))
    assert tris[:, :, 2].max() == pytest.approx(rep["relief"]["relief_mm"], abs=1e-3)
    assert rep["mesh"]["triangles_total"] == len(tris)
    assert rep["model"]["shape"] == "circle" and rep["model"]["diameter_mm"] == 56.0 and rep["model"]["vexag"] == 3.0
    assert rep["model"]["scale_mm_per_m"] == pytest.approx(0.2)
    assert rep["relief"]["relief_m"] == pytest.approx(80.0) and rep["relief"]["relief_true_mm"] == pytest.approx(16.0)
    assert rep["relief"]["suggested_vexag"] == 1.0
    assert rep["source"]["input"] == ["synthetique.tif"] and rep["source"]["used"] is None
    assert rep["argv"][-1] == "--report" and rep["schema"] == 1


def test_report_with_explicit_name_and_dry_run_writes_nothing(cone, monkeypatch, capsys, tmp_path):
    target = tmp_path / "sous" / "rapport.json"
    out = run(monkeypatch, capsys, tmp_path, "--dry-run", "--report", str(target))
    assert "rapport.json" in out and list(tmp_path.iterdir()) == []
    run(monkeypatch, capsys, tmp_path, "--report", str(target))
    assert json.loads(target.read_text(encoding="utf-8"))["mesh"]["triangles_total"] > 0


def test_unwritable_report_does_not_lose_the_model(cone, monkeypatch, capsys, tmp_path):
    def boom(path, data):
        raise OSError("disque plein")
    monkeypatch.setattr(circle_dem_to_stl.report, "write_report", boom)
    out = run(monkeypatch, capsys, tmp_path, "--report")
    assert "ATTENTION : rapport non écrit (disque plein)" in out and list(tmp_path.glob("*.stl"))
