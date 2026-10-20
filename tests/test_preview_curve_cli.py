# SPDX-License-Identifier: GPL-3.0-or-later
"""Combinaison --preview + --vexag-curve de circle_dem_to_stl.py : l'aperçu doit montrer le relief IMPRIMÉ, donc les
altitudes après la courbe, pas les altitudes brutes (cône synthétique de test_report_cli, lecture GDAL court-circuitée)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from test_preview import decode_png  # noqa: E402
from test_report_cli import cone, run  # noqa: E402,F401

WEST, EAST = (56, 44), (56, 68)          # même ligne que le centre, 12 pixels (30 m) de part et d'autre du sommet


def _gray(png_path):
    img = decode_png(Path(png_path).read_bytes())
    return int(img[WEST[0], WEST[1], 0]), int(img[EAST[0], EAST[1], 0])


def test_preview_shows_the_curved_relief_not_the_raw_one(cone, monkeypatch, capsys, tmp_path):
    raw_dir, curved_dir = tmp_path / "brut", tmp_path / "courbe"
    run(monkeypatch, capsys, raw_dir, "--preview")
    west_raw, east_raw = _gray(next(raw_dir.glob("*.png")))
    assert west_raw != east_raw                                        # versant ouest éclairé, versant est dans l'ombre

    # courbe qui écrête au-dessus de 40 m : tout le sommet (r < 67 m) devient un plateau, donc plat à l'aperçu
    run(monkeypatch, capsys, curved_dir, "--preview", "--vexag-curve", "0:0,40:40,80:40")
    west_curved, east_curved = _gray(next(curved_dir.glob("*.png")))
    assert west_curved == east_curved
    assert west_curved == pytest.approx(180, abs=1)                    # 255 x sin(45°) : sol plat, soleil à 45°
