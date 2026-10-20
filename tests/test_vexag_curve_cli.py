# SPDX-License-Identifier: GPL-3.0-or-later
"""Option --vexag-curve de circle_dem_to_stl.py, de bout en bout (main()), sur le cône synthétique de test_report_cli
(île de 80 m ; la lecture GDAL est court-circuitée : seul le câblage est testé)."""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
from test_report_cli import cone, run  # noqa: E402,F401


def test_curve_changes_the_printed_height_and_is_announced(cone, monkeypatch, capsys, tmp_path):
    out = run(monkeypatch, capsys, tmp_path, "--vexag-curve", "0:0,40:40,80:60", "--report")
    # 80 m réels -> 60 m équivalents ; échelle 0,2 mm/m ; exagération 3 : 60 x 0,2 x 3 = 36 mm
    assert "soit 36.0 mm au-dessus du socle" in out and "dénivelé équivalent 60 m" in out
    assert "ATTENTION : --vexag-curve" in out and "0 à 40 m : x3.00 | 40 à 80 m : x1.50" in out
    stl = next(tmp_path.glob("*.stl"))
    assert stl.name.endswith("_x3_courbe.stl")
    tris = check_stl.read_stl(str(stl))
    assert tris[:, :, 2].max() == pytest.approx(36.0, abs=1e-3)
    assert all(ok for _, ok, _ in check_stl.check(tris, 56.0, 3.0, tol_mm=5.0, topology=True))
    rep = json.loads(next(tmp_path.glob("*.json")).read_text(encoding="utf-8"))
    assert rep["model"]["vexag_curve"] == [[0.0, 0.0], [40.0, 40.0], [80.0, 60.0]]
    assert rep["relief"]["relief_m"] == pytest.approx(80.0) and rep["relief"]["relief_equivalent_m"] == pytest.approx(60.0)
    assert rep["relief"]["relief_mm"] == pytest.approx(36.0)


def test_without_the_option_nothing_changes(cone, monkeypatch, capsys, tmp_path):
    out = run(monkeypatch, capsys, tmp_path, "--report")
    assert "vexag-curve" not in out and "soit 48.0 mm au-dessus du socle" in out
    assert not next(tmp_path.glob("*.stl")).name.endswith("_courbe.stl")
    rep = json.loads(next(tmp_path.glob("*.json")).read_text(encoding="utf-8"))
    assert rep["model"]["vexag_curve"] is None and rep["relief"]["relief_equivalent_m"] is None


@pytest.mark.parametrize("bad", ["0:0", "0:0,10", "100:0,50:10", "0:0,100:50,200:20"])
def test_invalid_curve_stops_before_any_work_or_file(cone, monkeypatch, capsys, tmp_path, bad):
    with pytest.raises(SystemExit) as e:
        run(monkeypatch, capsys, tmp_path, "--vexag-curve", bad)
    assert "--vexag-curve" in str(e.value.code) and list(tmp_path.iterdir()) == []
