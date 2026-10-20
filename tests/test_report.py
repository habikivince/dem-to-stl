# SPDX-License-Identifier: GPL-3.0-or-later
"""report.py : hauteur de relief imprimé, repère d'exagération (indicatif) et rapport JSON.

N'a besoin ni de GDAL ni de fichiers (numpy et scipy). Le relief annoncé est confronté à la hauteur réelle du
maillage produit par mesh.py, qui est la seule vérité : un chiffre annoncé qui différerait du STL serait pire que
pas de chiffre."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import mesh  # noqa: E402
import report  # noqa: E402


def test_relief_stats_with_and_without_exaggeration():
    st = report.relief_stats(100.0, 2347.0, scale_mm_per_m=0.025, vexag=2.0, diameter_mm=150.0)
    assert st["relief_m"] == pytest.approx(2247.0)
    assert st["relief_true_mm"] == pytest.approx(2247.0 * 0.025)
    assert st["relief_mm"] == pytest.approx(2247.0 * 0.025 * 2.0)
    assert st["diameter_percent"] == pytest.approx(100.0 * 2247.0 * 0.025 * 2.0 / 150.0)
    assert report.relief_stats(50.0, 50.0, 0.1, 3.0, 100.0)["relief_mm"] == 0.0
    assert report.relief_stats(50.0, 40.0, 0.1, 3.0, 100.0)["relief_m"] == 0.0        # jamais négatif


def test_suggested_exaggeration_is_bounded_and_none_for_flat_relief():
    assert report.suggest_vexag(0.0, 100.0) is None and report.suggest_vexag(5.0, 0.0) is None
    assert report.suggest_vexag(100.0, 200.0) == 1.0                    # déjà grand (45 mm visés) : jamais sous 1
    assert report.suggest_vexag(0.001, 200.0) == 8.0                    # presque plat : plafonné à 8
    assert report.suggest_vexag(15.0, 200.0) == pytest.approx(3.0)      # 0,225 x 200 / 15


def test_description_gives_the_fact_and_a_hint_only_when_it_differs_clearly():
    plat = report.relief_stats(0.0, 10.0, 0.01, 1.0, 100.0)             # 0,1 mm : très plat à vexag 1
    lines = report.describe_relief(plat)
    assert lines[0].startswith("relief : altitudes de 0.0 à 10.0 m") and "0.1 mm au-dessus du socle" in lines[0]
    assert len(lines) == 2 and "--vexag 8.0" in lines[1] and "indicatif" in lines[1]
    ok = report.relief_stats(0.0, 1000.0, 0.0225, 1.0, 100.0)           # 22,5 mm = le repère : pas de conseil
    assert len(report.describe_relief(ok)) == 1
    assert report.describe_relief(report.relief_stats(7.0, 7.0, 0.1, 1.0, 100.0))[0].startswith("relief : nul")
    haut = report.relief_stats(0.0, 8000.0, 0.05, 1.0, 100.0)           # 400 mm à l'échelle vraie, 1 est déjà le minimum
    assert len(report.describe_relief(haut)) == 1


def test_touched_nodes_are_the_corners_of_included_cells():
    q = np.zeros((4, 5), dtype=bool)
    q[1, 2] = True
    t = report.touched_nodes(q)
    assert t.shape == (5, 6)
    assert {tuple(x) for x in np.argwhere(t)} == {(1, 2), (1, 3), (2, 2), (2, 3)}


@pytest.mark.parametrize("seed", range(8))
@pytest.mark.parametrize("vexag", [1.0, 3.0])
def test_announced_relief_equals_the_height_of_the_mesh(seed, vexag):
    noise = ndimage.gaussian_filter(np.random.default_rng(seed).normal(size=(44, 44)), 2.2)
    mask = noise > np.quantile(noise, 1 - (0.4 + 0.03 * (seed % 6)))
    elev = 500.0 + ndimage.gaussian_filter(np.random.default_rng(seed + 100).normal(size=mask.shape), 3) * 40.0
    xs = np.arange(mask.shape[1]) - (mask.shape[1] - 1) / 2.0
    ys = (mask.shape[0] - 1) / 2.0 - np.arange(mask.shape[0])
    clean, quad_ok = mesh.clean_mask(mask)
    scale, base = 0.37, 2.0
    min_e = elev[clean].min()
    tris = mesh.build_mesh_from_arrays(elev, quad_ok, xs, ys, scale, vexag, base, min_e)[0]
    st = report.relief_stats(min_e, elev[report.touched_nodes(quad_ok)].max(), scale, vexag, 100.0)
    assert tris[:, :, 2].max() == pytest.approx(st["relief_mm"], rel=1e-12)


def test_report_path():
    paths = {"stl": "/tmp/x/cervin_150mm.stl", "3mf": "/tmp/x/cervin_150mm.3mf"}
    assert report.report_path(paths, "auto") == "/tmp/x/cervin_150mm.json"
    assert report.report_path(paths, "/tmp/autre/rapport.json") == "/tmp/autre/rapport.json"


def test_report_is_valid_json_with_numpy_values_and_non_finite_numbers(tmp_path):
    relief = report.relief_stats(np.float32(1.5), np.float64(99.5), 0.1, 2.0, 100.0)
    data = report.build_report(
        outputs={"stl": str(tmp_path / "a.stl")}, model={"diameter_mm": np.float64(100.0), "grid": np.array([3, 4])},
        relief=relief, source={"used": None, "attribution": "© swisstopo"},
        mesh_info={"triangles_total": np.int64(12)}, location={"lat": float("nan"), "radius_m": 1000.0},
        options={"top": "grid"}, argv=["--lat", "46.0"])
    path = tmp_path / "a.json"
    report.write_report(path, data)
    back = json.loads(path.read_text(encoding="utf-8"))
    assert back["schema"] == report.SCHEMA_VERSION and back["argv"] == ["--lat", "46.0"]
    assert back["location"]["lat"] is None and back["model"]["grid"] == [3, 4]
    assert back["mesh"]["triangles_total"] == 12 and back["source"]["attribution"] == "© swisstopo"
    assert back["created_utc"].endswith("Z")
    assert back["relief"]["relief_m"] == pytest.approx(98.0)
    assert not list(tmp_path.glob("*.part"))


def test_failed_write_leaves_no_partial_file_and_keeps_the_previous_report(tmp_path):
    path = tmp_path / "r.json"
    report.write_report(path, {"ok": 1})
    with pytest.raises(TypeError):
        report.write_report(path, {"bad": object()})
    assert json.loads(path.read_text(encoding="utf-8")) == {"ok": 1}
    assert not list(tmp_path.glob("*.part"))
