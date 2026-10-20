# SPDX-License-Identifier: GPL-3.0-or-later
"""vexag_curve.py : courbe d'exagération non linéaire (--vexag-curve), PCHIP monotone. Sans GDAL : numpy et scipy.

Les propriétés vérifiées sont celles qui comptent pour un relief : monotonie (le plus haut reste le plus haut),
aucun dépassement entre deux points, continuité aux extrémités, passage exact par les points, puis intégration au
maillage de mesh.py (solide fermé, hauteur imprimée = altitude équivalente x échelle x exagération)."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import mesh  # noqa: E402
import report  # noqa: E402
import vexag_curve  # noqa: E402


# ------------------------------------------------------------------ lecture et validation

def test_parse_accepts_commas_semicolons_and_spaces():
    for text in ("0:0,500:400,2000:1100", "0:0; 500:400; 2000:1100", "0:0 500:400 2000:1100", " 0:0 , 500:400 ,2000:1100,"):
        assert vexag_curve.parse(text).points == [(0.0, 0.0), (500.0, 400.0), (2000.0, 1100.0)]
    assert vexag_curve.parse("-50:-50,1e3:900").points == [(-50.0, -50.0), (1000.0, 900.0)]


@pytest.mark.parametrize("text, needle", [
    ("", "au moins deux points"),
    ("0:0", "au moins deux points"),
    ("0:0,500", "illisible"),
    ("0:0,abc:5", "illisible"),
    ("0,5:0,10:5", "illisible"),                          # virgule décimale : ambiguë, refusée avec un message clair
    ("0:0,0:5", "strictement croissantes"),
    ("0:0,500:400,300:500", "strictement croissantes"),
    ("0:0,500:400,1000:300", "ne peuvent pas décroître"),
    ("0:100,500:100", "relief serait nul"),
    ("0:0,nan:5", "nombres finis"),
    ("0:0,inf:5", "nombres finis"),
])
def test_invalid_curves_are_refused_with_a_clear_message(text, needle):
    with pytest.raises(vexag_curve.CurveError) as e:
        vexag_curve.parse(text)
    assert needle in str(e.value)


# ------------------------------------------------------------------ propriétés de la courbe

def test_identity_curve_is_the_identity_everywhere_including_beyond_the_points():
    c = vexag_curve.parse("0:0,1000:1000")
    x = np.linspace(-500.0, 3000.0, 701)
    assert np.allclose(c(x), x)


def test_passes_exactly_through_its_points_and_is_monotone_without_overshoot():
    rng = np.random.default_rng(1)
    for _ in range(40):
        n = int(rng.integers(2, 8))
        xs = np.sort(rng.uniform(-200, 5000, n)) + np.arange(n)          # strictement croissantes
        ys = np.cumsum(np.abs(rng.normal(300, 200, n)))                  # croissantes
        flat = rng.random(n) < 0.2
        flat[0] = False
        ys[flat] = np.nan                                                # (rempli ci-dessous : paliers)
        for i in range(1, n):
            if np.isnan(ys[i]):
                ys[i] = ys[i - 1]
        if ys[-1] == ys[0]:
            ys[-1] += 1.0
        c = vexag_curve.AltitudeCurve(zip(xs, ys))
        assert np.allclose(c(xs), ys, atol=1e-9)
        fine = np.linspace(xs[0], xs[-1], 4001)
        v = c(fine)
        assert (np.diff(v) >= -1e-9).all(), "courbe non monotone"
        k = np.clip(np.searchsorted(xs, fine, side="right") - 1, 0, n - 2)
        assert (v >= ys[k] - 1e-9).all() and (v <= ys[k + 1] + 1e-9).all(), "dépassement entre deux points"


def test_extension_beyond_the_ends_is_linear_with_the_end_segment_slope_and_continuous():
    c = vexag_curve.parse("100:0,600:300,2100:900")                      # pentes 0,6 puis 0,4
    assert c(100.0) == pytest.approx(0.0) and c(2100.0) == pytest.approx(900.0)
    assert c(0.0) == pytest.approx(-60.0) and c(2600.0) == pytest.approx(900.0 + 0.4 * 500)
    eps = 1e-6
    assert abs(c(100.0 - eps) - c(100.0)) < 1e-5 and abs(c(2100.0 + eps) - c(2100.0)) < 1e-5


def test_flat_segment_clips_the_relief_exactly():
    c = vexag_curve.parse("0:0,1000:1000,2000:1000,3000:2000")
    v = c(np.linspace(1000.0, 2000.0, 101))
    assert np.allclose(v, 1000.0)
    assert c(3000.0) - c(2000.0) == pytest.approx(1000.0)


def test_nan_stays_nan_and_scalars_stay_scalars():
    c = vexag_curve.parse("0:0,100:50")
    out = c(np.array([np.nan, 50.0, np.nan]))
    assert np.isnan(out[0]) and np.isnan(out[2]) and out[1] == pytest.approx(25.0)
    assert isinstance(c(30.0), float) and np.isnan(c(float("nan")))


def test_segments_and_log_lines():
    c = vexag_curve.parse("0:0,500:500,2000:1100")
    assert c.segments() == [(0.0, 500.0, 1.0), (500.0, 2000.0, pytest.approx(0.4))]
    lines = c.describe(vexag=2.0)
    assert lines[0].startswith("ATTENTION : --vexag-curve") and "réplique fidèle" in lines[0]
    assert "0 à 500 m : x2.00" in lines[1] and "500 à 2000 m : x0.80" in lines[1]
    assert c.as_list() == [[0.0, 0.0], [500.0, 500.0], [2000.0, 1100.0]]


# ------------------------------------------------------------------ relief annoncé avec une courbe

def test_relief_stats_with_an_equivalent_relief():
    st = report.relief_stats(0.0, 4000.0, 0.05, 2.0, 100.0, equivalent_relief_m=1800.0)
    assert st["relief_m"] == 4000.0 and st["relief_equivalent_m"] == 1800.0
    assert st["relief_true_mm"] == pytest.approx(90.0) and st["relief_mm"] == pytest.approx(180.0)
    line = report.describe_relief(st)[0]
    assert "dénivelé 4000 m" in line and "dénivelé équivalent 1800 m" in line and "avec la courbe" in line
    assert report.relief_stats(0.0, 4000.0, 0.05, 2.0, 100.0)["relief_equivalent_m"] is None


# ------------------------------------------------------------------ intégration au maillage

def _solid(curve, vexag, seed):
    noise = ndimage.gaussian_filter(np.random.default_rng(seed).normal(size=(44, 44)), 2.2)
    mask = noise > np.quantile(noise, 0.45)
    elev = 200.0 + np.abs(ndimage.gaussian_filter(np.random.default_rng(seed + 50).normal(size=mask.shape), 3)) * 9000.0
    xs = np.arange(mask.shape[1]) - (mask.shape[1] - 1) / 2.0
    ys = (mask.shape[0] - 1) / 2.0 - np.arange(mask.shape[0])
    clean, quad_ok = mesh.clean_mask(mask)
    min_e = float(elev[clean].min())
    scale, base = 0.01, 2.0
    if curve is None:
        e, m = elev, min_e
    else:
        e, m = curve(elev), float(curve(min_e))
    tris = mesh.build_mesh_from_arrays(e, quad_ok, xs, ys, scale, vexag, base, m)[0]
    top = elev[report.touched_nodes(quad_ok)]
    return tris, elev, quad_ok, min_e, float(top.max()), scale


@pytest.mark.parametrize("seed", range(6))
def test_curved_solid_is_closed_oriented_and_has_the_announced_height(seed):
    curve = vexag_curve.parse("0:0,1000:1000,3000:1600,9000:2400")
    tris, elev, quad_ok, min_e, max_e, scale = _solid(curve, 2.0, seed)
    rep = check_stl.edge_report(tris)
    assert rep["boundary_edges"] == rep["nonmanifold_edges"] == rep["duplicate_directed_edges"] == 0
    assert check_stl.signed_volume(tris) > 0
    expected = (curve(max_e) - curve(min_e)) * scale * 2.0
    assert tris[:, :, 2].max() == pytest.approx(expected, rel=1e-12)
    st = report.relief_stats(min_e, max_e, scale, 2.0, 100.0, float(curve(max_e)) - float(curve(min_e)))
    assert st["relief_mm"] == pytest.approx(expected, rel=1e-12)


def test_curve_keeps_the_ranking_of_heights_and_compresses_the_summits():
    curve = vexag_curve.parse("0:0,1000:1000,9000:2000")                 # strictement croissante
    plain, elev, quad_ok, min_e, max_e, scale = _solid(None, 1.0, 3)
    curved = _solid(curve, 1.0, 3)[0]
    assert (np.argsort(elev.ravel(), kind="stable") == np.argsort(curve(elev).ravel(), kind="stable")).all()
    assert curved[:, :, 2].max() < plain[:, :, 2].max()                  # sommets comprimés
    assert not np.allclose(np.sort(plain[:, :, 2].ravel()), np.sort(curved[:, :, 2].ravel()))
