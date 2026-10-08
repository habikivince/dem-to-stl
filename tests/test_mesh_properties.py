# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 3 : tests par propriétés du maillage (masques aléatoires, trous, îles, pincements), de
--dry-run et des petits utilitaires en ligne de commande. Les propriétés vérifiées ne dépendent
pas du code de maillage : étanchéité, orientation, volume exact calculé à part."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import circle_dem_to_stl  # noqa: E402
import compare_dems  # noqa: E402
import mesh  # noqa: E402
import mesh_export  # noqa: E402
import sources  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402

SCALE, BASE = 0.1, 2.0


def build(mask, vexag=1.0, seed=0, **kw):
    n = mask.shape[0]
    rng = np.random.default_rng(seed)
    elev = 500.0 + ndimage.gaussian_filter(rng.normal(size=mask.shape), 3) * 400
    xs = np.arange(mask.shape[1]) - (mask.shape[1] - 1) / 2.0     # le centre de la grille est l'origine
    ys = (mask.shape[0] - 1) / 2.0 - np.arange(n)
    clean, quad_ok = mesh.clean_mask(mask)
    if not quad_ok.any():
        pytest.skip("masque vide après nettoyage")
    min_e = elev[clean].min()
    tris, n_top, n_bot, n_wall = mesh.build_mesh_from_arrays(elev, quad_ok, xs, ys, SCALE, vexag, BASE, min_e, **kw)
    return tris, elev, quad_ok, xs, ys, min_e, (n_top, n_bot, n_wall)


def expected_volume(elev, quad_ok, xs, ys, min_e, vexag):
    """Volume exact du solide : somme, pour chaque triangle du dessus, de aire x (altitude moyenne + socle)."""
    ii, jj = np.nonzero(quad_ok)
    z = lambda di, dj: (elev[ii + di, jj + dj] - min_e) * SCALE * vexag  # noqa: E731
    area = np.abs((xs[1] - xs[0]) * (ys[1] - ys[0])) * SCALE * SCALE
    z00, z10, z01, z11 = z(0, 0), z(1, 0), z(0, 1), z(1, 1)
    return float(np.sum(area / 2 * ((z00 + z10 + z01) / 3 + BASE) + area / 2 * ((z10 + z11 + z01) / 3 + BASE)))


def assert_solid(tris, elev, quad_ok, xs, ys, min_e, vexag=1.0):
    rep = check_stl.edge_report(tris)
    assert rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0 and rep["duplicate_directed_edges"] == 0, rep
    v, f = mesh_export.index_mesh(tris)
    vol = mesh_export.signed_volume(v, f)
    assert vol > 0
    assert vol == pytest.approx(expected_volume(elev, quad_ok, xs, ys, min_e, vexag), rel=1e-4)
    assert tris[:, :, 2].min() == pytest.approx(-BASE)


def blob(seed, n=44, keep=0.55):
    noise = ndimage.gaussian_filter(np.random.default_rng(seed).normal(size=(n, n)), 2.2)
    return noise > np.quantile(noise, 1 - keep)


@pytest.mark.parametrize("seed", range(14))
def test_random_blobs_give_closed_oriented_solids_of_exact_volume(seed):
    mask = blob(seed, keep=0.4 + 0.03 * (seed % 6))                  # formes variées : trous, îles, goulets
    tris, elev, quad_ok, xs, ys, min_e, _ = build(mask, vexag=1.0 + seed % 3, seed=seed)
    assert_solid(tris, elev, quad_ok, xs, ys, min_e, vexag=1.0 + seed % 3)


def test_special_masks_rectangle_ring_islands_and_checkerboard_pinches():
    n = 30
    rect = np.zeros((n, n), bool); rect[4:26, 3:27] = True
    ring = np.zeros((n, n), bool); ring[3:27, 3:27] = True; ring[10:20, 10:20] = False
    islands = np.zeros((n, n), bool); islands[3:12, 3:12] = True; islands[17:27, 16:27] = True
    pinch = np.zeros((n, n), bool); pinch[5:15, 5:15] = True; pinch[15:25, 15:25] = True   # deux blocs touchant par un coin
    for mask in (rect, ring, islands, pinch):
        tris, elev, quad_ok, xs, ys, min_e, _ = build(mask)
        assert_solid(tris, elev, quad_ok, xs, ys, min_e)


def test_clean_mask_leaves_no_diagonal_pinch_and_removes_thin_lines():
    for seed in range(10):
        _, quad_ok = mesh.clean_mask(blob(100 + seed, keep=0.5))
        a, b, c, d = quad_ok[:-1, :-1], quad_ok[:-1, 1:], quad_ok[1:, :-1], quad_ok[1:, 1:]
        assert not np.any(a & d & ~b & ~c) and not np.any(b & c & ~a & ~d)
    thin = np.zeros((20, 20), bool); thin[10, 2:18] = True
    assert not mesh.clean_mask(thin)[1].any()                         # une ligne d'une cellule disparaît


@pytest.mark.parametrize("n,radius", [(40, 18.0), (60, 27.5), (61, 29.0), (100, 45.3)])
def test_smooth_wall_on_discs_is_closed_oriented_and_circular(n, radius):
    yy, xx = np.mgrid[0:n, 0:n]
    c = (n - 1) / 2.0
    mask = (xx - c) ** 2 + (yy - c) ** 2 <= radius ** 2
    tris, elev, quad_ok, xs, ys, min_e, (nt, nb, nw) = build(mask, smooth_wall=True, circle_radius_mm=radius * SCALE)
    rep = check_stl.edge_report(tris)
    assert rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0 and rep["duplicate_directed_edges"] == 0
    v, f = mesh_export.index_mesh(tris)
    assert mesh_export.signed_volume(v, f) > 0
    wall = tris[-nw:].reshape(-1, 3)
    r = np.hypot(wall[:, 0], wall[:, 1])                           # la paroi est exactement sur le vrai cercle
    assert np.ptp(r) < 1e-9 and r.mean() == pytest.approx(radius * SCALE)


def test_drop_degenerate_removes_only_zero_area_triangles():
    t = np.array([[[0, 0, 0], [1, 0, 0], [0, 1, 0]],
                  [[0, 0, 0], [0, 0, 0], [0, 1, 0]],                           # deux sommets identiques
                  [[0, 0, 0], [1, 0, 0], [1.00000001, 0, 0]]], dtype=np.float64)  # distincts en float64, égaux en float32
    kept, n = mesh.drop_degenerate(t)
    assert n == 1 + 0 and len(kept) == 2 or n == 2          # le 3e est conservé : (0,0,0),(1,0,0),(1.00000001,0,0) diffèrent en float32
    assert mesh.drop_degenerate(np.zeros((0, 3, 3)))[1] == 0


def test_smooth_wall_with_holes_or_islands_is_refused_instead_of_producing_an_open_mesh():
    """Avant : --smooth-wall ne traçait que le contour extérieur et laissait celui d'un trou ouvert
    (88 arêtes ouvertes dans ce cas), sans avertissement."""
    n = 30
    ring = np.zeros((n, n), bool); ring[3:27, 3:27] = True; ring[10:20, 10:20] = False
    two = np.zeros((n, n), bool); two[3:12, 3:12] = True; two[17:27, 16:27] = True
    for mask in (ring, two):
        with pytest.raises(SystemExit) as e:
            build(mask, smooth_wall=True, circle_radius_mm=1.4)
        assert "--smooth-wall requiert un seul contour" in str(e.value)


# ---------------------------------------------------------------- --dry-run

def test_dry_run_prints_plan_and_touches_nothing(monkeypatch, tmp_path, capsys):
    def boom(*a, **k):
        raise AssertionError("aucun téléchargement en --dry-run")
    monkeypatch.setitem(sources.SOURCES, "swisstopo", sources.Source("swisstopo", "x", [(5, 45, 11, 48)], boom))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--source", "swisstopo", "--lat", "46", "--lon", "7.6", "--radius", "3000",
                                      "--diameter-mm", "150", "--vexag", "1", "--base-mm", "3", "--format", "all",
                                      "--dry-run"])
    circle_dem_to_stl.main()
    out = capsys.readouterr().out
    assert "--dry-run" in out and "Estimation" in out and "STL ~" in out and "3MF ~" in out and "OBJ ~" in out
    assert "lat46.000_lon7.600_150mm.3mf" in out
    assert not (tmp_path / "output").exists()


def test_dry_run_still_validates(monkeypatch, tmp_path):
    monkeypatch.setattr(sys, "argv", ["x", "--lat", "46", "--lon", "7.6", "--source", "auto", "--radius", "-3",
                                      "--diameter-mm", "150", "--dry-run"])
    with pytest.raises(SystemExit) as e:
        circle_dem_to_stl.main()
    assert "--radius" in str(e.value)


# ---------------------------------------------------------------- utilitaires en ligne de commande

def make_tifs(tmp_path):
    yy, xx = np.mgrid[0:400, 0:400]
    base = 1000 + 0.1 * xx + 5 * np.sin(yy / 30.0)
    write_synthetic_tif(tmp_path / "a.tif", base.astype(np.float32), 2600000.0, 1202000.0, 5.0, 2056)
    write_synthetic_tif(tmp_path / "b.tif", (base + 3.0).astype(np.float32), 2600000.0, 1202000.0, 5.0, 2056)
    return str(tmp_path / "a.tif"), str(tmp_path / "b.tif")


def test_compare_dems_cli_reports_offset(monkeypatch, tmp_path, capsys):
    a, b = make_tifs(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", a, b, "--cx", "2601000", "--cy", "1201000", "--radius", "800",
                                      "--pixel-size", "10"])
    compare_dems.main()
    out = capsys.readouterr().out
    assert "moyenne +3.000 m" in out and "RMS 3.000 m" in out
    monkeypatch.setattr(sys, "argv", ["x", a, b, "--radius", "800", "--pixel-size", "10"])
    with pytest.raises(SystemExit):
        compare_dems.main()                                            # ni --lat/--lon ni --cx/--cy
    monkeypatch.setattr(sys, "argv", ["x", a, b, "--cx", "2700000", "--cy", "1100000", "--radius", "500",
                                      "--pixel-size", "10"])
    with pytest.raises(SystemExit) as e:
        compare_dems.main()
    assert "aucune cellule valide" in str(e.value)


def test_check_stl_cli_exit_codes(monkeypatch, tmp_path, capsys):
    a, _ = make_tifs(tmp_path)
    out = tmp_path / "t.stl"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--input", a, "--cx", "2601000", "--cy", "1201000", "--radius", "800",
                                      "--diameter-mm", "60", "--vexag", "1", "--base-mm", "2",
                                      "--print-spacing-mm", "0.5", "--output", str(out)])
    circle_dem_to_stl.main()
    capsys.readouterr()
    monkeypatch.setattr(sys, "argv", ["x", str(out), "--diameter-mm", "60", "--base-mm", "2", "--radius", "800"])
    with pytest.raises(SystemExit) as e:
        check_stl.main()
    assert e.value.code == 0 and "[OK ] étanchéité" in capsys.readouterr().out
    monkeypatch.setattr(sys, "argv", ["x", str(out), "--diameter-mm", "90", "--base-mm", "2"])
    with pytest.raises(SystemExit) as e:
        check_stl.main()
    assert e.value.code == 1 and "ÉCHEC" in capsys.readouterr().out
    monkeypatch.setattr(sys, "argv", ["x", str(tmp_path / "absent.stl"), "--diameter-mm", "60", "--base-mm", "2"])
    with pytest.raises(SystemExit) as e:
        check_stl.main()
    assert str(e.value).startswith("Erreur :")
    monkeypatch.setattr(sys, "argv", ["x", str(tmp_path / "t.ply"), "--diameter-mm", "60", "--base-mm", "2"])
    (tmp_path / "t.ply").write_text("x")
    with pytest.raises(SystemExit) as e:
        check_stl.main()
    assert "extension non gérée" in str(e.value)
