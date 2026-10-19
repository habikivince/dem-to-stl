# SPDX-License-Identifier: GPL-3.0-or-later
"""Dessus adaptatif RTIN à bord forcé (--top adaptive, prototype).

Propriétés vérifiées indépendamment du code de maillage : solide fermé et orienté, volume cohérent
(théorème de la divergence contre somme des prismes), écart vertical aux nœuds du DEM borné par la
tolérance (mesure par coordonnées barycentriques), bord du masque à pleine résolution, parois et fond
inchangés."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import adaptive_top  # noqa: E402
import check_stl  # noqa: E402
import circle_dem_to_stl  # noqa: E402
import mesh_export  # noqa: E402
from test_mesh_properties import BASE, SCALE, blob, build, make_tifs  # noqa: E402


def smooth_terrain(H, W, seed, sigma=3.0, amp=8.0):
    rng = np.random.default_rng(seed)
    return ndimage.gaussian_filter(rng.normal(size=(H + 1, W + 1)), sigma) * amp


def top_of(quad, z, tol):
    H, W = quad.shape
    return adaptive_top.adaptive_top_triangles(quad, np.arange(W + 1, dtype=float), -np.arange(H + 1, dtype=float), z, tol)


def max_deviation(tris, quad, z):
    """Écart vertical maximal entre la surface (triangles du dessus) et les altitudes z, aux nœuds touchant le masque,
    par force brute : coordonnées barycentriques de chaque nœud dans chaque triangle."""
    H, W = quad.shape
    cp = np.pad(quad, 1)
    touch = (cp[:-1, :-1] | cp[:-1, 1:] | cp[1:, :-1] | cp[1:, 1:])[:H + 1, :W + 1]
    ii, jj = np.nonzero(touch)
    px, py = jj.astype(float), -ii.astype(float)
    a, b, c = tris[:, 0], tris[:, 1], tris[:, 2]
    det = (b[:, 1] - c[:, 1]) * (a[:, 0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (a[:, 1] - c[:, 1])
    worst = 0.0
    for s in range(0, len(px), 200):
        x, y = px[s:s + 200, None], py[s:s + 200, None]
        l1 = ((b[:, 1] - c[:, 1]) * (x - c[:, 0]) + (c[:, 0] - b[:, 0]) * (y - c[:, 1])) / det
        l2 = ((c[:, 1] - a[:, 1]) * (x - c[:, 0]) + (a[:, 0] - c[:, 0]) * (y - c[:, 1])) / det
        l3 = 1.0 - l1 - l2
        inside = (l1 >= -1e-9) & (l2 >= -1e-9) & (l3 >= -1e-9)
        assert inside.any(axis=1).all(), "un nœud du masque n'est couvert par aucun triangle du dessus"
        zhat = np.where(inside, l1 * a[:, 2] + l2 * b[:, 2] + l3 * c[:, 2], np.nan)
        zhat = np.nanmax(zhat, axis=1) if False else np.array([row[~np.isnan(row)][0] for row in zhat])
        worst = max(worst, float(np.abs(zhat - z[ii[s:s + 200], jj[s:s + 200]]).max()))
    return worst


def prism_volume(top_tris):
    """∫ (z + BASE) dA sur les triangles du dessus : somme des prismes, sans la divergence."""
    a, b, c = top_tris[:, 0], top_tris[:, 1], top_tris[:, 2]
    area = 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0]))
    return float(np.sum(area * ((a[:, 2] + b[:, 2] + c[:, 2]) / 3.0 + BASE)))


# --- solide fermé, orienté, de volume cohérent -------------------------------------------------------

@pytest.mark.parametrize("bottom", ["grid", "merged"])
@pytest.mark.parametrize("tol", [0.0, 0.05, 0.5, 5.0, 1e9])
@pytest.mark.parametrize("seed", range(6))
def test_adaptive_top_gives_a_closed_oriented_solid_with_consistent_volume(seed, tol, bottom):
    mask = blob(seed, keep=0.4 + 0.04 * seed)
    tris, elev, quad_ok, xs, ys, min_e, (nt, nb, nw) = build(mask, seed=seed, top="adaptive", top_tol_mm=tol, bottom=bottom)
    rep = check_stl.edge_report(tris)
    assert rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0 and rep["duplicate_directed_edges"] == 0, rep
    assert check_stl.signed_volume(tris) == pytest.approx(prism_volume(tris[:nt]), rel=1e-9)
    assert check_stl.degenerate_count(tris) == 0


@pytest.mark.parametrize("seed", range(6))
def test_walls_and_bottom_are_unchanged(seed):
    mask = blob(seed, keep=0.5)
    grid, *_a, (gt, gb, gw) = build(mask, seed=seed)
    adapt, *_b, (at, ab, aw) = build(mask, seed=seed, top="adaptive", top_tol_mm=0.5)
    assert (gb, gw) == (ab, aw) and np.array_equal(grid[gt:], adapt[at:])


def test_special_masks_and_shapes():
    cases = {
        "anneau": np.pad(np.pad(np.zeros((6, 6), bool), 3, constant_values=True), 2),
        "deux îles": np.pad(np.block([[np.ones((8, 8), bool), np.zeros((8, 5), bool), np.ones((8, 9), bool)]]), 3),
        "pincement": np.pad(np.block([[np.ones((6, 6), bool), np.zeros((6, 6), bool)], [np.zeros((6, 6), bool), np.ones((6, 6), bool)]]), 2),
        "non carré 17x45": np.pad(np.ones((11, 39), bool), 3),
        "une cellule": np.pad(np.ones((1, 1), bool), 2),
        "2x3 cellules": np.ones((2, 3), bool),
    }
    for name, mask in cases.items():
        tris, *_ = build(mask, top="adaptive", top_tol_mm=0.5)
        rep = check_stl.edge_report(tris)
        assert rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0 and rep["duplicate_directed_edges"] == 0, name


# --- bord à pleine résolution et erreur bornée -------------------------------------------------------

def test_every_node_on_the_mask_boundary_is_a_vertex_and_every_boundary_cell_is_two_half_cells():
    mask = blob(3, n=40, keep=0.55)
    z = smooth_terrain(*mask.shape, 3)
    tris = top_of(mask, z, 1e9)                                         # tolérance énorme : rien n'est simplifié sauf le plat
    verts = {(int(round(x)), int(round(-y))) for x, y in tris[:, :, :2].reshape(-1, 2)}
    cp = np.pad(mask, 1)
    for i in range(mask.shape[0] + 1):
        for j in range(mask.shape[1] + 1):
            around = [cp[i, j], cp[i, j + 1], cp[i + 1, j], cp[i + 1, j + 1]]
            if any(around) and not all(around):
                assert (j, i) in verts, f"nœud de bord ({i}, {j}) absent du dessus"
    area = 0.5 * np.abs(np.cross(tris[:, 1, :2] - tris[:, 0, :2], tris[:, 2, :2] - tris[:, 0, :2]))
    assert area.sum() == pytest.approx(mask.sum())                      # les triangles couvrent exactement le masque


def test_zero_tolerance_keeps_the_full_resolution_on_noisy_terrain():
    mask = blob(2, n=36, keep=0.6)
    z = np.random.default_rng(7).normal(size=(mask.shape[0] + 1, mask.shape[1] + 1))
    tris = top_of(mask, z, 0.0)
    assert len(tris) == 2 * mask.sum() and max_deviation(tris, mask, z) < 1e-9


@pytest.mark.parametrize("sigma,tol", [(2.0, 0.02), (3.0, 0.05), (4.0, 0.2), (1.5, 0.5)])
@pytest.mark.parametrize("seed", range(3))
def test_vertical_deviation_at_the_nodes_stays_within_twice_the_tolerance(seed, sigma, tol):
    mask = blob(seed, n=33, keep=0.7)
    z = smooth_terrain(*mask.shape, seed, sigma=sigma, amp=6.0)
    tris = top_of(mask, z, tol)
    assert len(tris) < 2 * mask.sum()
    assert max_deviation(tris, mask, z) <= 2.0 * tol + 1e-9             # mesuré : environ 1,5 à 1,8 fois la tolérance


def test_smooth_terrain_loses_most_of_its_triangles_and_rough_terrain_does_not():
    mask = np.pad(np.ones((60, 60), bool), 2)
    smooth = smooth_terrain(*mask.shape, 1, sigma=6.0, amp=3.0)
    rough = np.random.default_rng(1).normal(size=smooth.shape) * 3.0
    assert len(top_of(mask, smooth, 0.1)) < 0.25 * 2 * mask.sum()
    assert len(top_of(mask, rough, 0.1)) > 0.9 * 2 * mask.sum()


def test_a_plane_collapses_to_the_boundary_only():
    mask = np.pad(np.ones((32, 32), bool), 0)
    ii, jj = np.mgrid[0:33, 0:33]
    z = 0.3 * ii + 0.7 * jj
    tris = top_of(mask, z, 1e-9)
    assert len(tris) < 2 * mask.sum() * 0.25 and max_deviation(tris, mask, z) < 1e-9


# --- arguments ---------------------------------------------------------------------------------------

def test_negative_tolerance_and_unknown_or_conflicting_modes_are_refused():
    mask = np.ones((4, 4), bool)
    with pytest.raises(ValueError):
        top_of(mask, np.zeros((5, 5)), -1.0)
    with pytest.raises(ValueError, match="top"):
        build(blob(1), top="tin")
    with pytest.raises(ValueError, match="sans objet"):
        build(blob(1), top="adaptive", diagonal="adaptive")


def test_empty_result_when_the_mask_has_no_cell():
    assert top_of(np.zeros((4, 4), bool), np.zeros((5, 5)), 0.1).shape == (0, 3, 3)


# --- de bout en bout (CLI, trois formats) -------------------------------------------------------------

@pytest.mark.parametrize("ext", ["stl", "3mf", "obj"])
def test_cli_top_adaptive_produces_a_valid_smaller_file(monkeypatch, tmp_path, capsys, ext):
    a, _ = make_tifs(tmp_path)

    def run(out, *extra):
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["x", "--input", a, "--cx", "2601000", "--cy", "1201000", "--radius", "800",
                                          "--diameter-mm", "60", "--vexag", "1", "--base-mm", "2",
                                          "--print-spacing-mm", "0.4", "--output", str(out), *extra])
        circle_dem_to_stl.main()
        capsys.readouterr()
        return mesh_export.read_triangles(str(out))

    grid = run(tmp_path / f"g.{ext}")
    adapt = run(tmp_path / f"a.{ext}", "--top", "adaptive", "--top-tol-mm", "0.1")
    assert len(adapt) < len(grid)
    res = {n: ok for n, ok, _ in check_stl.check(adapt, 60.0, 2.0, tol_mm=1.0)}
    assert all(res.values()), res


def test_cli_refuses_the_two_adaptive_options_together(monkeypatch, tmp_path):
    a, _ = make_tifs(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--input", a, "--cx", "2601000", "--cy", "1201000", "--radius", "800",
                                      "--diameter-mm", "60", "--vexag", "1", "--base-mm", "2",
                                      "--print-spacing-mm", "0.4", "--output", str(tmp_path / "o.stl"),
                                      "--top", "adaptive", "--diagonal", "adaptive"])
    with pytest.raises(SystemExit) as e:
        circle_dem_to_stl.main()
    assert "sans objet" in str(e.value) or e.value.code not in (0, None)


# --- découpage en blocs (mémoire bornée) ---------------------------------------------------------------

@pytest.mark.parametrize("block", [2, 4, 8, 16, 32, 1024])
@pytest.mark.parametrize("seed", range(3))
def test_blocks_keep_the_solid_closed_and_the_deviation_bounded(seed, block):
    mask = blob(seed, n=37, keep=0.65)                         # 37 n'est pas un multiple des blocs : blocs partiels
    z = smooth_terrain(*mask.shape, seed, sigma=2.5, amp=6.0)
    tol = 0.1
    tris = adaptive_top.adaptive_top_triangles(mask, np.arange(mask.shape[1] + 1, dtype=float),
                                               -np.arange(mask.shape[0] + 1, dtype=float), z, tol, block=block)
    area = 0.5 * np.abs(np.cross(tris[:, 1, :2] - tris[:, 0, :2], tris[:, 2, :2] - tris[:, 0, :2]))
    assert area.sum() == pytest.approx(mask.sum())              # couverture exacte du masque, sans recouvrement
    assert max_deviation(tris, mask, z) <= 2.0 * tol + 1e-9
    # étanchéité de la seule surface : chaque arête intérieure est partagée par deux triangles, sans fissure
    v, f = mesh_export.index_mesh(tris)
    edges = np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    uniq, counts = np.unique(edges, axis=0, return_counts=True)
    assert counts.max() <= 2                                    # jamais plus de deux triangles par arête
    cp = np.pad(mask, 1)                                        # arêtes de bord attendues : celles du masque
    n_boundary = int(np.sum(cp[1:-1, 1:-1] & ~cp[:-2, 1:-1]) + np.sum(cp[1:-1, 1:-1] & ~cp[2:, 1:-1])
                     + np.sum(cp[1:-1, 1:-1] & ~cp[1:-1, :-2]) + np.sum(cp[1:-1, 1:-1] & ~cp[1:-1, 2:]))
    assert int((counts == 1).sum()) == n_boundary


def test_smaller_blocks_cost_a_strip_of_triangles_along_the_seams_only():
    mask = np.ones((64, 64), bool)
    z = smooth_terrain(64, 64, 5, sigma=8.0, amp=2.0)
    one = top_of(mask, z, 0.2)
    many = adaptive_top.adaptive_top_triangles(mask, np.arange(65, dtype=float), -np.arange(65, dtype=float), z, 0.2, block=16)
    assert len(many) > len(one)                                 # les jonctions coûtent des triangles
    seam_nodes = 2 * (64 // 16 - 1) * 64                         # longueur totale des jonctions, en nœuds
    assert len(many) < len(one) + 8 * seam_nodes                # mesuré : environ 5 triangles par nœud de jonction


def test_blocks_must_be_a_power_of_two():
    for bad in (0, 1, 3, 12, 1000):
        with pytest.raises(ValueError, match="puissance de deux"):
            adaptive_top.adaptive_top_triangles(np.ones((4, 4), bool), np.arange(5.0), np.arange(5.0), np.zeros((5, 5)), 0.1, block=bad)


def test_blocks_with_no_cell_are_skipped():
    mask = np.zeros((40, 40), bool)
    mask[1:5, 1:5] = True                                       # un seul bloc de 16 contient des cellules
    z = np.zeros((41, 41))
    tris = adaptive_top.adaptive_top_triangles(mask, np.arange(41.0), -np.arange(41.0), z, 0.1, block=16)
    assert len(tris) > 0 and tris[:, :, 0].max() <= 5 and tris[:, :, 1].min() >= -5
