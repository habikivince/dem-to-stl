# SPDX-License-Identifier: GPL-3.0-or-later
"""Fond fusionné (--bottom merged / merged-rows, prototype) : mêmes propriétés que le fond cellule par cellule.

Les propriétés vérifiées ne dépendent pas du code de maillage : étanchéité sans jonction en T, orientation,
volume exact calculé à part, et identité de la géométrie avec le fond « grid »."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import bottom as bottom_mod  # noqa: E402
import check_stl  # noqa: E402
import circle_dem_to_stl  # noqa: E402
import mesh_export  # noqa: E402
from test_mesh_properties import BASE, SCALE, assert_solid, blob, build, make_tifs  # noqa: E402

MERGED_MODES = ["merged", "merged-rows"]


def _cells(rects, shape):
    out = np.zeros(shape, int)
    for i0, i1, j0, j1 in rects:
        out[i0:i1 + 1, j0:j1 + 1] += 1
    return out


def _disc(n, radius):
    yy, xx = np.mgrid[0:n, 0:n]
    c = (n - 1) / 2.0
    return (xx - c) ** 2 + (yy - c) ** 2 <= radius ** 2


def _aspects(t):
    e = np.stack([np.linalg.norm(t[:, 1] - t[:, 0], axis=1), np.linalg.norm(t[:, 2] - t[:, 1], axis=1),
                  np.linalg.norm(t[:, 0] - t[:, 2], axis=1)])
    area = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2
    return e.max(axis=0) ** 2 / (2 * area)            # 1,15 : équilatéral ; 2 : demi-carré ; grand : aiguille


# --- couverture par rectangles gloutons (merged-rows) -----------------------------------------------

@pytest.mark.parametrize("seed", range(10))
@pytest.mark.parametrize("cap", [None, 1, 3, 16])
def test_rectangles_cover_the_mask_exactly_without_overlap_and_respect_the_cap(seed, cap):
    q = blob(seed, n=40, keep=0.35 + 0.05 * (seed % 6))
    rects = bottom_mod.cover_with_rectangles(q, cap)
    cover = _cells(rects, q.shape)
    assert cover.max() <= 1 and (cover == 1).tolist() == q.tolist()
    if cap is not None:
        assert all(i1 - i0 + 1 <= cap and j1 - j0 + 1 <= cap for i0, i1, j0, j1 in rects)


def test_cover_rejects_invalid_cap_and_handles_empty_mask():
    with pytest.raises(ValueError):
        bottom_mod.cover_with_rectangles(np.ones((3, 3), bool), 0)
    xs = np.arange(5.0)
    assert bottom_mod.cover_with_rectangles(np.zeros((4, 4), bool)) == []
    assert bottom_mod.cover_with_graded_squares(np.zeros((4, 4), bool)) == []
    for cover in bottom_mod.COVERS:
        assert bottom_mod.merged_bottom_triangles(np.zeros((4, 4), bool), xs, xs, -1.0, cover=cover).shape == (0, 3, 3)
    with pytest.raises(ValueError, match="cover"):
        bottom_mod.merged_bottom_triangles(np.ones((2, 2), bool), np.arange(3.0), np.arange(3.0), 0.0, cover="x")


def test_full_rectangle_is_a_single_rectangle_for_rows_cover():
    q = np.ones((12, 20), bool)
    assert bottom_mod.cover_with_rectangles(q) == [(0, 11, 0, 19)]
    xs, ys = np.arange(21.0), -np.arange(13.0)
    t = bottom_mod.merged_bottom_triangles(q, xs, ys, -2.0, cover="rows")
    assert len(t) == 2 * (12 + 20)          # contour : tous les nœuds du bord (parois), éventail de 64 triangles


# --- couverture par carrés équilibrés (merged) ------------------------------------------------------

@pytest.mark.parametrize("seed", range(10))
def test_graded_squares_cover_exactly_are_aligned_powers_of_two_and_balanced(seed):
    q = blob(seed, n=44, keep=0.35 + 0.05 * (seed % 6))
    rects = bottom_mod.cover_with_graded_squares(q)
    cover = _cells(rects, q.shape)
    assert cover.max() <= 1 and (cover == 1).tolist() == q.tolist()
    size = np.zeros(q.shape, int)
    for i0, i1, j0, j1 in rects:
        s = i1 - i0 + 1
        assert s == j1 - j0 + 1 and s & (s - 1) == 0 and i0 % s == 0 and j0 % s == 0      # carré 2^k aligné
        size[i0:i1 + 1, j0:j1 + 1] = s
    for a, b in ((size[:-1, :], size[1:, :]), (size[:, :-1], size[:, 1:])):                # voisins par un côté
        both = (a > 0) & (b > 0)
        if both.any():
            assert (np.maximum(a, b)[both] / np.minimum(a, b)[both]).max() <= 2.0


def test_graded_cover_of_a_full_square_is_a_single_leaf_and_boundary_cells_stay_unit_squares():
    full = np.ones((16, 16), bool)
    assert bottom_mod.cover_with_graded_squares(full) == [(0, 15, 0, 15)]
    q = np.zeros((16, 16), bool); q[1:15, 1:15] = True
    sizes = {i1 - i0 + 1 for i0, i1, _, _ in bottom_mod.cover_with_graded_squares(q)}
    assert 1 in sizes and max(sizes) >= 4


# --- solide fermé, orienté, de volume exact ---------------------------------------------------------

@pytest.mark.parametrize("bottom", MERGED_MODES)
@pytest.mark.parametrize("seed", range(14))
def test_random_blobs_with_merged_bottom_are_closed_oriented_solids_of_exact_volume(bottom, seed):
    mask = blob(seed, keep=0.4 + 0.03 * (seed % 6))
    vexag = 1.0 + seed % 3
    cap = 4 if bottom == "merged-rows" and seed % 2 else None
    tris, elev, quad_ok, xs, ys, min_e, _ = build(mask, vexag=vexag, seed=seed, bottom=bottom, bottom_max_cells=cap)
    assert_solid(tris, elev, quad_ok, xs, ys, min_e, vexag=vexag)


@pytest.mark.parametrize("bottom", MERGED_MODES)
def test_special_masks_with_merged_bottom(bottom):
    n = 30
    rect = np.zeros((n, n), bool); rect[4:26, 3:27] = True
    ring = np.zeros((n, n), bool); ring[3:27, 3:27] = True; ring[10:20, 10:20] = False
    islands = np.zeros((n, n), bool); islands[3:12, 3:12] = True; islands[17:27, 16:27] = True
    pinch = np.zeros((n, n), bool); pinch[5:15, 5:15] = True; pinch[15:25, 15:25] = True
    for mask in (rect, ring, islands, pinch):
        tris, elev, quad_ok, xs, ys, min_e, _ = build(mask, bottom=bottom)
        assert_solid(tris, elev, quad_ok, xs, ys, min_e)


@pytest.mark.parametrize("bottom", MERGED_MODES)
@pytest.mark.parametrize("seed", range(6))
def test_merged_and_grid_bottoms_have_the_same_geometry(bottom, seed):
    mask = blob(seed, keep=0.5)
    g, *_rest, (gt, gb, gw) = build(mask, seed=seed)
    m, *_rest2, (mt, mb, mw) = build(mask, seed=seed, bottom=bottom)
    assert gt == mt and gw == mw                                              # dessus et parois identiques
    assert np.array_equal(g[:gt], m[:mt]) and np.array_equal(g[gt + gb:], m[mt + mb:])
    gv, gf = mesh_export.index_mesh(g)
    mv, mf = mesh_export.index_mesh(m)
    assert mesh_export.signed_volume(mv, mf) == pytest.approx(mesh_export.signed_volume(gv, gf), rel=1e-12)
    assert mb < gb                                                            # moins de triangles au fond


@pytest.mark.parametrize("bottom", MERGED_MODES)
def test_merged_bottom_triangles_face_down_have_no_degenerate_and_lie_on_the_base_plane(bottom):
    mask = blob(3, keep=0.55)
    tris, *_rest, (nt, nb, nw) = build(mask, bottom=bottom)
    b = tris[nt:nt + nb]
    assert np.all(b[:, :, 2] == -BASE)
    cross = np.cross(b[:, 1] - b[:, 0], b[:, 2] - b[:, 0])
    assert np.all(cross[:, 2] < 0)                                             # normales vers -z
    assert check_stl.degenerate_count(b) == 0


@pytest.mark.parametrize("bottom", MERGED_MODES)
@pytest.mark.parametrize("n,radius", [(40, 18.0), (61, 29.0), (100, 45.3)])
def test_merged_bottom_with_smooth_wall_is_closed_oriented_and_same_volume_as_grid(bottom, n, radius):
    mask = _disc(n, radius)
    kw = dict(smooth_wall=True, circle_radius_mm=radius * SCALE)
    g, *_ = build(mask, **kw)
    m, *_rest, (nt, nb, nw) = build(mask, bottom=bottom, **kw)
    rep = check_stl.edge_report(m)
    assert rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0 and rep["duplicate_directed_edges"] == 0
    gv, gf = mesh_export.index_mesh(g)
    mv, mf = mesh_export.index_mesh(m)
    assert mesh_export.signed_volume(mv, mf) == pytest.approx(mesh_export.signed_volume(gv, gf), rel=1e-12)
    assert len(m) < (0.7 if n >= 60 else 1.0) * len(g)


# --- réduction et qualité des triangles sur un disque réaliste ---------------------------------------

def test_a_realistic_disc_loses_almost_the_whole_bottom_with_compact_triangles():
    mask = _disc(400, 190.0)                                                   # disque de 380 cellules de diamètre
    *_a, (gt, gb, gw) = build(mask)
    tg, *_b, (mt, mb, mw) = build(mask, bottom="merged")
    tr, *_c, (rt, rb, rw) = build(mask, bottom="merged-rows")
    assert mb < 0.05 * gb and rb < 0.02 * gb                                  # fond : -96 % (quadtree), -99 % (rectangles)
    assert (mt + mb + mw) < 0.55 * (gt + gb + gw)                              # total : environ -50 %
    assert _aspects(tg[mt:mt + mb]).max() <= 4.0 + 1e-9                        # quadtree : triangles compacts
    assert _aspects(tr[rt:rt + rb]).max() > 1000.0                             # rectangles gloutons : aiguilles (constat)


def test_rows_cap_trades_triangles_for_less_elongated_ones():
    mask = np.ones((60, 200), bool)
    xs, ys = np.arange(201.0), -np.arange(61.0)
    free = bottom_mod.merged_bottom_triangles(mask, xs, ys, 0.0, cover="rows")
    capped = bottom_mod.merged_bottom_triangles(mask, xs, ys, 0.0, max_cells=10, cover="rows")
    assert len(capped) > len(free) and _aspects(capped).max() < _aspects(free).max()


def test_unknown_bottom_mode_is_refused():
    with pytest.raises(ValueError, match="bottom"):
        build(blob(1), bottom="coarse")


# --- de bout en bout (CLI, trois formats) ----------------------------------------------------------

def _run(monkeypatch, tmp_path, capsys, out, *extra):
    a, _ = make_tifs(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--input", a, "--cx", "2601000", "--cy", "1201000", "--radius", "800",
                                      "--diameter-mm", "60", "--vexag", "1", "--base-mm", "2",
                                      "--print-spacing-mm", "0.5", "--output", str(out), *extra])
    circle_dem_to_stl.main()
    capsys.readouterr()
    return mesh_export.read_triangles(str(out))


@pytest.mark.parametrize("bottom", MERGED_MODES)
@pytest.mark.parametrize("ext", ["stl", "3mf", "obj"])
def test_cli_bottom_merged_produces_a_valid_file_in_all_formats(monkeypatch, tmp_path, capsys, bottom, ext):
    grid = _run(monkeypatch, tmp_path, capsys, tmp_path / f"g.{ext}")
    merged = _run(monkeypatch, tmp_path, capsys, tmp_path / f"m.{ext}", "--bottom", bottom)
    assert len(merged) < 0.65 * len(grid)
    res = {n: ok for n, ok, _ in check_stl.check(merged, 60.0, 2.0, tol_mm=1.0)}
    assert all(res.values()), res
    gv, gf = mesh_export.index_mesh(grid)
    mv, mf = mesh_export.index_mesh(merged)
    assert mesh_export.signed_volume(mv, mf) == pytest.approx(mesh_export.signed_volume(gv, gf), rel=1e-4)


def test_cli_default_is_still_the_grid_bottom(monkeypatch, tmp_path, capsys):
    default = _run(monkeypatch, tmp_path, capsys, tmp_path / "d.stl")
    explicit = _run(monkeypatch, tmp_path, capsys, tmp_path / "e.stl", "--bottom", "grid")
    assert np.array_equal(default, explicit)


# --- garde-fou de ressources --------------------------------------------------------------------------

def test_resource_estimate_halves_the_triangles_with_a_merged_bottom():
    grid = circle_dem_to_stl.estimate_resources(250.0, 0.2)
    for mode in ("merged", "merged-rows"):
        merged = circle_dem_to_stl.estimate_resources(250.0, 0.2, bottom=mode)
        assert merged["triangles"] == pytest.approx(grid["triangles"] / 2)
        assert merged["ram_bytes"] < grid["ram_bytes"]
    assert circle_dem_to_stl.estimate_resources(250.0, 0.2, bottom="grid") == grid


def test_a_machine_that_cannot_run_the_grid_bottom_can_run_the_merged_bottom():
    # Ø250 mm à 0,1 mm : mesuré 2,2 Go de pic avec le fond fusionné ; le fond cellule par cellule en demande plus de 4
    ram = lambda: 3_700_000_000  # noqa: E731
    quiet = lambda *a, **k: None  # noqa: E731
    with pytest.raises(SystemExit):
        circle_dem_to_stl.check_resources(250.0, 0.1, False, False, log=quiet, ram_available=ram, bottom="grid")
    circle_dem_to_stl.check_resources(250.0, 0.1, False, False, log=quiet, ram_available=ram, bottom="merged")
