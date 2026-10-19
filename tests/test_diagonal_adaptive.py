# SPDX-License-Identifier: GPL-3.0-or-later
"""Diagonale adaptative du dessus (--diagonal adaptive, prototype).

Seule la surface change : même nombre de triangles, mêmes sommets, mêmes parois et même fond que la
diagonale fixe ; le solide reste fermé, orienté, et son volume est celui que donne la règle appliquée
cellule par cellule (oracle écrit ici, indépendant de mesh.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import circle_dem_to_stl  # noqa: E402
import mesh  # noqa: E402
import mesh_export  # noqa: E402
from test_mesh_properties import BASE, SCALE, assert_solid, blob, build, make_tifs  # noqa: E402


def adaptive_expected_volume(elev, quad_ok, xs, ys, min_e, vexag):
    ii, jj = np.nonzero(quad_ok)
    z = lambda di, dj: (elev[ii + di, jj + dj] - min_e) * SCALE * vexag  # noqa: E731
    area = np.abs((xs[1] - xs[0]) * (ys[1] - ys[0])) * SCALE * SCALE
    z00, z10, z01, z11 = z(0, 0), z(1, 0), z(0, 1), z(1, 1)
    ac = np.abs(z00 - z11) <= np.abs(z10 - z01)
    v_ac = area / 6 * (2 * z00 + z10 + 2 * z11 + z01)          # diagonale c00-c11 : deux triangles de aire/2
    v_bd = area / 6 * (z00 + 2 * z10 + z11 + 2 * z01)          # diagonale c10-c01
    return float(np.sum(np.where(ac, v_ac, v_bd)) + area * len(ii) * BASE)


def _assert_adaptive_solid(tris, elev, quad_ok, xs, ys, min_e, vexag=1.0):
    rep = check_stl.edge_report(tris)
    assert rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0 and rep["duplicate_directed_edges"] == 0, rep
    v, f = mesh_export.index_mesh(tris)
    assert mesh_export.signed_volume(v, f) == pytest.approx(
        adaptive_expected_volume(elev, quad_ok, xs, ys, min_e, vexag), rel=1e-4)


@pytest.mark.parametrize("bottom", ["grid", "merged"])
@pytest.mark.parametrize("seed", range(14))
def test_adaptive_diagonal_gives_a_closed_oriented_solid_with_the_expected_volume(bottom, seed):
    mask = blob(seed, keep=0.4 + 0.03 * (seed % 6))
    vexag = 1.0 + seed % 3
    tris, elev, quad_ok, xs, ys, min_e, _ = build(mask, vexag=vexag, seed=seed, diagonal="adaptive", bottom=bottom)
    _assert_adaptive_solid(tris, elev, quad_ok, xs, ys, min_e, vexag)


@pytest.mark.parametrize("seed", range(6))
def test_only_the_top_surface_changes(seed):
    mask = blob(seed, keep=0.5)
    fixed, *_a, (ft, fb, fw) = build(mask, seed=seed)
    adapt, *_b, (at, ab, aw) = build(mask, seed=seed, diagonal="adaptive")
    assert (ft, fb, fw) == (at, ab, aw)                                       # mêmes nombres de triangles
    assert np.array_equal(fixed[ft:], adapt[at:])                             # fond et parois identiques
    assert not np.array_equal(fixed[:ft], adapt[:at])                         # la surface, elle, change
    key = lambda t: {tuple(p) for tri in t for p in tri}                     # noqa: E731
    assert key(fixed[:ft]) == key(adapt[:at])                                 # mêmes sommets


def test_fixed_is_still_the_default_and_unchanged():
    mask = blob(4, keep=0.5)
    default, *_ = build(mask, seed=4)
    explicit, *_ = build(mask, seed=4, diagonal="fixed")
    assert np.array_equal(default, explicit)


def _one_cell(z00, z10, z01, z11):
    elev = np.array([[z00, z01], [z10, z11]], dtype=float)
    mask = np.ones((1, 1), bool)
    tris, *_ = mesh.build_mesh_from_arrays(elev, mask, np.array([0.0, 1.0]), np.array([0.0, -1.0]), 1.0, 1.0, 1.0,
                                           0.0, diagonal="adaptive")             # min_elev = 0 : altitudes inchangées
    return tris[:2]                                                           # les deux triangles du dessus


def _has_edge(tris, p, q):
    inside = [{tuple(v) for v in t} for t in tris]
    return all(p in s and q in s for s in inside)


def test_the_diagonal_joins_the_two_corners_with_the_smaller_height_difference():
    # |z00 - z11| = 1 < |z10 - z01| = 8 : diagonale c00-c11, soit les nœuds (0, 0) et (1, -1)
    t = _one_cell(z00=10, z10=2, z01=18, z11=11)
    assert _has_edge(t, (0.0, 0.0, 10.0), (1.0, -1.0, 11.0))
    # |z00 - z11| = 8 > |z10 - z01| = 1 : diagonale c10-c01 (celle de l'ancien code)
    t = _one_cell(z00=10, z10=2, z01=3, z11=18)
    assert _has_edge(t, (0.0, -1.0, 2.0), (1.0, 0.0, 3.0))


def test_a_tie_keeps_the_c00_c11_diagonal_and_a_flat_cell_is_valid():
    t = _one_cell(z00=5, z10=5, z01=5, z11=5)
    assert _has_edge(t, (0.0, 0.0, 5.0), (1.0, -1.0, 5.0))
    assert np.all(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])[:, 2] > 0)    # normales vers +z


def test_adaptive_triangles_face_up():
    mask = blob(2, keep=0.55)
    tris, *_rest, (nt, nb, nw) = build(mask, diagonal="adaptive")
    top = tris[:nt]
    assert np.all(np.cross(top[:, 1] - top[:, 0], top[:, 2] - top[:, 0])[:, 2] > 0)


def test_unknown_diagonal_mode_is_refused():
    with pytest.raises(ValueError, match="diagonal"):
        build(blob(1), diagonal="zigzag")


def test_adaptive_diagonal_with_smooth_wall_is_closed():
    n = 61
    yy, xx = np.mgrid[0:n, 0:n]
    mask = (xx - 30) ** 2 + (yy - 30) ** 2 <= 29 ** 2
    tris, *_ = build(mask, diagonal="adaptive", smooth_wall=True, circle_radius_mm=29 * SCALE)
    rep = check_stl.edge_report(tris)
    assert rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0 and rep["duplicate_directed_edges"] == 0


@pytest.mark.parametrize("ext", ["stl", "3mf", "obj"])
def test_cli_adaptive_diagonal_produces_a_valid_file(monkeypatch, tmp_path, capsys, ext):
    a, _ = make_tifs(tmp_path)

    def run(out, *extra):
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["x", "--input", a, "--cx", "2601000", "--cy", "1201000", "--radius", "800",
                                          "--diameter-mm", "60", "--vexag", "1", "--base-mm", "2",
                                          "--print-spacing-mm", "0.5", "--output", str(out), *extra])
        circle_dem_to_stl.main()
        capsys.readouterr()
        return mesh_export.read_triangles(str(out))

    fixed = run(tmp_path / f"f.{ext}")
    adapt = run(tmp_path / f"a.{ext}", "--diagonal", "adaptive")
    assert len(adapt) == len(fixed)
    res = {n: ok for n, ok, _ in check_stl.check(adapt, 60.0, 2.0, tol_mm=1.0)}
    assert all(res.values()), res
    assert not np.array_equal(np.sort(fixed.reshape(len(fixed), -1), axis=0), np.sort(adapt.reshape(len(adapt), -1), axis=0))
