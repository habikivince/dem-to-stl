# SPDX-License-Identifier: GPL-3.0-or-later
"""Contrôles de topologie de check_stl.py : coques, caractéristique d'Euler V − E + F, triangles retournés.

Les solides de référence sont écrits à la main (cube, tore, cube ouvert, deux cubes) : leur topologie est connue
sans appeler le code de maillage. Les solides de mesh.py servent ensuite à vérifier l'absence de faux positifs
(îles, trous, fortes pentes) et qu'un sommet déplacé de l'autre côté de ses voisins est bien vu, alors que
l'étanchéité et le volume signé ne le voient pas. N'a besoin ni de GDAL ni de fichiers : numpy et scipy."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import mesh  # noqa: E402
import stl_io  # noqa: E402

TOPO = "topologie (coques, Euler)"
INV = "triangles retournés"


def cube(origin=(0.0, 0.0, 0.0), size=1.0):
    """12 triangles, faces tournées vers l'extérieur."""
    o = np.asarray(origin, dtype=float)
    v = o + size * np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                             [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
    quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    tris = []
    for a, b, c, d in quads:
        tris += [[v[a], v[b], v[c]], [v[a], v[c], v[d]]]
    return np.array(tris)


def torus(n=12, m=8, big=3.0, small=1.0):
    """Tore triangulé, faces tournées vers l'extérieur : V − E + F = 0."""
    u = np.linspace(0, 2 * np.pi, n, endpoint=False)
    w = np.linspace(0, 2 * np.pi, m, endpoint=False)

    def p(i, j):
        i, j = i % n, j % m
        r = big + small * np.cos(w[j])
        return np.array([r * np.cos(u[i]), r * np.sin(u[i]), small * np.sin(w[j])])

    tris = []
    for i in range(n):
        for j in range(m):
            a, b, c, d = p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1)
            tris += [[a, b, c], [a, c, d]]
    t = np.array(tris)
    if check_stl.signed_volume(t) < 0:
        t = t[:, [0, 2, 1]]
    return t


def results(tris, **kw):
    return {n: (ok, d) for n, ok, d in check_stl.check(tris, 10.0, 0.0, tol_mm=100.0, topology=True, **kw)}


# ------------------------------------------------------------------ solides de topologie connue

def test_cube_is_one_closed_ball():
    t = check_stl.topology_report(cube())
    assert (t["faces"], t["vertices"], t["edges"]) == (12, 8, 18)
    assert t["shells"] == 1 and t["shell_euler"] == [2] and t["euler"] == 2
    assert t["shared_vertices"] == 0 and t["inverted_faces"] == 0
    r = results(cube(), expect_shells=1)
    assert r[TOPO][0] and r[INV][0]


def test_two_separate_cubes_are_two_shells():
    two = np.concatenate([cube(), cube(origin=(5.0, 0.0, 0.0))])
    t = check_stl.topology_report(two)
    assert t["shells"] == 2 and t["shell_euler"] == [2, 2] and t["euler"] == 4 and t["shared_vertices"] == 0
    assert results(two)[TOPO][0]                                    # sans attente : décrit, ne refuse pas
    assert results(two, expect_shells=2)[TOPO][0]
    ok, detail = results(two, expect_shells=1)[TOPO]
    assert not ok and "2 coque(s)" in detail


def test_cubes_touching_at_a_single_vertex_are_refused():
    touching = np.concatenate([cube(), cube(origin=(1.0, 1.0, 1.0))])
    t = check_stl.topology_report(touching)
    assert t["shells"] == 2 and t["shared_vertices"] == 1
    ok, detail = results(touching)[TOPO]
    assert not ok and "partagé" in detail


def test_torus_has_a_handle():
    t = check_stl.topology_report(torus())
    assert t["shells"] == 1 and t["shell_euler"] == [0]
    ok, detail = results(torus())[TOPO]
    assert ok and "anse" in detail                                  # décrit : un tore peut être voulu
    assert not results(torus(), expect_shells=1)[TOPO][0]           # mais pas si on attend un disque


def test_open_surface_is_not_a_closed_ball():
    opened = cube()[2:]                                             # une face carrée (2 triangles) retirée
    t = check_stl.topology_report(opened)
    assert t["shells"] == 1 and t["shell_euler"] == [1]             # un disque : V − E + F = 8 − 17 + 10
    r = results(opened)
    assert not r["étanchéité"][0]
    assert not results(opened, expect_shells=1)[TOPO][0]


def test_topology_is_opt_in_so_existing_callers_are_unchanged():
    names = [n for n, _, _ in check_stl.check(cube(), 10.0, 0.0, tol_mm=100.0)]
    assert TOPO not in names and INV not in names
    assert check_stl.edge_report(cube()) == check_stl.edge_report(cube(), check_stl._weld(cube()))


def test_degenerate_triangles_are_ignored_not_counted():
    tris = np.concatenate([cube(), np.array([[[0, 0, 0], [0, 0, 0], [1, 0, 0]]], dtype=float)])
    t = check_stl.topology_report(tris)
    assert t["ignored_faces"] == 1 and t["shells"] == 1 and t["shell_euler"] == [2]


# ------------------------------------------------------------------ solides de mesh.py

def solid(mask, vexag=1.0, seed=0, amp=4.0, scale=1.0, **kw):
    """Même construction que tests/test_mesh_properties.build, sans GDAL ; amp est l'amplitude du relief en unités
    de cellule (4 : pente d'environ 20° à vexag 1)."""
    n = mask.shape[0]
    elev = 500.0 + ndimage.gaussian_filter(np.random.default_rng(seed).normal(size=mask.shape), 3) * amp
    xs = np.arange(mask.shape[1]) - (mask.shape[1] - 1) / 2.0
    ys = (mask.shape[0] - 1) / 2.0 - np.arange(n)
    clean, quad_ok = mesh.clean_mask(mask)
    tris = mesh.build_mesh_from_arrays(elev, quad_ok, xs, ys, scale, vexag, 2.0, elev[clean].min(), **kw)[0]
    return tris, quad_ok


def blob(seed, n=44, keep=0.55):
    noise = ndimage.gaussian_filter(np.random.default_rng(seed).normal(size=(n, n)), 2.2)
    return noise > np.quantile(noise, 1 - keep)


def disc(n=40, r=17):
    yy, xx = np.mgrid[:n, :n]
    return (yy - n / 2 + 0.5) ** 2 + (xx - n / 2 + 0.5) ** 2 <= r * r


def test_plain_disc_is_one_ball_and_ring_has_a_handle():
    t = check_stl.topology_report(solid(disc())[0])
    assert t["shells"] == 1 and t["shell_euler"] == [2] and t["inverted_faces"] == 0
    ring = disc(40, 17) & ~disc(40, 8)
    tr = check_stl.topology_report(solid(ring)[0])
    assert tr["shells"] == 1 and tr["shell_euler"] == [0]           # anneau = solide à une anse


def test_two_islands_give_two_shells():
    mask = np.zeros((40, 40), dtype=bool)
    mask[4:16, 4:16] = True
    mask[22:36, 20:36] = True
    t = check_stl.topology_report(solid(mask)[0])
    assert t["shells"] == 2 and t["shell_euler"] == [2, 2] and t["shared_vertices"] == 0


@pytest.mark.parametrize("seed", range(10))
@pytest.mark.parametrize("variant", [{}, {"diagonal": "adaptive"}, {"top": "adaptive", "top_tol_mm": 0.05},
                                     {"bottom": "merged"}])
def test_random_masks_meshed_by_mesh_py_are_clean(seed, variant):
    """Îles, trous et goulets (masques de test_mesh_properties) : une coque par composante du masque, aucune
    anse impaire (surfaces fermées orientables : V − E + F pair et ≤ 2), aucun contact par un sommet, aucun
    triangle retourné, y compris avec une forte exagération verticale."""
    mask = blob(seed, keep=0.4 + 0.03 * (seed % 6))
    for vexag in (1.0, 8.0):
        tris, quad_ok = solid(mask, vexag=vexag, seed=seed, amp=10.0, **variant)
        t = check_stl.topology_report(tris)
        assert t["shells"] == ndimage.label(quad_ok)[1]
        assert all(e <= 2 and e % 2 == 0 for e in t["shell_euler"])
        assert t["shared_vertices"] == 0 and t["inverted_faces"] == 0 and t["ignored_faces"] == 0


def test_steep_exaggerated_relief_is_not_reported_as_inverted():
    """Falaises : pente moyenne du dessus d'environ 80° (amp 40, vexag 8) ; une arête vive n'est pas un retournement."""
    for seed in range(6):
        tris, _ = solid(blob(seed), vexag=8.0, seed=seed, amp=40.0)
        assert check_stl.topology_report(tris)["inverted_faces"] == 0


def test_smooth_wall_disc_is_clean():
    n, r = 40, 17
    tris, _ = solid(disc(n, r), smooth_wall=True, circle_radius_mm=float(r))
    t = check_stl.topology_report(tris)
    assert t["shells"] == 1 and t["shell_euler"] == [2] and t["inverted_faces"] == 0
    assert check_stl.edge_report(tris)["boundary_edges"] == 0


# ------------------------------------------------------------------ retournement géométrique

def test_a_vertex_pushed_across_its_neighbours_is_seen_only_by_the_new_check():
    """Un sommet du dessus déplacé de l'autre côté d'une arête de son voisinage retourne un triangle sans toucher à la
    connectivité : l'étanchéité et le volume signé restent bons, seul le contrôle des triangles retournés le voit."""
    tris, _ = solid(disc(), seed=3, amp=4.0)
    assert check_stl.topology_report(tris)["inverted_faces"] == 0
    ids, _ = check_stl._weld(tris)
    count = np.bincount(ids.ravel())
    corners = tris.reshape(-1, 3)
    first = np.zeros(len(count), dtype=np.int64)
    first[ids.ravel()[::-1]] = np.arange(len(corners))[::-1]
    interior = [corners[first[v]] for v in np.flatnonzero(count == 6) if corners[first[v]][2] > 0.0][:6]
    assert interior, "aucun sommet intérieur du dessus trouvé"
    found = None
    for v in interior:
        sel = np.all(np.isclose(tris, v), axis=2)
        for ang in np.arange(0, 360, 45):
            for dist in (0.7, 0.9, 1.1, 1.4, 2.0):
                moved = tris.copy()
                moved[sel, 0] += dist * np.cos(np.radians(ang))
                moved[sel, 1] += dist * np.sin(np.radians(ang))
                t = check_stl.topology_report(moved)
                if t["inverted_faces"] >= 1:
                    found = (moved, t)
                    break
            if found:
                break
        if found:
            break
    assert found, "le déplacement d'un sommet n'a jamais retourné de triangle détecté"
    moved, t = found
    rep = check_stl.edge_report(moved)
    assert rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0 and rep["duplicate_directed_edges"] == 0
    assert check_stl.signed_volume(moved) > 0
    r = results(moved)
    assert not r[INV][0] and "près de" in r[INV][1]


# ------------------------------------------------------------------ ligne de commande

def test_cli_prints_topology_lines_and_expect_shells_fails_on_two_islands(tmp_path, monkeypatch, capsys):
    mask = np.zeros((40, 40), dtype=bool)
    mask[4:16, 4:16] = True
    mask[22:36, 20:36] = True
    tris, _ = solid(mask)
    path = tmp_path / "deux_iles.stl"
    stl_io.write_stl_binary(tris, str(path))
    zmin = float(tris[..., 2].min())
    argv = ["check_stl.py", str(path), "--diameter-mm", "40", "--base-mm", f"{-zmin:g}", "--tol-mm", "1000"]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit) as e:
        check_stl.main()
    out = capsys.readouterr().out
    assert e.value.code == 0 and "topologie (coques, Euler)" in out and "2 coque(s)" in out
    monkeypatch.setattr(sys, "argv", argv + ["--expect-shells", "1"])
    with pytest.raises(SystemExit) as e:
        check_stl.main()
    assert e.value.code == 1 and "[ÉCHEC] topologie" in capsys.readouterr().out
    monkeypatch.setattr(sys, "argv", argv + ["--no-topology"])
    with pytest.raises(SystemExit):
        check_stl.main()
    assert "topologie" not in capsys.readouterr().out
