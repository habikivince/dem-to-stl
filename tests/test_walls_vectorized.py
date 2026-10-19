# SPDX-License-Identifier: GPL-3.0-or-later
"""Parois sans boucle Python (mesh.build_plain_walls) et non-régression du maillage par défaut.

1. Les parois vectorisées sont comparées, tableau contre tableau (mêmes valeurs, même ordre), à l'ancienne
   double boucle sur tous les nœuds, conservée ici comme implémentation de référence.
2. Des empreintes SHA-256 de sorties calculées avec le code ORIGINAL de ``mesh.py`` (avant les chantiers de
   l'audit) sur des entrées entièrement déterministes (entiers et multiples exacts, sans générateur
   aléatoire) vérifient que le maillage par défaut (fond « grid ») n'a pas changé d'un seul bit. Si un
   changement de maillage est voulu un jour, ces empreintes sont à régénérer volontairement."""
import hashlib
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import mesh  # noqa: E402
from test_mesh_properties import blob  # noqa: E402


def reference_plain_walls(quad_ok, Xmm, Ymm, z_top_mm, zb):
    """L'ancienne construction (double boucle sur tous les nœuds), reprise telle quelle du code d'origine."""
    wall_tris = []

    def add_wall(p1t, p2t, outward):
        p1b = (p1t[0], p1t[1], zb)
        p2b = (p2t[0], p2t[1], zb)
        if outward:
            wall_tris.append((p1t, p2t, p1b))
            wall_tris.append((p2t, p2b, p1b))
        else:
            wall_tris.append((p2t, p1t, p1b))
            wall_tris.append((p2b, p2t, p1b))

    H1, W1 = quad_ok.shape
    for i in range(H1):
        for j in range(W1 + 1):
            left = quad_ok[i, j - 1] if j > 0 else False
            right = quad_ok[i, j] if j < W1 else False
            if left != right:
                p1 = (Xmm[j], Ymm[i], z_top_mm[i, j])
                p2 = (Xmm[j], Ymm[i + 1], z_top_mm[i + 1, j])
                add_wall(p1, p2, outward=bool(left))
    for j in range(W1):
        for i in range(H1 + 1):
            up = quad_ok[i - 1, j] if i > 0 else False
            down = quad_ok[i, j] if i < H1 else False
            if up != down:
                p1 = (Xmm[j], Ymm[i], z_top_mm[i, j])
                p2 = (Xmm[j + 1], Ymm[i], z_top_mm[i, j + 1])
                add_wall(p1, p2, outward=bool(down))
    return np.array(wall_tris, dtype=np.float64)


def _inputs(mask, seed=0):
    H, W = mask.shape
    rng = np.random.default_rng(seed)
    z = rng.normal(10.0, 3.0, (H + 1, W + 1))
    return z, np.arange(W + 1) * 0.1, -np.arange(H + 1) * 0.1


@pytest.mark.parametrize("seed", range(20))
def test_vectorized_walls_equal_the_original_loops_exactly(seed):
    mask = blob(seed, n=36, keep=0.35 + 0.04 * (seed % 7))
    z, X, Y = _inputs(mask, seed)
    got = mesh.build_plain_walls(mask, X, Y, z, -2.0)
    ref = reference_plain_walls(mask, X, Y, z, -2.0)
    assert got.shape == ref.shape and np.array_equal(got, ref)


def test_vectorized_walls_edge_cases():
    cases = {
        "une cellule": np.pad(np.ones((1, 1), bool), 3),
        "grille pleine (le masque touche les quatre bords)": np.ones((7, 9), bool),
        "bande d'une cellule contre le bord": np.pad(np.ones((1, 8), bool), ((0, 6), (0, 0))),
        "damier de 2x2 (le pincement est déjà traité en amont)": np.array([[1, 0], [0, 1]], bool),
        "trou central": np.pad(np.pad(np.zeros((3, 3), bool), 2, constant_values=True), 1),
    }
    for name, mask in cases.items():
        z, X, Y = _inputs(mask)
        got = mesh.build_plain_walls(mask, X, Y, z, -1.0)
        assert np.array_equal(got, reference_plain_walls(mask, X, Y, z, -1.0)), name


def test_vectorized_walls_return_an_empty_array_for_an_empty_mask():
    mask = np.zeros((5, 6), bool)
    z, X, Y = _inputs(mask)
    got = mesh.build_plain_walls(mask, X, Y, z, -1.0)
    assert got.shape == (0, 3, 3)


# --- non-régression du maillage par défaut (empreintes calculées avec le code d'origine) -----------------

def _golden_case(name):
    n = 41
    ii, jj = np.mgrid[0:n, 0:n]
    elev = (((ii * 7 + jj * 13) % 17) * 3.5 + (ii * jj % 5) * 1.25).astype(np.float64)    # valeurs exactes
    ci, cj = np.mgrid[0:n - 1, 0:n - 1]
    d2 = (ci - 19) ** 2 + (cj - 19) ** 2
    if name == "disc":
        mask = d2 <= 16 ** 2
    elif name == "ring":
        mask = (d2 <= 17 ** 2) & (d2 >= 6 ** 2)
    else:
        mask = np.zeros((n - 1, n - 1), bool)
        mask[3:35, 3:15] = True
        mask[25:37, 3:36] = True
    return elev, mask, jj[0].astype(np.float64), -ii[:, 0].astype(np.float64)


def _digest(t):
    return hashlib.sha256(np.ascontiguousarray(t, dtype=np.float64).tobytes()).hexdigest()


GOLDEN = {
    "disc": ("b67ab82a118848316918f9b5384c9aaea0b7d44919788d8173e4d43f4aca66db", (1594, 1594, 264)),
    "ring": ("1186a23b84c987f24dbb37c906e3aca33914af5fe3b3e181fa78ac9e63ec5f8e", (1584, 1584, 368)),
    "L": ("a207f228bf75f8975bc18d955db3cfc5d7502ca07cffd6403a950dd6acebaf53", (1320, 1320, 268)),
}


@pytest.mark.parametrize("name", sorted(GOLDEN))
def test_default_mesh_is_bit_identical_to_the_original_code(name):
    elev, mask, xs, ys = _golden_case(name)
    tris, n_top, n_bot, n_wall = mesh.build_mesh_from_arrays(elev, mask.copy(), xs, ys, 0.1, 1.5, 2.0, float(elev.min()))
    digest, counts = GOLDEN[name]
    assert (n_top, n_bot, n_wall) == counts and _digest(tris) == digest


def test_default_mesh_with_smooth_wall_is_bit_identical_to_the_original_code():
    elev, mask, xs, ys = _golden_case("disc")
    tris, n_top, n_bot, n_wall = mesh.build_mesh_from_arrays(
        elev, mask.copy(), xs, ys, 0.1, 1.5, 2.0, float(elev.min()), smooth_wall=True, circle_radius_mm=16.5 * 0.1)
    assert (n_top, n_bot, n_wall) == (1594, 2122, 264)
    assert _digest(tris) == "ab28b98842dff622c351a9dcd984bc306063f774f33278619aecc08a0ae8060b"
