"""
Test unitaire de fix_diagonal_pinches(), isolé du reste du pipeline —
pas besoin de GDAL ni de DEM, juste un tableau booléen 2D. Reproduit le
motif en damier (deux cellules incluses qui ne se touchent que par un
coin) qui créait un sommet non-manifold indétectable par le test
d'étanchéité par arêtes.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mesh import fix_diagonal_pinches, erode3x3, dilate3x3  # noqa: E402


def test_checkerboard_pattern_is_filled():
    quad_ok = np.array([
        [True, False, False],
        [False, True, False],
        [False, False, False],
    ], dtype=bool)

    fixed = fix_diagonal_pinches(quad_ok.copy())

    expected = np.array([
        [True, True, False],
        [True, True, False],
        [False, False, False],
    ], dtype=bool)
    assert np.array_equal(fixed, expected)


def test_second_diagonal_orientation_is_filled():
    quad_ok = np.array([
        [False, True],
        [True, False],
    ], dtype=bool)
    fixed = fix_diagonal_pinches(quad_ok.copy())
    assert fixed.all()


def test_no_pinch_leaves_mask_unchanged():
    quad_ok = np.array([
        [True, True],
        [True, True],
    ], dtype=bool)
    fixed = fix_diagonal_pinches(quad_ok.copy())
    assert np.array_equal(fixed, quad_ok)


def test_erode_then_dilate_removes_single_pixel_protrusion():
    mask = np.zeros((7, 7), dtype=bool)
    mask[2:5, 2:5] = True  # bloc plein 3x3
    mask[1, 3] = True      # languette d'un seul pixel qui dépasse

    opened = dilate3x3(erode3x3(mask))
    assert not opened[1, 3], "la languette d'un seul pixel aurait dû être supprimée"
    assert opened[2:5, 2:5].all(), "le bloc plein ne doit pas être entamé"


def _count_pinches(q):
    d1 = q[:-1, :-1] & q[1:, 1:] & ~q[:-1, 1:] & ~q[1:, :-1]
    d2 = q[:-1, 1:] & q[1:, :-1] & ~q[:-1, :-1] & ~q[1:, 1:]
    return int(d1.sum() + d2.sum())


def test_fix_runs_until_stable_not_a_fixed_number_of_passes():
    # Masque aléatoire (16x16, graine 87) où l'ancien plafond de 4 passes
    # laissait un pincement : le comblement d'un pincement peut en créer un voisin.
    q = np.random.default_rng(87).random((16, 16)) < 0.3
    assert _count_pinches(fix_diagonal_pinches(q.copy(), max_iterations=4)) > 0
    assert _count_pinches(fix_diagonal_pinches(q.copy())) == 0


def test_fix_leaves_no_pinch_on_random_masks_and_only_adds_cells():
    for density in (0.3, 0.5, 0.7):
        for seed in range(25):
            q = np.random.default_rng(seed).random((30, 30)) < density
            before = q.copy()
            out = fix_diagonal_pinches(q)
            assert _count_pinches(out) == 0
            assert (out | before).sum() == out.sum()      # monotone : on n'enlève jamais de matière
            assert (fix_diagonal_pinches(out.copy()) == out).all()   # idempotent
