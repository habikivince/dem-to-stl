# SPDX-License-Identifier: GPL-3.0-or-later
"""
mesh.py
-------
Logique géométrique pure (numpy/scipy uniquement, aucun accès disque ni
GDAL) : nettoyage du masque d'inclusion et construction du maillage
dessus/dessous/parois. Séparé de raster.py pour être testable sans
ouvrir le moindre fichier.
"""
import sys

import numpy as np
from scipy import ndimage


def erode3x3(m):
    """Érosion binaire 3x3 (tous les voisins doivent être True pour survivre)."""
    p = np.pad(m, 1, constant_values=False)
    out = m.copy()
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            out &= p[1 + di:1 + di + m.shape[0], 1 + dj:1 + dj + m.shape[1]]
    return out


def dilate3x3(m):
    """Dilatation binaire 3x3 (un seul voisin True suffit)."""
    p = np.pad(m, 1, constant_values=False)
    out = np.zeros_like(m)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            out |= p[1 + di:1 + di + m.shape[0], 1 + dj:1 + dj + m.shape[1]]
    return out


def fix_diagonal_pinches(quad_ok, max_iterations=4):
    """Élimine les pincements en diagonale : deux cellules incluses qui ne se
    touchent que par un coin (motif "damier" sur un bloc 2x2) créent un
    sommet où deux nappes de surface se rejoignent en un point unique —
    épaisseur quasi nulle, indétectable par un simple test d'étanchéité par
    arêtes. On comble les deux cases manquantes plutôt que de retirer de
    la matière. Modifie quad_ok sur place et le retourne."""
    for _ in range(max_iterations):
        diag1 = quad_ok[:-1, :-1] & quad_ok[1:, 1:] & ~quad_ok[:-1, 1:] & ~quad_ok[1:, :-1]
        diag2 = quad_ok[:-1, 1:] & quad_ok[1:, :-1] & ~quad_ok[:-1, :-1] & ~quad_ok[1:, 1:]
        if not (diag1.any() or diag2.any()):
            break
        fix = diag1 | diag2
        quad_ok[:-1, 1:] |= fix
        quad_ok[1:, :-1] |= fix
        quad_ok[:-1, :-1] |= diag2
        quad_ok[1:, 1:] |= diag2
    return quad_ok


def clean_mask(mask):
    """Ouverture morphologique + anti-pincement diagonal, à appliquer à tout
    masque d'inclusion (cercle ou île) avant triangulation."""
    mask = dilate3x3(erode3x3(mask))
    quad_ok = mask[:-1, :-1] & mask[1:, :-1] & mask[:-1, 1:] & mask[1:, 1:]
    quad_ok = fix_diagonal_pinches(quad_ok)
    return mask, quad_ok


def largest_connected_component(seed_mask):
    """Étiquette les composantes connexes de seed_mask (connectivité 8) et
    retourne (masque_composante_principale, nb_composantes, tailles)."""
    labeled, n_components = ndimage.label(seed_mask, structure=np.ones((3, 3)))
    if n_components == 0:
        return np.zeros_like(seed_mask), 0, []
    sizes = ndimage.sum(seed_mask, labeled, index=range(1, n_components + 1))
    largest_label = 1 + int(np.argmax(sizes))
    return labeled == largest_label, n_components, sizes


def build_mesh_from_arrays(elev, quad_ok, xs, ys, scale, vexag, base_mm, min_elev, progress=None):
    """Construit le maillage (triangles dessus/dessous/parois) à partir d'un
    tableau d'altitudes déjà découpé et d'un masque de cellules (quad_ok)
    déjà nettoyé. Ne touche jamais au disque ni à GDAL — testable avec de
    simples tableaux numpy. quad_ok doit avoir la forme (H-1, W-1) par
    rapport à elev. xs/ys doivent déjà être recentrées sur le centre du
    cercle (en unités du CRS source, pas encore mises à l'échelle)."""
    def _progress(i, n, label):
        if progress is not None:
            progress(i, n, label)

    z_top_mm = (elev - min_elev) * scale * vexag
    Xmm = xs * scale
    Ymm = ys * scale

    if not quad_ok.any():
        sys.exit("Erreur : masque vide, rien à mailler.")

    ii, jj = np.nonzero(quad_ok)

    def corner(di, dj):
        i, j = ii + di, jj + dj
        return np.stack([Xmm[j], Ymm[i], z_top_mm[i, j]], axis=1)

    c00, c10, c01, c11 = corner(0, 0), corner(1, 0), corner(0, 1), corner(1, 1)

    top_tris = np.concatenate([
        np.stack([c00, c10, c01], axis=1),
        np.stack([c10, c11, c01], axis=1),
    ])

    zb = -base_mm
    c00b, c10b, c01b, c11b = c00.copy(), c10.copy(), c01.copy(), c11.copy()
    for c in (c00b, c10b, c01b, c11b):
        c[:, 2] = zb

    bot_tris = np.concatenate([
        np.stack([c00b, c01b, c10b], axis=1),
        np.stack([c10b, c01b, c11b], axis=1),
    ])

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
        _progress(i, H1, "Parois (bord vertical)")
        for j in range(W1 + 1):
            left = quad_ok[i, j - 1] if j > 0 else False
            right = quad_ok[i, j] if j < W1 else False
            if left != right:
                p1 = (Xmm[j], Ymm[i], z_top_mm[i, j])
                p2 = (Xmm[j], Ymm[i + 1], z_top_mm[i + 1, j])
                add_wall(p1, p2, outward=bool(left))
    for j in range(W1):
        _progress(j, W1, "Parois (bord horizontal)")
        for i in range(H1 + 1):
            up = quad_ok[i - 1, j] if i > 0 else False
            down = quad_ok[i, j] if i < H1 else False
            if up != down:
                p1 = (Xmm[j], Ymm[i], z_top_mm[i, j])
                p2 = (Xmm[j + 1], Ymm[i], z_top_mm[i, j + 1])
                add_wall(p1, p2, outward=bool(down))

    wall_tris = np.array(wall_tris, dtype=np.float64)
    all_tris = np.concatenate([top_tris, bot_tris, wall_tris])
    return all_tris, len(top_tris), len(bot_tris), len(wall_tris)
