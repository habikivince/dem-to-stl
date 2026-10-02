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


def trace_boundary_loop(quad_ok):
    """Trace la boucle fermée simple formant le contour du masque de
    cellules incluses (quad_ok), en suivant les arêtes de bord — les
    mêmes arêtes que celles où build_mesh_from_arrays construit
    normalement la paroi verticale. Retourne la liste ordonnée des
    points de grille (i, j) formant le contour.

    Suppose une région simplement connexe à un seul contour — lève une
    ValueError explicite sinon (chaque point de bord doit avoir
    exactement 2 arêtes de bord incidentes) plutôt que de produire un
    résultat silencieusement incorrect. C'est le cas normal après
    clean_mask() pour un disque ou une île simple ; une topologie plus
    complexe (plusieurs îles séparées, un contour qui se referme sur
    lui-même) n'est pas supportée par --smooth-wall."""
    H1, W1 = quad_ok.shape
    edges = []
    for i in range(H1):
        for j in range(W1 + 1):
            left = quad_ok[i, j - 1] if j > 0 else False
            right = quad_ok[i, j] if j < W1 else False
            if left != right:
                edges.append(((i, j), (i + 1, j)))
    for j in range(W1):
        for i in range(H1 + 1):
            up = quad_ok[i - 1, j] if i > 0 else False
            down = quad_ok[i, j] if i < H1 else False
            if up != down:
                edges.append(((i, j), (i, j + 1)))

    adjacency = {}
    for a, b in edges:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)

    for v, nbrs in adjacency.items():
        if len(nbrs) != 2:
            raise ValueError(
                f"--smooth-wall requiert un contour unique et simple : le point de grille {v} "
                f"a {len(nbrs)} arête(s) de bord au lieu de 2 (topologie trop complexe : "
                "plusieurs composantes, ou contour non simple)."
            )

    start = next(iter(adjacency))
    loop = [start]
    prev = None
    current = start
    max_len = 4 * (H1 + 1) * (W1 + 1)
    while True:
        nbrs = adjacency[current]
        nxt = nbrs[0] if nbrs[0] != prev else nbrs[1]
        if nxt == start:
            break
        loop.append(nxt)
        prev, current = current, nxt
        if len(loop) > max_len:
            raise ValueError("--smooth-wall : le contour ne s'est pas refermé (masque probablement "
                              "non simplement connexe).")
    return loop


def build_smooth_rim(loop_ij, Xmm, Ymm, z_top_mm, radius_mm, base_mm):
    """Construit le raccord lisse entre le contour en escalier du masque
    nettoyé et le vrai cercle : chaque point du contour est projeté
    radialement (même hauteur) vers le cercle de rayon radius_mm, puis
    relié au socle. Remplace la paroi verticale habituelle (qui partait
    directement du contour en escalier) — ne pas construire les deux
    pour les mêmes arêtes.

    Retourne (band_tris, wall_tris, bottom_band_tris) — trois listes de
    triangles (chaque triangle = 3 tuples (x, y, z)) :
      - band_tris   : bande dessus, contour en escalier -> cercle vrai (hauteur réelle)
      - wall_tris   : paroi verticale, cercle vrai (hauteur réelle) -> socle
      - bottom_band_tris : bande dessous, plate, contour en escalier -> cercle vrai
    """
    loop_xyz = [(Xmm[j], Ymm[i], z_top_mm[i, j]) for i, j in loop_ij]
    n = len(loop_xyz)

    projected = []
    for x, y, z in loop_xyz:
        r = np.hypot(x, y)
        if r < 1e-9:
            projected.append((0.0, 0.0, z))
        else:
            s = radius_mm / r
            projected.append((x * s, y * s, z))

    zb = -base_mm
    band_tris, wall_tris, bottom_band_tris = [], [], []
    for k in range(n):
        k2 = (k + 1) % n
        a, b = loop_xyz[k], loop_xyz[k2]
        pa, pb = projected[k], projected[k2]

        band_tris.append((a, b, pa))
        band_tris.append((b, pb, pa))

        pa_base = (pa[0], pa[1], zb)
        pb_base = (pb[0], pb[1], zb)
        wall_tris.append((pa, pb, pa_base))
        wall_tris.append((pb, pb_base, pa_base))

        a_base = (a[0], a[1], zb)
        b_base = (b[0], b[1], zb)
        bottom_band_tris.append((b_base, a_base, pa_base))
        bottom_band_tris.append((pb_base, b_base, pa_base))

    return band_tris, wall_tris, bottom_band_tris


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


def build_mesh_from_arrays(elev, quad_ok, xs, ys, scale, vexag, base_mm, min_elev,
                            progress=None, smooth_wall=False, circle_radius_mm=None):
    """Construit le maillage (triangles dessus/dessous/parois) à partir d'un
    tableau d'altitudes déjà découpé et d'un masque de cellules (quad_ok)
    déjà nettoyé. Ne touche jamais au disque ni à GDAL — testable avec de
    simples tableaux numpy. quad_ok doit avoir la forme (H-1, W-1) par
    rapport à elev. xs/ys doivent déjà être recentrées sur le centre du
    cercle (en unités du CRS source, pas encore mises à l'échelle).

    smooth_wall=True : au lieu de suivre le contour en escalier du masque,
    chaque point du contour est projeté radialement vers le vrai cercle
    (même altitude), formant un bord parfaitement circulaire. Requiert un
    masque à un seul contour simple (voir trace_boundary_loop)."""
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

    if smooth_wall:
        loop_ij = trace_boundary_loop(quad_ok)
        if circle_radius_mm is None:
            raise ValueError("smooth_wall=True nécessite circle_radius_mm (le vrai rayon, en mm).")
        radius_mm = circle_radius_mm
        band_tris, wall_tris, bottom_band_tris = build_smooth_rim(
            loop_ij, Xmm, Ymm, z_top_mm, radius_mm, base_mm)
        extra_tris = np.array(band_tris + bottom_band_tris, dtype=np.float64)
        wall_tris = np.array(wall_tris, dtype=np.float64)
        all_tris = np.concatenate([top_tris, bot_tris, extra_tris, wall_tris])
        return all_tris, len(top_tris), len(bot_tris) + len(bottom_band_tris) + len(band_tris), len(wall_tris)

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
