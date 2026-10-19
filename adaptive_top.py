# SPDX-License-Identifier: GPL-3.0-or-later
"""
adaptive_top.py
---------------
Dessus adaptatif (PROTOTYPE, option --top adaptive) : maillage RTIN (Right-Triangulated Irregular Network,
Evans et al. 2001 ; formulation « Martini » de Mapbox) à erreur bornée, avec bord du masque forcé.

Principe. La grille de nœuds est complétée jusqu'à une taille 2^k + 1. Deux grands triangles rectangles la
couvrent ; chacun se coupe en deux par la médiane de son hypoténuse, récursivement. Pour chaque milieu
d'hypoténuse m, l'erreur est l'écart entre l'altitude du nœud et l'interpolation de l'hypoténuse, majoré par
l'erreur des descendants. Un triangle est conservé tel quel si l'erreur de son milieu ne dépasse pas la
tolérance ; sinon il est coupé. Le maillage obtenu est sans fissure : deux triangles voisins partagent leurs
sommets.

Adaptations au solide :
  * les nœuds du bord du masque (cellules incluses contre cellules exclues ou hors masque) ont une erreur
    infinie : ils sont toujours des sommets et les triangles qui les touchent sont des demi-cellules, donc
    les parois verticales, construites arête de cellule par arête de cellule, s'y raccordent exactement ;
  * les nœuds sans aucune cellule voisine dans le masque ont une erreur nulle (ils ne forcent rien) ;
  * un triangle feuille n'est conservé que si sa cellule centrale est dans le masque.

Garantie : l'écart entre l'altitude d'un nœud et la surface est borné par une petite multiple de la tolérance
(l'erreur de chaque milieu est mesurée par rapport à l'interpolation de ses voisins hiérarchiques, qui
peuvent eux-mêmes être omis ; le cumul est borné par la série géométrique, voir les tests). La fidélité n'est
donc plus exacte aux nœuds, contrairement à la grille complète : c'est le prix de la réduction du nombre de
triangles.

Mémoire bornée : la grille est traitée par blocs de ``block`` cellules (puissance de deux, 1 024 par défaut),
chacun avec son propre RTIN. Les cellules voisines situées hors du bloc comptent comme « hors masque » : les
nœuds de jonction entre deux blocs sont donc forcés exactement comme le bord du masque, et les deux blocs y
partagent les mêmes sommets (aucune fissure). Surcoût : une bande de demi-cellules le long des jonctions.

Tout est vectorisé niveau par niveau (aucune récursion Python).
"""
import numpy as np


def _children(ax, ay, bx, by, cx, cy):
    """Triangles enfants (c, a, m) puis (b, c, m), m = milieu de l'hypoténuse ab."""
    mx, my = (ax + bx) >> 1, (ay + by) >> 1
    return (np.concatenate([cx, bx]), np.concatenate([cy, by]),
            np.concatenate([ax, cx]), np.concatenate([ay, cy]),
            np.concatenate([mx, mx]), np.concatenate([my, my]))


def _rtin_block(quad_ok, Xmm, Ymm, z_top_mm, tol_mm):
    """RTIN d'un bloc : cellules ``quad_ok`` (H, W), nœuds (H+1, W+1). Voir adaptive_top_triangles."""
    H, W = quad_ok.shape
    k = max(1, int(np.ceil(np.log2(max(H, W, 2)))))
    T = 1 << k                                   # taille de la grille complétée, en cellules
    N = T + 1                                    # en nœuds
    z = np.pad(np.asarray(z_top_mm, dtype=np.float64), ((0, N - H - 1), (0, N - W - 1)), mode="edge")
    X = np.pad(np.asarray(Xmm, dtype=np.float64), (0, N - W - 1), mode="edge")
    Y = np.pad(np.asarray(Ymm, dtype=np.float64), (0, N - H - 1), mode="edge")
    cells = np.zeros((T, T), dtype=bool)
    cells[:H, :W] = quad_ok

    cp = np.pad(cells, 1)                        # cellule (i, j) en cp[i + 1, j + 1]
    around = np.stack([cp[0:N, 0:N], cp[0:N, 1:N + 1], cp[1:N + 1, 0:N], cp[1:N + 1, 1:N + 1]])
    any_in, all_in = around.any(axis=0), around.all(axis=0)
    forced = any_in & ~all_in                    # nœuds du bord du masque
    outside = ~any_in                            # nœuds sans cellule voisine dans le masque

    # --- niveaux, du plus grand au plus petit triangle : racines (a, b, c), c = angle droit -------------
    i32 = np.int32
    lv = (np.array([T, 0], i32), np.array([T, 0], i32), np.array([0, T], i32),
          np.array([0, T], i32), np.array([T, 0], i32), np.array([0, T], i32))
    levels = [lv]
    for _ in range(2 * k - 1):                   # niveaux 0 .. 2k-1 : triangles dont l'hypoténuse a un milieu entier
        lv = _children(*lv)
        levels.append(lv)
    # le niveau 2k (triangles d'une demi-cellule) n'a pas de milieu entier : il est généré à la demande

    # --- erreurs, du plus petit au plus grand ------------------------------------------------------------
    E = np.zeros((N, N), dtype=np.float64)
    for depth in range(2 * k - 1, -1, -1):
        ax, ay, bx, by, cx, cy = levels[depth]
        mx, my = (ax + bx) >> 1, (ay + by) >> 1
        own = np.abs(z[my, mx] - (z[ay, ax] + z[by, bx]) * 0.5)
        own = np.where(forced[my, mx], np.inf, np.where(outside[my, mx], 0.0, own))
        if depth < 2 * k - 1:                    # les enfants ont un milieu entier : on reprend leur erreur
            child = np.maximum(E[(ay + cy) >> 1, (ax + cx) >> 1], E[(by + cy) >> 1, (bx + cx) >> 1])
            own = np.maximum(own, child)
        np.maximum.at(E, (my, mx), own)          # deux triangles partagent chaque milieu : on garde le maximum
    root = levels[0]
    del levels                                    # les niveaux ne servent plus : on libère la mémoire avant l'extraction

    # --- extraction des triangles feuilles, du plus grand au plus petit -----------------------------------
    leaves = []
    cur = root
    for depth in range(2 * k + 1):
        ax, ay, bx, by, cx, cy = cur
        if depth == 2 * k:
            keep, split = np.ones(len(ax), dtype=bool), np.zeros(len(ax), dtype=bool)
        else:
            mx, my = (ax + bx) >> 1, (ay + by) >> 1
            split = E[my, mx] > tol_mm
            keep = ~split
        if keep.any():
            lax, lay, lbx, lby, lcx, lcy = (a[keep] for a in cur)
            inside = cells[(lay + lby + lcy) // 3, (lax + lbx + lcx) // 3]          # cellule du centre de gravité
            if inside.any():
                leaves.append(tuple(a[inside] for a in (lax, lay, lbx, lby, lcx, lcy)))
        if depth < 2 * k and split.any():
            cur = _children(*(a[split] for a in cur))
        elif depth < 2 * k:
            break
    if not leaves:
        return np.empty((0, 3, 3), dtype=np.float64)

    ax, ay, bx, by, cx, cy = (np.concatenate([lf[n] for lf in leaves]) for n in range(6))
    # orientation : normale vers +z (échange b et c quand le produit vectoriel est négatif)
    cross = (X[bx] - X[ax]) * (Y[cy] - Y[ay]) - (Y[by] - Y[ay]) * (X[cx] - X[ax])
    flip = cross < 0
    bx2, by2 = np.where(flip, cx, bx), np.where(flip, cy, by)
    cx2, cy2 = np.where(flip, bx, cx), np.where(flip, by, cy)
    tris = np.empty((len(ax), 3, 3), dtype=np.float64)
    for v, (px, py) in enumerate(((ax, ay), (bx2, by2), (cx2, cy2))):
        tris[:, v, 0], tris[:, v, 1], tris[:, v, 2] = X[px], Y[py], z[py, px]
    return tris


def adaptive_top_triangles(quad_ok, Xmm, Ymm, z_top_mm, tol_mm, block=1024):
    """Triangles (N, 3, 3) du dessus, normales vers +z, pour les cellules de ``quad_ok`` (H, W).

    Xmm[j], Ymm[i] : coordonnées des nœuds ; z_top_mm[i, j] : altitudes des nœuds (H+1, W+1), en mm ;
    tol_mm : tolérance verticale en mm (>= 0) ; block : côté des blocs traités séparément, en cellules
    (puissance de deux). Voir le docstring du module."""
    if tol_mm < 0:
        raise ValueError("tol_mm doit être >= 0")
    if block < 2 or block & (block - 1):
        raise ValueError("block doit être une puissance de deux >= 2")
    quad_ok = np.asarray(quad_ok, dtype=bool)
    Xmm, Ymm = np.asarray(Xmm, dtype=np.float64), np.asarray(Ymm, dtype=np.float64)
    z = np.asarray(z_top_mm, dtype=np.float64)
    H, W = quad_ok.shape
    out = []
    for i0 in range(0, H, block):
        for j0 in range(0, W, block):
            sub = quad_ok[i0:i0 + block, j0:j0 + block]
            if not sub.any():
                continue
            h, w = sub.shape
            out.append(_rtin_block(sub, Xmm[j0:j0 + w + 1], Ymm[i0:i0 + h + 1], z[i0:i0 + h + 1, j0:j0 + w + 1], tol_mm))
    return np.concatenate(out) if out else np.empty((0, 3, 3), dtype=np.float64)
