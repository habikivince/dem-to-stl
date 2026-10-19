# SPDX-License-Identifier: GPL-3.0-or-later
"""
bottom.py
---------
Fond plat « fusionné » (PROTOTYPE, option --bottom merged).

Le fond du socle est un plan : le trianguler cellule par cellule (2 triangles par cellule, comme le
dessus) double le nombre de triangles sans rien apporter. Ici les cellules du masque sont couvertes
par des rectangles (couverture gloutonne), et chaque rectangle est triangulé en éventail depuis son
centre en ne gardant que les nœuds *requis* de son contour :

  * les nœuds du bord du masque (les parois verticales y sont accrochées, une par arête de cellule) ;
  * les quatre coins de chaque rectangle (les rectangles voisins s'y accrochent).

Deux rectangles voisins voient donc exactement les mêmes nœuds le long de leur arête commune : aucune
jonction en T, le solide reste étanche (chaque arête partagée par exactement deux triangles) et la
géométrie est identique à celle du fond cellule par cellule (même plan, même contour, même volume).

Deux couvertures sont disponibles :

  * ``graded`` (défaut) : carrés d'un quadtree équilibré (deux carrés voisins diffèrent d'un facteur 2
    au plus ; les cellules du bord restent des carrés d'une cellule). Les triangles restent compacts.
  * ``rows`` : rectangles gloutons ligne par ligne, c'est la méthode du dépôt my3dmaps. Très peu de
    triangles, mais des rectangles d'une cellule de haut sur des centaines de large : l'éventail depuis
    le centre y produit des triangles-aiguilles (mesuré : allongement jusqu'à 10^5 et angle minimal
    de l'ordre de 10^-3 degré sur un disque de 600 cellules). Conservée pour comparaison ;
    ``max_cells`` borne alors la taille des rectangles.

Idée inspirée du dépôt my3dmaps (licence MIT), ré-implémentée ici de façon indépendante.
"""
import numpy as np


def cover_with_rectangles(quad_ok, max_cells=None):
    """Couvre exactement les cellules vraies de ``quad_ok`` (H, W) par des rectangles disjoints.

    Balayage ligne par ligne : chaque plage horizontale de cellules libres est étendue vers le bas tant
    que les lignes suivantes ont les mêmes cellules libres sur toute la plage. ``max_cells`` (entier >= 1)
    borne le côté des rectangles (None : pas de borne). Renvoie une liste de (i0, i1, j0, j1), bornes
    incluses, en indices de cellules."""
    if max_cells is not None and max_cells < 1:
        raise ValueError("max_cells doit être >= 1 (ou None)")
    todo = np.array(quad_ok, dtype=bool, copy=True)
    H, W = todo.shape
    rects = []
    for i in range(H):
        row = todo[i]
        if not row.any():
            continue
        edges = np.diff(np.concatenate(([0], row.astype(np.int8), [0])))
        starts, ends = np.nonzero(edges == 1)[0], np.nonzero(edges == -1)[0] - 1
        for j0, j1 in zip(starts, ends):
            segs = [(j0, j1)] if max_cells is None else [
                (a, min(a + max_cells - 1, j1)) for a in range(j0, j1 + 1, max_cells)]
            for a, b in segs:
                i1 = i
                while i1 + 1 < H and (max_cells is None or i1 - i + 1 < max_cells) and todo[i1 + 1, a:b + 1].all():
                    i1 += 1
                todo[i:i1 + 1, a:b + 1] = False
                rects.append((i, i1, int(a), int(b)))
    return rects


def cover_with_graded_squares(quad_ok):
    """Couvre exactement les cellules vraies de ``quad_ok`` par des carrés de côté 2^k alignés (quadtree),
    avec la contrainte d'équilibre : deux carrés qui se touchent par un côté ont des tailles qui diffèrent
    d'un facteur 2 au plus. Un carré n'est conservé que s'il est entièrement dans le masque : les cellules du
    bord restent des carrés d'une cellule. Renvoie une liste de (i0, i1, j0, j1), bornes incluses."""
    quad_ok = np.asarray(quad_ok, dtype=bool)
    H, W = quad_ok.shape
    N = 1
    while N < max(H, W):
        N *= 2
    q = np.zeros((N, N), dtype=bool)
    q[:H, :W] = quad_ok
    sat = np.zeros((N + 1, N + 1), dtype=np.int64)          # table des sommes cumulées : comptage O(1) par bloc
    sat[1:, 1:] = q.cumsum(0).cumsum(1)

    def count(i, j, s):
        return sat[i + s, j + s] - sat[i, j + s] - sat[i + s, j] + sat[i, j]

    leaves, stack = [], [(0, 0, N)]
    while stack:
        i, j, s = stack.pop()
        c = count(i, j, s)
        if c == 0:
            continue
        if c == s * s:
            leaves.append((i, j, s))
        else:                                               # s >= 2 ici : un carré d'une cellule est plein ou vide
            h = s // 2
            stack += [(i, j, h), (i, j + h, h), (i + h, j, h), (i + h, j + h, h)]

    size_map = np.zeros((N, N), dtype=np.int32)
    for i, j, s in leaves:
        size_map[i:i + s, j:j + s] = s
    big = N * 4
    changed = True
    while changed:                                          # équilibrage : on coupe un carré dont un voisin est plus de 2 fois plus petit
        changed, nxt = False, []
        for i, j, s in leaves:
            if s == 1:
                nxt.append((i, j, s))
                continue
            sides = []
            if i > 0:
                sides.append(size_map[i - 1, j:j + s])
            if i + s < N:
                sides.append(size_map[i + s, j:j + s])
            if j > 0:
                sides.append(size_map[i:i + s, j - 1])
            if j + s < N:
                sides.append(size_map[i:i + s, j + s])
            smallest = min((int(np.where(v > 0, v, big).min()) for v in sides), default=big)
            if smallest * 2 < s:
                h = s // 2
                size_map[i:i + s, j:j + s] = h
                nxt += [(i, j, h), (i, j + h, h), (i + h, j, h), (i + h, j + h, h)]
                changed = True
            else:
                nxt.append((i, j, s))
        leaves = nxt
    return [(i, i + s - 1, j, j + s - 1) for i, j, s in leaves]


def required_nodes(quad_ok, rects):
    """Masque booléen (H+1, W+1) des nœuds à conserver : nœuds dont les quatre cellules voisines ne sont
    pas toutes dans le masque (bord, trous), plus les quatre coins de chaque rectangle."""
    q = np.pad(np.asarray(quad_ok, dtype=bool), 1, constant_values=False)
    req = ~(q[1:, 1:] & q[1:, :-1] & q[:-1, :-1] & q[:-1, 1:])
    for i0, i1, j0, j1 in rects:
        req[i0, j0] = req[i0, j1 + 1] = req[i1 + 1, j1 + 1] = req[i1 + 1, j0] = True
    return req


COVERS = ("graded", "rows")


def merged_bottom_triangles(quad_ok, Xmm, Ymm, z_mm, max_cells=None, cover="graded"):
    """Triangles (N, 3, 3) du fond plat à l'altitude ``z_mm``, normales vers -z, couvrant exactement les
    cellules de ``quad_ok``. Xmm[j] et Ymm[i] sont les coordonnées des nœuds (colonne j, ligne i), comme dans
    mesh.build_mesh_from_arrays. ``cover`` : "graded" (défaut, carrés équilibrés) ou "rows" (rectangles gloutons,
    ``max_cells`` borne alors leur taille). Voir le docstring du module."""
    if cover == "graded":
        rects = cover_with_graded_squares(quad_ok)
    elif cover == "rows":
        rects = cover_with_rectangles(quad_ok, max_cells)
    else:
        raise ValueError(f"cover : valeur inconnue {cover!r} (choix : {', '.join(COVERS)})")
    req = required_nodes(quad_ok, rects)
    tris = []
    for i0, i1, j0, j1 in rects:
        # contour dans le sens horaire vu d'en haut (x vers la droite, y vers le haut ; i croît vers le bas) :
        # haut de gauche à droite, droite de haut en bas, bas de droite à gauche, gauche de bas en haut.
        top = [(i0, j) for j in range(j0, j1 + 1) if req[i0, j]]
        right = [(i, j1 + 1) for i in range(i0, i1 + 1) if req[i, j1 + 1]]
        bottom = [(i1 + 1, j) for j in range(j1 + 1, j0, -1) if req[i1 + 1, j]]
        left = [(i, j0) for i in range(i1 + 1, i0, -1) if req[i, j0]]
        ring = top + right + bottom + left
        pts = np.array([(Xmm[j], Ymm[i], z_mm) for i, j in ring], dtype=np.float64)
        if len(ring) == 4:                       # coins seulement : deux triangles suffisent
            tris.append(np.stack([pts[[0, 1, 2]], pts[[0, 2, 3]]]))
            continue
        centre = np.array([(Xmm[j0] + Xmm[j1 + 1]) / 2.0, (Ymm[i0] + Ymm[i1 + 1]) / 2.0, z_mm])
        nxt = np.roll(pts, -1, axis=0)
        tris.append(np.stack([np.broadcast_to(centre, pts.shape), pts, nxt], axis=1))
    if not tris:
        return np.empty((0, 3, 3), dtype=np.float64)
    return np.concatenate(tris)
