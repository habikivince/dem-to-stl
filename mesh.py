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

from adaptive_top import adaptive_top_triangles
from bottom import merged_bottom_triangles


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


def drop_degenerate(tris):
    """Retire les triangles dont deux sommets sont identiques une fois arrondis en float32 (ce que
    voient les fichiers de sortie). Ils n'ont ni aire ni volume ; conservés, ils font apparaître des
    arêtes non-manifold dès qu'on fusionne les sommets (3MF, OBJ). Renvoie (triangles, nombre retiré)."""
    tris = np.asarray(tris, dtype=np.float64).reshape(-1, 3, 3)
    t = tris.astype(np.float32)
    same = (np.all(t[:, 0] == t[:, 1], axis=1) | np.all(t[:, 1] == t[:, 2], axis=1)
            | np.all(t[:, 0] == t[:, 2], axis=1))
    return tris[~same], int(same.sum())


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
    # Le sens de parcours du contour (trace_boundary_loop) est arbitraire, or l'orientation des
    # triangles de la bande et de la paroi en dépend : sans cette normalisation, ils sortaient
    # retournés (normales vers l'intérieur, volume signé négatif, arêtes mal orientées). Les
    # triangles ci-dessous sont orientés vers l'extérieur quand la boucle est parcourue dans le
    # sens horaire en (x, y), vérifié par check_stl.py et lib3mf (aire signée négative).
    area2 = sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(loop_xyz, loop_xyz[1:] + loop_xyz[:1]))
    if area2 > 0:
        loop_xyz.reverse()
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


def split_collinear_slivers(tris, rel_tol=1e-12):
    """Retire les triangles d'aire nulle à trois sommets DISTINCTS alignés (u, m, v avec m entre u et v) sans ouvrir
    le maillage : l'arête longue (u, v) est partagée avec un voisin (v, u, w), qui est scindé en (v, m, w) et
    (m, u, w). Même surface, même étanchéité, plus de triangle plat (ce que check_stl.py signale comme dégénéré).
    ``tris`` : liste de triangles, chacun trois tuples (x, y, z) ; les sommets partagés doivent être strictement
    identiques. Un triangle dont le voisin est absent de la liste ou lui-même plat est conservé tel quel.
    Les triangles à sommets confondus ne sont pas concernés (voir drop_degenerate)."""
    if not tris:
        return tris
    arr = np.asarray(tris, dtype=np.float64)
    e = [arr[:, (i + 1) % 3] - arr[:, i] for i in range(3)]          # e[i] : sommet i -> sommet i+1
    len2 = np.stack([np.einsum("ij,ij->i", v, v) for v in e], axis=1)
    longest = len2.max(axis=1)
    area2 = np.linalg.norm(np.cross(e[0], -e[2]), axis=1)
    flat = (area2 <= rel_tol * longest) & (longest > 0) & (len2.min(axis=1) > 0)
    if not flat.any():
        return tris
    owner = {}
    for idx, (p0, p1, p2) in enumerate(tris):
        owner[(p0, p1)] = idx
        owner[(p1, p2)] = idx
        owner[(p2, p0)] = idx
    drop, replace = set(), {}
    for i in np.nonzero(flat)[0]:
        t = tris[i]
        k = int(np.argmax(len2[i]))                    # arête la plus longue : de t[k] à t[k + 1], le milieu est t[k + 2]
        u, v, m = t[k], t[(k + 1) % 3], t[(k + 2) % 3]
        j = owner.get((v, u))
        if j is None or flat[j] or j in replace or j in drop:
            continue
        a, b, c = tris[j]
        s = [(a, b, c), (b, c, a), (c, a, b)][[(a, b), (b, c), (c, a)].index((v, u))]   # (v, u, w)
        w = s[2]
        drop.add(int(i))
        replace[j] = [(v, m, w), (m, u, w)]
    out = []
    for idx, t in enumerate(tris):
        if idx in drop:
            continue
        out.extend(replace[idx] if idx in replace else [t])
    return out


def build_outline_rim(loop_ij, Xmm, Ymm, z_top_mm, outline, base_mm, merge_tol_mm=1e-4):
    """Généralisation de build_smooth_rim à un contour convexe exact quelconque (polygone régulier, éventuellement
    à sommets arrondis : voir regular_polygon.py) au lieu du cercle. Même structure de retour :
    (band_tris, wall_tris, bottom_band_tris).

    ``outline`` fournit ``project(x, y)`` (projection sur le point du contour exact le plus proche : le contour est
    convexe et contient l'origine) et ``points`` (N, 2) : ses sommets exacts (angles vifs ou points des arcs
    d'arrondi), en mm comme Xmm/Ymm.

    Comme pour le cercle, chaque point du contour en escalier est projeté sur le contour, à la même hauteur. Mais
    contrairement au cercle, la projection doit être perpendiculaire et non radiale : le long d'un côté oblique, la
    projection radiale n'est pas monotone (voir RegularPolygon.project) et la paroi repasserait sur ses pas. Et une
    simple projection ne suffit pas pour un polygone : un angle vif tomberait entre deux points projetés (coin
    coupé en biais) et un petit arc pourrait disparaître. Tous les sommets exacts du contour sont donc insérés dans
    la chaîne projetée, à leur place angulaire (vue du centre, la chaîne projetée est monotone), avec une hauteur
    interpolée entre les deux points projetés qui l'encadrent. Un sommet exact à moins de ``merge_tol_mm`` du
    segment qui joint ses deux voisins projetés n'est pas inséré (évite des triangles quasi nuls ou alignés). Deux
    nœuds voisins qui se projettent au même point sont fusionnés. Les triangles de la bande, de la paroi et de la
    bande du fond sont découpés selon la
    chaîne complète : un éventail issu du point b du contour en escalier relie b à chaque segment de la chaîne
    entre les deux points projetés qui lui correspondent, sans jonction en T, donc étanche."""
    loop_xyz = [(Xmm[j], Ymm[i], z_top_mm[i, j]) for i, j in loop_ij]
    # sens horaire, comme build_smooth_rim (les triangles ci-dessous sont orientés vers l'extérieur dans ce sens)
    area2 = sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(loop_xyz, loop_xyz[1:] + loop_xyz[:1]))
    if area2 > 0:
        loop_xyz.reverse()
    n = len(loop_xyz)
    L = np.asarray(loop_xyz, dtype=np.float64)
    px, py = outline.project(L[:, 0], L[:, 1])
    P = np.stack([px, py, L[:, 2]], axis=1)
    on_outline = np.hypot(px - L[:, 0], py - L[:, 1]) <= 1e-9     # nœud déjà sur le contour : pas de bruit flottant
    P[on_outline, :2] = L[on_outline, :2]
    # Deux nœuds voisins qui se projettent au même point (marche perpendiculaire au côté) mais à des hauteurs
    # différentes donneraient une « ailette » verticale d'aire nulle dans la paroi : on leur donne la hauteur moyenne
    # du groupe et des coordonnées strictement identiques ; la paroi n'a plus alors de triangle sans aire (retiré
    # ensuite par drop_degenerate, comme pour le cercle).
    same_next = np.hypot(P[:, 0] - np.roll(P[:, 0], -1), P[:, 1] - np.roll(P[:, 1], -1)) <= 1e-9
    if same_next.any() and not same_next.all():
        start = int(np.nonzero(~same_next)[0][0]) + 1               # premier indice qui ouvre un groupe
        order = [(start + i) % n for i in range(n)]
        i = 0
        while i < n:
            j = i
            while j < n - 1 and same_next[order[j]]:
                j += 1
            if j > i:
                idx = order[i:j + 1]
                P[idx, 2] = P[idx, 2].mean()
                P[idx, :2] = P[idx[0], :2]          # coordonnées strictement identiques (pas de quasi-doublon)
            i = j + 1

    # angle polaire cumulé (sens horaire) des points projetés, ramené à [0, 2 pi] : sert à placer les sommets exacts
    two_pi = 2.0 * np.pi
    psi = np.arctan2(P[:, 1], P[:, 0])
    step = -((np.diff(np.concatenate([psi, psi[:1]])) + np.pi) % two_pi - np.pi)
    step = np.maximum(step, 0.0)                      # une marche d'escalier qui « recule » n'avance pas
    total = float(step.sum())
    if total < 1e-9:
        raise ValueError("contour exact : la boucle du masque ne fait pas le tour du centre de la forme.")
    s = np.concatenate([[0.0], np.cumsum(step[:-1])]) * (two_pi / total)

    E = np.asarray(outline.points, dtype=np.float64)
    t_e = (psi[0] - np.arctan2(E[:, 1], E[:, 0])) % two_pi
    gap = np.clip(np.searchsorted(s, t_e, side="right") - 1, 0, n - 1)
    nxt = (gap + 1) % n
    s_next = np.where(nxt == 0, two_pi, s[nxt])
    span = s_next - s[gap]
    frac = np.divide(t_e - s[gap], span, out=np.zeros_like(span), where=span > 0)
    z_e = P[gap, 2] + (P[nxt, 2] - P[gap, 2]) * frac
    # un sommet exact déjà sur le segment qui joint ses deux voisins projetés (à merge_tol_mm près) n'apporte rien :
    # l'insérer créerait un triangle aligné d'aire nulle
    A, B = P[gap, :2], P[nxt, :2]
    AB = B - A
    len2 = np.einsum("ij,ij->i", AB, AB)
    tt = np.divide(np.einsum("ij,ij->i", E - A, AB), len2, out=np.zeros_like(len2), where=len2 > 0)
    foot = A + np.clip(tt, 0.0, 1.0)[:, None] * AB
    near = np.hypot(E[:, 0] - foot[:, 0], E[:, 1] - foot[:, 1]) <= merge_tol_mm
    extra = {}
    for idx in np.lexsort((t_e, gap)):
        if not near[idx]:
            extra.setdefault(int(gap[idx]), []).append((float(E[idx, 0]), float(E[idx, 1]), float(z_e[idx])))

    zb = -base_mm
    Lt = [tuple(row) for row in L.tolist()]
    Pt = [tuple(row) for row in P.tolist()]
    band_tris, wall_tris, bottom_band_tris = [], [], []
    for k in range(n):
        k2 = (k + 1) % n
        a, b = Lt[k], Lt[k2]
        chain = [Pt[k]] + extra.get(k, []) + [Pt[k2]]
        a_base, b_base = (a[0], a[1], zb), (b[0], b[1], zb)
        band_tris.append((a, b, chain[0]))
        bottom_band_tris.append((b_base, a_base, (chain[0][0], chain[0][1], zb)))
        for q0, q1 in zip(chain[:-1], chain[1:]):
            q0_base, q1_base = (q0[0], q0[1], zb), (q1[0], q1[1], zb)
            band_tris.append((b, q1, q0))
            wall_tris.append((q0, q1, q0_base))
            wall_tris.append((q1, q1_base, q0_base))
            bottom_band_tris.append((q1_base, b_base, q0_base))

    # un nœud qui se projette dans l'alignement du côté de l'escalier donne, au fond (plan z = -base), un triangle
    # plat : on le remplace en scindant son voisin
    return (split_collinear_slivers(band_tris), split_collinear_slivers(wall_tris),
            split_collinear_slivers(bottom_band_tris))


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


def fix_diagonal_pinches(quad_ok, max_iterations=None):
    """Élimine les pincements en diagonale : deux cellules incluses qui ne se
    touchent que par un coin (motif "damier" sur un bloc 2x2) créent un
    sommet où deux nappes de surface se rejoignent en un point unique —
    épaisseur quasi nulle, indétectable par un simple test d'étanchéité par
    arêtes. On comble les deux cases manquantes plutôt que de retirer de
    la matière. Modifie quad_ok sur place et le retourne.

    Le comblement peut créer de nouveaux pincements voisins : on itère donc
    jusqu'à stabilité. Il ne fait qu'ajouter des cellules (monotone), donc la
    boucle termine en au plus quad_ok.size passes. ``max_iterations`` reste
    un plafond optionnel (None = jusqu'à stabilité)."""
    limit = quad_ok.size if max_iterations is None else max_iterations
    for _ in range(limit):
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


def build_plain_walls(quad_ok, Xmm, Ymm, z_top_mm, zb):
    """Parois verticales sur tout le bord du masque (deux triangles par arête de cellule), sans boucle Python.

    Même résultat, dans le même ordre, que l'ancienne double boucle sur tous les nœuds : d'abord les arêtes
    verticales (ligne par ligne), puis les arêtes horizontales (colonne par colonne). Une arête de bord sépare
    une cellule incluse d'une cellule exclue ou hors grille ; la paroi est orientée vers l'extérieur."""
    H1, W1 = quad_ok.shape
    cols = np.pad(quad_ok, ((0, 0), (1, 1)))                     # (H1, W1 + 2)
    left, right = cols[:, :-1], cols[:, 1:]                        # côtés gauche et droit de chaque arête verticale
    iv, jv = np.nonzero(left != right)                             # ordre : ligne puis colonne, comme la boucle
    rows = np.pad(quad_ok, ((1, 1), (0, 0)))                     # (H1 + 2, W1)
    up, down = rows[:-1, :], rows[1:, :]                           # côtés haut et bas de chaque arête horizontale
    jh, ih = np.nonzero((up != down).T)                            # ordre : colonne puis ligne, comme la boucle

    def walls(p1, p2, outward):
        n = len(p1)
        p1b, p2b = p1.copy(), p2.copy()
        p1b[:, 2] = zb
        p2b[:, 2] = zb
        out = np.empty((n, 2, 3, 3), dtype=np.float64)
        o = outward[:, None, None]
        out[:, 0] = np.where(o, np.stack([p1, p2, p1b], axis=1), np.stack([p2, p1, p1b], axis=1))
        out[:, 1] = np.where(o, np.stack([p2, p2b, p1b], axis=1), np.stack([p2b, p2, p1b], axis=1))
        return out.reshape(-1, 3, 3)

    p1v = np.stack([Xmm[jv], Ymm[iv], z_top_mm[iv, jv]], axis=1)
    p2v = np.stack([Xmm[jv], Ymm[iv + 1], z_top_mm[iv + 1, jv]], axis=1)
    p1h = np.stack([Xmm[jh], Ymm[ih], z_top_mm[ih, jh]], axis=1)
    p2h = np.stack([Xmm[jh + 1], Ymm[ih], z_top_mm[ih, jh + 1]], axis=1)
    return np.concatenate([walls(p1v, p2v, left[iv, jv]), walls(p1h, p2h, down[ih, jh])])


BOTTOM_MODES = ("grid", "merged", "merged-rows")
DIAGONAL_MODES = ("fixed", "adaptive")
TOP_MODES = ("grid", "adaptive")


def build_mesh_from_arrays(elev, quad_ok, xs, ys, scale, vexag, base_mm, min_elev,
                            progress=None, smooth_wall=False, circle_radius_mm=None,
                            bottom="grid", bottom_max_cells=None, diagonal="fixed",
                            top="grid", top_tol_mm=0.05, outline=None):
    """Construit le maillage (triangles dessus/dessous/parois) à partir d'un
    tableau d'altitudes déjà découpé et d'un masque de cellules (quad_ok)
    déjà nettoyé. Ne touche jamais au disque ni à GDAL — testable avec de
    simples tableaux numpy. quad_ok doit avoir la forme (H-1, W-1) par
    rapport à elev. xs/ys doivent déjà être recentrées sur le centre du
    cercle (en unités du CRS source, pas encore mises à l'échelle).

    smooth_wall=True : au lieu de suivre le contour en escalier du masque,
    chaque point du contour est projeté radialement vers le vrai cercle
    (même altitude), formant un bord parfaitement circulaire. Requiert un
    masque à un seul contour simple (voir trace_boundary_loop).

    bottom="grid" (défaut) : fond triangulé cellule par cellule (comme le dessus).
    bottom="merged" (PROTOTYPE) : fond plat couvert par des carrés d'un quadtree équilibré, triangulés en
    éventail sans jonction en T (voir bottom.py) ; même géométrie, environ moitié moins de triangles au total,
    triangles compacts. bottom="merged-rows" : couverture par rectangles gloutons (très peu de triangles mais
    triangles-aiguilles) ; bottom_max_cells borne alors la taille des rectangles (None : pas de borne).

    diagonal="fixed" (défaut) : chaque cellule du dessus est coupée selon la même diagonale (c10-c01).
    diagonal="adaptive" (PROTOTYPE) : la diagonale qui relie les deux coins de plus faible écart d'altitude
    (c00-c11 si |z00 - z11| <= |z10 - z01|, sinon c10-c01). Même nombre de triangles, mêmes sommets, mêmes
    parois et même fond ; seule la surface change (mesure : au centre des cellules, erreur moyenne réduite
    de 4 à 7 % sur six relevés réels, au pire +2,5 %, voir rapport 9).

    top="grid" (défaut) : un ou deux triangles par cellule, exact aux nœuds.
    top="adaptive" (PROTOTYPE) : dessus en maillage RTIN à erreur bornée (voir adaptive_top.py) : bien moins de
    triangles, à ``top_tol_mm`` mm de tolérance verticale ; le bord du masque reste à pleine résolution, donc
    les parois et le fond sont inchangés. Incompatible avec diagonal="adaptive" (sans objet).

    outline (None par défaut) : contour exact convexe (polygone régulier, voir regular_polygon.py et
    build_outline_rim), en mm. Même principe que smooth_wall, mais le contour en escalier est projeté sur ce
    contour au lieu du cercle, et ses sommets exacts (angles, arcs d'arrondi) sont insérés dans la paroi. Implique
    un bord lisse : smooth_wall et circle_radius_mm sont alors ignorés. Sans outline, rien ne change."""
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

    if top not in TOP_MODES:
        raise ValueError(f"top : valeur inconnue {top!r} (choix : {', '.join(TOP_MODES)})")
    if top == "adaptive" and diagonal != "fixed":
        raise ValueError("diagonal='adaptive' est sans objet avec top='adaptive'")
    if top == "grid" or bottom == "grid":
        c00, c10, c01, c11 = corner(0, 0), corner(1, 0), corner(0, 1), corner(1, 1)

    if top == "adaptive":
        top_tris = adaptive_top_triangles(quad_ok, Xmm, Ymm, z_top_mm, top_tol_mm)
    elif diagonal == "fixed":
        top_tris = np.concatenate([
            np.stack([c00, c10, c01], axis=1),
            np.stack([c10, c11, c01], axis=1),
        ])
    elif diagonal == "adaptive":
        # même structure (premier triangle de chaque cellule, puis second) ; diagonale c00-c11 quand elle relie
        # les deux coins de plus faible écart d'altitude (orientation vers +z dans les deux cas)
        ac = np.abs(c00[:, 2] - c11[:, 2]) <= np.abs(c10[:, 2] - c01[:, 2])
        first = np.stack([c00, c10, np.where(ac[:, None], c11, c01)], axis=1)
        second = np.where(ac[:, None, None], np.stack([c00, c11, c01], axis=1), np.stack([c10, c11, c01], axis=1))
        top_tris = np.concatenate([first, second])
    else:
        raise ValueError(f"diagonal : valeur inconnue {diagonal!r} (choix : {', '.join(DIAGONAL_MODES)})")

    zb = -base_mm
    if bottom == "merged":
        bot_tris = merged_bottom_triangles(quad_ok, Xmm, Ymm, zb, cover="graded")
    elif bottom == "merged-rows":
        bot_tris = merged_bottom_triangles(quad_ok, Xmm, Ymm, zb, bottom_max_cells, cover="rows")
    elif bottom == "grid":
        c00b, c10b, c01b, c11b = c00.copy(), c10.copy(), c01.copy(), c11.copy()
        for c in (c00b, c10b, c01b, c11b):
            c[:, 2] = zb

        bot_tris = np.concatenate([
            np.stack([c00b, c01b, c10b], axis=1),
            np.stack([c10b, c01b, c11b], axis=1),
        ])
    else:
        raise ValueError(f"bottom : valeur inconnue {bottom!r} (choix : {', '.join(BOTTOM_MODES)})")

    if smooth_wall or outline is not None:
        loop_ij = trace_boundary_loop(quad_ok)
        padded = np.pad(quad_ok, 1, constant_values=False)
        n_boundary = int(sum((quad_ok & ~padded[a, b]).sum() for a, b in (
            (slice(None, -2), slice(1, -1)), (slice(2, None), slice(1, -1)),
            (slice(1, -1), slice(None, -2)), (slice(1, -1), slice(2, None)))))
        if n_boundary != len(loop_ij):
            if outline is not None:
                sys.exit("Erreur : la forme (--sides) requiert un seul contour : le masque a un trou (zone sans "
                         "donnée à l'intérieur de la forme) ou plusieurs parties. Relance avec --sea-level si le "
                         "trou est de la mer, ou avec --fill-with-copernicus.")
            sys.exit("Erreur : --smooth-wall requiert un seul contour : le masque a un trou (zone sans donnée "
                     "à l'intérieur du disque) ou plusieurs parties. Relance sans --smooth-wall, ou avec "
                     "--sea-level si le trou est de la mer.")
        if outline is not None:
            band_tris, wall_tris, bottom_band_tris = build_outline_rim(
                loop_ij, Xmm, Ymm, z_top_mm, outline, base_mm)
        else:
            if circle_radius_mm is None:
                raise ValueError("smooth_wall=True nécessite circle_radius_mm (le vrai rayon, en mm).")
            radius_mm = circle_radius_mm
            band_tris, wall_tris, bottom_band_tris = build_smooth_rim(
                loop_ij, Xmm, Ymm, z_top_mm, radius_mm, base_mm)
        # un sommet du contour déjà exactement sur le cercle donne des triangles sans aire : on les retire
        extra_tris, _ = drop_degenerate(band_tris + bottom_band_tris)
        wall_tris, _ = drop_degenerate(wall_tris)
        all_tris = np.concatenate([top_tris, bot_tris, extra_tris, wall_tris])
        return all_tris, len(top_tris), len(bot_tris) + len(bottom_band_tris) + len(band_tris), len(wall_tris)

    H1, W1 = quad_ok.shape
    _progress(0, 2, "Parois (bord vertical)")
    wall_tris = build_plain_walls(quad_ok, Xmm, Ymm, z_top_mm, zb)
    _progress(1, 2, "Parois (bord horizontal)")
    all_tris = np.concatenate([top_tris, bot_tris, wall_tris])
    return all_tris, len(top_tris), len(bot_tris), len(wall_tris)
