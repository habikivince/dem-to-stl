#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
check_stl.py
------------
Contrôle indépendant d'un STL binaire (ou d'un 3MF / OBJ) produit par circle_dem_to_stl.py : dimensions,
épaisseur du socle, relief implicite, étanchéité du maillage, orientation (volume signé), triangles
dégénérés, topologie (nombre de coques, caractéristique d'Euler V − E + F) et triangles retournés.

Le contrôle ne réutilise pas le code de maillage : il relit le fichier et compare ce
qu'il contient à ce que les paramètres imposent.

Exemple (Cervin, 150 mm, rayon 3 km, socle 3 mm) :
  python check_stl.py cervin.stl --diameter-mm 150 --base-mm 3 --radius 3000 --relief-m 2247

Sans --radius (mode polygone) l'échelle n'est pas connue : seuls les contrôles de
dimension, de socle et d'étanchéité sont faits.

Pièce générée avec --sides (polygone régulier) : --diameter-mm reste le diamètre du cercle de référence et les
dimensions attendues sont celles de l'encombrement exact du polygone (regular_polygon.py, géométrie seule : le
code de maillage n'est toujours pas réutilisé) :
  python check_stl.py hexagone.stl --diameter-mm 60 --base-mm 3 --radius 300 \
      --sides 6 --circle-mode circumscribed --corner-radius 2
"""
import argparse
import sys

import numpy as np

STL_DTYPE = np.dtype([("n", "<f4", (3,)), ("v", "<f4", (3, 3)), ("a", "<u2")])


def read_stl(path):
    with open(path, "rb") as f:
        data = f.read()
    if len(data) < 84:
        raise ValueError("fichier trop court pour un STL binaire")
    n = int(np.frombuffer(data, dtype="<u4", count=1, offset=80)[0])
    if len(data) != 84 + 50 * n:
        raise ValueError(f"taille incohérente avec un STL binaire ({n} triangles annoncés, "
                         f"{len(data)} octets) — STL ASCII ou fichier tronqué ?")
    return np.frombuffer(data, dtype=STL_DTYPE, count=n, offset=84)["v"].astype(np.float64)


def _weld(tris):
    """Soude les sommets identiques une fois arrondis en float32 (ce que voient les fichiers de sortie).
    Retourne (ids, n_vert) : ids (T, 3) = indice de sommet soudé de chaque coin de triangle."""
    flat = tris.reshape(-1, 3).astype(np.float32)
    order = np.lexsort((flat[:, 2], flat[:, 1], flat[:, 0]))
    sorted_pts = flat[order]
    new = np.ones(len(flat), dtype=bool)
    new[1:] = np.any(sorted_pts[1:] != sorted_pts[:-1], axis=1)
    ids = np.empty(len(flat), dtype=np.int64)
    ids[order] = np.cumsum(new) - 1
    return ids.reshape(-1, 3), int(ids.max()) + 1


def edge_report(tris, welded=None):
    """Étanchéité : nombre d'arêtes utilisées une seule fois (trous), plus de deux fois
    (non-manifold), et d'arêtes orientées en double (faces mal orientées).
    welded : résultat de _weld(tris) si l'appelant l'a déjà calculé (évite de souder deux fois)."""
    ids, n_vert = _weld(tris) if welded is None else welded
    a = np.concatenate([ids[:, 0], ids[:, 1], ids[:, 2]])
    b = np.concatenate([ids[:, 1], ids[:, 2], ids[:, 0]])
    directed = a * n_vert + b
    undirected = np.minimum(a, b) * n_vert + np.maximum(a, b)
    _, counts = np.unique(undirected, return_counts=True)
    _, dcounts = np.unique(directed, return_counts=True)
    return {"vertices": n_vert, "boundary_edges": int((counts == 1).sum()),
            "nonmanifold_edges": int((counts > 2).sum()), "duplicate_directed_edges": int((dcounts > 1).sum())}


def topology_report(tris, welded=None, min_quality=1e-3, oppose_dot=-0.5):
    """Topologie du maillage soudé, en une seule passe sur les arêtes.

    Coques : composantes connexes par ARÊTE commune (deux solides qui ne se touchent que par un sommet sont
    deux coques). Pour chaque coque, V − E + F : 2 pour une surface fermée sans anse (2 − 2g pour g anses ; les
    arêtes ouvertes ou non-manifold la faussent, mais l'étanchéité les signale déjà). `euler` est la somme
    des caractéristiques par coque ; `shared_vertices` compte les sommets partagés entre deux coques (contact
    par un point). Les triangles dont deux coins soudés sont identiques sont ignorés (`ignored_faces`).

    Triangles retournés : un triangle est « retourné » si ses trois arêtes sont manifold et orientées de façon
    cohérente avec le voisin, et si sa normale fait avec celles des TROIS voisins un produit scalaire < oppose_dot
    (-0,5 : plus de 120°). Un retournement géométrique (sommet déplacé de l'autre côté de ses voisins, par ex.
    par une projection de contour) garde la cohérence des arêtes, donc échappe à l'étanchéité et au volume
    signé. Une arête vive isolée (crête, falaise au bord, forte exagération verticale) n'est PAS un
    retournement : au moins un voisin reste dans le même sens. Mesure sur 14 masques aléatoires (îles, trous,
    goulets) maillés par mesh.py, 4 variantes de maillage chacun : aucun signalement jusqu'à une pente
    moyenne du dessus de 82° (90e centile 87°, vexag 8 inclus) ; 1 solide sur 14 à 89° de pente moyenne
    (terrain quasi vertical partout). Un retournement n'est vu que si l'angle réel entre le triangle et
    ses voisins est inférieur à 60°. Les triangles-aiguilles (2·aire / plus grand côté² <
    min_quality) ont une normale peu fiable : ni jugés, ni utilisés comme voisins.
    `inverted_examples` : centres (x, y, z) des 5 premiers triangles retournés.

    Ne détecte pas les auto-intersections (il faudrait un arbre de boîtes englobantes)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    ids, n_vert = _weld(tris) if welded is None else welded
    valid = (ids[:, 0] != ids[:, 1]) & (ids[:, 1] != ids[:, 2]) & (ids[:, 0] != ids[:, 2])
    ignored = int((~valid).sum())
    ids = ids[valid]
    n_faces = len(ids)
    empty = {"faces": 0, "vertices": 0, "edges": 0, "shells": 0, "shell_euler": [], "euler": 0,
             "shared_vertices": 0, "ignored_faces": ignored, "inverted_faces": 0, "inverted_examples": []}
    if n_faces == 0:
        return empty
    t = np.asarray(tris, dtype=np.float64)[valid]

    a = np.concatenate([ids[:, 0], ids[:, 1], ids[:, 2]])
    b = np.concatenate([ids[:, 1], ids[:, 2], ids[:, 0]])
    face = np.tile(np.arange(n_faces), 3)
    key = np.minimum(a, b) * n_vert + np.maximum(a, b)
    order = np.argsort(key, kind="stable")
    ks, fs, fwd = key[order], face[order], (a < b)[order]
    same = ks[1:] == ks[:-1]                            # entrée i et i+1 : même arête
    first = np.concatenate([[True], ~same])             # première entrée de chaque arête

    # coques : graphe face-face par arête commune (une chaîne par arête suffit même au-delà de 2 faces)
    graph = coo_matrix((np.ones(int(same.sum()), dtype=np.int8), (fs[:-1][same], fs[1:][same])),
                       shape=(n_faces, n_faces))
    n_shells, labels = connected_components(graph, directed=False)
    faces_per = np.bincount(labels, minlength=n_shells)
    edges_per = np.bincount(labels[fs[first]], minlength=n_shells)
    pair = np.unique((labels[:, None].astype(np.int64) * n_vert + ids).ravel())
    verts_per = np.bincount(pair // n_vert, minlength=n_shells)
    shell_euler = (verts_per - edges_per + faces_per).astype(int)
    n_used = len(np.unique(ids))

    # triangles retournés
    nvec = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    twice_area = np.linalg.norm(nvec, axis=1)
    longest2 = np.maximum.reduce([np.einsum("ij,ij->i", d, d)
                                  for d in (t[:, 1] - t[:, 0], t[:, 2] - t[:, 1], t[:, 0] - t[:, 2])])
    reliable = twice_area >= min_quality * np.maximum(longest2, 1e-300)
    unit = np.divide(nvec, twice_area[:, None], out=np.zeros_like(nvec), where=twice_area[:, None] > 0)
    run_start = np.flatnonzero(first)
    run_len = np.diff(np.concatenate([run_start, [len(ks)]]))
    two = run_start[run_len == 2]                       # arêtes manifold (exactement deux faces)
    f1, f2 = fs[two], fs[two + 1]
    usable = (fwd[two] != fwd[two + 1]) & reliable[f1] & reliable[f2]     # orientées de façon cohérente
    f1, f2 = f1[usable], f2[usable]
    opposed = np.einsum("ij,ij->i", unit[f1], unit[f2]) < oppose_dot
    n_ok = np.bincount(f1, minlength=n_faces) + np.bincount(f2, minlength=n_faces)
    n_opp = np.bincount(f1[opposed], minlength=n_faces) + np.bincount(f2[opposed], minlength=n_faces)
    inverted = np.flatnonzero((n_ok == 3) & (n_opp == 3) & reliable)
    return {"faces": n_faces, "vertices": n_used, "edges": int(first.sum()), "shells": int(n_shells),
            "shell_euler": [int(x) for x in shell_euler[:20]], "euler": int(shell_euler.sum()),
            "shared_vertices": int(verts_per.sum()) - n_used, "ignored_faces": ignored,
            "inverted_faces": int(len(inverted)),
            "inverted_examples": [tuple(float(x) for x in t[i].mean(axis=0)) for i in inverted[:5]]}


def signed_volume(tris, block=1_000_000):
    """Volume signé (mm³) d'un maillage de triangles (N, 3, 3), par la somme des produits mixtes
    v0 · (v1 x v2) / 6. Positif quand les faces sont tournées vers l'extérieur du solide. Un maillage
    entièrement retourné reste « étanche » et cohérent arête par arête : seul ce signe le révèle."""
    total = 0.0
    for i in range(0, len(tris), block):
        t = tris[i:i + block]
        total += float(np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum())
    return total / 6.0


def degenerate_count(tris, block=1_000_000):
    """Nombre de triangles d'aire nulle : deux sommets identiques une fois arrondis en float32 (ce
    que voient les fichiers de sortie) ou trois sommets exactement alignés (produit vectoriel nul)."""
    n = 0
    for i in range(0, len(tris), block):
        t = tris[i:i + block]
        f = t.astype(np.float32)
        same = (np.all(f[:, 0] == f[:, 1], axis=1) | np.all(f[:, 1] == f[:, 2], axis=1)
                | np.all(f[:, 0] == f[:, 2], axis=1))
        cross = np.cross(f[:, 1].astype(np.float64) - f[:, 0], f[:, 2].astype(np.float64) - f[:, 0])
        n += int((same | ~np.any(cross != 0.0, axis=1)).sum())
    return n


def check(tris, diameter_mm, base_mm, radius_m=None, vexag=1.0, relief_m=None, tol_mm=None, watertight=True,
          shape_mm=None, topology=False, expect_shells=None):
    """Liste de (nom, ok, détail). tol_mm : tolérance sur les dimensions (défaut 1 % du diamètre).
    shape_mm : RegularPolygon (regular_polygon.py) exprimé en mm imprimés, pour une pièce générée avec --sides :
    les dimensions attendues sont alors celles de son encombrement exact (et non plus le diamètre).
    topology : ajoute (avec l'étanchéité) le contrôle des coques / de la caractéristique d'Euler et celui des
    triangles retournés (voir topology_report). Désactivé par défaut pour ne rien changer aux appelants existants ;
    la ligne de commande l'active. expect_shells : nombre de coques attendu ; chaque coque doit alors avoir
    V − E + F = 2 (surface fermée sans anse). Sans expect_shells, la topologie est seulement décrite : plusieurs
    îles ou un trou dans le masque donnent légitimement plusieurs coques ou une anse ; seul un sommet partagé
    entre deux coques (contact par un point) est refusé."""
    tol_mm = 0.01 * diameter_mm if tol_mm is None else tol_mm
    mins, maxs = tris.reshape(-1, 3).min(axis=0), tris.reshape(-1, 3).max(axis=0)
    dx, dy = maxs[0] - mins[0], maxs[1] - mins[1]
    zmax = maxs[2]
    results = []
    if shape_mm is None:
        results.append(("plus grande dimension horizontale",
                        abs(max(dx, dy) - diameter_mm) <= tol_mm,
                        f"{max(dx, dy):.2f} mm (attendu {diameter_mm:g} ± {tol_mm:.2f})"))
    else:
        x0, x1, y0, y1 = shape_mm.bbox()
        w, h = x1 - x0, y1 - y0
        results.append(("dimensions de la pièce (forme --sides)",
                        abs(dx - w) <= tol_mm and abs(dy - h) <= tol_mm,
                        f"{dx:.2f} x {dy:.2f} mm (attendu {w:.2f} x {h:.2f} ± {tol_mm:.2f})"))
    results.append(("socle", abs(mins[2] + base_mm) <= 0.01,
                    f"fond à z = {mins[2]:.3f} mm (attendu {-base_mm:g})"))
    results.append(("hauteur totale", True, f"{zmax - mins[2]:.2f} mm (sommet z = {zmax:.2f} mm au-dessus de l'altitude minimale)"))
    if radius_m:
        scale = diameter_mm / (2.0 * radius_m)
        implied = zmax / (scale * vexag)
        ok = True if relief_m is None else abs(implied - relief_m) <= max(0.01 * relief_m, 1.0)
        detail = f"relief implicite {implied:.1f} m (échelle {scale:.5f} mm/m, exagération {vexag:g})"
        if relief_m is not None:
            detail += f", attendu {relief_m:g} m"
        results.append(("relief vertical", ok, detail))
    if watertight:
        welded = _weld(tris)
        r = edge_report(tris, welded)
        ok = r["boundary_edges"] == 0 and r["nonmanifold_edges"] == 0 and r["duplicate_directed_edges"] == 0
        results.append(("étanchéité", ok,
                        f"{r['vertices']} sommets, {r['boundary_edges']} arête(s) ouverte(s), "
                        f"{r['nonmanifold_edges']} non-manifold, {r['duplicate_directed_edges']} mal orientée(s)"))
        vol = signed_volume(tris)
        results.append(("orientation (volume signé)", vol > 0,
                        f"{vol:.1f} mm³" + (" (faces tournées vers l'extérieur)" if vol > 0
                                            else " : NÉGATIF ou nul, faces tournées vers l'intérieur")))
        nd = degenerate_count(tris)
        results.append(("triangles dégénérés", nd == 0, f"{nd} triangle(s) d'aire nulle"))
        if topology:
            topo = topology_report(tris, welded)
            per = topo["shell_euler"]
            detail = f"{topo['shells']} coque(s), V − E + F par coque = {per}"
            if topo["shared_vertices"]:
                detail += f", {topo['shared_vertices']} sommet(s) partagé(s) entre coques (contact par un point)"
            if expect_shells is None:
                ok = topo["shared_vertices"] == 0
                if all(e == 2 for e in per):
                    detail += " (surfaces fermées sans anse)"
                else:
                    detail += " (2 = surface fermée sans anse ; moins = anse(s) ou surface ouverte)"
            else:
                ok = topo["shells"] == expect_shells and all(e == 2 for e in per) and topo["shared_vertices"] == 0
                detail += f" (attendu {expect_shells} coque(s), chacune à 2)"
            results.append(("topologie (coques, Euler)", ok, detail))
            ex = topo["inverted_examples"]
            where = "" if not ex else " ; premier près de (" + ", ".join(f"{c:.2f}" for c in ex[0]) + ") mm"
            results.append(("triangles retournés", topo["inverted_faces"] == 0,
                            f"{topo['inverted_faces']} triangle(s) opposé(s) à leurs trois voisins{where}"))
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("stl", help="fichier .stl (binaire), .3mf ou .obj")
    p.add_argument("--diameter-mm", type=float, required=True, help="valeur donnée à circle_dem_to_stl.py")
    p.add_argument("--base-mm", type=float, required=True)
    p.add_argument("--radius", type=float, default=None, help="rayon en mètres (mode cercle) : active le contrôle du relief")
    p.add_argument("--vexag", type=float, default=1.0)
    p.add_argument("--relief-m", type=float, default=None,
                   help="dénivelé attendu en mètres (ex. maximum - minimum de gdalinfo -stats), tolérance 1 %%")
    p.add_argument("--tol-mm", type=float, default=None, help="tolérance sur les dimensions (défaut 1 %% du diamètre)")
    p.add_argument("--no-watertight", action="store_true", help="saute le test d'étanchéité (gros fichiers)")
    p.add_argument("--no-topology", action="store_true",
                   help="saute le contrôle des coques / d'Euler et des triangles retournés (inclus par défaut avec l'étanchéité)")
    p.add_argument("--expect-shells", type=int, default=None,
                   help="nombre de coques attendu (chacune doit avoir V - E + F = 2, sans anse) ; sans cette option "
                        "la topologie est seulement décrite, plusieurs îles ou un trou étant légitimes")
    p.add_argument("--sides", type=int, default=None,
                   help="pièce générée avec --sides : contrôle l'encombrement exact du polygone (avec --circle-mode, "
                        "--corner-radius en mm, --rotation ; --diameter-mm reste le diamètre du cercle de référence)")
    p.add_argument("--circle-mode", choices=("inscribed", "circumscribed"), default=None)
    p.add_argument("--corner-radius", type=float, default=0.0, help="rayon d'arrondi en mm imprimés")
    p.add_argument("--rotation", type=float, default=0.0)
    args = p.parse_args()
    shape_mm = None
    if args.sides is not None:
        if args.circle_mode is None:
            sys.exit("Erreur : --sides requiert --circle-mode.")
        import regular_polygon
        try:
            shape_mm = regular_polygon.RegularPolygon(args.sides, args.circle_mode, args.diameter_mm / 2.0,
                                                      args.corner_radius, args.rotation)
        except regular_polygon.ShapeError as e:
            sys.exit(f"Erreur : {e}")
    try:
        import mesh_export
        tris = mesh_export.read_triangles(args.stl)
    except (OSError, ValueError, KeyError) as e:
        sys.exit(f"Erreur : {e}")
    print(f"{args.stl} : {len(tris)} triangles")
    results = check(tris, args.diameter_mm, args.base_mm, args.radius, args.vexag, args.relief_m,
                    args.tol_mm, not args.no_watertight, shape_mm, topology=not args.no_topology,
                    expect_shells=args.expect_shells)
    for name, ok, detail in results:
        print(f"[{'OK ' if ok else 'ÉCHEC'}] {name} : {detail}")
    sys.exit(0 if all(ok for _, ok, _ in results) else 1)


if __name__ == "__main__":
    main()
