#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
check_stl.py
------------
Contrôle indépendant d'un STL binaire (ou d'un 3MF / OBJ) produit par circle_dem_to_stl.py : dimensions,
épaisseur du socle, relief implicite, étanchéité du maillage, orientation (volume signé) et triangles
dégénérés.

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


def edge_report(tris):
    """Étanchéité : nombre d'arêtes utilisées une seule fois (trous), plus de deux fois
    (non-manifold), et d'arêtes orientées en double (faces mal orientées)."""
    flat = tris.reshape(-1, 3).astype(np.float32)
    order = np.lexsort((flat[:, 2], flat[:, 1], flat[:, 0]))
    sorted_pts = flat[order]
    new = np.ones(len(flat), dtype=bool)
    new[1:] = np.any(sorted_pts[1:] != sorted_pts[:-1], axis=1)
    ids = np.empty(len(flat), dtype=np.int64)
    ids[order] = np.cumsum(new) - 1
    n_vert = int(ids.max()) + 1
    ids = ids.reshape(-1, 3)
    a = np.concatenate([ids[:, 0], ids[:, 1], ids[:, 2]])
    b = np.concatenate([ids[:, 1], ids[:, 2], ids[:, 0]])
    directed = a * n_vert + b
    undirected = np.minimum(a, b) * n_vert + np.maximum(a, b)
    _, counts = np.unique(undirected, return_counts=True)
    _, dcounts = np.unique(directed, return_counts=True)
    return {"vertices": n_vert, "boundary_edges": int((counts == 1).sum()),
            "nonmanifold_edges": int((counts > 2).sum()), "duplicate_directed_edges": int((dcounts > 1).sum())}


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
          shape_mm=None):
    """Liste de (nom, ok, détail). tol_mm : tolérance sur les dimensions (défaut 1 % du diamètre).
    shape_mm : RegularPolygon (regular_polygon.py) exprimé en mm imprimés, pour une pièce générée avec --sides :
    les dimensions attendues sont alors celles de son encombrement exact (et non plus le diamètre)."""
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
        r = edge_report(tris)
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
                    args.tol_mm, not args.no_watertight, shape_mm)
    for name, ok, detail in results:
        print(f"[{'OK ' if ok else 'ÉCHEC'}] {name} : {detail}")
    sys.exit(0 if all(ok for _, ok, _ in results) else 1)


if __name__ == "__main__":
    main()
