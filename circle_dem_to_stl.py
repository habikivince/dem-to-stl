#!/usr/bin/env python3
"""
circle_dem_to_stl.py
---------------------
Découpe un DEM (GeoTIFF, CRS projeté en mètres) selon un cercle et exporte
un STL solide et étanche (surface du relief + paroi verticale + fond plat),
sans trou, prêt pour impression 3D.

Licence : GNU GPLv3 (voir LICENSE). Ce programme est un logiciel libre :
vous pouvez le redistribuer et/ou le modifier selon les termes de la
Licence Publique Générale GNU publiée par la Free Software Foundation,
version 3. Distribué SANS AUCUNE GARANTIE, pas même implicite de
qualité marchande ou d'adéquation à un usage particulier.

Usage :
    python circle_dem_to_stl.py --input mont_fuji_dem1a.tif --output fuji_disc.stl \
        --cx 293526 --cy 3915399 --radius 9300 --diameter-mm 250 --vexag 1.0
"""
import argparse
import struct
import sys

import numpy as np
from osgeo import gdal
from scipy import ndimage

gdal.UseExceptions()


def log(msg):
    print(f"[DEM2STL] {msg}", flush=True)


def progress(i, n, label):
    """Barre de progression légère, sans dépendance externe, pour les
    boucles Python pures (parois, écriture STL) qui ne bénéficient pas
    du vectorisme numpy utilisé pour le dessus/dessous du disque."""
    step = max(1, n // 200)  # ~200 mises à jour max, pour ne pas ralentir par excès d'I/O
    if i % step == 0 or i == n - 1:
        pct = 100.0 * (i + 1) / n
        print(f"\r[DEM2STL] {label} : {pct:5.1f}%", end="", flush=True)
        if i == n - 1:
            print()


def ask_float(question, default):
    """Demande une valeur numérique à l'utilisateur, avec une valeur par
    défaut si Entrée est pressée sans rien saisir."""
    raw = input(f"[DEM2STL] {question} [{default}] : ").strip()
    if raw == "":
        return default
    try:
        return float(raw)
    except ValueError:
        print(f"[DEM2STL]   valeur invalide, on garde {default}")
        return default


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


def build_mesh_from_arrays(elev, quad_ok, xs, ys, scale, vexag, base_mm, min_elev):
    """Construit le maillage (triangles dessus/dessous/parois) à partir d'un
    tableau d'altitudes déjà découpé et d'un masque de cellules (quad_ok)
    déjà nettoyé (ouverture morphologique + anti-pincement diagonal déjà
    appliqués par l'appelant). Séparée de build_mesh() pour être testable
    sans GDAL. quad_ok doit avoir la forme (H-1, W-1) par rapport à elev."""
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
        progress(i, H1, "Parois (bord vertical)")
        for j in range(W1 + 1):
            left = quad_ok[i, j - 1] if j > 0 else False
            right = quad_ok[i, j] if j < W1 else False
            if left != right:
                p1 = (Xmm[j], Ymm[i], z_top_mm[i, j])
                p2 = (Xmm[j], Ymm[i + 1], z_top_mm[i + 1, j])
                add_wall(p1, p2, outward=bool(left))
    for j in range(W1):
        progress(j, W1, "Parois (bord horizontal)")
        for i in range(H1 + 1):
            up = quad_ok[i - 1, j] if i > 0 else False
            down = quad_ok[i, j] if i < H1 else False
            if up != down:
                p1 = (Xmm[j], Ymm[i], z_top_mm[i, j])
                p2 = (Xmm[j + 1], Ymm[i], z_top_mm[i, j + 1])
                add_wall(p1, p2, outward=bool(down))

    wall_tris = np.array(wall_tris, dtype=np.float64)
    all_tris = np.concatenate([top_tris, bot_tris, wall_tris])
    log(f"{len(top_tris)} triangles (dessus) + {len(bot_tris)} (dessous) + {len(wall_tris)} (paroi) "
        f"= {len(all_tris)} triangles au total")
    return all_tris


def build_mesh(input_tif, cx, cy, radius, diameter_mm, vexag, base_mm, print_spacing_mm,
                sea_level=None, min_elevation=-50.0, land_threshold=1.0):
    ds = gdal.Open(input_tif)
    nodata = ds.GetRasterBand(1).GetNoDataValue()

    scale = (diameter_mm / 2.0) / radius  # mm par mètre réel
    pixel_size_m = print_spacing_mm / scale
    log(f"échelle : {scale:.6f} mm/m | résolution de ré-échantillonnage : {pixel_size_m:.3f} m/pixel")

    xmin, xmax = cx - radius, cx + radius
    ymin, ymax = cy - radius, cy + radius
    warped = gdal.Warp(
        "", input_tif, format="MEM",
        outputBounds=(xmin, ymin, xmax, ymax),
        xRes=pixel_size_m, yRes=pixel_size_m,
        resampleAlg="average",
        dstNodata=nodata,
        callback=gdal.TermProgress_nocb,
    )
    band = warped.GetRasterBand(1)
    gdal.FillNodata(band, None, maxSearchDist=50, smoothingIterations=0,
                     callback=gdal.TermProgress_nocb)
    elev = band.ReadAsArray().astype(np.float64)
    wgt = warped.GetGeoTransform()
    H, W = elev.shape
    log(f"grille ré-échantillonnée : {W} x {H} points")

    xs = wgt[0] + (np.arange(W) + 0.5) * wgt[1]
    ys = wgt[3] + (np.arange(H) + 0.5) * wgt[5]

    dist2 = (xs[None, :] - cx) ** 2 + (ys[:, None] - cy) ** 2
    valid = (elev != nodata) if nodata is not None else np.ones_like(elev, dtype=bool)
    plausible = valid & (elev > min_elevation)

    if sea_level is not None:
        # Le bruit de surface de mer (vagues/reflets classés par erreur
        # comme "sol" par le LiDAR) reste presque toujours très proche du
        # niveau de la mer, ET il touche souvent directement le vrai
        # littoral — la connectivité seule ne peut donc pas le séparer de
        # la vraie terre. On identifie la terre "sûre" par un critère
        # d'altitude strict (nettement au-dessus du niveau de la mer),
        # puis seule sa plus grande composante connexe est retenue comme
        # terre ; tout le reste (y compris ce bruit proche de zéro) est
        # aplati avec la mer, même si sa valeur passait le seuil
        # --min-elevation (bien plus permissif, pensé pour les gouffres
        # d'interpolation, pas le bruit de surface d'eau).
        land_seed = plausible & (elev > sea_level + land_threshold)
        labeled, n_components = ndimage.label(land_seed, structure=np.ones((3, 3)))
        if n_components == 0:
            sys.exit("Erreur : aucune terre détectée au-dessus de --land-threshold.")
        sizes = ndimage.sum(land_seed, labeled, index=range(1, n_components + 1))
        largest_label = 1 + int(np.argmax(sizes))
        is_land = labeled == largest_label
        if n_components > 1:
            n_other = int(land_seed.sum() - sizes[largest_label - 1])
            log(f"{n_components} composantes au-dessus de --land-threshold, {n_other} pixel(s) "
                f"hors composante principale traité(s) comme artefacts")

        n_replaced = int((~is_land).sum())
        elev = np.where(is_land, elev, sea_level)
        log(f"{n_replaced} pixels remplacés par le niveau de la mer ({sea_level} m)")
        mask = dist2 <= radius ** 2
    else:
        mask = (dist2 <= radius ** 2) & plausible

    if not mask.any():
        sys.exit("Erreur : aucun pixel valide à l'intérieur du cercle demandé.")

    # Ouverture morphologique (érosion puis dilatation) pour éliminer les
    # languettes d'un seul pixel de large que peut créer le masquage
    # circulaire naïf (parois trop fines pour être imprimées en FDM).
    mask = dilate3x3(erode3x3(mask))
    if not mask.any():
        sys.exit("Erreur : le masque est vide après nettoyage morphologique (cercle trop petit ?).")

    min_elev = sea_level if sea_level is not None else elev[mask].min()

    quad_ok = mask[:-1, :-1] & mask[1:, :-1] & mask[:-1, 1:] & mask[1:, 1:]
    quad_ok = fix_diagonal_pinches(quad_ok)
    if not quad_ok.any():
        sys.exit("Erreur : cercle trop petit par rapport à la résolution de ré-échantillonnage.")
    log(f"{int(quad_ok.sum())} cellules incluses dans le disque")

    # xs/ys sont recentrées sur (cx, cy) AVANT mise à l'échelle, pour que le
    # disque produit soit centré sur (0, 0) en mm plutôt que sur les
    # coordonnées absolues (souvent énormes) du CRS source.
    return build_mesh_from_arrays(elev, quad_ok, xs - cx, ys - cy, scale, vexag, base_mm, min_elev)


def write_stl_binary(tris, path):
    n = len(tris)
    with open(path, "wb") as f:
        f.write(b"\x00" * 80)
        f.write(struct.pack("<I", n))
        for idx, t in enumerate(tris):
            progress(idx, n, "Écriture STL")
            v0, v1, v2 = t[0], t[1], t[2]
            nvec = np.cross(v1 - v0, v2 - v0)
            norm = np.linalg.norm(nvec)
            if norm > 0:
                nvec = nvec / norm
            f.write(struct.pack("<3f", *nvec))
            for v in (v0, v1, v2):
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))


def main():
    p = argparse.ArgumentParser(description="Découpe circulaire d'un DEM + export STL solide.")
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--cx", type=float, default=None,
                    help="Centre X du cercle (demandé interactivement si omis)")
    p.add_argument("--cy", type=float, default=None,
                    help="Centre Y du cercle (demandé interactivement si omis)")
    p.add_argument("--radius", type=float, required=True, help="Rayon du cercle, en mètres")
    p.add_argument("--diameter-mm", type=float, default=None,
                    help="Diamètre final imprimé, en mm (demandé interactivement si omis)")
    p.add_argument("--vexag", type=float, default=None,
                    help="Exagération verticale (demandée interactivement si omise)")
    p.add_argument("--base-mm", type=float, default=None,
                    help="Épaisseur du socle plat, en mm (demandé interactivement si omis)")
    p.add_argument("--print-spacing-mm", type=float, default=0.2,
                    help="Résolution cible du maillage, en mm (~0.2 recommandé pour buse 0.4mm)")
    p.add_argument("--sea-level", type=float, default=None,
                    help="Mode île : remplace les pixels invalides/aberrants par un plat à cette altitude "
                         "(m) au lieu de les exclure — le cercle est alors toujours plein.")
    p.add_argument("--min-elevation", type=float, default=-50.0,
                    help="Utilisé avec --sea-level : seuil sous lequel une valeur (hors nodata déclaré) "
                         "est traitée comme un artefact et remplacée, en m")
    p.add_argument("--land-threshold", type=float, default=1.0,
                    help="Utilisé avec --sea-level : élévation au-dessus du niveau de la mer à partir de "
                         "laquelle un pixel est considéré comme terre 'sûre' pour identifier la composante "
                         "principale — le bruit de surface d'eau (vagues classées comme sol par le LiDAR) "
                         "reste en général sous ce seuil, en m au-dessus de --sea-level")
    args = p.parse_args()

    if args.cx is None:
        args.cx = ask_float("Centre X du cercle", 0.0)
    if args.cy is None:
        args.cy = ask_float("Centre Y du cercle", 0.0)
    if args.diameter_mm is None:
        args.diameter_mm = ask_float("Diamètre final imprimé (mm)", 250.0)
    if args.base_mm is None:
        args.base_mm = ask_float("Épaisseur du socle plat (mm)", 3.0)
    if args.vexag is None:
        args.vexag = ask_float("Exagération verticale", 1.0)

    tris = build_mesh(args.input, args.cx, args.cy, args.radius, args.diameter_mm,
                       args.vexag, args.base_mm, args.print_spacing_mm,
                       sea_level=args.sea_level, min_elevation=args.min_elevation,
                       land_threshold=args.land_threshold)
    write_stl_binary(tris, args.output)
    log(f"Terminé : {args.output}")


if __name__ == "__main__":
    main()
