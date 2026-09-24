#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
circle_dem_to_stl.py
---------------------
Découpe un DEM (GeoTIFF ou VRT, CRS géographique ou projeté) selon un
cercle et exporte un STL solide et étanche (surface du relief + paroi
verticale + fond plat), sans trou, prêt pour impression 3D.

Point d'entrée fin : toute la logique vit dans raster.py (GDAL),
mesh.py (numpy/scipy pur) et stl_io.py (écriture STL).

Licence : GNU GPLv3 (voir LICENSE). Ce programme est un logiciel libre :
vous pouvez le redistribuer et/ou le modifier selon les termes de la
Licence Publique Générale GNU publiée par la Free Software Foundation,
version 3. Distribué SANS AUCUNE GARANTIE, pas même implicite de
qualité marchande ou d'adéquation à un usage particulier.

Usage :
    python circle_dem_to_stl.py --input mont_fuji_dem1a.tif --output fuji_disc.stl \
        --cx 293526 --cy 3915399 --radius 9300 --diameter-mm 250 --vexag 1.0

    # ou directement en coordonnées géographiques, source projetée ou non :
    python circle_dem_to_stl.py --input mont_fuji_dem1a.tif --output fuji_disc.stl \
        --lat 35.3606 --lon 138.7274 --radius 9300 --diameter-mm 250
"""
import argparse
import sys

import numpy as np

import download
import mesh
import raster
import stl_io

log = raster.log


def progress(i, n, label):
    """Barre de progression légère, sans dépendance externe, pour les
    boucles Python pures (parois, écriture STL) qui ne bénéficient pas
    du vectorisme numpy utilisé pour le dessus/dessous du disque."""
    step = max(1, n // 200)
    if i % step == 0 or i == n - 1:
        pct = 100.0 * (i + 1) / n
        print(f"\r[DEM2STL] {label} : {pct:5.1f}%", end="", flush=True)
        if i == n - 1:
            print()


def ask_float(question, default):
    raw = input(f"[DEM2STL] {question} [{default}] : ").strip()
    if raw == "":
        return default
    try:
        return float(raw)
    except ValueError:
        print(f"[DEM2STL]   valeur invalide, on garde {default}")
        return default


def build_mesh(input_path, lat, lon, cx, cy, radius, diameter_mm, vexag, base_mm,
                print_spacing_mm, sea_level=None, min_elevation=-50.0, land_threshold=1.0):
    dst_srs_wkt, cx, cy = raster.resolve_target_crs_and_center(input_path, lat, lon, cx, cy)

    scale = (diameter_mm / 2.0) / radius
    pixel_size_m = print_spacing_mm / scale
    log(f"échelle : {scale:.6f} mm/m | résolution de ré-échantillonnage : {pixel_size_m:.3f} m/pixel")

    elev, xs, ys, nodata = raster.warp_window(input_path, dst_srs_wkt, cx, cy, radius, pixel_size_m)
    H, W = elev.shape
    log(f"grille ré-échantillonnée : {W} x {H} points")

    dist2 = (xs[None, :] - cx) ** 2 + (ys[:, None] - cy) ** 2
    valid = (elev != nodata) if nodata is not None else np.ones_like(elev, dtype=bool)
    plausible = valid & (elev > min_elevation)

    if sea_level is not None:
        # Le bruit de surface de mer (vagues/reflets classés par erreur
        # comme "sol" par le LiDAR) reste presque toujours très proche du
        # niveau de la mer, ET il touche souvent directement le vrai
        # littoral — la connectivité seule ne peut donc pas le séparer de
        # la vraie terre. On identifie la terre "sûre" par un critère
        # d'altitude strict, puis seule sa plus grande composante connexe
        # est retenue ; tout le reste est aplati avec la mer.
        land_seed = plausible & (elev > sea_level + land_threshold)
        is_land, n_components, sizes = mesh.largest_connected_component(land_seed)
        if n_components == 0:
            sys.exit("Erreur : aucune terre détectée au-dessus de --land-threshold.")
        if n_components > 1:
            n_other = int(land_seed.sum() - sizes[int(np.argmax(sizes))])
            log(f"{n_components} composantes au-dessus de --land-threshold, {n_other} pixel(s) "
                f"hors composante principale traité(s) comme artefacts")
            if len(sizes) >= 2:
                sorted_sizes = sorted(sizes, reverse=True)
                if sorted_sizes[1] > 0.3 * sorted_sizes[0]:
                    log(f"  avertissement : une composante secondaire représente "
                        f"{100 * sorted_sizes[1] / sorted_sizes[0]:.0f}% de la principale — "
                        f"possible deuxième île supprimée à tort, vérifie visuellement")

        n_replaced = int((~is_land).sum())
        elev = np.where(is_land, elev, sea_level)
        log(f"{n_replaced} pixels remplacés par le niveau de la mer ({sea_level} m)")
        raw_mask = dist2 <= radius ** 2
    else:
        raw_mask = (dist2 <= radius ** 2) & plausible

    if not raw_mask.any():
        sys.exit("Erreur : aucun pixel valide à l'intérieur du cercle demandé.")

    clean_mask, quad_ok = mesh.clean_mask(raw_mask)
    if not quad_ok.any():
        sys.exit("Erreur : cercle trop petit par rapport à la résolution de ré-échantillonnage "
                  "(ou vide après nettoyage morphologique).")
    log(f"{int(quad_ok.sum())} cellules incluses dans le disque")

    min_elev = sea_level if sea_level is not None else elev[clean_mask].min()

    tris, n_top, n_bot, n_wall = mesh.build_mesh_from_arrays(
        elev, quad_ok, xs - cx, ys - cy, scale, vexag, base_mm, min_elev, progress=progress)
    log(f"{n_top} triangles (dessus) + {n_bot} (dessous) + {n_wall} (paroi) = {len(tris)} triangles au total")
    return tris


def main():
    p = argparse.ArgumentParser(description="Découpe circulaire d'un DEM + export STL solide.")
    p.add_argument("--input", default=None,
                    help="GeoTIFF ou VRT source ; omis si --download-copernicus est utilisé")
    p.add_argument("--output", required=True)
    p.add_argument("--cx", type=float, default=None,
                    help="Centre X du cercle, dans le CRS du fichier source")
    p.add_argument("--cy", type=float, default=None,
                    help="Centre Y du cercle, dans le CRS du fichier source")
    p.add_argument("--lat", type=float, default=None,
                    help="Latitude du centre (WGS84) — alternative à --cx/--cy")
    p.add_argument("--lon", type=float, default=None,
                    help="Longitude du centre (WGS84) — alternative à --cx/--cy")
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
                    help="Utilisé avec --sea-level : élévation au-dessus de --sea-level à partir de "
                         "laquelle un pixel est considéré comme terre 'sûre', en m")
    p.add_argument("--download-copernicus", action="store_true",
                    help="Télécharge automatiquement les tuiles Copernicus DEM GLO-30 nécessaires "
                         "(requiert --lat/--lon) au lieu de fournir --input")
    p.add_argument("--copernicus-product", type=int, choices=[30, 90], default=30,
                    help="Résolution Copernicus DEM à télécharger : 30 (GLO-30) ou 90 (GLO-90)")
    p.add_argument("--cache-dir", default="./copernicus_cache",
                    help="Dossier de cache pour les tuiles Copernicus téléchargées")
    args = p.parse_args()

    have_latlon = args.lat is not None and args.lon is not None
    have_xy = args.cx is not None and args.cy is not None
    if have_latlon and have_xy:
        sys.exit("Erreur : fournis --cx/--cy OU --lat/--lon, pas les deux.")
    if not have_latlon and not have_xy:
        sys.exit("Erreur : fournis --cx/--cy (CRS du fichier source) ou --lat/--lon (WGS84).")

    if args.download_copernicus:
        if not have_latlon:
            sys.exit("Erreur : --download-copernicus requiert --lat/--lon (coordonnées géographiques "
                      "nécessaires pour déterminer les tuiles à récupérer).")
        if args.input is not None:
            sys.exit("Erreur : ne fournis pas --input en même temps que --download-copernicus.")
        tile_paths = download.download_copernicus(
            args.lat, args.lon, args.radius, args.cache_dir,
            product_m=args.copernicus_product, log=raster.log)
        vrt_path = f"{args.cache_dir.rstrip('/')}/_mosaic.vrt"
        from osgeo import gdal
        gdal.BuildVRT(vrt_path, tile_paths)
        args.input = vrt_path
        raster.log(f"Mosaïque de {len(tile_paths)} tuile(s) prête : {vrt_path}")
    elif args.input is None:
        sys.exit("Erreur : fournis --input, ou utilise --download-copernicus.")

    if args.diameter_mm is None:
        args.diameter_mm = ask_float("Diamètre final imprimé (mm)", 250.0)
    if args.base_mm is None:
        args.base_mm = ask_float("Épaisseur du socle plat (mm)", 3.0)
    if args.vexag is None:
        args.vexag = ask_float("Exagération verticale", 1.0)

    tris = build_mesh(args.input, args.lat, args.lon, args.cx, args.cy, args.radius,
                       args.diameter_mm, args.vexag, args.base_mm, args.print_spacing_mm,
                       sea_level=args.sea_level, min_elevation=args.min_elevation,
                       land_threshold=args.land_threshold)
    stl_io.write_stl_binary(tris, args.output, progress=progress)
    log(f"Terminé : {args.output}")


if __name__ == "__main__":
    main()
