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
import math
import os
import re
import sys
import unicodedata

import numpy as np

import mesh
import mesh_export
import geocode
import polygon
import raster
import sources
import stl_io

log = raster.log


RAM_BYTES_PER_TRIANGLE = 220    # mesuré : 1,05 Go de pic pour 4,9 M de triangles
RAM_BASE_BYTES = 300_000_000
RAM_SAFETY = 0.8
# Taille des fichiers par rapport au STL, mesurée sur 4,9 M de triangles (STL 245 Mo, 3MF 61 Mo, OBJ 192 Mo) ;
# varie avec le relief (nombre de chiffres des coordonnées)
FORMAT_SIZE_RATIO = {"stl": 1.0, "3mf": 0.25, "obj": 0.78}
MIN_GRID_CELLS = 8


def fail(msg):
    sys.exit(f"Erreur : {msg}")


def _pos(v):
    return v is not None and math.isfinite(v) and v > 0


def available_ram_bytes():
    """RAM disponible (Linux, /proc/meminfo) ou None si inconnue."""
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except (OSError, ValueError):
        pass
    return None


def estimate_resources(diameter_mm, spacing_mm, polygon_mode=False):
    """Estimation (grille n x n, triangles, taille STL, RAM) d'après le pas d'impression.
    Cercle : pi*n^2 triangles ; contour : borne haute 4*n^2. Constantes mesurées sur un cas
    de 4,9 M de triangles ; au-delà de ce cas, c'est une extrapolation."""
    n = diameter_mm / spacing_mm
    triangles = (4.0 if polygon_mode else math.pi) * n * n
    return {"grid": n, "triangles": triangles, "stl_bytes": 84 + 50 * triangles,
            "ram_bytes": RAM_BASE_BYTES + RAM_BYTES_PER_TRIANGLE * triangles}


def check_resources(diameter_mm, spacing_mm, polygon_mode, force, log=print, ram_available=available_ram_bytes,
                    formats=("stl",)):
    est = estimate_resources(diameter_mm, spacing_mm, polygon_mode)
    avail = ram_available()
    sizes = ", ".join(f"{f.upper()} ~{est['stl_bytes'] * FORMAT_SIZE_RATIO[f] / 1e6:.0f} Mo" for f in formats)
    log(f"Estimation : grille ~{est['grid']:.0f}², ~{est['triangles'] / 1e6:.1f} M triangles, "
        f"{sizes}, RAM ~{est['ram_bytes'] / 1e9:.1f} Go"
        + (f" (disponible : {avail / 1e9:.1f} Go)" if avail else ""))
    if avail and est["ram_bytes"] > RAM_SAFETY * avail and not force:
        fail(f"cette taille de grille demanderait environ {est['ram_bytes'] / 1e9:.1f} Go de RAM pour "
             f"{avail / 1e9:.1f} Go disponibles. Augmente --print-spacing-mm ou réduis --diameter-mm, "
             f"ou ajoute --force pour tenter quand même.")
    return est


DEFAULT_OUTPUT_DIR = "output"
OUTPUT_FORMATS = mesh_export.FORMATS          # ("stl", "3mf", "obj")

# Mentions de source (voir « Sources de données, licences et mentions » du README)
SOURCE_ATTRIBUTION = {
    "swisstopo": "© swisstopo",
    "ign": "Source : IGN (Licence Ouverte Etalab 2.0)",
    "kartverket": "© Kartverket (CC BY 4.0)",
    "usgs": "U.S. Geological Survey, 3D Elevation Program",
    "gsi": "出典:国土地理院ウェブサイト 地理院タイル(標高タイル)を加工して作成",
    "copernicus": ("produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and "
                   "Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved"),
    "osm": "© OpenStreetMap contributors (ODbL)",
}


def slugify(text):
    text = unicodedata.normalize("NFKD", str(text)).encode("ascii", "ignore").decode("ascii")
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-") or "relief"


def default_stem(args):
    if args.place:
        base = slugify(args.place) + ("-contour" if args.clip_to_place else "")
    elif args.polygon:
        base = slugify(os.path.splitext(os.path.basename(args.polygon))[0])
    elif args.lat is not None and args.lon is not None:
        base = f"lat{args.lat:.3f}_lon{args.lon:.3f}".replace("-", "m")
    else:
        base = f"x{args.cx:.0f}_y{args.cy:.0f}".replace("-", "m")
    name = f"{base}_{args.diameter_mm:g}mm"
    return name if args.vexag == 1 else f"{name}_x{args.vexag:g}"


def parse_formats(args):
    if args.format:
        items = [t.strip().lower() for t in args.format.split(",") if t.strip()]
        out = []
        for t in items:
            for f in (OUTPUT_FORMATS if t == "all" else (t,)):
                if f not in OUTPUT_FORMATS:
                    fail(f"format inconnu : {t} (attendu : {', '.join(OUTPUT_FORMATS)} ou all)")
                if f not in out:
                    out.append(f)
        return out
    if args.output:
        ext = os.path.splitext(args.output)[1].lower().lstrip(".")
        if ext in OUTPUT_FORMATS:
            return [ext]
    return ["stl"]


def resolve_outputs(args):
    """{format: chemin}. --output sans dossier va dans --output-dir (défaut « output ») ; avec un
    dossier, il est respecté tel quel. Sans --output, le nom est déduit du lieu/contour et du
    diamètre, et ne remplace jamais un fichier existant (suffixe _2, _3...)."""
    formats = parse_formats(args)
    if args.output:
        known_ext = os.path.splitext(args.output)[1].lower().lstrip(".") in OUTPUT_FORMATS
        target = args.output if os.path.dirname(args.output) else os.path.join(args.output_dir, args.output)
        if not known_ext and not args.format:
            return {"stl": target}            # comportement historique : nom donné repris tel quel (STL)
        stem = os.path.splitext(target)[0] if known_ext else target
        return {f: f"{stem}.{f}" for f in formats}
    base = os.path.join(args.output_dir, default_stem(args))
    n, stem = 1, base
    while any(os.path.exists(f"{stem}.{f}") for f in formats):
        n += 1
        stem = f"{base}_{n}"
    return {f: f"{stem}.{f}" for f in formats}


def prepare_outputs(paths):
    """Crée les dossiers parents et vérifie qu'on peut écrire, AVANT le calcul."""
    for path in paths.values():
        if os.path.isdir(path):
            fail(f"la sortie est un dossier ({path}) : donne un nom de fichier.")
        parent = os.path.dirname(os.path.abspath(path))
        probe = path + ".part"
        try:
            os.makedirs(parent, exist_ok=True)
            open(probe, "wb").close()
            os.remove(probe)
        except OSError as e:
            fail(f"impossible d'écrire dans {parent} ({e.strerror or e})")


def write_outputs(tris, paths, title, description, attribution, progress=None):
    """Écrit chaque format demandé à partir du même maillage."""
    if "stl" in paths:
        stl_io.write_stl_binary(tris, paths["stl"], progress=progress)
        log(f"écrit : {os.path.abspath(paths['stl'])}")
    if "3mf" in paths or "obj" in paths:
        vertices, faces = mesh_export.index_mesh(tris)
        log(f"maillage indexé : {len(vertices)} sommets, {len(faces)} triangles")
        if "3mf" in paths:
            mesh_export.write_3mf(paths["3mf"], vertices, faces, title=title, description=description,
                                  copyright_text=attribution)
            log(f"écrit : {os.path.abspath(paths['3mf'])}")
        if "obj" in paths:
            mesh_export.write_obj(paths["obj"], vertices, faces, title=title,
                                  comment="; ".join(x for x in (description, attribution) if x))
            log(f"écrit : {os.path.abspath(paths['obj'])}")


def validate_geometry(args, from_contour):
    if not from_contour:
        if not _pos(args.radius):
            fail("--radius doit être un nombre strictement positif (en mètres).")
    if args.lat is not None and not (math.isfinite(args.lat) and -90 <= args.lat <= 90):
        fail(f"latitude hors limites : {args.lat} (attendu entre -90 et 90).")
    if args.lon is not None and not (math.isfinite(args.lon) and -180 <= args.lon <= 180):
        fail(f"longitude hors limites : {args.lon} (attendu entre -180 et 180).")
    for name, v in (("--cx", args.cx), ("--cy", args.cy)):
        if v is not None and not math.isfinite(v):
            fail(f"{name} doit être un nombre fini.")
    if not _pos(args.print_spacing_mm):
        fail("--print-spacing-mm doit être un nombre strictement positif (en mm).")
    if args.place_pick < 1:
        fail("--place-pick commence à 1.")


def validate_print(args):
    if not _pos(args.diameter_mm):
        fail("--diameter-mm doit être un nombre strictement positif (en mm).")
    if not _pos(args.base_mm):
        fail("--base-mm doit être strictement positif : un socle d'épaisseur nulle ou négative "
             "donne un maillage dégénéré (arêtes non-manifold constatées pour 0).")
    if not math.isfinite(args.vexag) or args.vexag < 0:
        fail("--vexag doit être positif : un relief inversé passe sous le socle.")
    if args.vexag == 0:
        log("ATTENTION : --vexag 0 donne un disque plat, sans relief.")
    if args.diameter_mm / args.print_spacing_mm < MIN_GRID_CELLS:
        fail(f"--print-spacing-mm ({args.print_spacing_mm:g}) est trop grand pour --diameter-mm "
             f"({args.diameter_mm:g}) : il faut au moins {MIN_GRID_CELLS} cellules sur le diamètre.")


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
                print_spacing_mm, sea_level=None, min_elevation=-50.0, land_threshold=1.0,
                smooth_wall=False, clip=None):
    dst_srs_wkt, cx, cy = raster.resolve_target_crs_and_center(input_path, lat, lon, cx, cy)

    clip_rings = None
    if clip is not None:
        # mode polygone : fenêtre carrée englobant le contour, « rayon » = demi-plus grande dimension
        clip_geom = clip.to_srs(dst_srs_wkt)
        x_min, x_max, y_min, y_max = clip_geom.GetEnvelope()
        cx, cy = (x_min + x_max) / 2.0, (y_min + y_max) / 2.0
        radius = max(x_max - x_min, y_max - y_min) / 2.0
        clip_rings = polygon.rings_of(clip_geom)
        log(f"contour : emprise {x_max - x_min:.0f} x {y_max - y_min:.0f} m ; plus grande dimension "
            f"imprimée = {diameter_mm:g} mm")
    inside_fn = (lambda xs_, ys_: polygon.polygon_mask(clip_rings, xs_, ys_)) if clip is not None else None

    scale = (diameter_mm / 2.0) / radius
    pixel_size_m = print_spacing_mm / scale
    log(f"échelle : {scale:.6f} mm/m | résolution de ré-échantillonnage : {pixel_size_m:.3f} m/pixel")

    elev, xs, ys, nodata = raster.warp_window(input_path, dst_srs_wkt, cx, cy, radius, pixel_size_m,
                                              inside_fn=inside_fn)
    H, W = elev.shape
    log(f"grille ré-échantillonnée : {W} x {H} points")

    inside = (inside_fn(xs, ys) if inside_fn is not None
              else (xs[None, :] - cx) ** 2 + (ys[:, None] - cy) ** 2 <= radius ** 2)
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
        raw_mask = inside
    else:
        raw_mask = inside & plausible

    if not raw_mask.any():
        sys.exit("Erreur : aucun pixel valide à l'intérieur de la zone demandée (cercle ou contour).")

    clean_mask, quad_ok = mesh.clean_mask(raw_mask)
    if not quad_ok.any():
        sys.exit("Erreur : zone trop petite par rapport à la résolution de ré-échantillonnage "
                  "(ou vide après nettoyage morphologique).")
    log(f"{int(quad_ok.sum())} cellules incluses dans {'le contour' if clip is not None else 'le disque'}")

    min_elev = sea_level if sea_level is not None else elev[clean_mask].min()

    tris, n_top, n_bot, n_wall = mesh.build_mesh_from_arrays(
        elev, quad_ok, xs - cx, ys - cy, scale, vexag, base_mm, min_elev, progress=progress,
        smooth_wall=smooth_wall, circle_radius_mm=(radius * scale if smooth_wall else None))
    log(f"{n_top} triangles (dessus) + {n_bot} (dessous) + {n_wall} (paroi) = {len(tris)} triangles au total")
    return tris


def main():
    p = argparse.ArgumentParser(description="Découpe circulaire d'un DEM + export STL solide.")
    p.add_argument("--input", default=None,
                    help="GeoTIFF ou VRT source ; omis si --source / --download-copernicus est utilisé")
    p.add_argument("--output", default=None,
                    help="Fichier de sortie. Sans dossier, il est placé dans --output-dir ; avec un dossier, il est "
                         "respecté. L'extension (.stl, .3mf, .obj) choisit le format. Facultatif : le nom est alors "
                         "déduit du lieu et du diamètre sans écraser de fichier existant.")
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                    help="Dossier des fichiers générés (défaut : %(default)s, créé au besoin)")
    p.add_argument("--format", default=None,
                    help="Format(s) de sortie : stl (défaut), 3mf, obj, une liste séparée par des virgules, ou all. "
                         "Même géométrie, mêmes dimensions (millimètres) dans les trois.")
    p.add_argument("--force", action="store_true",
                    help="Ignore le garde-fou de taille (estimation de RAM supérieure à 80 %% de la RAM disponible)")
    p.add_argument("--debug", action="store_true", help="Affiche la trace complète des erreurs")
    p.add_argument("--dry-run", action="store_true",
                    help="Valide les paramètres, affiche l'estimation (grille, triangles, taille des fichiers, RAM) "
                         "et les fichiers prévus, puis s'arrête sans rien télécharger ni écrire")
    p.add_argument("--cx", type=float, default=None,
                    help="Centre X du cercle, dans le CRS du fichier source")
    p.add_argument("--cy", type=float, default=None,
                    help="Centre Y du cercle, dans le CRS du fichier source")
    p.add_argument("--lat", type=float, default=None,
                    help="Latitude du centre (WGS84) — alternative à --cx/--cy")
    p.add_argument("--lon", type=float, default=None,
                    help="Longitude du centre (WGS84) — alternative à --cx/--cy")
    p.add_argument("--full-res", action="store_true",
                    help="Avec --source : télécharge la résolution native maximale de la source (swisstopo "
                         "0,5 m, USGS 1 m, Kartverket 1 m, IGN 0,5 m...) au lieu de celle suffisante pour le "
                         "maillage demandé. Le STL obtenu est quasi identique (le maillage reste limité par "
                         "--print-spacing-mm) mais les téléchargements sont beaucoup plus lourds.")
    p.add_argument("--fill-with-copernicus", action="store_true",
                    help="Comble les parties du disque sans donnée (frontière, couverture incomplète) avec "
                         "Copernicus GLO-30, recalé en altitude sur la source principale d'après leur zone "
                         "commune. Copernicus n'est téléchargé que s'il y a des trous. Requiert --lat/--lon "
                         "(ou --place) ; fonctionne avec --source ou --input.")
    p.add_argument("--place", default=None,
                    help="Nom de lieu à la place de --lat/--lon (ex. \"Mont Fuji\", \"Mount Everest\"), "
                         "résolu via Nominatim/OpenStreetMap. Sans --input ni --source, la source est "
                         "choisie automatiquement (--source auto).")
    p.add_argument("--clip-to-place", action="store_true",
                    help="Avec --place : découpe selon le contour OpenStreetMap du lieu (île, commune, pays) "
                         "au lieu d'un cercle ; --diameter-mm devient la plus grande dimension imprimée. "
                         "Incompatible avec --radius, --polygon, --smooth-wall.")
    p.add_argument("--place-pick", type=int, default=1,
                    help="Avec --place : numéro du résultat à retenir si le premier n'est pas le bon (1-5)")
    p.add_argument("--polygon", default=None,
                    help="Contour (GeoJSON, GPKG, shapefile...) à la place du cercle : seule la zone "
                         "intérieure est modélisée. --diameter-mm devient la plus grande dimension imprimée. "
                         "Incompatible avec --radius, --lat/--lon, --cx/--cy, --place, --smooth-wall.")
    p.add_argument("--polygon-where", default=None,
                    help="Filtre d'attributs OGR pour choisir des entités du contour, ex. \"NAME='Corse'\"")
    p.add_argument("--radius", type=float, default=None,
                    help="Rayon du cercle, en mètres (requis sans --polygon)")
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
    p.add_argument("--smooth-wall", action="store_true",
                    help="Bord parfaitement circulaire (contour en escalier projeté radialement "
                         "vers le vrai cercle) au lieu de suivre la grille de ré-échantillonnage. "
                         "Requiert un masque à un seul contour simple.")
    p.add_argument("--download-copernicus", action="store_true",
                    help="Équivalent de --source copernicus (conservé pour compatibilité)")
    p.add_argument("--source", choices=["auto", "copernicus"] + sorted(sources.SOURCES), default=None,
                    help="Télécharge automatiquement le DEM (requiert --lat/--lon) au lieu de fournir "
                         "--input : une source nationale (swisstopo, kartverket, usgs, gsi), copernicus "
                         "(mondial), ou auto (source nationale couvrant le point, sinon Copernicus)")
    p.add_argument("--copernicus-product", type=int, choices=[30, 90], default=30,
                    help="Résolution Copernicus DEM à télécharger : 30 (GLO-30) ou 90 (GLO-90)")
    p.add_argument("--cache-dir", default="./copernicus_cache",
                    help="Dossier de cache des données téléchargées (Copernicus et sources nationales)")
    args = p.parse_args()

    clip = None
    if args.polygon:
        used = [n for n, v in (("--radius", args.radius), ("--lat/--lon", args.lat if args.lat is not None else args.lon),
                               ("--cx/--cy", args.cx if args.cx is not None else args.cy),
                               ("--place", args.place), ("--smooth-wall", args.smooth_wall or None)) if v is not None]
        if used:
            sys.exit("Erreur : --polygon est incompatible avec " + ", ".join(used)
                     + " (le centre et la taille viennent du contour).")
        try:
            clip = polygon.load_polygon(args.polygon, args.polygon_where, log=raster.log)
        except polygon.PolygonError as e:
            sys.exit(f"Erreur : {e}")
        args.lat, args.lon = clip.center_latlon()
        args.radius = clip.approx_radius_m() * 1.02   # sert seulement à choisir/télécharger les données
        raster.log(f"Contour : centre lat={args.lat:.5f} lon={args.lon:.5f}")
        if args.input is None and args.source is None and not args.download_copernicus:
            args.source = "auto"
    elif args.clip_to_place:
        bad = [n for n, v in (("--radius", args.radius), ("--smooth-wall", args.smooth_wall or None)) if v is not None]
        if not args.place:
            sys.exit("Erreur : --clip-to-place requiert --place.")
        if bad:
            sys.exit("Erreur : --clip-to-place est incompatible avec " + ", ".join(bad)
                     + " (la taille vient du contour).")
    elif args.radius is None:
        sys.exit("Erreur : --radius est requis (ou --polygon, ou --place avec --clip-to-place).")

    if args.place:
        if any(v is not None for v in (args.lat, args.lon, args.cx, args.cy)):
            sys.exit("Erreur : --place est incompatible avec --lat/--lon et --cx/--cy.")
        try:
            if args.clip_to_place:
                clip = geocode.resolve_polygon(args.place, args.place_pick, args.cache_dir, log=raster.log)
                args.lat, args.lon = clip.center_latlon()
                args.radius = clip.approx_radius_m() * 1.02   # sert seulement à choisir/télécharger les données
            else:
                args.lat, args.lon = geocode.resolve(args.place, args.place_pick, args.cache_dir, log=raster.log)
        except geocode.GeocodeError as e:
            sys.exit(f"Erreur : {e}")
        if args.input is None and args.source is None and not args.download_copernicus:
            args.source = "auto"

    validate_geometry(args, from_contour=clip is not None)
    have_latlon = args.lat is not None and args.lon is not None
    have_xy = args.cx is not None and args.cy is not None
    if have_latlon and have_xy:
        sys.exit("Erreur : fournis --cx/--cy OU --lat/--lon, pas les deux.")
    if not have_latlon and not have_xy:
        sys.exit("Erreur : fournis --cx/--cy (CRS du fichier source), --lat/--lon (WGS84) ou --place.")

    source = args.source or ("copernicus" if args.download_copernicus else None)
    if args.source and args.download_copernicus and args.source != "copernicus":
        sys.exit("Erreur : --download-copernicus est incompatible avec --source " + args.source + ".")
    if source:
        if not have_latlon:
            sys.exit("Erreur : --source / --download-copernicus requiert --lat/--lon (coordonnées "
                      "géographiques nécessaires pour déterminer les données à récupérer).")
        if args.input is not None:
            sys.exit("Erreur : ne fournis pas --input en même temps que --source / --download-copernicus.")
    elif args.input is None:
        sys.exit("Erreur : fournis --input, ou utilise --source (ou --download-copernicus).")

    if args.fill_with_copernicus and not have_latlon:
        sys.exit("Erreur : --fill-with-copernicus requiert --lat/--lon (ou --place).")

    # Paramètres d'impression demandés AVANT le téléchargement : la résolution de
    # travail (pixel_size_m) en découle et sert à choisir la résolution à récupérer.
    if args.diameter_mm is None:
        args.diameter_mm = ask_float("Diamètre final imprimé (mm)", 250.0)
    if args.base_mm is None:
        args.base_mm = ask_float("Épaisseur du socle plat (mm)", 3.0)
    if args.vexag is None:
        args.vexag = ask_float("Exagération verticale", 1.0)

    validate_print(args)
    paths = resolve_outputs(args)
    check_resources(args.diameter_mm, args.print_spacing_mm, clip is not None, args.force, log=raster.log,
                    formats=tuple(paths))
    if args.dry_run:
        log("--dry-run : aucun téléchargement, calcul ni écriture. Fichiers prévus : "
            + ", ".join(os.path.abspath(p) for p in paths.values()))
        return
    prepare_outputs(paths)

    used = None
    if source:
        pixel_size_m = args.print_spacing_mm / ((args.diameter_mm / 2.0) / args.radius)
        if args.full_res:
            raster.log(f"--full-res : résolution native maximale demandée (le maillage reste à "
                       f"{pixel_size_m:.2f} m/pixel d'après --print-spacing-mm)")
            pixel_size_m = sources.FULL_RES_PIXEL_M
        try:
            used, args.input = sources.acquire(
                source, args.lat, args.lon, args.radius, args.cache_dir, pixel_size_m,
                copernicus_product=args.copernicus_product, log=raster.log)
        except (sources.SourceUnavailable, RuntimeError) as e:
            if args.debug:
                raise
            fail(f"source de données inutilisable : {e}")
        raster.log(f"source utilisée : {used}")

    if args.fill_with_copernicus:
        if used == "copernicus":
            raster.log("--fill-with-copernicus ignoré : la source principale est déjà Copernicus")
        else:
            args.input = sources.with_copernicus_fill(
                args.input, args.lat, args.lon, args.radius, args.cache_dir,
                product_m=args.copernicus_product, log=raster.log)

    try:
        raster.check_readable(args.input)
        tris = build_mesh(args.input, args.lat, args.lon, args.cx, args.cy, args.radius,
                          args.diameter_mm, args.vexag, args.base_mm, args.print_spacing_mm,
                          sea_level=args.sea_level, min_elevation=args.min_elevation,
                          land_threshold=args.land_threshold, smooth_wall=args.smooth_wall, clip=clip)
    except raster.InputError as e:
        fail(str(e))
    except RuntimeError as e:
        if args.debug:
            raise
        fail(f"traitement GDAL impossible : {str(e).strip().splitlines()[0] if str(e).strip() else e} "
             f"(--debug pour la trace complète)")
    stem = os.path.splitext(os.path.basename(next(iter(paths.values()))))[0]
    srcs = [used] if used else []
    if args.fill_with_copernicus:
        srcs.append("copernicus")
    if args.clip_to_place:
        srcs.append("osm")
    attribution = "; ".join(SOURCE_ATTRIBUTION[x] for x in dict.fromkeys(srcs) if x in SOURCE_ATTRIBUTION)
    description = (f"Relief DEM : {'plus grande dimension' if clip is not None else 'diamètre'} {args.diameter_mm:g} mm, "
                   f"exagération verticale x{args.vexag:g}, socle {args.base_mm:g} mm"
                   + (f", source {used}" if used else ""))
    write_outputs(tris, paths, stem, description, attribution, progress=progress)
    log("Terminé : " + ", ".join(os.path.abspath(p) for p in paths.values()))


if __name__ == "__main__":
    main()
