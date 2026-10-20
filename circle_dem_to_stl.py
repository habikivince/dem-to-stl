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

import formats
import mesh
import mesh_export
import geocode
import polygon
import preview
import raster
import regular_polygon
import report
import sources
import stl_io
import vexag_curve

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


def estimate_resources(diameter_mm, spacing_mm, polygon_mode=False, bottom="grid", area_mm2=None, size_mm=None):
    """Estimation (grille n x n, triangles, taille STL, RAM) d'après le pas d'impression.
    Cercle : pi*n^2 triangles ; contour : borne haute 4*n^2. Constantes mesurées sur un cas
    de 4,9 M de triangles ; au-delà de ce cas, c'est une extrapolation. Avec un fond fusionné
    (--bottom merged ou merged-rows), le fond ne compte plus : environ la moitié des triangles
    (mesuré à Ø250 mm : pic de 645 Mo pour 840 Mo estimés). Avec --top adaptive, c'est une borne
    haute (le nombre réel de triangles dépend du relief). Avec une forme (--sides) : area_mm2 est son aire
    imprimée et size_mm sa plus grande dimension ; la borne est alors 4 x aire / pas² (même formule que
    le disque : pour Ø D, 4 x (pi D²/4) / pas² = pi n²)."""
    n = (diameter_mm if size_mm is None else size_mm) / spacing_mm
    triangles = (4.0 if polygon_mode else math.pi) * n * n
    if area_mm2 is not None:
        triangles = 4.0 * area_mm2 / (spacing_mm * spacing_mm)
    if bottom != "grid":
        triangles /= 2
    return {"grid": n, "triangles": triangles, "stl_bytes": 84 + 50 * triangles,
            "ram_bytes": RAM_BASE_BYTES + RAM_BYTES_PER_TRIANGLE * triangles}


def check_resources(diameter_mm, spacing_mm, polygon_mode, force, log=print, ram_available=available_ram_bytes,
                    formats=("stl",), bottom="grid", area_mm2=None, size_mm=None):
    est = estimate_resources(diameter_mm, spacing_mm, polygon_mode, bottom, area_mm2, size_mm)
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
    "austria": "© BEV – Bundesamt für Eich- und Vermessungswesen (CC BY 4.0)",
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
    sides = getattr(args, "sides", None)
    if sides is not None:
        name += f"_{sides}cotes-{'inscrit' if args.circle_mode == 'inscribed' else 'circonscrit'}"
        if getattr(args, "corner_radius", None):
            name += f"_r{args.corner_radius:g}mm"
    name = name if args.vexag == 1 else f"{name}_x{args.vexag:g}"
    return name + "_courbe" if getattr(args, "vexag_curve", None) else name


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


def source_options(args):
    return formats.SourceOptions(
        layer=args.layer, band=args.band, source_crs=args.source_crs, vertical_unit=args.vertical_unit,
        lidar_select=args.lidar_select, lidar_cell=args.lidar_cell, lidar_interp=args.lidar_interp,
        lidar_max_gap=args.lidar_max_gap, lidar_spike_m=args.lidar_spike_m)


def validate_source_options(args):
    opts = source_options(args)
    given = opts.raster_options() + opts.lidar_options()
    if given and args.input is None:
        fail(f"{', '.join(given)} ne s'applique qu'avec --input (fichier fourni par l'utilisateur).")
    if args.band is not None and args.band < 1:
        fail("--band commence à 1.")
    if args.lidar_cell is not None and not _pos(args.lidar_cell):
        fail("--lidar-cell doit être un nombre strictement positif (mètres).")
    for name, value in (("--lidar-max-gap", args.lidar_max_gap), ("--lidar-spike-m", args.lidar_spike_m)):
        if value is not None and (not math.isfinite(value) or value < 0):
            fail(f"{name} doit être un nombre positif ou nul (mètres).")


DEFAULT_CORNER_TOL_MM = 0.005     # flèche maximale admise sur les arcs d'arrondi, en mm imprimés
MIN_CORNER_TOL_MM = 1e-6          # en dessous, la précision du STL (float32) ne suit plus


def build_shape(args, log=log):
    """Polygone régulier (en mètres de terrain) décrit par --sides, --circle-mode, --corner-radius, --rotation.
    Le cercle de référence est celui de --radius (mètres de terrain) et de --diameter-mm (son diamètre imprimé) :
    même échelle que pour un disque. --corner-radius est donné en mm imprimés. Échoue avec un message clair ;
    renseigne args.corner_tol_mm (défaut DEFAULT_CORNER_TOL_MM)."""
    if not regular_polygon.MIN_SIDES <= args.sides <= regular_polygon.MAX_SIDES:
        fail(f"--sides doit être compris entre {regular_polygon.MIN_SIDES} et {regular_polygon.MAX_SIDES} "
             f"(reçu : {args.sides}).")
    scale = (args.diameter_mm / 2.0) / args.radius                  # mm imprimés par mètre de terrain
    corner_mm = 0.0 if args.corner_radius is None else args.corner_radius
    rotation = 0.0 if args.rotation is None else args.rotation
    if not (math.isfinite(corner_mm) and corner_mm >= 0):
        fail("--corner-radius doit être un nombre positif ou nul (en mm imprimés ; 0 = angles vifs).")
    if not math.isfinite(rotation):
        fail("--rotation doit être un nombre fini (en degrés, sens antihoraire).")
    if args.corner_tol_mm is None:
        args.corner_tol_mm = DEFAULT_CORNER_TOL_MM
    if not (math.isfinite(args.corner_tol_mm) and args.corner_tol_mm >= MIN_CORNER_TOL_MM):
        fail(f"--corner-tol-mm doit être un nombre >= {MIN_CORNER_TOL_MM:g} (en mm imprimés).")
    if args.corner_segments is not None and args.corner_segments < 1:
        fail("--corner-segments doit être un entier >= 1.")
    try:
        sharp = regular_polygon.RegularPolygon(args.sides, args.circle_mode, args.radius, 0.0, rotation)
    except regular_polygon.ShapeError as e:
        fail(str(e))
    apothem_mm = sharp.apothem * scale
    if corner_mm > apothem_mm * (1.0 + 1e-9):
        fail(f"--corner-radius ({corner_mm:g} mm) dépasse l'apothème de la forme ({apothem_mm:.4g} mm) : les "
             f"arrondis de deux sommets voisins se chevaucheraient. Maximum : {apothem_mm:.4g} mm "
             f"(à cette valeur, le contour devient un cercle de rayon {apothem_mm:.4g} mm).")
    shape = regular_polygon.RegularPolygon(args.sides, args.circle_mode, args.radius, min(corner_mm / scale, sharp.apothem),
                                           rotation)
    shape_mm = shape.scaled(scale)
    try:
        k = shape_mm.arc_segments(args.corner_tol_mm, args.corner_segments)
    except regular_polygon.ShapeError as e:
        fail(f"{e} Augmente --corner-tol-mm ou fixe --corner-segments.")
    xmin, xmax, ymin, ymax = shape_mm.bbox()
    log("forme : " + shape_mm.describe("mm", 2))
    log(f"pièce : {xmax - xmin:.2f} x {ymax - ymin:.2f} mm (cercle de référence : Ø {2 * shape_mm.ref_radius:.2f} mm "
        f"= --diameter-mm), aire {shape_mm.area():.1f} mm²")
    if corner_mm > 0:
        sag = corner_mm * (1.0 - math.cos(shape.half_angle / k))
        log(f"arrondi : {k} segment(s) par arc, flèche maximale {sag:.5f} mm "
            + ("(imposé par --corner-segments)" if args.corner_segments is not None
               else f"(tolérance demandée {args.corner_tol_mm:g} mm)"))
        if corner_mm < args.print_spacing_mm:
            log(f"ATTENTION : rayon d'arrondi ({corner_mm:g} mm) inférieur au pas d'impression "
                f"({args.print_spacing_mm:g} mm) : le maillage du dessus ne le résout pas, seul le bord exact le suit.")
    return shape


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


def check_shape_coverage(shape, inside, kept, xs, ys, depth_m, scale):
    """Refuse un masque de données qui ne couvre pas la forme. Le bord exact d'un polygone projette les nœuds de
    bord du masque sur le contour théorique : si une partie de la forme est sans donnée (DEM trop petit, nodata,
    valeurs sous --min-elevation), le bord du masque est loin du contour et la projection étirerait le maillage
    sur des dizaines de pas. Seules comptent les cellules à plus de `depth_m` (mètres) du bord de la forme ;
    l'ouverture morphologique n'en retire jamais (elle ne mange qu'une cellule)."""
    missing = inside & ~kept
    if not missing.any():
        return
    ii, jj = np.nonzero(missing)
    gx, gy = xs[jj], ys[ii]
    px, py = shape.project(gx, gy)
    depth = np.hypot(px - gx, py - gy)
    deep = depth > depth_m
    if deep.any():
        fail(f"{int(deep.sum())} point(s) de la grille, à l'intérieur de la forme et jusqu'à "
             f"{depth[deep].max() * scale:.1f} mm imprimés de son bord, n'ont pas de donnée valide : le DEM ne couvre "
             f"pas toute la forme (rappel : un polygone circonscrit dépasse le cercle de référence), ou contient des "
             f"valeurs sans donnée ou inférieures à --min-elevation. Un bord exact ne peut pas suivre un masque "
             f"incomplet. Fournis un DEM plus grand, ou utilise --sea-level (zone d'eau) ou --fill-with-copernicus.")


def write_report_file(path, args, paths, info, used, attribution, clip, shape, preview_file=None):
    """Écrit le rapport JSON (--report). Un échec d'écriture est signalé sans faire échouer la génération : les
    modèles sont déjà écrits."""
    if shape is not None:
        kind = "regular-polygon"
    elif args.clip_to_place:
        kind = "place-contour"
    elif clip is not None:
        kind = "polygon"
    else:
        kind = "circle"
    inputs = args.input if isinstance(args.input, (list, tuple)) else [args.input]
    data = report.build_report(
        outputs={**paths, **({"preview": preview_file} if preview_file else {})},
        model={"shape": kind, "diameter_mm": args.diameter_mm, "base_mm": args.base_mm, "vexag": args.vexag,
               "scale_mm_per_m": info.get("scale_mm_per_m"), "pixel_size_m": info.get("pixel_size_m"),
               "print_spacing_mm": args.print_spacing_mm, "grid": info.get("grid"), "cells": info.get("cells"),
               "vexag_curve": info.get("vexag_curve")},
        relief=info.get("relief"),
        source={"used": used, "input": [os.fspath(x) for x in inputs if x is not None],
                "fill_with_copernicus": bool(args.fill_with_copernicus), "attribution": attribution or None},
        mesh_info={k: info.get(k) for k in ("triangles_top", "triangles_bottom", "triangles_wall", "triangles_total")},
        location={"lat": args.lat, "lon": args.lon, "cx": args.cx, "cy": args.cy, "radius_m": info.get("radius_m"),
                  "place": args.place, "polygon": args.polygon},
        options={"top": args.top, "top_tol_mm": args.top_tol_mm, "bottom": args.bottom, "diagonal": args.diagonal,
                 "resample": args.resample, "smooth_wall": bool(args.smooth_wall), "sea_level": args.sea_level,
                 "min_elevation": args.min_elevation, "land_threshold": args.land_threshold})
    try:
        report.write_report(path, data)
        log(f"écrit : {os.path.abspath(path)}")
    except (OSError, ValueError, TypeError) as e:
        log(f"ATTENTION : rapport non écrit ({e}) ; les modèles, eux, sont écrits.")


def build_mesh(input_path, lat, lon, cx, cy, radius, diameter_mm, vexag, base_mm,
                print_spacing_mm, sea_level=None, min_elevation=-50.0, land_threshold=1.0,
                smooth_wall=False, clip=None, resample="auto", bottom="grid", bottom_max_cells=None,
                diagonal="fixed", top="grid", top_tol_mm=0.05, shape=None, corner_tol_mm=DEFAULT_CORNER_TOL_MM,
                corner_segments=None, info=None, preview_path=None, vexag_curve=None):
    """info : dict optionnel, rempli avec ce que --report écrit (échelle, relief, nombres de triangles...).
    preview_path : si donné, aperçu ombré PNG de la pièce, écrit avant le maillage (voir preview.py).
    vexag_curve : AltitudeCurve (vexag_curve.py) ou None : correspondance altitude -> altitude équivalente appliquée
    avant l'échelle et l'exagération."""
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

    window_radius, outline = radius, None
    if shape is not None:
        # polygone régulier (--sides) : `radius` reste celui du cercle de référence (il fixe l'échelle) ; la fenêtre
        # du DEM enveloppe la forme, le masque est le test exact d'appartenance, le bord est le contour exact
        inside_fn = lambda xs_, ys_: shape.contains(xs_[None, :] - cx, ys_[:, None] - cy)  # noqa: E731
        window_radius = shape.half_extent()
        outline = shape.scaled(scale).discretize(corner_tol_mm, corner_segments)

    elev, xs, ys, nodata = raster.warp_window(input_path, dst_srs_wkt, cx, cy, window_radius, pixel_size_m,
                                              inside_fn=inside_fn, resample=resample)
    H, W = elev.shape
    log(f"grille ré-échantillonnée : {W} x {H} points")

    inside = (inside_fn(xs, ys) if inside_fn is not None
              else (xs[None, :] - cx) ** 2 + (ys[:, None] - cy) ** 2 <= radius ** 2)
    valid = np.isfinite(elev) & ((elev != nodata) if nodata is not None else True)   # NaN/inf ne sont jamais des altitudes
    plausible = valid & (elev > min_elevation)
    n_zone = int(inside.sum())
    n_no_data = int((inside & ~valid).sum())
    n_low = int((inside & valid & ~plausible).sum())
    if n_no_data:
        log(f"{n_no_data} cellule(s) sans donnée valide dans la zone ({100.0 * n_no_data / max(n_zone, 1):.1f} %) : "
            + ("remplacées par le niveau de la mer" if sea_level is not None else
               "exclues du modèle (trous), le relief n'est jamais inventé"))
    if n_low:
        log(f"{n_low} cellule(s) valide(s) de la zone sous --min-elevation ({min_elevation:g} m) "
            + ("traitée(s) comme artefacts" if sea_level is not None else
               "EXCLUE(S) : abaisse --min-elevation pour conserver un relief situé sous cette altitude"))

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
    if shape is not None:
        check_shape_coverage(shape, inside, clean_mask, xs - cx, ys - cy, 3.0 * pixel_size_m, scale)
    log(f"{int(quad_ok.sum())} cellules incluses dans "
        f"{'le contour' if clip is not None else 'la forme' if shape is not None else 'le disque'}")

    min_elev = sea_level if sea_level is not None else elev[clean_mask].min()

    # hauteur du relief imprimé, affichée AVANT le maillage (étape longue) : on peut l'interrompre si elle surprend
    top_alt = elev[report.touched_nodes(quad_ok)]
    top_alt = top_alt[np.isfinite(top_alt)]
    max_alt = float(top_alt.max()) if top_alt.size else min_elev
    elev_mesh, min_mesh = elev, min_elev
    equivalent_relief = None
    if vexag_curve is not None:
        elev_mesh, min_mesh = vexag_curve(elev), float(vexag_curve(min_elev))
        equivalent_relief = float(vexag_curve(max_alt)) - min_mesh
    relief = report.relief_stats(min_elev, max_alt, scale, vexag, diameter_mm, equivalent_relief)
    for line in report.describe_relief(relief):
        log(line)
    if vexag_curve is not None:
        for line in vexag_curve.describe(vexag):
            log(line)
    if preview_path:
        try:
            # aperçu du relief IMPRIMÉ : altitudes après la courbe éventuelle (elev_mesh), pas les altitudes brutes
            w_px, h_px = preview.write_preview(preview_path, elev_mesh, quad_ok, pixel_size_m, vexag)
            log(f"aperçu : {os.path.abspath(preview_path)} ({w_px} x {h_px} px, nord en haut, ombrage du relief imprimé)")
        except (OSError, ValueError) as e:
            log(f"ATTENTION : aperçu non écrit ({e})")

    tris, n_top, n_bot, n_wall = mesh.build_mesh_from_arrays(
        elev_mesh, quad_ok, xs - cx, ys - cy, scale, vexag, base_mm, min_mesh, progress=progress,
        smooth_wall=smooth_wall, circle_radius_mm=(radius * scale if smooth_wall else None),
        bottom=bottom, bottom_max_cells=bottom_max_cells, diagonal=diagonal,
        top=top, top_tol_mm=top_tol_mm, outline=outline)
    log(f"{n_top} triangles (dessus) + {n_bot} (dessous) + {n_wall} (paroi) = {len(tris)} triangles au total")
    if info is not None:
        info.update(scale_mm_per_m=scale, pixel_size_m=pixel_size_m, radius_m=radius, grid=[W, H],
                    cells=int(quad_ok.sum()), relief=relief, triangles_top=n_top, triangles_bottom=n_bot,
                    triangles_wall=n_wall, triangles_total=len(tris),
                    vexag_curve=None if vexag_curve is None else vexag_curve.as_list())
    return tris


def main():
    p = argparse.ArgumentParser(description="Découpe circulaire d'un DEM + export STL solide.")
    p.add_argument("--input", nargs="+", default=None, metavar="FICHIER",
                    help="Fichier(s) source : GeoTIFF (.tif), ESRI ASCII Grid (.asc), SRTM (.hgt), GeoPackage "
                         "raster (.gpkg), NetCDF (.nc), HDF5 géoréférencé (.h5), VRT, tout raster lisible par "
                         "GDAL, ou nuage de points LAS/LAZ ; plusieurs fichiers = mosaïque (tuiles voisines). "
                         "Omis si --source / --download-copernicus est utilisé")
    g = p.add_argument_group("Formats d'entrée (voir README, section « Formats d'entrée »)")
    g.add_argument("--layer", default=None, metavar="NOM",
                   help="Table raster GeoPackage, variable NetCDF ou dataset HDF5 à utiliser (obligatoire s'il "
                        "y en a plusieurs et qu'une seule ne se déclare pas comme une altitude)")
    g.add_argument("--band", type=int, default=None, metavar="N",
                   help="Bande à utiliser (à partir de 1) d'un fichier multibande ou d'une variable à dimension "
                        "supplémentaire (temps...) : jamais choisie à ta place")
    g.add_argument("--source-crs", default=None, metavar="CRS",
                   help="CRS à attribuer à un fichier qui n'en déclare aucun (ex. EPSG:2056) ; refusé si le "
                        "fichier en a déjà un")
    g.add_argument("--vertical-unit", choices=sorted(formats.VERTICAL_UNIT_CHOICES), default=None,
                   help="Unité des altitudes si le fichier ne la déclare pas (sinon supposée en mètres, avec "
                        "avertissement) ou pour corriger celle qu'il déclare")
    g.add_argument("--lidar-select", choices=["ground", "lowest", "highest"], default=None,
                   help="LAS/LAZ — points retenus : ground = classe 2 (défaut, exige un nuage classifié), "
                        "lowest = minimum par cellule (nuage non classifié : approximation du sol), "
                        "highest = maximum par cellule (surface avec végétation et bâti)")
    g.add_argument("--lidar-cell", type=float, default=None, metavar="M",
                   help="LAS/LAZ — taille de la cellule de la grille, en mètres (défaut : pixel de travail, "
                        "dérivé de --print-spacing-mm)")
    g.add_argument("--lidar-interp", choices=["linear", "nearest", "none"], default=None,
                   help="LAS/LAZ — comblement des cellules sans point : linear = triangulation (défaut), "
                        "nearest = valeur de la cellule remplie la plus proche, none = aucun")
    g.add_argument("--lidar-max-gap", type=float, default=None, metavar="M",
                   help="LAS/LAZ — distance maximale de comblement, en mètres (défaut : 50 cellules) ; au-delà, "
                        "les cellules restent sans donnée")
    g.add_argument("--lidar-spike-m", type=float, default=None, metavar="M",
                   help="LAS/LAZ — une cellule isolée dépassant tous ses voisins de plus de cette marge "
                        "(en plus de leur étendue) est écartée comme aberrante (défaut 10 ; 0 = désactivé)")
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
    p.add_argument("--report", nargs="?", const="auto", default=None, metavar="FICHIER",
                    help="Écrit un rapport JSON (source et mention de licence, échelle, hauteur de relief, exagération, "
                         "triangles, ligne de commande). Sans nom : à côté du modèle, même nom en .json")
    p.add_argument("--preview", nargs="?", const="auto", default=None, metavar="FICHIER",
                    help="Écrit un aperçu ombré PNG de la zone découpée, avant le maillage (relief imprimé, nord en "
                         "haut). Sans nom : à côté du modèle, même nom en .png")
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
                    help="Rayon du cercle, en mètres (requis sans --polygon). Avec --sides : rayon du cercle de "
                         "référence du polygone")
    p.add_argument("--diameter-mm", type=float, default=None,
                    help="Diamètre final imprimé, en mm (demandé interactivement si omis)")
    p.add_argument("--vexag", type=float, default=None,
                    help="Exagération verticale (demandée interactivement si omise)")
    p.add_argument("--vexag-curve", default=None, metavar="ALT:EQ,...",
                    help="Courbe d'exagération non linéaire : points « altitude:altitude_équivalente » en mètres "
                         "(ex. 0:0,500:500,2000:1100,4000:1800), interpolés en PCHIP monotone, appliqués avant l'échelle "
                         "et --vexag. Déforme le relief : la pièce n'est pas une réplique fidèle")
    p.add_argument("--base-mm", type=float, default=None,
                    help="Épaisseur du socle plat, en mm (demandé interactivement si omis)")
    p.add_argument("--print-spacing-mm", type=float, default=0.2,
                    help="Résolution cible du maillage, en mm (~0.2 recommandé pour buse 0.4mm)")
    p.add_argument("--resample", choices=raster.RESAMPLE_CHOICES, default="auto",
                    help="Ré-échantillonnage du DEM vers la grille de travail. auto (défaut) : moyenne si "
                         "la grille est plus grossière que la source, bilinéaire si elle est plus fine "
                         "(la moyenne de GDAL donnerait alors des paliers plats). cubic dépasse près des "
                         "falaises (≈ 7 %% de la marche) ; cubicspline lisse sans dépasser.")
    p.add_argument("--bottom", choices=mesh.BOTTOM_MODES, default="grid",
                    help="Triangulation du fond plat. grid (défaut) : cellule par cellule. merged "
                         "(PROTOTYPE) : carrés d'un quadtree équilibré en éventail, même géométrie, environ "
                         "moitié moins de triangles au total, triangles compacts. merged-rows : rectangles "
                         "gloutons, encore moins de triangles mais très allongés (comparaison).")
    p.add_argument("--bottom-max-cells", type=int, default=0,
                    help="Avec --bottom merged-rows : taille maximale, en cellules, d'un côté de rectangle "
                         "(0 = pas de limite).")
    p.add_argument("--diagonal", choices=mesh.DIAGONAL_MODES, default="fixed",
                    help="Diagonale de coupe des cellules du dessus. fixed (défaut) : toujours la même. "
                         "adaptive (PROTOTYPE) : celle qui relie les deux coins de plus faible écart "
                         "d'altitude ; même nombre de triangles, erreur au centre des cellules réduite de 4 à "
                         "7 %% en moyenne sur six relevés réels.")
    p.add_argument("--top", choices=mesh.TOP_MODES, default="grid",
                    help="Triangulation du dessus. grid (défaut) : un ou deux triangles par cellule, exacte aux "
                         "nœuds. adaptive (PROTOTYPE) : maillage RTIN à erreur bornée, très peu de triangles là où "
                         "le relief est lisse, bord et parois à pleine résolution.")
    p.add_argument("--top-tol-mm", type=float, default=0.05,
                    help="Avec --top adaptive : tolérance verticale, en mm d'impression (défaut 0,05 mm). L'écart "
                         "entre un nœud du DEM et la surface reste de l'ordre de cette valeur (voir adaptive_top.py).")
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
    p.add_argument("--sides", type=int, default=None,
                    help=f"Polygone régulier à N côtés ({regular_polygon.MIN_SIDES} à {regular_polygon.MAX_SIDES}) "
                         "à la place du disque, associé au cercle de référence (--radius, --diameter-mm). "
                         "Requiert --circle-mode. Incompatible avec --polygon, --clip-to-place, --smooth-wall. "
                         "Orientation par défaut : un sommet au nord pour N impair, un côté au nord pour N pair.")
    p.add_argument("--circle-mode", choices=regular_polygon.CIRCLE_MODES, default=None,
                    help="Avec --sides : inscribed = les sommets sont sur le cercle de référence ; circumscribed = "
                         "les côtés sont tangents au cercle de référence (la pièce dépasse alors de --diameter-mm). "
                         "Le cercle de référence décrit toujours le polygone à angles vifs, même arrondi.")
    p.add_argument("--corner-radius", type=float, default=None,
                    help="Avec --sides : rayon d'arrondi des sommets, en mm imprimés (raccordement tangent exact "
                         "aux côtés ; 0 ou omis = angles vifs). Maximum : l'apothème de la pièce (à ce maximum, le "
                         "contour devient un cercle).")
    p.add_argument("--rotation", type=float, default=None,
                    help="Avec --sides : rotation du polygone autour du centre, en degrés, sens antihoraire "
                         "(défaut 0 = orientation par défaut)")
    p.add_argument("--corner-tol-mm", type=float, default=None,
                    help=f"Avec --sides et --corner-radius : flèche maximale admise entre un arc d'arrondi et ses "
                         f"cordes, en mm imprimés (défaut {DEFAULT_CORNER_TOL_MM:g}). Les sommets du maillage sont "
                         f"sur l'arc exact : l'écart est d'un seul côté (cordes à l'intérieur) et borné.")
    p.add_argument("--corner-segments", type=int, default=None,
                    help="Avec --sides et --corner-radius : nombre de segments par arc d'arrondi, imposé "
                         "(remplace --corner-tol-mm)")
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

    if args.sides is None:
        stray = [n for n, v in (("--circle-mode", args.circle_mode), ("--corner-radius", args.corner_radius),
                                ("--rotation", args.rotation), ("--corner-tol-mm", args.corner_tol_mm),
                                ("--corner-segments", args.corner_segments)) if v is not None]
        if stray:
            sys.exit("Erreur : " + ", ".join(stray) + " requiert --sides (polygone régulier).")
    else:
        bad = [n for n, v in (("--polygon", args.polygon), ("--clip-to-place", args.clip_to_place or None),
                              ("--smooth-wall", args.smooth_wall or None)) if v is not None]
        if bad:
            sys.exit("Erreur : --sides est incompatible avec " + ", ".join(bad)
                     + " (le contour exact du polygone remplace déjà le bord lissé ; la forme vient de --sides).")
        if args.circle_mode is None:
            sys.exit("Erreur : --sides requiert --circle-mode inscribed|circumscribed (polygone inscrit dans le "
                     "cercle de référence, ou circonscrit à lui).")

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
    validate_source_options(args)

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
    curve = None
    if args.vexag_curve is not None:
        try:
            curve = vexag_curve.parse(args.vexag_curve)
        except vexag_curve.CurveError as e:
            fail(f"--vexag-curve : {e}")
    shape = build_shape(args, log=raster.log) if args.sides is not None else None
    # rayon du disque dont les données sont à récupérer : le cercle de référence, ou (forme) le disque qui contient
    # toute la pièce (rayon extrême du contour : un polygone circonscrit dépasse le cercle de référence). Les
    # sources ajoutent leurs propres marges ; la fenêtre de travail (demi-côté du carré englobant) est incluse.
    data_radius = args.radius if shape is None else shape.extreme_radius
    paths = resolve_outputs(args)
    report_file = report.report_path(paths, args.report) if args.report else None
    preview_file = preview.preview_path(paths, args.preview) if args.preview else None
    est_kw = {}
    if shape is not None:
        scale = (args.diameter_mm / 2.0) / args.radius
        shape_mm = shape.scaled(scale)
        xmin, xmax, ymin, ymax = shape_mm.bbox()
        est_kw = {"area_mm2": shape_mm.area(), "size_mm": max(xmax - xmin, ymax - ymin)}
    check_resources(args.diameter_mm, args.print_spacing_mm, clip is not None, args.force, log=raster.log,
                    formats=tuple(paths), bottom=args.bottom, **est_kw)
    if args.dry_run:
        log("--dry-run : aucun téléchargement, calcul ni écriture. Fichiers prévus : "
            + ", ".join(os.path.abspath(p) for p in list(paths.values()) + [f for f in (preview_file, report_file) if f]))
        return
    input_kind = "raster"
    if args.input is not None:
        # détection du format -> lecteur adapté -> normalisation (voir formats.py) ; échoue ici, avant d'écrire
        try:
            input_kind, args.input = formats.prepare(args.input, source_options(args), log=raster.log)
        except raster.InputError as e:
            fail(str(e))
        if input_kind == "lidar" and args.fill_with_copernicus:
            fail("--fill-with-copernicus est incompatible avec une entrée LAS/LAZ.")
    prepare_outputs({**paths, **({"png": preview_file} if preview_file else {}),
                     **({"json": report_file} if report_file else {})})

    used = None
    if source:
        pixel_size_m = args.print_spacing_mm / ((args.diameter_mm / 2.0) / args.radius)
        if args.full_res:
            raster.log(f"--full-res : résolution native maximale demandée (le maillage reste à "
                       f"{pixel_size_m:.2f} m/pixel d'après --print-spacing-mm)")
            pixel_size_m = sources.FULL_RES_PIXEL_M
        try:
            used, args.input = sources.acquire(
                source, args.lat, args.lon, data_radius, args.cache_dir, pixel_size_m,
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
                args.input, args.lat, args.lon, data_radius, args.cache_dir,
                product_m=args.copernicus_product, log=raster.log)

    if args.top == "adaptive" and args.diagonal != "fixed":
        fail("--diagonal adaptive est sans objet avec --top adaptive")
    if args.top_tol_mm < 0:
        fail("--top-tol-mm doit être >= 0")
    info = {}
    try:
        raster.check_readable(args.input)
        tris = build_mesh(args.input, args.lat, args.lon, args.cx, args.cy, args.radius,
                          args.diameter_mm, args.vexag, args.base_mm, args.print_spacing_mm,
                          sea_level=args.sea_level, min_elevation=args.min_elevation,
                          land_threshold=args.land_threshold, smooth_wall=args.smooth_wall, clip=clip,
                          resample=args.resample, bottom=args.bottom,
                          bottom_max_cells=args.bottom_max_cells or None, diagonal=args.diagonal,
                          top=args.top, top_tol_mm=args.top_tol_mm, shape=shape,
                          corner_tol_mm=args.corner_tol_mm if shape is not None else DEFAULT_CORNER_TOL_MM,
                          corner_segments=args.corner_segments, info=info, preview_path=preview_file,
                          vexag_curve=curve)
    except raster.InputError as e:
        fail(str(e))
    except RuntimeError as e:
        if args.debug:
            raise
        fail(f"traitement GDAL impossible : {str(e).strip().splitlines()[0] if str(e).strip() else e} "
             f"(--debug pour la trace complète)")
    preview_done = preview_file if preview_file and os.path.exists(preview_file) else None
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
    if shape is not None:
        description = (f"Relief DEM : {shape.name()} {args.circle_mode} (cercle de référence Ø {args.diameter_mm:g} mm"
                       + (f", sommets arrondis r = {args.corner_radius:g} mm" if args.corner_radius else "")
                       + f"), exagération verticale x{args.vexag:g}, socle {args.base_mm:g} mm"
                       + (f", source {used}" if used else ""))
    write_outputs(tris, paths, stem, description, attribution, progress=progress)
    if report_file:
        write_report_file(report_file, args, paths, info, used, attribution, clip, shape, preview_done)
    log("Terminé : " + ", ".join(os.path.abspath(p) for p in list(paths.values()) + [f for f in (preview_done, report_file) if f]))


if __name__ == "__main__":
    main()
