# SPDX-License-Identifier: GPL-3.0-or-later
"""
raster.py
---------
Tout ce qui touche GDAL : ouverture du fichier source, détection et
reprojection automatique du CRS, ré-échantillonnage fenêtré autour du
centre choisi. Aucune logique géométrique ici (voir mesh.py).
"""
import math
import sys

import numpy as np
from osgeo import gdal, osr
from scipy import ndimage

gdal.UseExceptions()
osr.UseExceptions()


def log(msg):
    print(f"[DEM2STL] {msg}", flush=True)


class InputError(ValueError):
    """Raster d'entrée absent ou illisible."""


def check_readable(input_path):
    """Ouvre le raster principal pour échouer tôt avec un message clair (au lieu d'un
    traceback GDAL au milieu du traitement)."""
    path = _primary(input_path)
    if getattr(path, "is_lidar", False):
        return   # nuage de points : déjà ouvert et validé par formats.prepare
    try:
        ds = gdal.Open(path)
    except RuntimeError as e:
        reason = str(e).strip().splitlines()[0] if str(e).strip() else "format non reconnu"
        raise InputError(f"raster illisible : {path} ({reason})")
    if ds is None or ds.RasterCount < 1:
        raise InputError(f"raster sans bande : {path}")


def _primary(input_path):
    """input_path peut être un chemin unique ou une liste de rasters (sources
    nationales, éventuellement dans des CRS différents) : on utilise le
    premier pour lire CRS et nodata, gdal.Warp se charge du reste."""
    return input_path[0] if isinstance(input_path, (list, tuple)) else input_path


def _source_srs(input_path):
    """CRS de la source principale (raster, ou nuage LAS/LAZ qui expose `srs_wkt`). Une source sans CRS est
    refusée : sans lui, rien ne dit si les coordonnées sont des degrés ou des mètres."""
    primary = _primary(input_path)
    wkt = primary.srs_wkt if getattr(primary, "is_lidar", False) else gdal.Open(primary).GetProjection()
    if not wkt:
        raise InputError("la source n'a pas de système de coordonnées (CRS) : degrés ou mètres, impossible de le "
                         "savoir. Indique-le avec --source-crs (ex. EPSG:2056) ou géoréférence le fichier.")
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return srs


def _is_web_mercator(srs):
    """EPSG:3857 est « projeté » mais ses mètres sont dilatés de 1/cos(lat) :
    inutilisable comme CRS de travail métrique."""
    ref = osr.SpatialReference()
    ref.ImportFromEPSG(3857)
    ref.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return bool(srs.IsSame(ref))


def resolve_target_crs_and_center(input_path, lat, lon, cx, cy):
    """
    Détermine le CRS de travail (toujours projeté, en mètres) et le
    centre du cercle dans ce CRS, quelle que soit la combinaison de
    paramètres fournie :
      - --cx/--cy déjà exprimés dans le CRS du fichier source (comportement
        historique, toujours supporté) ;
      - --lat/--lon en WGS84 (EPSG:4326), converti automatiquement.

    Si le fichier source est en CRS géographique (degrés), une zone UTM
    est calculée automatiquement à partir du point de centre et la
    reprojection se fera à la volée pendant le ré-échantillonnage —
    plus besoin de gdalwarp manuel au préalable.

    Retourne (dst_srs_wkt, cx, cy) où cx/cy sont dans dst_srs_wkt.
    """
    src_srs = _source_srs(input_path)

    have_latlon = lat is not None and lon is not None
    have_xy = cx is not None and cy is not None
    if not have_latlon and not have_xy:
        sys.exit("Erreur interne : ni --cx/--cy ni --lat/--lon résolus avant resolve_target_crs_and_center.")

    web_mercator = (not src_srs.IsGeographic()) and _is_web_mercator(src_srs)
    # CRS projeté dont l'unité n'est pas le mètre (pieds, pieds US...) : le pipeline raisonne en mètres, on passe
    # donc aussi par une zone UTM plutôt que de lire des pieds comme des mètres
    non_metric = (src_srs.IsProjected() and not web_mercator and abs(src_srs.GetLinearUnits() - 1.0) > 1e-9)
    if src_srs.IsGeographic() or web_mercator or non_metric:
        # Point de référence pour déterminer la zone UTM : --lat/--lon en
        # priorité, sinon --cx/--cy réinterprétés comme lon/lat (ils sont
        # nécessairement dans ce système si la source l'est ; en Web Mercator
        # ils sont d'abord reconvertis en lon/lat).
        if have_latlon:
            ref_lon, ref_lat = lon, lat
        elif web_mercator or non_metric:
            to_geo = osr.CoordinateTransformation(src_srs, _wgs84())
            ref_lon, ref_lat, _ = to_geo.TransformPoint(cx, cy)
        else:
            ref_lon, ref_lat = cx, cy

        utm_zone = int((ref_lon + 180) / 6) + 1
        dst_epsg = (32600 if ref_lat >= 0 else 32700) + utm_zone
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(dst_epsg)
        dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        kind = ("Web Mercator (EPSG:3857, mètres non métriques)" if web_mercator
                else f"projeté en {src_srs.GetLinearUnitsName()} (1 unité = {src_srs.GetLinearUnits():.7g} m)"
                if non_metric else "géographique")
        log(f"CRS source {kind} détecté -> reprojection à la volée vers "
            f"UTM {utm_zone}{'N' if ref_lat >= 0 else 'S'} (EPSG:{dst_epsg})")

        wgs84 = osr.SpatialReference()
        wgs84.ImportFromEPSG(4326)
        wgs84.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        ct = osr.CoordinateTransformation(wgs84, dst_srs)
        x, y, _ = ct.TransformPoint(ref_lon, ref_lat)
        return dst_srs.ExportToWkt(), x, y

    # Source déjà projetée.
    if have_latlon:
        wgs84 = osr.SpatialReference()
        wgs84.ImportFromEPSG(4326)
        wgs84.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        ct = osr.CoordinateTransformation(wgs84, src_srs)
        x, y, _ = ct.TransformPoint(lon, lat)
        log(f"--lat/--lon converti vers le CRS projeté de la source : ({x:.1f}, {y:.1f})")
        return src_srs.ExportToWkt(), x, y

    return src_srs.ExportToWkt(), cx, cy


def _wgs84():
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return srs


FILL_NODATA = -9999.0       # nodata imposé quand un complément est possible et que la source n'en déclare pas
FILL_MIN_DISC_COVERAGE = 0.999
FILL_MAX_OFFSET_M = 50.0
FILL_MIN_OVERLAP_CELLS = 100


def _complete_with_secondary(warped, loader, dst_srs_wkt, bounds, pixel_size_m, nodata, cx, cy, radius,
                             inside_fn=None):
    """Comble les cellules du disque sans donnée avec une source secondaire grossière
    (Copernicus), recalée en altitude sur la source principale d'après leur zone commune.
    Le complément est facultatif : toute erreur laisse le disque tel quel."""
    band = warped.GetRasterBand(1)
    main = band.ReadAsArray().astype(np.float64)
    gt = warped.GetGeoTransform()
    xs = gt[0] + (np.arange(main.shape[1]) + 0.5) * gt[1]
    ys = gt[3] + (np.arange(main.shape[0]) + 0.5) * gt[5]
    inside = (inside_fn(xs, ys) if inside_fn is not None
              else (xs[None, :] - cx) ** 2 + (ys[:, None] - cy) ** 2 <= radius ** 2)
    valid_main = np.isfinite(main) & (main != nodata)
    coverage = float(valid_main[inside].mean())
    if coverage >= FILL_MIN_DISC_COVERAGE:
        return
    log(f"Complément : {100 * (1 - coverage):.1f} % du disque sans donnée dans la source principale")
    try:
        secondary_path = loader()
        sec_ds = gdal.Warp(
            "", secondary_path, format="MEM", dstSRS=dst_srs_wkt, outputBounds=bounds,
            xRes=pixel_size_m, yRes=pixel_size_m,
            resampleAlg="bilinear" if pixel_size_m < 25.0 else "average", dstNodata=nodata)
        sec = sec_ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    except (RuntimeError, OSError) as e:
        log(f"Complément indisponible ({str(e).splitlines()[0][:200]}) : disque laissé partiel")
        return
    if sec.shape != main.shape:
        log("Complément ignoré : grilles incompatibles")
        return
    valid_sec = np.isfinite(sec) & (sec != nodata)
    need = inside & ~valid_main & valid_sec
    if not need.any():
        log("Complément : la source secondaire n'a rien à ajouter ici")
        return
    both = inside & valid_main & valid_sec
    offset = 0.0
    if both.sum() >= FILL_MIN_OVERLAP_CELLS:
        diff = main[both] - sec[both]
        offset = float(np.median(diff))
        iqr = float(np.percentile(diff, 75) - np.percentile(diff, 25))
        if abs(offset) > FILL_MAX_OFFSET_M:
            log(f"Complément abandonné : décalage d'altitude de {offset:+.1f} m entre les deux sources "
                f"(> {FILL_MAX_OFFSET_M:g} m), référence d'altitude incompatible ?")
            return
        log(f"Complément : décalage d'altitude principal - secondaire {offset:+.2f} m "
            f"(médiane sur {int(both.sum())} cellules communes, écart interquartile {iqr:.1f} m) appliqué")
    else:
        log(f"Complément : zone commune trop petite ({int(both.sum())} cellules) pour recaler les altitudes, "
            f"aucun décalage appliqué")
    main[need] = sec[need] + offset
    band.WriteArray(main)
    log(f"Complément : {int(need.sum())} cellules ({100 * need.sum() / inside.sum():.1f} % du disque) "
        f"ajoutées depuis la source secondaire, à plus faible résolution : raccord visible à la limite")


FILL_MAX_DISTANCE_PX = 50   # rayon de recherche de gdal.FillNodata, et taille maximale d'un trou comblé


def closed_gaps(missing, max_distance_px=FILL_MAX_DISTANCE_PX):
    """Parmi les cellules sans donnée (`missing`, booléen), celles d'un TROU : composante connexe (8 voisins) dont
    chaque cellule est à `max_distance_px` pixels au plus d'une cellule avec donnée, et qui, soit ne touche pas le
    bord de la grille (elle est entourée de données), soit le touche en deux endroits distincts au moins (elle
    traverse la grille de part en part : rivière, joint entre deux dalles). Interpoler un trou ne fait que
    raccorder les données qui l'entourent.

    Tout le reste n'est pas un trou : la zone au-delà de l'emprise de la source (elle touche le bord de la fenêtre
    d'un seul côté), la mer d'une source qui la laisse sans donnée, toute zone vide trop grande. Les interpoler
    reviendrait à EXTRAPOLER les altitudes du bord des données sur `max_distance_px` pixels, soit un relief inventé :
    le README annonce que ce qui est plus grand que 50 pixels « is left out of the model ».
    Retourne un masque booléen de même forme que `missing`."""
    out = np.zeros(missing.shape, dtype=bool)
    if not missing.any() or missing.all():
        return out
    labels, n = ndimage.label(missing, structure=np.ones((3, 3), dtype=bool))
    dist = ndimage.distance_transform_edt(missing)
    farthest = np.asarray(ndimage.maximum(dist, labels, index=np.arange(1, n + 1)))
    ok = farthest <= max_distance_px
    h, w = labels.shape
    if min(h, w) < 2:
        ring = labels.ravel()      # grille dégénérée : tout est bord
    else:
        # le pourtour dans l'ordre (haut, droite, bas, gauche) : deux cellules consécutives sont voisines, donc de la
        # même composante si toutes deux sont sans donnée
        ring = np.concatenate([labels[0, :], labels[1:, -1], labels[-1, -2::-1], labels[-2:0:-1, 0]])
    starts = (ring > 0) & (ring != np.roll(ring, 1))                 # début de chaque tronçon sans donnée du pourtour
    runs = np.bincount(ring[starts], minlength=n + 1)[1:]            # nombre de tronçons par composante
    touches = np.zeros(n + 1, dtype=bool)
    touches[np.unique(ring)] = True
    ok &= ~touches[1:] | (runs >= 2)
    return np.concatenate(([False], ok))[labels]


RESAMPLE_CHOICES = ("auto", "average", "bilinear", "cubic", "cubicspline")
# Au-delà de ce rapport (pixel de la source / pixel de travail), on considère que la source est
# AGRANDIE : `average` de GDAL se comporte alors comme `nearest` (blocs constants, terrasses).
UPSAMPLE_RATIO = 1.05


def source_pixel_size_m(input_path, dst_srs_wkt, cx, cy, max_files=8, reduce="min"):
    """Taille du pixel de la source, en mètres dans le CRS cible, au voisinage de (cx, cy), parmi au plus
    ``max_files`` fichiers de la liste : la plus fine (``reduce="min"``, défaut) ou la plus grossière
    (``reduce="max"``) ; sur chaque fichier, l'axe le plus grossier (un pixel en degrés est plus étroit
    en x qu'en y). None si elle ne peut pas être déterminée (l'appelant retombe alors sur `average`,
    comportement historique)."""
    if reduce not in ("min", "max"):
        raise ValueError("reduce doit valoir 'min' ou 'max'")
    paths = list(input_path) if isinstance(input_path, (list, tuple)) else [input_path]
    if any(getattr(p, "is_lidar", False) for p in paths):
        return None   # nuage de points : pas de pixel natif (voir lidar.py)
    if len(paths) > max_files:
        step = len(paths) / max_files
        paths = [paths[int(i * step)] for i in range(max_files)]
    dst = osr.SpatialReference()
    dst.ImportFromWkt(dst_srs_wkt)
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    best = None
    for path in paths:
        try:
            ds = gdal.Open(path)
            gt, wkt = ds.GetGeoTransform(), ds.GetProjection()
            if not wkt:
                continue
            src = osr.SpatialReference()
            src.ImportFromWkt(wkt)
            src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            fwd = osr.CoordinateTransformation(dst, src)
            inv = osr.CoordinateTransformation(src, dst)
            sx, sy, _ = fwd.TransformPoint(cx, cy)
            x0, y0, _ = inv.TransformPoint(sx, sy)
            xx, xy, _ = inv.TransformPoint(sx + gt[1], sy + gt[4])   # un pixel vers la droite
            yx, yy, _ = inv.TransformPoint(sx + gt[2], sy + gt[5])   # un pixel vers le bas
            size = max(math.hypot(xx - x0, xy - y0), math.hypot(yx - x0, yy - y0))
        except (RuntimeError, ValueError, TypeError):
            continue
        if size > 0 and math.isfinite(size):
            best = size if best is None else (min(best, size) if reduce == "min" else max(best, size))
    return best


def choose_resample_alg(mode, src_pixel_m, pixel_size_m):
    """Algorithme de gdal.Warp. ``auto`` : `average` quand la grille de travail est plus grossière que
    la source (réduction : moyenne = anticrénelage), `bilinear` quand elle est plus fine (agrandissement :
    `average` donnerait des paliers plats, voir README). Sans connaître le pixel source : `average`."""
    if mode != "auto":
        return mode
    if src_pixel_m is not None and src_pixel_m > UPSAMPLE_RATIO * pixel_size_m:
        return "bilinear"
    return "average"


def _working_type(gdal_type):
    """Type de la grille de travail : Float32, ou Float64 pour les sources à 32 bits entiers ou plus. Une source
    entière ne doit jamais être ré-échantillonnée en entiers (la moyenne serait arrondie au mètre)."""
    wide = (gdal.GDT_Int32, gdal.GDT_UInt32, gdal.GDT_Float64, getattr(gdal, "GDT_Int64", -1),
            getattr(gdal, "GDT_UInt64", -1))
    return gdal.GDT_Float64 if gdal_type in wide else gdal.GDT_Float32


def warp_window(input_path, dst_srs_wkt, cx, cy, radius, pixel_size_m, progress_cb=None, fill=True,
                inside_fn=None, resample="auto", extrapolate_edges=True, stats=None):
    """Ré-échantillonne une fenêtre carrée de 2*radius de côté autour de
    (cx, cy), dans dst_srs_wkt, à pixel_size_m/pixel. Reprojette à la
    volée si dst_srs_wkt diffère du CRS de la source — gdal.Warp gère
    reprojection et découpe en une seule passe, donc même une source
    volumineuse en CRS géographique n'est jamais traitée en entier.

    Si input_path porte un attribut `fill_loader` (voir sources.with_copernicus_fill),
    les trous du disque sont comblés par cette source secondaire, recalée en altitude.
    fill=False : valeurs brutes, sans comblement ni complément (compare_dems).
    inside_fn(xs, ys) -> masque : zone à compléter si elle n'est pas le disque (mode polygone).
    resample : "auto" (défaut), "average", "bilinear", "cubic" ou "cubicspline". En "auto", un DEM plus
    grossier que la grille de travail est interpolé (bilinear) au lieu d'être répliqué par blocs ; avec
    plusieurs sources de résolutions différentes, la plus grossière décide.

    Valeurs sans donnée : le nodata déclaré par la source, sinon (non déclaré ou NaN) FILL_NODATA ; la zone
    hors de l'emprise de la source est toujours marquée sans donnée, jamais remplie de zéros. Une source
    entière est ré-échantillonnée en virgule flottante (pas d'arrondi des moyennes). Un nuage de points
    (lidar.LidarInput) n'est pas ré-échantillonné : il est rastérisé directement sur cette grille.

    Comblement (fill=True) : par défaut (extrapolate_edges=True), interpolation sur FILL_MAX_DISTANCE_PX pixels
    autour de toute donnée, donc aussi extrapolation au bord de la source. extrapolate_edges=False
    (--fill-holes-only) : seuls les TROUS (closed_gaps : entourés de données, à FILL_MAX_DISTANCE_PX pixels au plus
    d'une donnée) sont interpolés ; les cellules sans donnée hors de l'emprise de la source ou dans une zone vide
    plus grande restent sans donnée : le relief n'y est pas inventé. stats (dict, facultatif) reçoit "edge_cells" :
    masque des cellules qu'une interpolation sur FILL_MAX_DISTANCE_PX pixels remplirait sans qu'elles soient dans un
    trou (remplies si extrapolate_edges, laissées sans donnée sinon)."""
    if resample not in RESAMPLE_CHOICES:
        raise InputError(f"--resample : valeur inconnue {resample!r} (choix : {', '.join(RESAMPLE_CHOICES)})")
    primary = _primary(input_path)
    if getattr(primary, "is_lidar", False):
        return primary.grid_window(dst_srs_wkt, cx, cy, radius, pixel_size_m, progress_cb=progress_cb,
                                   inside_fn=inside_fn, fill=fill, resample=resample)
    ds = gdal.Open(primary)
    src_band = ds.GetRasterBand(1)
    declared = src_band.GetNoDataValue()
    loader = getattr(input_path, "fill_loader", None) if fill else None
    nodata = float(declared) if declared is not None and math.isfinite(declared) else FILL_NODATA
    # NaN dans une grille flottante sans nodata déclaré : traité comme « sans donnée », jamais comme altitude
    src_nodata = "nan" if declared is None and src_band.DataType in (gdal.GDT_Float32, gdal.GDT_Float64) else None

    xmin, xmax = cx - radius, cx + radius
    ymin, ymax = cy - radius, cy + radius

    # La source la plus GROSSIÈRE décide : dès qu'une source est agrandie, `average` la répliquerait par blocs
    # (cas d'une mosaïque qui mêle une couche fine et un repli grossier, ex. GSI DEM1A + DEM5A).
    src_px = source_pixel_size_m(input_path, dst_srs_wkt, cx, cy, reduce="max") if resample == "auto" else None
    alg = choose_resample_alg(resample, src_px, pixel_size_m)
    if alg != "average":
        why = (f"pixel de la source ≈ {src_px:.1f} m > pixel de travail {pixel_size_m:.2f} m : "
               "la source est agrandie" if resample == "auto" else "demandé par --resample")
        log(f"ré-échantillonnage {alg} ({why})")
    elif src_px is not None:
        log(f"résolution : pixel natif de la source ≈ {src_px:.3g} m, grille de travail {pixel_size_m:.3g} m "
            f"({'moyenne des pixels sources' if src_px < pixel_size_m / UPSAMPLE_RATIO else 'résolutions voisines'}"
            f", ré-échantillonnage average)")

    warped = gdal.Warp(
        "", input_path, format="MEM",
        dstSRS=dst_srs_wkt,
        outputBounds=(xmin, ymin, xmax, ymax),
        xRes=pixel_size_m, yRes=pixel_size_m,
        resampleAlg=alg,
        outputType=_working_type(src_band.DataType),
        srcNodata=src_nodata,
        dstNodata=nodata,
        callback=progress_cb or gdal.TermProgress_nocb,
    )
    band = warped.GetRasterBand(1)
    if loader is not None:
        _complete_with_secondary(warped, loader, dst_srs_wkt, (xmin, ymin, xmax, ymax),
                                 pixel_size_m, nodata, cx, cy, radius, inside_fn)
    if fill:  # fill=False : valeurs brutes, pour mesurer sans créer de données (compare_dems)
        before = band.ReadAsArray()
        missing = ~np.isfinite(before) | (before == nodata)
        n_missing = int(missing.sum())
        gdal.FillNodata(band, None, maxSearchDist=FILL_MAX_DISTANCE_PX, smoothingIterations=0,
                         callback=progress_cb or gdal.TermProgress_nocb)
        if n_missing:
            after = band.ReadAsArray()
            filled = missing & np.isfinite(after) & (after != nodata)
            edge = filled & ~closed_gaps(missing)    # remplies par FillNodata sans être dans un trou
            if stats is not None:
                stats["edge_cells"] = edge
            if edge.any() and not extrapolate_edges:
                after[edge] = nodata                  # le relief hors des données n'est pas inventé
                band.WriteArray(after)
                filled &= ~edge
            n_filled = int(filled.sum())
            if n_filled > 0:
                log(f"{n_filled} cellule(s) sans donnée comblée(s) par interpolation pondérée par l'inverse de la "
                    f"distance (gdal.FillNodata, rayon ≤ {FILL_MAX_DISTANCE_PX} pixels) ; "
                    f"{n_missing - n_filled} restent sans donnée")
    elev = band.ReadAsArray().astype(np.float64)
    wgt = warped.GetGeoTransform()
    H, W = elev.shape

    xs = wgt[0] + (np.arange(W) + 0.5) * wgt[1]
    ys = wgt[3] + (np.arange(H) + 0.5) * wgt[5]
    return elev, xs, ys, nodata
