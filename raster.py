# SPDX-License-Identifier: GPL-3.0-or-later
"""
raster.py
---------
Tout ce qui touche GDAL : ouverture du fichier source, détection et
reprojection automatique du CRS, ré-échantillonnage fenêtré autour du
centre choisi. Aucune logique géométrique ici (voir mesh.py).
"""
import sys

import numpy as np
from osgeo import gdal, osr

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
    ds = gdal.Open(_primary(input_path))
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(ds.GetProjection())
    src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    have_latlon = lat is not None and lon is not None
    have_xy = cx is not None and cy is not None
    if not have_latlon and not have_xy:
        sys.exit("Erreur interne : ni --cx/--cy ni --lat/--lon résolus avant resolve_target_crs_and_center.")

    web_mercator = (not src_srs.IsGeographic()) and _is_web_mercator(src_srs)
    if src_srs.IsGeographic() or web_mercator:
        # Point de référence pour déterminer la zone UTM : --lat/--lon en
        # priorité, sinon --cx/--cy réinterprétés comme lon/lat (ils sont
        # nécessairement dans ce système si la source l'est ; en Web Mercator
        # ils sont d'abord reconvertis en lon/lat).
        if have_latlon:
            ref_lon, ref_lat = lon, lat
        elif web_mercator:
            wm_to_geo = osr.CoordinateTransformation(src_srs, _wgs84())
            ref_lon, ref_lat, _ = wm_to_geo.TransformPoint(cx, cy)
        else:
            ref_lon, ref_lat = cx, cy

        utm_zone = int((ref_lon + 180) / 6) + 1
        dst_epsg = (32600 if ref_lat >= 0 else 32700) + utm_zone
        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(dst_epsg)
        dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        kind = "Web Mercator (EPSG:3857, mètres non métriques)" if web_mercator else "géographique"
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


def warp_window(input_path, dst_srs_wkt, cx, cy, radius, pixel_size_m, progress_cb=None, fill=True,
                inside_fn=None):
    """Ré-échantillonne une fenêtre carrée de 2*radius de côté autour de
    (cx, cy), dans dst_srs_wkt, à pixel_size_m/pixel. Reprojette à la
    volée si dst_srs_wkt diffère du CRS de la source — gdal.Warp gère
    reprojection et découpe en une seule passe, donc même une source
    volumineuse en CRS géographique n'est jamais traitée en entier.

    Si input_path porte un attribut `fill_loader` (voir sources.with_copernicus_fill),
    les trous du disque sont comblés par cette source secondaire, recalée en altitude.
    fill=False : valeurs brutes, sans comblement ni complément (compare_dems).
    inside_fn(xs, ys) -> masque : zone à compléter si elle n'est pas le disque (mode polygone)."""
    ds = gdal.Open(_primary(input_path))
    nodata = ds.GetRasterBand(1).GetNoDataValue()
    loader = getattr(input_path, "fill_loader", None) if fill else None
    if loader is not None and nodata is None:
        nodata = FILL_NODATA

    xmin, xmax = cx - radius, cx + radius
    ymin, ymax = cy - radius, cy + radius

    warped = gdal.Warp(
        "", input_path, format="MEM",
        dstSRS=dst_srs_wkt,
        outputBounds=(xmin, ymin, xmax, ymax),
        xRes=pixel_size_m, yRes=pixel_size_m,
        resampleAlg="average",
        dstNodata=nodata,
        callback=progress_cb or gdal.TermProgress_nocb,
    )
    band = warped.GetRasterBand(1)
    if loader is not None:
        _complete_with_secondary(warped, loader, dst_srs_wkt, (xmin, ymin, xmax, ymax),
                                 pixel_size_m, nodata, cx, cy, radius, inside_fn)
    if fill:  # fill=False : valeurs brutes, pour mesurer sans créer de données (compare_dems)
        gdal.FillNodata(band, None, maxSearchDist=50, smoothingIterations=0,
                         callback=progress_cb or gdal.TermProgress_nocb)
    elev = band.ReadAsArray().astype(np.float64)
    wgt = warped.GetGeoTransform()
    H, W = elev.shape

    xs = wgt[0] + (np.arange(W) + 0.5) * wgt[1]
    ys = wgt[3] + (np.arange(H) + 0.5) * wgt[5]
    return elev, xs, ys, nodata
