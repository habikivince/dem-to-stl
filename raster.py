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


def warp_window(input_path, dst_srs_wkt, cx, cy, radius, pixel_size_m, progress_cb=None):
    """Ré-échantillonne une fenêtre carrée de 2*radius de côté autour de
    (cx, cy), dans dst_srs_wkt, à pixel_size_m/pixel. Reprojette à la
    volée si dst_srs_wkt diffère du CRS de la source — gdal.Warp gère
    reprojection et découpe en une seule passe, donc même une source
    volumineuse en CRS géographique n'est jamais traitée en entier."""
    ds = gdal.Open(_primary(input_path))
    nodata = ds.GetRasterBand(1).GetNoDataValue()

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
    gdal.FillNodata(band, None, maxSearchDist=50, smoothingIterations=0,
                     callback=progress_cb or gdal.TermProgress_nocb)
    elev = band.ReadAsArray().astype(np.float64)
    wgt = warped.GetGeoTransform()
    H, W = elev.shape

    xs = wgt[0] + (np.arange(W) + 0.5) * wgt[1]
    ys = wgt[3] + (np.arange(H) + 0.5) * wgt[5]
    return elev, xs, ys, nodata
