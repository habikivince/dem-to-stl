# SPDX-License-Identifier: GPL-3.0-or-later
"""
polygon.py
----------
Découpe par polygone (frontière, île, commune...) à la place du disque.

Un fichier vectoriel (GeoJSON, GPKG, shapefile... tout ce qu'OGR lit) est chargé, ses
polygones sont fusionnés, puis reprojetés dans le CRS de travail du pipeline. Le masque
d'inclusion est obtenu par remplissage de lignes (règle pair-impair) sur les centres de
cellules : trous et polygones disjoints (archipel) sont gérés, sans dépendance au-delà
de GDAL/OGR et numpy. Les polygones qui traversent l'antiméridien ne sont pas gérés.
"""
import json

import numpy as np
from osgeo import ogr, osr

ogr.UseExceptions()
osr.UseExceptions()

METERS_PER_DEGREE = 111_320.0


class PolygonError(ValueError):
    """Fichier de contour illisible, vide ou sans polygone."""


def _srs(epsg=None, wkt=None):
    s = osr.SpatialReference()
    if epsg is not None:
        s.ImportFromEPSG(epsg)
    else:
        s.ImportFromWkt(wkt)
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return s


def _as_multipolygon(geom):
    """Valide la géométrie et ne garde que des polygones ; None si aucun."""
    if not geom.IsValid():
        try:
            geom = geom.MakeValid()
        except RuntimeError:
            geom = geom.Buffer(0)
    flat = ogr.GT_Flatten(geom.GetGeometryType())
    if flat in (ogr.wkbPolygon, ogr.wkbMultiPolygon, ogr.wkbGeometryCollection):
        multi = ogr.ForceToMultiPolygon(geom)
        if multi is not None and not multi.IsEmpty():
            return multi
    return None


class PolygonClip:
    def __init__(self, geom_wgs84, n_features):
        self.geom_wgs84 = geom_wgs84      # (Multi)Polygon, EPSG:4326, ordre lon/lat
        self.n_features = n_features

    def center_latlon(self):
        lon_min, lon_max, lat_min, lat_max = self.geom_wgs84.GetEnvelope()
        return (lat_min + lat_max) / 2.0, (lon_min + lon_max) / 2.0

    def approx_radius_m(self):
        """Demi-plus grande dimension de l'emprise, en mètres (approximation sphérique)."""
        lon_min, lon_max, lat_min, lat_max = self.geom_wgs84.GetEnvelope()
        lat_c = np.radians((lat_min + lat_max) / 2.0)
        width = (lon_max - lon_min) * METERS_PER_DEGREE * max(np.cos(lat_c), 0.01)
        height = (lat_max - lat_min) * METERS_PER_DEGREE
        return max(width, height) / 2.0

    def to_srs(self, dst_wkt):
        """Copie de la géométrie reprojetée dans le CRS dst_wkt."""
        geom = self.geom_wgs84.Clone()
        geom.Transform(osr.CoordinateTransformation(_srs(epsg=4326), _srs(wkt=dst_wkt)))
        return geom


def from_geojson(geometry):
    """PolygonClip depuis une géométrie GeoJSON (dict) en lon/lat WGS84 (ex. Nominatim)."""
    try:
        geom = ogr.CreateGeometryFromJson(json.dumps(geometry))
    except (RuntimeError, TypeError, ValueError) as e:
        raise PolygonError(f"géométrie GeoJSON invalide ({e})")
    geom = _as_multipolygon(geom) if geom is not None else None
    if geom is None:
        raise PolygonError("la géométrie n'est pas un polygone")
    return PolygonClip(geom, 1)


def load_polygon(path, where=None, log=print):
    try:
        ds = ogr.Open(path)
    except RuntimeError as e:
        raise PolygonError(f"contour illisible : {path} ({e})")
    if ds is None or ds.GetLayerCount() == 0:
        raise PolygonError(f"contour illisible ou sans couche : {path}")
    layer = ds.GetLayer(0)
    if where:
        try:
            layer.SetAttributeFilter(where)
        except RuntimeError as e:
            raise PolygonError(f"filtre --polygon-where invalide ({e})")
    src = layer.GetSpatialRef()
    if src is None:
        log("Contour : aucun CRS déclaré, EPSG:4326 (lon/lat) supposé")
        src = _srs(epsg=4326)
    else:
        src = src.Clone()
        src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    union, n_used, n_skipped = None, 0, 0
    for feat in layer:
        g = feat.GetGeometryRef()
        g = _as_multipolygon(g.Clone()) if g is not None else None
        if g is None:
            n_skipped += 1
            continue
        union = g if union is None else union.Union(g)
        n_used += 1
    if union is None:
        raise PolygonError("aucun polygone dans le contour"
                           + (" après filtre" if where else "") + " (seuls les polygones sont gérés)")
    # densification avant reprojection : les arêtes droites d'un CRS ne le sont pas dans l'autre
    union.Segmentize(0.01 if src.IsGeographic() else 1000.0)
    union.Transform(osr.CoordinateTransformation(src, _srs(epsg=4326)))
    union = _as_multipolygon(union)
    if union is None:
        raise PolygonError("contour vide après reprojection")
    log(f"Contour : {n_used} entité(s) fusionnée(s)" + (f", {n_skipped} ignorée(s) (non polygonales)" if n_skipped else ""))
    return PolygonClip(union, n_used)


def rings_of(geom):
    """Anneaux (tableaux N x 2) de toutes les parties d'un (Multi)Polygon."""
    out = []
    polys = [geom.GetGeometryRef(i) for i in range(geom.GetGeometryCount())] \
        if ogr.GT_Flatten(geom.GetGeometryType()) == ogr.wkbMultiPolygon else [geom]
    for poly in polys:
        for k in range(poly.GetGeometryCount()):
            pts = poly.GetGeometryRef(k).GetPoints()
            if pts:
                out.append(np.asarray(pts, dtype=np.float64)[:, :2])
    return out


def polygon_mask(rings, xs, ys):
    """Masque (len(ys), len(xs)) des centres de cellules intérieurs au polygone (règle
    pair-impair : les trous et polygones disjoints sont gérés). xs croissant, ys dans
    n'importe quel ordre."""
    xs = np.asarray(xs, dtype=np.float64)
    ys = np.asarray(ys, dtype=np.float64)
    H, W = len(ys), len(xs)
    ex0, ey0, ex1, ey1 = [], [], [], []
    for p in rings:
        q = np.roll(p, -1, axis=0)
        ex0.append(p[:, 0]); ey0.append(p[:, 1]); ex1.append(q[:, 0]); ey1.append(q[:, 1])
    if not ex0:
        return np.zeros((H, W), dtype=bool)
    x0, y0, x1, y1 = (np.concatenate(a) for a in (ex0, ey0, ex1, ey1))
    keep = y0 != y1
    x0, y0, x1, y1 = x0[keep], y0[keep], x1[keep], y1[keep]
    swap = y0 > y1
    xa = np.where(swap, x1, x0); ya = np.where(swap, y1, y0)
    xb = np.where(swap, x0, x1); yb = np.where(swap, y0, y1)
    slope = (xb - xa) / (yb - ya)
    diff = np.zeros((H, W + 1), dtype=np.int32)
    for i, y in enumerate(ys):
        sel = (ya <= y) & (y < yb)
        if not sel.any():
            continue
        xi = np.sort(xa[sel] + (y - ya[sel]) * slope[sel])
        n = len(xi) - (len(xi) % 2)
        if n == 0:
            continue
        starts = np.searchsorted(xs, xi[0:n:2], side="left")
        ends = np.searchsorted(xs, xi[1:n:2], side="left")
        np.add.at(diff[i], starts, 1)
        np.add.at(diff[i], ends, -1)
    return np.cumsum(diff[:, :W], axis=1) > 0
