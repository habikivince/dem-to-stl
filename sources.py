# SPDX-License-Identifier: GPL-3.0-or-later
"""
sources.py
----------
Sources DEM nationales à haute résolution, récupérées automatiquement
pour la zone demandée (cercle centre/rayon). Chaque source expose une
fonction `fetch_<nom>(lat, lon, radius_m, cache_dir, pixel_size_m, log)`
qui renvoie une LISTE de rasters lisibles par GDAL (chemins locaux ou
/vsicurl/...). Aucune mosaïque ni reprojection ici : raster.warp_window
(gdal.Warp) accepte une liste et gère les CRS différents en une passe.

`pixel_size_m` est la résolution de travail du pipeline (pas du maillage
ramené au terrain) : chaque source choisit la résolution native la plus
proche plutôt que de télécharger du 1 m inutile.

État de vérification (octobre 2026) :
  - swisstopo : TESTÉ EN RÉEL (Cervin) : STAC, tuiles 2 m EPSG:2056.
  - Kartverket : GetCapabilities lu en réel (formats image/GeoTIFF, WCS 1.0.0 à
    1.1.2) ; syntaxe GetCoverage non confirmée, variantes essayées au runtime.
  - USGS : écrit de mémoire ; TESTÉ EN RÉEL pour le jeu 1/3 arc-seconde
    (« National Elevation Dataset (NED) 1/3 arc-second Current »). Le jeu
    1 m (utilisé si le pixel de travail < 5 m) reste non vérifié.
  - GSI : TESTÉ EN RÉEL (Mont Fuji) : la couche `dem1a_png` répond bien.
  - IGN : WMS-R Géoplateforme (annexe officielle lue), couche LiDAR HD brute
    non confirmée sur serveur ; repli RGE ALTI.
"""
import json
import math
import os
import shutil
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from dataclasses import dataclass
from typing import Callable, List, Tuple

import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()
osr.UseExceptions()

USER_AGENT = os.environ.get("DEM2STL_USER_AGENT", "dem-to-stl/1.0")
HTTP_TIMEOUT = 30
HTTP_RETRIES = 3
FULL_RES_PIXEL_M = 0.1   # « pixel » factice plus fin que toute source : chaque source prend son maximum
USGS_1M_BELOW_PIXEL_M = 12.0
NODATA = -9999.0


class SourceUnavailable(RuntimeError):
    """Source injoignable ou en erreur (réseau, service, réponse invalide)."""


class NoCoverage(SourceUnavailable):
    """La source répond mais n'a pas de données pour cette zone."""


# ---------------------------------------------------------------------------
# Utilitaires réseau / géométrie
# ---------------------------------------------------------------------------

def _http_get(url, timeout=None, retries=None):
    """Renvoie les octets de la réponse, ou None pour un 404. Réessaie avec
    attente exponentielle sur erreurs transitoires ; lève SourceUnavailable
    sinon (réseau, timeout, 4xx/5xx)."""
    timeout = HTTP_TIMEOUT if timeout is None else timeout
    retries = HTTP_RETRIES if retries is None else retries
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            try:
                body = e.read()[:200].decode("utf-8", "replace").replace("\n", " ")
            except Exception:
                body = ""
            last = f"{e}{' — ' + body if body else ''}"
            if e.code not in (429, 500, 502, 503, 504):
                break
        except OSError as e:  # URLError, timeout, connexion coupée
            last = e
        if attempt < retries - 1:
            time.sleep(min(2 ** attempt, 8))
    raise SourceUnavailable(f"{url} : {last}")


def _http_json(url, params=None):
    if params:
        url = f"{url}{'&' if '?' in url else '?'}{urllib.parse.urlencode(params)}"
    data = _http_get(url)
    if data is None:
        raise SourceUnavailable(f"{url} : HTTP 404")
    try:
        return json.loads(data.decode("utf-8"))
    except ValueError as e:
        raise SourceUnavailable(f"{url} : réponse JSON invalide ({e})")


def _download(url, dest):
    """Télécharge vers dest (écriture atomique) sauf si déjà en cache."""
    if os.path.exists(dest):
        return dest
    data = _http_get(url)
    if data is None:
        raise SourceUnavailable(f"{url} : HTTP 404")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    with open(tmp, "wb") as f:
        f.write(data)
    os.replace(tmp, dest)
    return dest


def circle_bbox_wgs84(lat, lon, radius_m, margin_m=2000.0):
    """(lon_min, lat_min, lon_max, lat_max) englobant le cercle + marge."""
    r = radius_m + margin_m
    dlat = r / 111_320.0
    dlon = r / (111_320.0 * max(math.cos(math.radians(lat)), 0.01))
    return (lon - dlon, lat - dlat, lon + dlon, lat + dlat)


def _srs(epsg):
    s = osr.SpatialReference()
    s.ImportFromEPSG(epsg)
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return s


def _to_epsg(epsg, lat, lon):
    ct = osr.CoordinateTransformation(_srs(4326), _srs(epsg))
    x, y, _ = ct.TransformPoint(lon, lat)
    return x, y


# ---------------------------------------------------------------------------
# swisstopo — swissALTI3D (STAC)
# ---------------------------------------------------------------------------

SWISS_STAC_ITEMS = "https://data.geo.admin.ch/api/stac/v0.9/collections/ch.swisstopo.swissalti3d/items"


def swisstopo_gsd(pixel_size_m):
    return 2.0 if pixel_size_m >= 2.0 else 0.5


def _stac_items(url, params, max_pages=50):
    items = []
    next_url = f"{url}?{urllib.parse.urlencode(params)}"
    for _ in range(max_pages):
        page = _http_json(next_url)
        items.extend(page.get("features", []))
        nxt = [l.get("href") for l in page.get("links", []) if l.get("rel") == "next"]
        if not nxt or not nxt[0]:
            return items
        next_url = nxt[0]
    raise SourceUnavailable("STAC swisstopo : trop de pages de résultats")


def select_swisstopo_assets(items, gsd):
    """Une seule URL par tuile (version la plus récente), au pas le plus
    proche de `gsd`. L'id d'item est de la forme swissalti3d_<année>_<tuile>."""
    best = {}
    for it in items:
        parts = str(it.get("id", "")).split("_", 2)
        tile = parts[2] if len(parts) == 3 else str(it.get("id"))
        stamp = (it.get("properties") or {}).get("datetime") or (parts[1] if len(parts) == 3 else "")
        cands = []
        for a in (it.get("assets") or {}).values():
            href = a.get("href", "")
            if not href.lower().endswith(".tif"):
                continue
            a_gsd = a.get("eo:gsd")
            cands.append((abs(a_gsd - gsd) if a_gsd is not None else 99.0, href))
        if not cands:
            continue
        href = min(cands)[1]
        if tile not in best or stamp > best[tile][0]:
            best[tile] = (stamp, href)
    return sorted(h for _, h in best.values())


def fetch_swisstopo(lat, lon, radius_m, cache_dir, pixel_size_m, log=print):
    bbox = circle_bbox_wgs84(lat, lon, radius_m)
    gsd = swisstopo_gsd(pixel_size_m)
    items = _stac_items(SWISS_STAC_ITEMS, {"bbox": ",".join(f"{v:.6f}" for v in bbox), "limit": 100})
    urls = select_swisstopo_assets(items, gsd)
    if not urls:
        raise NoCoverage("swissALTI3D : aucune tuile pour cette zone")
    log(f"swissALTI3D : {len(urls)} tuile(s), pas demandé {gsd} m (pixel de travail {pixel_size_m:.2f} m)")
    out = []
    for u in urls:
        dest = os.path.join(cache_dir, "swisstopo", os.path.basename(urllib.parse.urlparse(u).path))
        log(f"  {os.path.basename(dest)} : {'en cache' if os.path.exists(dest) else 'téléchargement...'}")
        out.append(_download(u, dest))
    return out


# ---------------------------------------------------------------------------
# Kartverket — Høydedata DTM (WCS ArcGIS ImageServer)
# ---------------------------------------------------------------------------
# Constats (GetCapabilities lu en réel) : serveur ArcGIS ImageServer ; formats
# acceptés image/GeoTIFF, image/HDF, image/JPEG2000, image/NetCDF — PAS image/tiff,
# ce qui explique les HTTP 400 obtenus avec le « nativeFormat » de DescribeCoverage
# (requêtes GDAL et requêtes 2.0.1 précédentes). Versions de GetCoverage annoncées :
# 1.0.0 à 1.1.2. Une application tierce utilise WCS 1.0.0 avec format=GeoTIFF sur ce
# service. La syntaxe exacte n'étant pas confirmée, plusieurs variantes sont essayées
# par la sonde du centre ; la première qui renvoie un TIFF de la bonne taille est
# conservée et indiquée dans les logs.

KARTVERKET_WCS_DEFAULT = "https://wms.geonorge.no/skwms1/wcs.hoyde-dtm-nhm-topobathy-25833"
KARTVERKET_COVERAGE_DEFAULT = "nhm_dtm_topobathy_25833"
KARTVERKET_EPSG = 25833
KARTVERKET_NATIVE_M = 1.0
KARTVERKET_CHUNK_PX = 2000      # côté max (en pixels de sortie) d'une requête
KARTVERKET_MAX_CHUNKS = 400
KARTVERKET_VARIANTS = ["v100_geotiff", "v100_image", "v201_image"]
TIFF_MAGICS = (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")


def kartverket_bounds(lat, lon, radius_m, res, margin_m=None):
    margin_m = max(50.0, 4 * res) if margin_m is None else margin_m
    x, y = _to_epsg(KARTVERKET_EPSG, lat, lon)
    half = radius_m + margin_m
    return (math.floor((x - half) / res) * res, math.floor((y - half) / res) * res,
            math.ceil((x + half) / res) * res, math.ceil((y + half) / res) * res)


def _fmt(v):
    return f"{v:.10g}"


def wcs_request_url(base, coverage, variant, b, size):
    """URL GetCoverage. b = (xmin, ymin, xmax, ymax) en EPSG:25833 ; size = (largeur,
    hauteur) en pixels de sortie. Variantes : v100_geotiff / v100_image (WCS 1.0.0,
    FORMAT=GeoTIFF ou image/GeoTIFF), v201_image (WCS 2.0.1, SUBSET + SCALESIZE)."""
    w, h = size
    if variant in ("v100_geotiff", "v100_image"):
        parts = [("SERVICE", "WCS"), ("VERSION", "1.0.0"), ("REQUEST", "GetCoverage"),
                 ("COVERAGE", coverage), ("CRS", f"EPSG:{KARTVERKET_EPSG}"),
                 ("BBOX", ",".join(_fmt(v) for v in b)), ("WIDTH", str(w)), ("HEIGHT", str(h)),
                 ("FORMAT", "GeoTIFF" if variant == "v100_geotiff" else "image/GeoTIFF")]
    elif variant == "v201_image":
        parts = [("SERVICE", "WCS"), ("VERSION", "2.0.1"), ("REQUEST", "GetCoverage"),
                 ("COVERAGEID", coverage),
                 ("SUBSET", f"x({_fmt(b[0])},{_fmt(b[2])})"),
                 ("SUBSET", f"y({_fmt(b[1])},{_fmt(b[3])})"),
                 ("SCALESIZE", f"x({w}),y({h})"), ("FORMAT", "image/GeoTIFF")]
    else:
        raise ValueError(variant)
    query = "&".join(f"{k}={urllib.parse.quote(v, safe='(),:/')}" for k, v in parts)
    return f"{base.split('?')[0]}?{query}"


class _RasterService:
    """Client HTTP d'un service raster (WCS/WMS) dont la syntaxe exacte n'est pas
    connue d'avance : `build_url(variante, fenêtre, taille)` est essayé pour chaque
    variante ; la première qui renvoie un GeoTIFF de la bonne taille, en valeurs
    flottantes (donc de vraies altitudes, pas une image ombrée), est conservée."""

    def __init__(self, label, build_url, variants, epsg, work_dir, log):
        self.label, self.build_url, self.epsg = label, build_url, epsg
        self.work_dir, self.log = work_dir, log
        self.variants = list(variants)
        self.variant = None
        os.makedirs(work_dir, exist_ok=True)

    def get(self, b, size):
        """Chemin d'un GeoTIFF local (taille exacte `size`) pour la fenêtre b."""
        errors = []
        for variant in list(self.variants):
            url = self.build_url(variant, b, size)
            try:
                data = _http_get(url)
                if data is None:
                    raise SourceUnavailable(f"{url} : HTTP 404")
                path = self._save(data, b, size, url)
            except SourceUnavailable as e:
                errors.append(f"[{variant}] {e}")
                continue
            if self.variant != variant:
                self.variant = variant
                self.variants.remove(variant)
                self.variants.insert(0, variant)
                self.log(f"{self.label} : variante {variant} acceptée par le serveur")
            return path
        raise SourceUnavailable(f"{self.label} : aucune variante de requête acceptée.\n  "
                                + "\n  ".join(errors))

    def _save(self, data, b, size, url):
        if data[:4] not in TIFF_MAGICS:
            snippet = data[:300].decode("utf-8", "replace").replace("\n", " ")
            raise SourceUnavailable(f"réponse non TIFF pour {url} : {snippet}")
        path = os.path.join(self.work_dir, f"c_{uuid.uuid4().hex}.tif")
        with open(path, "wb") as f:
            f.write(data)
        try:
            ds = gdal.Open(path)
        except RuntimeError as e:
            os.remove(path)
            raise SourceUnavailable(f"TIFF illisible pour {url} : {e}")
        dims = (ds.RasterXSize, ds.RasterYSize)
        dtype = ds.GetRasterBand(1).DataType
        georef = bool(ds.GetProjection()) and ds.GetGeoTransform() != (0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
        ds = None
        if dims != tuple(size):
            os.remove(path)
            raise SourceUnavailable(f"taille reçue {dims} différente de {tuple(size)} pour {url}")
        if dtype not in (gdal.GDT_Float32, gdal.GDT_Float64):
            os.remove(path)
            raise SourceUnavailable(f"valeurs de type {gdal.GetDataTypeName(dtype)} et non flottantes "
                                    f"(image rendue, pas des altitudes ?) pour {url}")
        if not georef:  # TIFF sans géoréférencement : on le déduit de la fenêtre demandée
            geo = path + ".geo.tif"
            gdal.Translate(geo, path, outputSRS=f"EPSG:{self.epsg}",
                           outputBounds=[b[0], b[3], b[2], b[1]])
            os.replace(geo, path)
        # valeurs sentinelles (-99999, -3.4e38, NaN...) -> nodata explicite
        ds = gdal.Open(path, gdal.GA_Update)
        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray().astype(np.float32)
        bad = ~np.isfinite(arr) | (arr < -9000) | (arr > 1e5)
        if bad.any():
            arr[bad] = NODATA
            band.WriteArray(arr)
        band.SetNoDataValue(NODATA)
        ds = None
        return path


def _valid_stats(path):
    ds = gdal.Open(path)
    arr = ds.GetRasterBand(1).ReadAsArray()
    valid = arr != NODATA
    return bool(valid.any()), bool(valid.any() and np.all(arr[valid] == 0))


def _kartverket_chunks(nb, out_res, max_px):
    """Découpe nb (entiers, mètres) en fenêtres d'au plus max_px pixels de sortie ;
    renvoie [(fenêtre, (largeur_px, hauteur_px))]."""
    step = max(1, int(max_px * out_res))
    xs = list(range(int(nb[0]), int(nb[2]), step)) or [int(nb[0])]
    ys = list(range(int(nb[1]), int(nb[3]), step)) or [int(nb[1])]
    out = []
    for y0 in ys:
        for x0 in xs:
            w = (x0, y0, min(x0 + step, int(nb[2])), min(y0 + step, int(nb[3])))
            out.append((w, (max(1, round((w[2] - w[0]) / out_res)),
                            max(1, round((w[3] - w[1]) / out_res)))))
    return out


def fetch_kartverket(lat, lon, radius_m, cache_dir, pixel_size_m, log=print):
    base = os.environ.get("DEM2STL_KARTVERKET_WCS", KARTVERKET_WCS_DEFAULT)
    coverage = os.environ.get("DEM2STL_KARTVERKET_COVERAGE", KARTVERKET_COVERAGE_DEFAULT)
    res = max(KARTVERKET_NATIVE_M, float(pixel_size_m))
    # le serveur ré-échantillonne au plus proche : on demande 2x plus fin que la
    # résolution de travail puis on moyenne (warp 'average') pour limiter l'aliasing
    out_res = max(KARTVERKET_NATIVE_M, res / 2.0)
    b = kartverket_bounds(lat, lon, radius_m, res)
    kdir = os.path.join(cache_dir, "kartverket")
    dest = os.path.join(kdir, f"dtm_{res:g}m_{int(b[0])}_{int(b[1])}_{int(b[2])}_{int(b[3])}.tif")
    if os.path.exists(dest):
        log(f"Kartverket : {os.path.basename(dest)} déjà en cache")
        return [dest]
    work = os.path.join(kdir, f"_tmp_{uuid.uuid4().hex[:8]}")
    wcs = _RasterService("Kartverket WCS",
                         lambda v, bb, sz: wcs_request_url(base, coverage, v, bb, sz),
                         KARTVERKET_VARIANTS, KARTVERKET_EPSG, work, log)
    try:
        # 1. Sonde de 400 m au centre : choisit la variante de requête, vérifie les données
        cx, cy = _to_epsg(KARTVERKET_EPSG, lat, lon)
        x0, y0 = math.floor(cx) - 200, math.floor(cy) - 200
        n = max(2, round(400 / out_res))
        probe = wcs.get((x0, y0, x0 + 400, y0 + 400), (n, n))
        has_data, all_zero = _valid_stats(probe)
        if not has_data:
            raise NoCoverage("Kartverket : aucune donnée au centre de la zone")
        if all_zero:
            log("Kartverket : ATTENTION, altitude nulle partout au centre (mer ou hors couverture ?)")
        # 2. Fenêtres
        nb = (math.floor(b[0]), math.floor(b[1]), math.ceil(b[2]), math.ceil(b[3]))
        chunks = _kartverket_chunks(nb, out_res, KARTVERKET_CHUNK_PX)
        if len(chunks) > KARTVERKET_MAX_CHUNKS:
            raise SourceUnavailable(f"Kartverket : {len(chunks)} requêtes nécessaires "
                                    f"(max {KARTVERKET_MAX_CHUNKS}) ; réduis --radius ou utilise --input")
        log(f"Kartverket DTM : {len(chunks)} requête(s) WCS à {out_res:g} m/pixel, "
            f"moyennées à {res:g} m/pixel (EPSG:{KARTVERKET_EPSG})")
        paths = [wcs.get(w, size) for w, size in chunks]
        tmp = dest + ".part"
        os.makedirs(kdir, exist_ok=True)
        gdal.Warp(tmp, paths, format="GTiff", dstSRS=f"EPSG:{KARTVERKET_EPSG}", outputBounds=b,
                  xRes=res, yRes=res, resampleAlg="average", srcNodata=NODATA, dstNodata=NODATA,
                  creationOptions=["COMPRESS=DEFLATE", "TILED=YES"])
        os.replace(tmp, dest)
    except SourceUnavailable:
        raise
    except RuntimeError as e:
        raise SourceUnavailable(f"Kartverket WCS : {e}")
    finally:
        shutil.rmtree(work, ignore_errors=True)
        if os.path.exists(dest + ".part"):
            os.remove(dest + ".part")
    return [dest]


# ---------------------------------------------------------------------------
# IGN — MNT LiDAR HD (WMS-R Géoplateforme), repli RGE ALTI
# ---------------------------------------------------------------------------
# Constats (annexe officielle altimetrie.xml, lue en réel) : WMS 1.3.0 sans clé sur
# data.geopf.fr/wms-r/wms, GetMap en image/geotiff | image/tiff | image/x-bil;bits=32,
# 5010 px max par côté, CRS dont EPSG:2154. Les couches « ...SHADOW » sont des
# estompages ; les couches ...ELEVATIONGRIDCOVERAGE.LAMB93 / .RGR92UTM40S (MNT dérivé
# du LiDAR HD, CRS natif) sont citées par le GetCapabilities complet. La doc IGN
# indique « TIFF + style normal pour des valeurs brutes » (couche MNS). Non confirmé
# sur serveur : valeurs réellement brutes, nodata hors couverture. Le client refuse
# donc toute réponse non flottante (image rendue) et essaie plusieurs variantes.

IGN_WMS_DEFAULT = "https://data.geopf.fr/wms-r/wms"
IGN_LIDAR_LAYERS = {
    "fxx": (2154, "IGNF_LIDAR-HD_MNT_ELEVATION.ELEVATIONGRIDCOVERAGE.LAMB93"),
    "reu": (2975, "IGNF_LIDAR-HD_MNT_ELEVATION.ELEVATIONGRIDCOVERAGE.RGR92UTM40S"),
}
IGN_REGION_BOXES = {"fxx": (-5.6, 41.0, 10.0, 51.6), "reu": (55.1, -21.5, 55.9, -20.8)}
IGN_RGEALTI_LAYER = "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES"
IGN_VARIANTS = ["geotiff_normal", "geotiff_default", "tiff_normal"]
IGN_CHUNK_PX = 4000
IGN_MAX_CHUNKS = 400
IGN_MIN_COVERAGE = 0.98


def ign_getmap_url(base, layer, variant, epsg, b, size):
    fmt = {"geotiff": "image/geotiff", "tiff": "image/tiff"}[variant.split("_")[0]]
    parts = [("SERVICE", "WMS"), ("VERSION", "1.3.0"), ("REQUEST", "GetMap"), ("LAYERS", layer),
             ("STYLES", "normal" if variant.endswith("normal") else ""),
             ("CRS", f"EPSG:{epsg}"), ("BBOX", ",".join(_fmt(v) for v in b)),
             ("WIDTH", str(size[0])), ("HEIGHT", str(size[1])), ("FORMAT", fmt)]
    query = "&".join(f"{k}={urllib.parse.quote(v, safe=',:/')}" for k, v in parts)
    return f"{base.split('?')[0]}?{query}"


def _aligned_window(epsg, lat, lon, radius_m, res, margin_m):
    x, y = _to_epsg(epsg, lat, lon)
    half = radius_m + margin_m
    return (math.floor((x - half) / res) * res, math.floor((y - half) / res) * res,
            math.ceil((x + half) / res) * res, math.ceil((y + half) / res) * res), (x, y)


def _disc_coverage(path, center, radius_m):
    ds = gdal.Open(path)
    gt = ds.GetGeoTransform()
    arr = ds.GetRasterBand(1).ReadAsArray()
    xs = gt[0] + (np.arange(arr.shape[1]) + 0.5) * gt[1]
    ys = gt[3] + (np.arange(arr.shape[0]) + 0.5) * gt[5]
    inside = (xs[None, :] - center[0]) ** 2 + (ys[:, None] - center[1]) ** 2 <= radius_m ** 2
    return float(np.mean(arr[inside] != NODATA)) if inside.any() else 0.0


def fetch_ign(lat, lon, radius_m, cache_dir, pixel_size_m, log=print):
    base = os.environ.get("DEM2STL_IGN_WMS", IGN_WMS_DEFAULT)
    region = next((k for k, bx in IGN_REGION_BOXES.items()
                   if bx[0] <= lon <= bx[2] and bx[1] <= lat <= bx[3]), None)
    if region is None:
        raise NoCoverage("IGN : point hors France métropolitaine / La Réunion")
    epsg, lidar_layer = IGN_LIDAR_LAYERS[region]
    res = max(0.5, float(pixel_size_m))
    kdir = os.path.join(cache_dir, "ign")
    last_error = None
    for tag, layer, native in (("lidarhd", lidar_layer, 0.5), ("rgealti", IGN_RGEALTI_LAYER, 1.0)):
        res_l = max(native, res)
        out_res = max(native, res_l / 2.0)
        b, center = _aligned_window(epsg, lat, lon, radius_m, res_l, max(50.0, 4 * res_l))
        dest = os.path.join(kdir, f"{tag}_{region}_{res_l:g}m_{int(b[0])}_{int(b[1])}_{int(b[2])}_{int(b[3])}.tif")
        if os.path.exists(dest):
            log(f"IGN : {os.path.basename(dest)} déjà en cache")
            return [dest]
        work = os.path.join(kdir, f"_tmp_{uuid.uuid4().hex[:8]}")
        svc = _RasterService(f"IGN {tag}", lambda v, bb, sz, _l=layer: ign_getmap_url(base, _l, v, epsg, bb, sz),
                             IGN_VARIANTS, epsg, work, log)
        try:
            n = max(2, round(400 / out_res))
            x0, y0 = math.floor(center[0]) - 200, math.floor(center[1]) - 200
            probe = svc.get((x0, y0, x0 + 400, y0 + 400), (n, n))
            if not _valid_stats(probe)[0]:
                last_error = NoCoverage(f"IGN {tag} : aucune donnée au centre de la zone")
                log(f"IGN : pas de {tag} au centre de la zone" + (" -> repli RGE ALTI" if tag == "lidarhd" else ""))
                continue
            nb = (math.floor(b[0]), math.floor(b[1]), math.ceil(b[2]), math.ceil(b[3]))
            chunks = _kartverket_chunks(nb, out_res, IGN_CHUNK_PX)
            if len(chunks) > IGN_MAX_CHUNKS:
                raise SourceUnavailable(f"IGN : {len(chunks)} requêtes nécessaires (max {IGN_MAX_CHUNKS}) ; "
                                        f"réduis --radius ou utilise --input")
            log(f"IGN {tag} ({layer}) : {len(chunks)} requête(s) WMS à {out_res:g} m/pixel, "
                f"moyennées à {res_l:g} m/pixel (EPSG:{epsg})")
            paths = [svc.get(w, size) for w, size in chunks]
            tmp = dest + ".part"
            os.makedirs(kdir, exist_ok=True)
            gdal.Warp(tmp, paths, format="GTiff", dstSRS=f"EPSG:{epsg}", outputBounds=b,
                      xRes=res_l, yRes=res_l, resampleAlg="average", srcNodata=NODATA, dstNodata=NODATA,
                      creationOptions=["COMPRESS=DEFLATE", "TILED=YES"])
            os.replace(tmp, dest)
            cov = _disc_coverage(dest, center, radius_m)
            if cov < IGN_MIN_COVERAGE:
                log(f"IGN : ATTENTION, {100 * (1 - cov):.1f} % du disque sans donnée dans {tag} "
                    f"(couverture LiDAR HD incomplète ?) ; vérifie le résultat ou complète avec une autre source")
            return [dest]
        except SourceUnavailable as e:
            last_error = e
            log(f"IGN : {tag} inutilisable ({str(e).splitlines()[0][:200]})" + (" -> repli RGE ALTI" if tag == "lidarhd" else ""))
        except RuntimeError as e:
            last_error = SourceUnavailable(f"IGN {tag} : {e}")
        finally:
            shutil.rmtree(work, ignore_errors=True)
            if os.path.exists(dest + ".part"):
                os.remove(dest + ".part")
    raise last_error or SourceUnavailable("IGN : aucune couche exploitable")


# ---------------------------------------------------------------------------
# USGS 3DEP — TNM Access API (non vérifiée)
# ---------------------------------------------------------------------------

TNM_PRODUCTS = "https://tnmaccess.nationalmap.gov/api/v1/products"
USGS_DATASETS_1M = ["Digital Elevation Model (DEM) 1 meter"]
USGS_DATASETS_13 = ["National Elevation Dataset (NED) 1/3 arc-second Current",
                    "National Elevation Dataset (NED) 1/3 arc-second - Current"]
USGS_MAX_PRODUCTS = 400


def _tnm_query(datasets, bbox):
    for name in datasets:
        data = _http_json(TNM_PRODUCTS, {
            "datasets": name, "bbox": ",".join(f"{v:.6f}" for v in bbox),
            "prodFormats": "GeoTIFF", "outputFormat": "JSON", "max": USGS_MAX_PRODUCTS})
        items = [i for i in data.get("items", []) if str(i.get("downloadURL", "")).lower().endswith(".tif")]
        if items:
            return name, items, int(data.get("total", len(items)))
    return None, [], 0


def fetch_usgs(lat, lon, radius_m, cache_dir, pixel_size_m, log=print):
    bbox = circle_bbox_wgs84(lat, lon, radius_m)
    attempts = []
    if pixel_size_m < USGS_1M_BELOW_PIXEL_M:
        attempts.append(USGS_DATASETS_1M)
    attempts.append(USGS_DATASETS_13)
    for datasets in attempts:
        name, items, total = _tnm_query(datasets, bbox)
        if not items:
            log(f"USGS : aucun produit « {datasets[0]} » pour cette zone")
            continue
        if total > USGS_MAX_PRODUCTS:
            log(f"USGS : {total} produits « {name} » (> {USGS_MAX_PRODUCTS}), produit trop fragmenté — essai plus grossier")
            continue
        # gdal.Warp : en cas de recouvrement la dernière source gagne -> plus récent en dernier
        items.sort(key=lambda i: str(i.get("publicationDate") or i.get("lastUpdated") or ""))
        log(f"USGS 3DEP « {name} » : {len(items)} produit(s), lecture partielle via /vsicurl/")
        return ["/vsicurl/" + i["downloadURL"] for i in items]
    raise NoCoverage("USGS 3DEP : aucun produit exploitable pour cette zone")


# ---------------------------------------------------------------------------
# GSI Japon — tuiles d'élévation PNG (Web Mercator)
# ---------------------------------------------------------------------------

GSI_URL = "https://cyberjapandata.gsi.go.jp/xyz/{layer}/{z}/{x}/{y}.png"
# (couche, zoom max) par ordre de qualité décroissante. Zoom max : page officielle
# « standard des tuiles d'élévation » (DEM1A 17, DEM5A/5B/5C 15, DEM10B 14).
GSI_LAYERS = [("dem1a_png", 17), ("dem5a_png", 15), ("dem5b_png", 15), ("dem_png", 14)]
GSI_MAX_TILES = 600
WM_ORIGIN = 20037508.342789244


def decode_gsi_rgb(rgb):
    """Décodage officiel : x = 2^16 R + 2^8 G + B ; x < 2^23 -> x*0.01 ;
    x == 2^23 -> invalide ; x > 2^23 -> (x - 2^24)*0.01. rgb : (3, H, W) uint8."""
    r, g, b = (rgb[i].astype(np.int64) for i in range(3))
    x = r * 65536 + g * 256 + b
    h = np.where(x < 2 ** 23, x, x - 2 ** 24) * 0.01
    h = np.where(x == 2 ** 23, np.nan, h)
    return h.astype(np.float32)


def gsi_zoom_for_pixel(lat, pixel_size_m, zmax=17):
    ground0 = (2 * WM_ORIGIN / 256) * math.cos(math.radians(lat))
    z = math.ceil(math.log2(ground0 / pixel_size_m)) if pixel_size_m > 0 else zmax
    return max(0, min(zmax, z))


def gsi_tile_xy(lat, lon, z):
    n = 2 ** z
    x = int(math.floor((lon + 180.0) / 360.0 * n))
    y = int(math.floor((1.0 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2.0 * n))
    return min(max(x, 0), n - 1), min(max(y, 0), n - 1)


def gsi_tile_range(bbox, z):
    x0, y0 = gsi_tile_xy(bbox[3], bbox[0], z)  # coin NW
    x1, y1 = gsi_tile_xy(bbox[1], bbox[2], z)  # coin SE
    return x0, y0, x1, y1


def _read_png_rgb(data):
    name = f"/vsimem/gsi_{uuid.uuid4().hex}.png"
    gdal.FileFromMemBuffer(name, data)
    try:
        ds = gdal.Open(name)
        return ds.ReadAsArray()[:3]
    finally:
        gdal.Unlink(name)


def _gsi_tile(z, x, y, cache_dir):
    """(altitudes 256x256, couche) ou (None, None) si aucune couche n'a cette tuile."""
    for layer, zmax in GSI_LAYERS:
        if z > zmax:
            continue
        path = os.path.join(cache_dir, "gsi", layer, str(z), f"{x}_{y}.png")
        if os.path.exists(path):
            with open(path, "rb") as f:
                data = f.read()
        else:
            data = _http_get(GSI_URL.format(layer=layer, z=z, x=x, y=y))
            if data is None:
                continue
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path + ".part", "wb") as f:
                f.write(data)
            os.replace(path + ".part", path)
        try:
            return decode_gsi_rgb(_read_png_rgb(data)), layer
        except RuntimeError:
            raise SourceUnavailable(f"GSI : tuile {layer}/{z}/{x}/{y} illisible")
    return None, None


def fetch_gsi(lat, lon, radius_m, cache_dir, pixel_size_m, log=print):
    bbox = circle_bbox_wgs84(lat, lon, radius_m)
    z = gsi_zoom_for_pixel(lat, pixel_size_m)
    x0, y0, x1, y1 = gsi_tile_range(bbox, z)
    while (x1 - x0 + 1) * (y1 - y0 + 1) > GSI_MAX_TILES and z > 0:
        z -= 1
        x0, y0, x1, y1 = gsi_tile_range(bbox, z)
        log(f"GSI : trop de tuiles, zoom abaissé à {z}")
    nx, ny = x1 - x0 + 1, y1 - y0 + 1
    log(f"GSI : zoom {z}, {nx * ny} tuile(s) (pixel de travail {pixel_size_m:.2f} m)")

    mosaic = np.full((ny * 256, nx * 256), NODATA, dtype=np.float32)
    used, missing = {}, 0
    for j in range(ny):
        for i in range(nx):
            h, layer = _gsi_tile(z, x0 + i, y0 + j, cache_dir)
            if h is None:
                missing += 1
                continue
            used[layer] = used.get(layer, 0) + 1
            mosaic[j * 256:(j + 1) * 256, i * 256:(i + 1) * 256] = np.where(np.isnan(h), NODATA, h)
    if not used:
        raise NoCoverage("GSI : aucune tuile d'élévation pour cette zone")
    log("GSI : couches utilisées " + ", ".join(f"{k}={v}" for k, v in used.items())
        + (f" ; {missing} tuile(s) absente(s)" if missing else ""))

    tile_m = 2 * WM_ORIGIN / 2 ** z
    px = tile_m / 256
    dest = os.path.join(cache_dir, "gsi", f"mosaic_z{z}_{x0}_{y0}_{nx}x{ny}.tif")
    ds = gdal.GetDriverByName("GTiff").Create(dest, nx * 256, ny * 256, 1, gdal.GDT_Float32,
                                              ["COMPRESS=DEFLATE", "TILED=YES"])
    ds.SetGeoTransform((-WM_ORIGIN + x0 * tile_m, px, 0, WM_ORIGIN - y0 * tile_m, 0, -px))
    ds.SetProjection(_srs(3857).ExportToWkt())
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(NODATA)
    band.WriteArray(mosaic)
    ds = None
    return [dest]


# ---------------------------------------------------------------------------
# Registre, sélection automatique
# ---------------------------------------------------------------------------

@dataclass
class Source:
    name: str
    label: str
    boxes: List[Tuple[float, float, float, float]]  # emprises grossières (lon_min, lat_min, lon_max, lat_max)
    fetch: Callable

    def covers(self, lat, lon):
        return any(b[0] <= lon <= b[2] and b[1] <= lat <= b[3] for b in self.boxes)


SOURCES = {s.name: s for s in [
    Source("swisstopo", "swissALTI3D (Suisse)", [(5.9, 45.7, 10.6, 47.9)], fetch_swisstopo),
    Source("kartverket", "Kartverket Høydedata DTM1 (Norvège)", [(4.0, 57.8, 31.5, 71.5)], fetch_kartverket),
    Source("usgs", "USGS 3DEP (États-Unis)",
           [(-125.0, 24.3, -66.8, 49.6), (-170.0, 51.0, -129.0, 72.0),
            (-161.0, 18.8, -154.7, 22.4), (-67.4, 17.8, -65.1, 18.6)], fetch_usgs),
    Source("gsi", "GSI tuiles d'élévation (Japon)", [(122.9, 24.0, 146.0, 45.6)], fetch_gsi),
    Source("ign", "IGN MNT LiDAR HD / RGE ALTI (France)", list(IGN_REGION_BOXES.values()), fetch_ign),
]}
AUTO_ORDER = ["swisstopo", "kartverket", "usgs", "gsi", "ign"]


class LayeredInputs(list):
    """Rasters principaux (liste de chemins GDAL) + `fill_loader`, appelable qui fournit
    à la demande une source secondaire pour combler les trous (voir raster.warp_window)."""
    fill_loader = None


def with_copernicus_fill(inputs, lat, lon, radius_m, cache_dir, product_m=30, log=print):
    """Enveloppe `inputs` (chemin ou liste) : Copernicus GLO-30/90 ne sera téléchargé que
    si le disque a des trous une fois la source principale ré-échantillonnée."""
    layered = LayeredInputs(list(inputs) if isinstance(inputs, (list, tuple)) else [inputs])

    def loader():
        import download
        log("Complément : téléchargement de Copernicus pour combler le disque")
        tiles = download.download_copernicus(lat, lon, radius_m, cache_dir, product_m=product_m, log=log)
        os.makedirs(cache_dir, exist_ok=True)
        vrt = f"{cache_dir.rstrip('/')}/_fill_mosaic.vrt"
        gdal.BuildVRT(vrt, tiles)
        return vrt

    layered.fill_loader = loader
    return layered


def acquire(source, lat, lon, radius_m, cache_dir, pixel_size_m, copernicus_product=30, log=print):
    """Renvoie (nom_source_utilisée, entrée GDAL : chemin ou liste de chemins).
    source : nom d'une source, 'copernicus' ou 'auto'. En 'auto', les sources
    nationales dont l'emprise contient le point sont essayées dans l'ordre ;
    en cas d'échec ou d'absence de données, repli sur Copernicus GLO-30."""
    if source == "auto":
        candidates = [n for n in AUTO_ORDER if SOURCES[n].covers(lat, lon)]
        if not candidates:
            log(f"auto : ({lat:.4f}, {lon:.4f}) hors des emprises des sources nationales -> Copernicus GLO-{copernicus_product}")
        for name in candidates:
            log(f"auto : ({lat:.4f}, {lon:.4f}) dans l'emprise de {SOURCES[name].label}, essai")
            try:
                return name, SOURCES[name].fetch(lat, lon, radius_m, cache_dir, pixel_size_m, log)
            except SourceUnavailable as e:
                log(f"auto : {name} inutilisable ({e}) -> source suivante")
        if candidates:
            log(f"auto : repli sur Copernicus GLO-{copernicus_product}")
        source = "copernicus"
    if source == "copernicus":
        import download
        tiles = download.download_copernicus(lat, lon, radius_m, cache_dir,
                                             product_m=copernicus_product, log=log)
        vrt = f"{cache_dir.rstrip('/')}/_mosaic.vrt"
        gdal.BuildVRT(vrt, tiles)
        log(f"Mosaïque de {len(tiles)} tuile(s) prête : {vrt}")
        return "copernicus", vrt
    if source not in SOURCES:
        raise ValueError(f"source inconnue : {source}")
    log(f"source demandée : {SOURCES[source].label}")
    return source, SOURCES[source].fetch(lat, lon, radius_m, cache_dir, pixel_size_m, log)
