# SPDX-License-Identifier: GPL-3.0-or-later
"""
Tests de sources.py et de son intégration (raster.py, circle_dem_to_stl.py).
Aucun accès réseau : les appels HTTP et le WCS sont simulés. Ils vérifient la
logique (sélection, décodage, CRS, cache, repli, erreurs), PAS la conformité
des vrais serveurs, qui n'a pas pu être testée ici.
"""
import math
import os
import re
import sys
import urllib.error
import urllib.parse
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal, osr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import circle_dem_to_stl  # noqa: E402
import download  # noqa: E402
import raster  # noqa: E402
import sources  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402

gdal.UseExceptions()
quiet = lambda *a, **k: None  # noqa: E731


# ---------------------------------------------------------------- GSI

def encode_gsi(h):
    h = np.asarray(h, dtype=np.float64)
    v = np.rint(np.nan_to_num(h) * 100).astype(np.int64)
    x = np.where(v >= 0, v, v + 2 ** 24)
    x = np.where(np.isnan(h), 2 ** 23, x)
    return np.stack([(x >> 16) & 255, (x >> 8) & 255, x & 255]).astype(np.uint8)


def make_png(h):
    rgb = encode_gsi(h)
    mem = gdal.GetDriverByName("MEM").Create("", 256, 256, 3, gdal.GDT_Byte)
    for i in range(3):
        mem.GetRasterBand(i + 1).WriteArray(rgb[i])
    name = "/vsimem/test_tile.png"
    gdal.GetDriverByName("PNG").CreateCopy(name, mem)
    f = gdal.VSIFOpenL(name, "rb")
    gdal.VSIFSeekL(f, 0, 2)
    size = gdal.VSIFTellL(f)
    gdal.VSIFSeekL(f, 0, 0)
    data = gdal.VSIFReadL(1, size, f)
    gdal.VSIFCloseL(f)
    gdal.Unlink(name)
    return data


def test_gsi_decode_matches_official_formula():
    h = np.full((256, 256), np.nan)
    h[0, 0], h[0, 1], h[0, 2], h[0, 3] = 0.0, 5.28, -1.23, 3776.12
    out = sources.decode_gsi_rgb(encode_gsi(h))
    assert out[0, 0] == pytest.approx(0.0)
    assert out[0, 1] == pytest.approx(5.28, abs=1e-4)
    assert out[0, 2] == pytest.approx(-1.23, abs=1e-4)
    assert out[0, 3] == pytest.approx(3776.12, abs=1e-3)
    assert np.isnan(out[1, 1])  # (128, 0, 0) = invalide


def test_gsi_zoom_gives_ground_resolution_at_most_pixel_size():
    lat = 35.36
    for pixel in (0.8, 3.0, 15.0, 80.0):
        z = sources.gsi_zoom_for_pixel(lat, pixel, zmax=19)
        res = lambda zz: (2 * sources.WM_ORIGIN / 256) * math.cos(math.radians(lat)) / 2 ** zz  # noqa: E731
        assert res(z) <= pixel < res(z - 1)


def test_gsi_tile_contains_point():
    lat, lon, z = 35.3606, 138.7274, 14
    x, y = sources.gsi_tile_xy(lat, lon, z)
    tile_m = 2 * sources.WM_ORIGIN / 2 ** z
    ct = osr.CoordinateTransformation(sources._srs(4326), sources._srs(3857))
    mx, my, _ = ct.TransformPoint(lon, lat)
    assert -sources.WM_ORIGIN + x * tile_m <= mx < -sources.WM_ORIGIN + (x + 1) * tile_m
    assert sources.WM_ORIGIN - (y + 1) * tile_m < my <= sources.WM_ORIGIN - y * tile_m


@pytest.fixture
def gsi_fake(monkeypatch):
    calls = []

    def fake_get(url, timeout=None, retries=None):
        calls.append(url)
        if "/dem1a_png/" in url:
            return None  # pas de DEM1A ici -> repli sur DEM5A
        return make_png(np.full((256, 256), 123.45))

    monkeypatch.setattr(sources, "_http_get", fake_get)
    return calls


def test_gsi_mosaic_fallback_cache_and_pipeline_crs(tmp_path, gsi_fake):
    lat, lon = 35.3606, 138.7274
    paths = sources.fetch_gsi(lat, lon, 500, str(tmp_path), 15.0, log=quiet)
    assert len(paths) == 1
    n_calls = len(gsi_fake)
    assert any("/dem5a_png/" in u for u in gsi_fake)  # repli DEM1A -> DEM5A

    ds = gdal.Open(paths[0])
    assert ds.RasterXSize % 256 == 0
    arr = ds.GetRasterBand(1).ReadAsArray()
    assert np.allclose(arr[arr != sources.NODATA], 123.45, atol=1e-3)

    # 2e appel : tout vient du cache disque
    sources.fetch_gsi(lat, lon, 500, str(tmp_path), 15.0, log=quiet)
    assert len(gsi_fake) == n_calls + sum(1 for u in gsi_fake[:n_calls] if "/dem1a_png/" in u)

    # Le pipeline doit traiter ce Web Mercator comme à reprojeter vers UTM 54N
    wkt, cx, cy = raster.resolve_target_crs_and_center(paths[0], lat, lon, None, None)
    dst = osr.SpatialReference()
    dst.ImportFromWkt(wkt)
    assert dst.GetAuthorityCode(None) == "32654"
    ex, ey = sources._to_epsg(32654, lat, lon)
    assert (cx, cy) == pytest.approx((ex, ey), abs=0.01)
    elev, xs, ys, nodata = raster.warp_window(paths[0], wkt, cx, cy, 500, 15.0, progress_cb=quiet)
    assert elev[elev.shape[0] // 2, elev.shape[1] // 2] == pytest.approx(123.45, abs=1e-2)


def test_gsi_web_mercator_with_cx_cy_uses_utm(tmp_path, gsi_fake):
    lat, lon = 35.3606, 138.7274
    paths = sources.fetch_gsi(lat, lon, 500, str(tmp_path), 15.0, log=quiet)
    mx, my = sources._to_epsg(3857, lat, lon)
    wkt, cx, cy = raster.resolve_target_crs_and_center(paths[0], None, None, mx, my)
    ex, ey = sources._to_epsg(32654, lat, lon)
    assert (cx, cy) == pytest.approx((ex, ey), abs=0.01)


def test_gsi_no_tile_anywhere_is_no_coverage(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", lambda *a, **k: None)
    with pytest.raises(sources.NoCoverage):
        sources.fetch_gsi(35.36, 138.72, 500, str(tmp_path), 15.0, log=quiet)


# ---------------------------------------------------------------- swisstopo

def _item(year, tile, assets):
    return {"id": f"swissalti3d_{year}_{tile}",
            "properties": {"datetime": f"{year}-01-01T00:00:00Z"},
            "assets": {f"a{i}": {"href": h, "eo:gsd": g} for i, (h, g) in enumerate(assets)}}


def test_swisstopo_selects_latest_version_and_closest_gsd():
    items = [
        _item(2019, "2600-1199", [("https://x/a_2019_0.5.tif", 0.5), ("https://x/a_2019_2.tif", 2.0)]),
        _item(2023, "2600-1199", [("https://x/a_2023_0.5.tif", 0.5), ("https://x/a_2023_2.tif", 2.0),
                                   ("https://x/a_2023.xyz.zip", 0.5)]),
        _item(2021, "2601-1199", [("https://x/b_2021_2.tif", 2.0)]),
    ]
    assert sources.select_swisstopo_assets(items, 2.0) == ["https://x/a_2023_2.tif", "https://x/b_2021_2.tif"]
    assert sources.select_swisstopo_assets(items, 0.5) == ["https://x/a_2023_0.5.tif", "https://x/b_2021_2.tif"]
    assert sources.swisstopo_gsd(15.0) == 2.0 and sources.swisstopo_gsd(0.8) == 0.5


def test_swisstopo_fetch_paginates_downloads_and_caches(tmp_path, monkeypatch):
    tif = tmp_path / "src.tif"
    write_synthetic_tif(tif, np.full((10, 10), 500.0), 2600000.0, 1200000.0, 2.0, 2056)
    payload = tif.read_bytes()
    pages = {
        "p1": {"features": [_item(2023, "2600-1199", [("https://x/t1_2.tif", 2.0)])],
               "links": [{"rel": "next", "href": "p2"}]},
        "p2": {"features": [_item(2023, "2601-1199", [("https://x/t2_2.tif", 2.0)])], "links": []},
    }
    seen_params = {}

    def fake_json(url, params=None):
        if url.startswith(sources.SWISS_STAC_ITEMS):
            seen_params["url"] = url
            return pages["p1"]
        return pages[url]

    gets = []
    monkeypatch.setattr(sources, "_http_json", fake_json)
    monkeypatch.setattr(sources, "_http_get", lambda url, **k: gets.append(url) or payload)

    out = sources.fetch_swisstopo(46.0, 7.6, 1000, str(tmp_path / "c"), 10.0, log=quiet)
    assert [os.path.basename(p) for p in out] == ["t1_2.tif", "t2_2.tif"]
    assert "bbox=" in seen_params["url"]
    sources.fetch_swisstopo(46.0, 7.6, 1000, str(tmp_path / "c"), 10.0, log=quiet)
    assert len(gets) == 2  # 2e appel entièrement servi par le cache


def test_swisstopo_no_items_is_no_coverage(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_json", lambda *a, **k: {"features": [], "links": []})
    with pytest.raises(sources.NoCoverage):
        sources.fetch_swisstopo(46.0, 7.6, 1000, str(tmp_path), 10.0, log=quiet)


# ---------------------------------------------------------------- Kartverket

class FakeWcs:
    """Serveur WCS ArcGIS simulé. `accept` = variantes acceptées ; toute autre
    requête -> HTTP 400. La valeur renvoyée vaut 600 + 0,01·x (x = abscisse EPSG:25833)
    pour vérifier le placement des morceaux."""

    def __init__(self, accept=("v100_geotiff",), empty=False, html=False, wrong_size=False, sentinel=False):
        self.urls, self.accept = [], accept
        self.empty, self.html, self.wrong_size, self.sentinel = empty, html, wrong_size, sentinel

    @staticmethod
    def variant_of(url):
        q = dict(urllib.parse.parse_qsl(urllib.parse.urlparse(url).query))
        if q.get("VERSION") == "1.0.0":
            return {"GeoTIFF": "v100_geotiff", "image/GeoTIFF": "v100_image"}.get(q.get("FORMAT"))
        if q.get("VERSION") == "2.0.1" and q.get("FORMAT") == "image/GeoTIFF":
            return "v201_image"
        return None

    def __call__(self, url, timeout=None, retries=None):
        self.urls.append(url)
        variant = self.variant_of(url)
        if variant is None or variant not in self.accept:
            raise sources.SourceUnavailable(f"{url} : HTTP 400")
        if self.html:
            return b"<html>ArcGIS Server Error</html>"
        q = dict(urllib.parse.parse_qsl(urllib.parse.urlparse(url).query))
        if variant.startswith("v100"):
            assert q["CRS"] == "EPSG:25833" and q["COVERAGE"] == sources.KARTVERKET_COVERAGE_DEFAULT
            x0, y0, x1, y1 = [float(v) for v in q["BBOX"].split(",")]
            w, h = int(q["WIDTH"]), int(q["HEIGHT"])
        else:
            subs = [v for k, v in urllib.parse.parse_qsl(urllib.parse.urlparse(url).query) if k == "SUBSET"]
            x0, x1 = [float(v) for v in subs[0][2:-1].split(",")]
            y0, y1 = [float(v) for v in subs[1][2:-1].split(",")]
            w, h = [int(v) for v in re.findall(r"\((\d+)\)", q["SCALESIZE"])]
        if self.wrong_size:
            w, h = w + 1, h
        xs = x0 + (np.arange(w) + 0.5) * (x1 - x0) / w
        arr = np.tile((600.0 + 0.01 * xs).astype(np.float32), (h, 1))
        if self.empty:
            arr[:] = -3.4028235e38
        elif self.sentinel:
            arr[:, :3] = -3.4028235e38
        mem = gdal.GetDriverByName("MEM").Create("", w, h, 1, gdal.GDT_Float32)
        mem.SetGeoTransform((x0, (x1 - x0) / w, 0, y1, 0, -(y1 - y0) / h))
        mem.SetProjection(sources._srs(25833).ExportToWkt())
        mem.GetRasterBand(1).WriteArray(arr)
        name = f"/vsimem/fakewcs_{len(self.urls)}.tif"
        gdal.GetDriverByName("GTiff").CreateCopy(name, mem)
        f = gdal.VSIFOpenL(name, "rb")
        gdal.VSIFSeekL(f, 0, 2)
        size = gdal.VSIFTellL(f)
        gdal.VSIFSeekL(f, 0, 0)
        data = gdal.VSIFReadL(1, size, f)
        gdal.VSIFCloseL(f)
        gdal.Unlink(name)
        return data


LAT_K, LON_K = 58.9863, 6.1904


def _check_kartverket_output(path, res, lat=LAT_K, lon=LON_K):
    ds = gdal.Open(path)
    assert ds.GetGeoTransform()[1] == pytest.approx(res)
    arr = ds.GetRasterBand(1).ReadAsArray()
    x, _ = sources._to_epsg(25833, lat, lon)
    gt = ds.GetGeoTransform()
    col = int((x - gt[0]) / res)
    assert arr[arr.shape[0] // 2, col] == pytest.approx(600.0 + 0.01 * x, abs=0.2 + 0.01 * res)


def test_wcs_request_urls_use_image_geotiff_not_image_tiff():
    b, size = (-5, 10, 15, 30), (20, 20)
    u100 = sources.wcs_request_url("https://h/s?x=1", "cov", "v100_geotiff", b, size)
    u100i = sources.wcs_request_url("https://h/s", "cov", "v100_image", b, size)
    u201 = sources.wcs_request_url("https://h/s", "cov", "v201_image", b, size)
    assert u100.startswith("https://h/s?SERVICE=WCS&VERSION=1.0.0&REQUEST=GetCoverage&COVERAGE=cov")
    assert "CRS=EPSG:25833&BBOX=-5,10,15,30&WIDTH=20&HEIGHT=20&FORMAT=GeoTIFF" in u100
    assert "FORMAT=image/GeoTIFF" in u100i and "FORMAT=image/GeoTIFF" in u201
    assert "SUBSET=x(-5,15)&SUBSET=y(10,30)&SCALESIZE=x(20),y(20)" in u201
    assert "image/tiff" not in u100 + u100i + u201


@pytest.mark.parametrize("accepted", [("v100_geotiff",), ("v100_image",), ("v201_image",)])
def test_kartverket_discovers_accepted_variant_and_caches(tmp_path, monkeypatch, accepted):
    fake = FakeWcs(accept=accepted)
    monkeypatch.setattr(sources, "_http_get", fake)
    out = sources.fetch_kartverket(LAT_K, LON_K, 3000, str(tmp_path), 8.0, log=quiet)
    _check_kartverket_output(out[0], 8.0)
    # une fois la bonne variante trouvée, elle est réutilisée (pas de nouvel essai des autres)
    good = [u for u in fake.urls if FakeWcs.variant_of(u) == accepted[0]]
    assert len(fake.urls) - len(good) <= 2 and len(good) >= 2
    assert not list((tmp_path / "kartverket").glob("_tmp_*"))
    monkeypatch.setattr(sources, "_http_get", lambda *a, **k: pytest.fail("réseau malgré le cache"))
    assert sources.fetch_kartverket(LAT_K, LON_K, 3000, str(tmp_path), 8.0, log=quiet) == out


def test_kartverket_requests_twice_the_working_resolution(tmp_path, monkeypatch):
    fake = FakeWcs()
    monkeypatch.setattr(sources, "_http_get", fake)
    sources.fetch_kartverket(LAT_K, LON_K, 1000, str(tmp_path), 8.0, log=quiet)
    probe = dict(urllib.parse.parse_qsl(urllib.parse.urlparse(fake.urls[0]).query))
    assert (probe["WIDTH"], probe["HEIGHT"]) == ("100", "100")  # 400 m à 4 m/pixel


def test_kartverket_multiple_chunks_are_mosaicked(tmp_path, monkeypatch):
    fake = FakeWcs()
    monkeypatch.setattr(sources, "_http_get", fake)
    monkeypatch.setattr(sources, "KARTVERKET_CHUNK_PX", 100)
    out = sources.fetch_kartverket(LAT_K, LON_K, 1000, str(tmp_path), 10.0, log=quiet)
    _check_kartverket_output(out[0], 10.0)
    assert len(fake.urls) > 10  # sonde + de nombreux morceaux


def test_kartverket_native_resolution_floor(tmp_path, monkeypatch):
    fake = FakeWcs()
    monkeypatch.setattr(sources, "_http_get", fake)
    out = sources.fetch_kartverket(LAT_K, LON_K, 300, str(tmp_path), 0.8, log=quiet)  # pixel < 1 m -> 1 m
    _check_kartverket_output(out[0], 1.0)


def test_kartverket_sentinel_values_become_nodata(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", FakeWcs(sentinel=True))
    out = sources.fetch_kartverket(LAT_K, LON_K, 1000, str(tmp_path), 10.0, log=quiet)
    ds = gdal.Open(out[0])
    arr = ds.GetRasterBand(1).ReadAsArray()
    assert ds.GetRasterBand(1).GetNoDataValue() == sources.NODATA
    assert not np.any((arr < -1e5) & (arr != sources.NODATA))  # aucune valeur sentinelle brute
    assert (arr == sources.NODATA).any()                        # remplacées par le nodata explicite


def test_kartverket_empty_area_is_no_coverage(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", FakeWcs(empty=True))
    with pytest.raises(sources.NoCoverage):
        sources.fetch_kartverket(LAT_K, LON_K, 1000, str(tmp_path), 10.0, log=quiet)
    assert not list((tmp_path / "kartverket").glob("*"))


def test_kartverket_non_tiff_answer_is_source_unavailable(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", FakeWcs(html=True))
    with pytest.raises(sources.SourceUnavailable, match="non TIFF"):
        sources.fetch_kartverket(LAT_K, LON_K, 1000, str(tmp_path), 10.0, log=quiet)


def test_kartverket_wrong_size_answer_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", FakeWcs(wrong_size=True))
    with pytest.raises(sources.SourceUnavailable, match="taille reçue"):
        sources.fetch_kartverket(LAT_K, LON_K, 1000, str(tmp_path), 10.0, log=quiet)


def test_kartverket_all_variants_refused_lists_every_attempt(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", FakeWcs(accept=()))
    with pytest.raises(sources.SourceUnavailable) as exc:
        sources.fetch_kartverket(LAT_K, LON_K, 1000, str(tmp_path), 10.0, log=quiet)
    msg = str(exc.value)
    assert all(v in msg for v in sources.KARTVERKET_VARIANTS)


def test_kartverket_too_many_chunks_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", FakeWcs())
    monkeypatch.setattr(sources, "KARTVERKET_MAX_CHUNKS", 2)
    monkeypatch.setattr(sources, "KARTVERKET_CHUNK_PX", 100)
    with pytest.raises(sources.SourceUnavailable, match="requêtes"):
        sources.fetch_kartverket(LAT_K, LON_K, 3000, str(tmp_path), 10.0, log=quiet)


# ---------------------------------------------------------------- IGN

class FakeWms:
    """WMS-R Géoplateforme simulé. `accept` = variantes acceptées ; `lidar` = la
    couche LiDAR HD a des données ; `rendered` = renvoie une image 8 bits (estompage).
    Valeur = 300 + 0,01·x (EPSG:2154) pour contrôler le placement."""

    def __init__(self, accept=("geotiff_normal",), lidar=True, rgealti=True, rendered=False, partial=False):
        self.urls, self.accept = [], accept
        self.lidar, self.rgealti, self.rendered, self.partial = lidar, rgealti, rendered, partial

    @staticmethod
    def variant_of(q):
        fmt = {"image/geotiff": "geotiff", "image/tiff": "tiff"}.get(q.get("FORMAT"))
        if fmt is None:
            return None
        return f"{fmt}_{'normal' if q.get('STYLES') == 'normal' else 'default'}"

    def __call__(self, url, timeout=None, retries=None):
        self.urls.append(url)
        q = dict(urllib.parse.parse_qsl(urllib.parse.urlparse(url).query, keep_blank_values=True))
        assert q["SERVICE"] == "WMS" and q["VERSION"] == "1.3.0" and q["REQUEST"] == "GetMap"
        if self.variant_of(q) not in self.accept:
            raise sources.SourceUnavailable(f"{url} : HTTP 400")
        x0, y0, x1, y1 = [float(v) for v in q["BBOX"].split(",")]
        w, h = int(q["WIDTH"]), int(q["HEIGHT"])
        is_lidar = "LIDAR-HD" in q["LAYERS"]
        has_data = self.lidar if is_lidar else self.rgealti
        xs = x0 + (np.arange(w) + 0.5) * (x1 - x0) / w
        arr = np.tile((300.0 + 0.01 * xs).astype(np.float32), (h, 1))
        if not has_data:
            arr[:] = -99999.0
        elif self.partial and is_lidar:
            arr[:, : w // 3] = -99999.0
        if self.rendered:
            mem = gdal.GetDriverByName("MEM").Create("", w, h, 3, gdal.GDT_Byte)
            for i in range(3):
                mem.GetRasterBand(i + 1).WriteArray(np.full((h, w), 128, dtype=np.uint8))
        else:
            mem = gdal.GetDriverByName("MEM").Create("", w, h, 1, gdal.GDT_Float32)
            mem.GetRasterBand(1).WriteArray(arr)
        mem.SetGeoTransform((x0, (x1 - x0) / w, 0, y1, 0, -(y1 - y0) / h))
        mem.SetProjection(sources._srs(int(q["CRS"].split(":")[1])).ExportToWkt())
        name = f"/vsimem/fakewms_{len(self.urls)}.tif"
        gdal.GetDriverByName("GTiff").CreateCopy(name, mem)
        f = gdal.VSIFOpenL(name, "rb")
        gdal.VSIFSeekL(f, 0, 2)
        size = gdal.VSIFTellL(f)
        gdal.VSIFSeekL(f, 0, 0)
        data = gdal.VSIFReadL(1, size, f)
        gdal.VSIFCloseL(f)
        gdal.Unlink(name)
        return data


LAT_I, LON_I = 45.1885, 5.7245  # Grenoble


def _layers_requested(fake):
    return {dict(urllib.parse.parse_qsl(urllib.parse.urlparse(u).query))["LAYERS"] for u in fake.urls}


def test_ign_url_uses_native_crs_and_raw_style():
    u = sources.ign_getmap_url("https://h/wms?SERVICE=WMS&", "LAYER", "geotiff_normal", 2154, (1, 2, 3, 4), (10, 20))
    assert u.startswith("https://h/wms?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&LAYERS=LAYER&STYLES=normal")
    assert "CRS=EPSG:2154&BBOX=1,2,3,4&WIDTH=10&HEIGHT=20&FORMAT=image/geotiff" in u
    assert "STYLES=&" in sources.ign_getmap_url("https://h/wms", "L", "tiff_default", 2154, (1, 2, 3, 4), (10, 20))


@pytest.mark.parametrize("accepted", [("geotiff_normal",), ("geotiff_default",), ("tiff_normal",)])
def test_ign_lidar_variant_discovery_cache_and_values(tmp_path, monkeypatch, accepted):
    fake = FakeWms(accept=accepted)
    monkeypatch.setattr(sources, "_http_get", fake)
    out = sources.fetch_ign(LAT_I, LON_I, 2000, str(tmp_path), 5.0, log=quiet)
    ds = gdal.Open(out[0])
    assert ds.GetGeoTransform()[1] == pytest.approx(5.0)
    assert os.path.basename(out[0]).startswith("lidarhd_fxx_")
    assert all("LIDAR-HD_MNT" in l for l in _layers_requested(fake))  # pas de repli inutile
    x, _ = sources._to_epsg(2154, LAT_I, LON_I)
    arr = ds.GetRasterBand(1).ReadAsArray()
    col = int((x - ds.GetGeoTransform()[0]) / 5.0)
    assert arr[arr.shape[0] // 2, col] == pytest.approx(300.0 + 0.01 * x, abs=0.2)
    assert not list((tmp_path / "ign").glob("_tmp_*"))
    monkeypatch.setattr(sources, "_http_get", lambda *a, **k: pytest.fail("réseau malgré le cache"))
    assert sources.fetch_ign(LAT_I, LON_I, 2000, str(tmp_path), 5.0, log=quiet) == out


def test_ign_falls_back_to_rgealti_when_no_lidar_hd(tmp_path, monkeypatch):
    fake = FakeWms(lidar=False)
    monkeypatch.setattr(sources, "_http_get", fake)
    out = sources.fetch_ign(LAT_I, LON_I, 1000, str(tmp_path), 5.0, log=quiet)
    assert os.path.basename(out[0]).startswith("rgealti_fxx_")
    assert sources.IGN_RGEALTI_LAYER in _layers_requested(fake)


def test_ign_rejects_rendered_images_instead_of_using_them_as_elevation(tmp_path, monkeypatch):
    fake = FakeWms(rendered=True)
    monkeypatch.setattr(sources, "_http_get", fake)
    with pytest.raises(sources.SourceUnavailable, match="non flottantes"):
        sources.fetch_ign(LAT_I, LON_I, 1000, str(tmp_path), 5.0, log=quiet)


def test_ign_nothing_anywhere_is_no_coverage(tmp_path, monkeypatch):
    monkeypatch.setattr(sources, "_http_get", FakeWms(lidar=False, rgealti=False))
    with pytest.raises(sources.NoCoverage):
        sources.fetch_ign(LAT_I, LON_I, 1000, str(tmp_path), 5.0, log=quiet)


def test_ign_partial_lidar_coverage_is_reported(tmp_path, monkeypatch):
    messages = []
    monkeypatch.setattr(sources, "_http_get", FakeWms(partial=True))
    sources.fetch_ign(LAT_I, LON_I, 1000, str(tmp_path), 5.0, log=messages.append)
    assert any("sans donnée" in m for m in messages)


def test_ign_outside_france_and_reunion_crs(tmp_path, monkeypatch):
    with pytest.raises(sources.NoCoverage):
        sources.fetch_ign(0.0, -30.0, 1000, str(tmp_path), 5.0, log=quiet)
    fake = FakeWms()
    monkeypatch.setattr(sources, "_http_get", fake)
    sources.fetch_ign(-21.2433, 55.7087, 1000, str(tmp_path), 5.0, log=quiet)
    q = dict(urllib.parse.parse_qsl(urllib.parse.urlparse(fake.urls[0]).query))
    assert q["CRS"] == "EPSG:2975" and q["LAYERS"].endswith("RGR92UTM40S")


def test_ign_in_auto_selection_and_cli_choices():
    assert "ign" in sources.SOURCES and sources.SOURCES["ign"].covers(LAT_I, LON_I)
    assert not sources.SOURCES["ign"].covers(35.0, 138.0)
    assert sources.AUTO_ORDER[-1] == "ign"


# ---------------------------------------------------------------- USGS

def _prod(url, date):
    return {"downloadURL": url, "publicationDate": date}


def test_usgs_prefers_1m_then_falls_back_and_orders_by_date(monkeypatch):
    queried = []

    def fake_json(url, params=None):
        queried.append(params["datasets"])
        if params["datasets"] == sources.USGS_DATASETS_1M[0]:
            return {"total": 2, "items": [_prod("https://s/b.tif", "2022-01-01"), _prod("https://s/a.tif", "2019-01-01")]}
        return {"total": 0, "items": []}

    monkeypatch.setattr(sources, "_http_json", fake_json)
    out = sources.fetch_usgs(36.97, -110.1, 1000, "unused", 1.0, log=quiet)
    assert out == ["/vsicurl/https://s/a.tif", "/vsicurl/https://s/b.tif"]  # le plus récent en dernier

    queried.clear()
    monkeypatch.setattr(sources, "_http_json", lambda url, params=None: (queried.append(params["datasets"]) or
                        {"total": 1, "items": [_prod("https://s/n.tif", "2020")]}))
    sources.fetch_usgs(36.97, -110.1, 9000, "unused", 15.0, log=quiet)
    assert sources.USGS_DATASETS_1M[0] not in queried  # pixel de 15 m : pas de 1 m inutile


def test_usgs_nothing_found_is_no_coverage(monkeypatch):
    monkeypatch.setattr(sources, "_http_json", lambda *a, **k: {"total": 0, "items": []})
    with pytest.raises(sources.NoCoverage):
        sources.fetch_usgs(36.97, -110.1, 1000, "unused", 1.0, log=quiet)


# ---------------------------------------------------------------- sélection auto

@pytest.fixture
def fake_copernicus(tmp_path, monkeypatch):
    tif = tmp_path / "cop.tif"
    write_synthetic_tif(tif, np.full((10, 10), 10.0), 7.0, 47.0, 0.001, 4326)
    monkeypatch.setattr(download, "download_copernicus", lambda *a, **k: [str(tif)])
    return tif


def test_auto_uses_national_source_when_it_covers(tmp_path, monkeypatch):
    monkeypatch.setitem(sources.SOURCES, "swisstopo",
                        sources.Source("swisstopo", "x", [(5, 45, 11, 48)], lambda *a, **k: ["/tmp/ch.tif"]))
    name, inp = sources.acquire("auto", 46.0, 7.6, 1000, str(tmp_path), 10.0, log=quiet)
    assert (name, inp) == ("swisstopo", ["/tmp/ch.tif"])


def test_auto_falls_back_to_copernicus_on_failure_and_outside_coverage(tmp_path, monkeypatch, fake_copernicus):
    def boom(*a, **k):
        raise sources.NoCoverage("rien")
    monkeypatch.setitem(sources.SOURCES, "swisstopo", sources.Source("swisstopo", "x", [(5, 45, 11, 48)], boom))
    name, inp = sources.acquire("auto", 46.0, 7.6, 1000, str(tmp_path), 10.0, log=quiet)
    assert name == "copernicus" and inp.endswith("_mosaic.vrt")
    name, _ = sources.acquire("auto", 0.0, -30.0, 1000, str(tmp_path), 10.0, log=quiet)  # plein Atlantique
    assert name == "copernicus"


def test_explicit_source_failure_is_not_silently_replaced(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise sources.SourceUnavailable("panne")
    monkeypatch.setitem(sources.SOURCES, "usgs", sources.Source("usgs", "x", [(-180, -90, 180, 90)], boom))
    with pytest.raises(sources.SourceUnavailable):
        sources.acquire("usgs", 36.9, -110.1, 1000, str(tmp_path), 10.0, log=quiet)


# ---------------------------------------------------------------- réseau

def test_http_get_retries_transient_errors_and_maps_404(monkeypatch):
    monkeypatch.setattr(sources.time, "sleep", lambda s: None)
    state = {"n": 0}

    class Resp:
        def __enter__(self): return self
        def __exit__(self, *a): return False
        def read(self): return b"ok"

    def flaky(req, timeout=None):
        state["n"] += 1
        if state["n"] < 3:
            raise urllib.error.URLError("connexion coupée")
        return Resp()

    monkeypatch.setattr(sources.urllib.request, "urlopen", flaky)
    assert sources._http_get("https://x/y") == b"ok" and state["n"] == 3

    def not_found(req, timeout=None):
        raise urllib.error.HTTPError("u", 404, "nf", {}, None)
    monkeypatch.setattr(sources.urllib.request, "urlopen", not_found)
    assert sources._http_get("https://x/y") is None

    def forbidden(req, timeout=None):
        state["n"] += 1
        raise urllib.error.HTTPError("u", 403, "no", {}, None)
    state["n"] = 0
    monkeypatch.setattr(sources.urllib.request, "urlopen", forbidden)
    with pytest.raises(sources.SourceUnavailable):
        sources._http_get("https://x/y")
    assert state["n"] == 1  # 403 : pas de nouvel essai


def test_http_get_gives_up_after_retries(monkeypatch):
    monkeypatch.setattr(sources.time, "sleep", lambda s: None)
    monkeypatch.setattr(sources.urllib.request, "urlopen",
                        lambda req, timeout=None: (_ for _ in ()).throw(TimeoutError("timeout")))
    with pytest.raises(sources.SourceUnavailable):
        sources._http_get("https://x/y")


def test_copernicus_network_failure_is_explicit_and_leaves_no_partial_file(tmp_path, monkeypatch):
    monkeypatch.setattr(download.urllib.request, "urlopen",
                        lambda req, timeout=None: (_ for _ in ()).throw(urllib.error.URLError("hors ligne")))
    with pytest.raises(RuntimeError, match="incomplet"):
        download.download_copernicus(35.3606, 138.7274, 9300, str(tmp_path), log=quiet)
    assert not list(tmp_path.glob("*"))


# ---------------------------------------------------------------- intégration CLI

def test_cli_with_national_source_returning_two_adjacent_tiles(tmp_path, monkeypatch):
    """Source projetée (EPSG:2056) fournie en deux dalles : le pipeline complet
    (acquire -> warp_window sur liste -> maillage -> STL) doit fonctionner."""
    lat, lon = 46.0, 7.6
    x, y = sources._to_epsg(2056, lat, lon)
    xx, yy = np.meshgrid(np.arange(1500), np.arange(1500))
    dem = 1000.0 + 0.05 * xx + 0.02 * yy
    left, right = tmp_path / "l.tif", tmp_path / "r.tif"
    ox, oy = x - 1500.0, y + 1500.0  # 3000 m de côté, pixel de 2 m
    write_synthetic_tif(left, dem[:, :750], ox, oy, 2.0, 2056)
    write_synthetic_tif(right, dem[:, 750:], ox + 1500.0, oy, 2.0, 2056)
    seen = {}

    def fake_fetch(lat_, lon_, radius_m, cache_dir, pixel_size_m, log=print):
        seen["pixel"] = pixel_size_m
        return [str(left), str(right)]

    monkeypatch.setitem(sources.SOURCES, "swisstopo", sources.Source("swisstopo", "x", [(5, 45, 11, 48)], fake_fetch))
    out = tmp_path / "out.stl"
    monkeypatch.setattr(sys, "argv", [
        "circle_dem_to_stl.py", "--source", "swisstopo", "--lat", str(lat), "--lon", str(lon),
        "--radius", "500", "--diameter-mm", "50", "--vexag", "1", "--base-mm", "2",
        "--print-spacing-mm", "0.5", "--output", str(out), "--cache-dir", str(tmp_path / "c")])
    circle_dem_to_stl.main()
    assert out.stat().st_size > 84
    # pixel de travail = 0,5 mm / (25 mm / 500 m) = 10 m, transmis à la source AVANT téléchargement
    assert seen["pixel"] == pytest.approx(10.0)


def test_cli_rejects_source_with_input(monkeypatch, tmp_path):
    monkeypatch.setattr(sys, "argv", ["x", "--source", "gsi", "--input", "a.tif", "--lat", "35", "--lon", "138",
                                      "--radius", "500", "--output", str(tmp_path / "o.stl")])
    with pytest.raises(SystemExit):
        circle_dem_to_stl.main()


# ---------------------------------------------------------------- résolution / qualité

def test_usgs_tries_1m_whenever_working_pixel_is_below_threshold(monkeypatch):
    queried = []

    def fake_json(url, params=None):
        queried.append(params["datasets"])
        return {"total": 1, "items": [_prod("https://s/a.tif", "2022-01-01")]}

    monkeypatch.setattr(sources, "_http_json", fake_json)
    sources.fetch_usgs(37.7, -119.6, 3000, "unused", 8.0, log=quiet)   # 8 m de travail : 1 m disponible -> 1 m
    assert queried[0] == sources.USGS_DATASETS_1M[0]
    queried.clear()
    sources.fetch_usgs(37.7, -119.6, 9000, "unused", 15.0, log=quiet)  # 15 m : le ~10 m suffit
    assert queried[0] == sources.USGS_DATASETS_13[0]


def test_full_res_asks_every_source_for_its_finest_resolution(tmp_path, monkeypatch):
    lat, lon = 46.0, 7.6
    x, y = sources._to_epsg(2056, lat, lon)
    tif = tmp_path / "ch.tif"
    write_synthetic_tif(tif, np.full((1500, 1500), 1500.0), x - 1500.0, y + 1500.0, 2.0, 2056)
    seen = []

    def fake_fetch(lat_, lon_, radius_m, cache_dir, pixel_size_m, log=print):
        seen.append(pixel_size_m)
        return [str(tif)]

    monkeypatch.setitem(sources.SOURCES, "swisstopo", sources.Source("swisstopo", "x", [(5, 45, 11, 48)], fake_fetch))
    base = ["x", "--source", "swisstopo", "--lat", str(lat), "--lon", str(lon), "--radius", "500",
            "--diameter-mm", "50", "--vexag", "1", "--base-mm", "2", "--print-spacing-mm", "0.5",
            "--cache-dir", str(tmp_path / "c")]
    monkeypatch.setattr(sys, "argv", base + ["--output", str(tmp_path / "a.stl")])
    circle_dem_to_stl.main()
    monkeypatch.setattr(sys, "argv", base + ["--full-res", "--output", str(tmp_path / "b.stl")])
    circle_dem_to_stl.main()
    assert seen[0] == pytest.approx(10.0) and seen[1] == sources.FULL_RES_PIXEL_M
    # le maillage, lui, reste identique : même nombre de triangles
    assert (tmp_path / "a.stl").stat().st_size == (tmp_path / "b.stl").stat().st_size
    assert sources.swisstopo_gsd(sources.FULL_RES_PIXEL_M) == 0.5


def test_compare_dems_reports_offset_and_small_scatter(tmp_path):
    import compare_dems

    def terrain(x, y):
        return 1000.0 + 0.05 * (x - 2600000) + 0.02 * (y - 1200000) + 15.0 * np.sin((x - 2600000) / 300.0)

    def write(path, pixel, offset):
        n = int(4000 / pixel)
        xs = 2600000.0 + (np.arange(n) + 0.5) * pixel
        ys = 1204000.0 - (np.arange(n) + 0.5) * pixel
        xx, yy = np.meshgrid(xs, ys)
        write_synthetic_tif(path, (terrain(xx, yy) + offset).astype(np.float32), 2600000.0, 1204000.0, pixel, 2056)

    a, b = tmp_path / "fine.tif", tmp_path / "coarse.tif"
    write(a, 0.5, 0.0)
    write(b, 2.0, 3.0)
    s = compare_dems.compare_dems(str(a), str(b), None, None, 2602000.0, 1202000.0, 1500.0, 8.0)
    assert s["mean"] == pytest.approx(3.0, abs=0.05)
    assert s["std"] < 0.2 and s["cells_both"] > 0.95 * s["cells_disc"]
    assert s["relief_a_m"] > 50


def test_compare_dems_without_overlap_is_an_error(tmp_path):
    import compare_dems
    a, b = tmp_path / "a.tif", tmp_path / "b.tif"
    write_synthetic_tif(a, np.full((100, 100), 10.0), 2600000.0, 1200100.0, 2.0, 2056)
    write_synthetic_tif(b, np.full((100, 100), 10.0), 2650000.0, 1250100.0, 2.0, 2056)
    with pytest.raises(ValueError):
        compare_dems.compare_dems(str(a), str(b), None, None, 2600100.0, 1200000.0, 90.0, 8.0)


def test_compare_dems_ignores_filled_halo_of_the_smaller_raster(tmp_path):
    """A ne couvre qu'une partie du disque : la comparaison ne doit porter que sur les
    cellules réellement mesurées (pas sur les 50 pixels que FillNodata inventerait)."""
    import compare_dems

    def terrain(x, y):
        return 1000.0 + 0.05 * (x - 2600000) + 0.02 * (y - 1200000)

    def write(path, x_min, x_max, pixel):
        n_x, n_y = int((x_max - x_min) / pixel), int(4000 / pixel)
        xs = x_min + (np.arange(n_x) + 0.5) * pixel
        ys = 1204000.0 - (np.arange(n_y) + 0.5) * pixel
        xx, yy = np.meshgrid(xs, ys)
        write_synthetic_tif(path, terrain(xx, yy).astype(np.float32), x_min, 1204000.0, pixel, 2056)

    a, b = tmp_path / "half.tif", tmp_path / "full.tif"
    write(a, 2600000.0, 2602000.0, 2.0)        # moitié ouest seulement
    write(b, 2600000.0, 2604000.0, 2.0)
    s = compare_dems.compare_dems(str(a), str(b), None, None, 2602000.0, 1202000.0, 1500.0, 8.0)
    assert s["cells_only_b"] > 0.4 * s["cells_disc"]   # la moitié est du disque n'existe que dans B
    assert s["cells_only_a"] == 0
    assert s["max_abs"] < 0.5                           # aucune différence inventée par le comblement
