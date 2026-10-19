# SPDX-License-Identifier: GPL-3.0-or-later
"""Repli de résolution GSI : une tuile absente d'une couche fine est remplacée par la tuile de la couche
suivante (zoom 15 puis 14), interpolée par le ré-échantillonnage, au lieu de laisser un trou ou de lever NoCoverage."""
import re
import sys
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal, osr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import raster  # noqa: E402
import sources  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402
from test_sources import make_png  # noqa: E402

LAT, LON = 35.36, 138.72
CX, CY, _ = osr.CoordinateTransformation(sources._srs(4326), sources._srs(3857)).TransformPoint(LON, LAT)
URL_RE = re.compile(r"/xyz/([a-z0-9_]+)/(\d+)/(\d+)/(\d+)\.png")


def tile_geo(z, x, y):
    tile_m = 2 * sources.WM_ORIGIN / 2 ** z
    return -sources.WM_ORIGIN + x * tile_m, sources.WM_ORIGIN - y * tile_m, tile_m / 256


def ramp(z, x, y, offset=0.0):
    """Altitudes 256x256 d'un plan incliné défini en mètres Web Mercator, identique à tout zoom."""
    x_o, y_o, px = tile_geo(z, x, y)
    X, Y = np.meshgrid(x_o + (np.arange(256) + 0.5) * px, y_o - (np.arange(256) + 0.5) * px)
    return 1000.0 + offset + 0.3 * (X - CX) + 0.5 * (Y - CY)


@pytest.fixture
def fake_gsi(monkeypatch):
    """Réseau GSI simulé. `avail` : {couche: prédicat(z, x, y)} ; `offset` : {couche: décalage d'altitude}."""
    state = {"avail": {}, "offset": {}, "calls": []}

    def fake_get(url, timeout=None, retries=None):
        layer, z, x, y = URL_RE.search(url).groups()
        z, x, y = int(z), int(x), int(y)
        state["calls"].append((layer, z, x, y))
        pred = state["avail"].get(layer)
        if pred is None or not pred(z, x, y):
            return None
        return make_png(ramp(z, x, y, state["offset"].get(layer, 0.0)))

    monkeypatch.setattr(sources, "_http_get", fake_get)
    return state


def wm_center():
    return CX, CY


def pixel_size(path):
    return gdal.Open(path).GetGeoTransform()[1]


def warp_wm(paths, cx, cy, radius, px, **kw):
    wkt = sources._srs(3857).ExportToWkt()
    elev = raster.warp_window(paths, wkt, cx, cy, radius, px, progress_cb=lambda *a: 1, fill=False, **kw)
    return elev[0]


def truth(cx, cy, radius, px):
    n = int(round(2 * radius / px))
    X, Y = np.meshgrid(cx - radius + (np.arange(n) + 0.5) * px, cy + radius - (np.arange(n) + 0.5) * px)
    return 1000.0 + 0.3 * (X - CX) + 0.5 * (Y - CY)


quiet = lambda *a, **k: None  # noqa: E731


def test_without_fine_layer_the_coarser_layer_is_used_instead_of_no_coverage(fake_gsi, tmp_path):
    fake_gsi["avail"] = {"dem5a_png": lambda z, x, y: z == 15}
    paths = sources.fetch_gsi(LAT, LON, 200, str(tmp_path), 1.0, log=quiet)          # pixel de 1 m : zoom 17
    assert len(paths) == 1 and pixel_size(paths[0]) == pytest.approx(2 * sources.WM_ORIGIN / 2 ** 15 / 256)


def test_fallback_is_interpolated_so_a_linear_surface_is_reproduced(fake_gsi, tmp_path):
    fake_gsi["avail"] = {"dem5a_png": lambda z, x, y: z == 15}
    cx, cy = wm_center()
    paths = sources.fetch_gsi(LAT, LON, 150, str(tmp_path), 1.0, log=quiet)
    auto = warp_wm(paths, cx, cy, 100.0, 2.0)                                       # auto : bilinéaire
    avg = warp_wm(paths, cx, cy, 100.0, 2.0, resample="average")                     # ancien comportement : par blocs
    t = truth(cx, cy, 100.0, 2.0)
    assert np.abs(auto - t).max() < 0.05                                             # un plan est reproduit exactement
    assert np.abs(avg - t).max() > 0.5                                               # terrasses de 4,8 m sur une pente de 0,3 à 0,5


def test_fine_tiles_override_the_fallback_where_they_exist(fake_gsi, tmp_path):
    fake_gsi["avail"] = {"dem5a_png": lambda z, x, y: z == 15,
                         "dem1a_png": lambda z, x, y: z == 17 and x % 2 == 0}
    fake_gsi["offset"] = {"dem1a_png": 100.0}                                        # décalage pour distinguer les couches
    paths = sources.fetch_gsi(LAT, LON, 300, str(tmp_path), 1.0, log=quiet)
    assert len(paths) == 2 and pixel_size(paths[0]) > pixel_size(paths[1])           # grossier d'abord, fin en dernier
    out = warp_wm(paths, CX, CY, 280.0, 1.0, resample="bilinear")                   # fenêtre de 560 m : plusieurs tuiles de ~250 m
    diff = out - truth(CX, CY, 280.0, 1.0)
    fine = np.abs(diff - 100.0) < 0.3
    coarse = np.abs(diff) < 0.3
    assert fine.mean() > 0.15 and coarse.mean() > 0.15                               # les deux couches sont visibles
    assert (fine | coarse).mean() > 0.97                                             # et rien d'autre (pas de zone aberrante)


def test_no_coarse_request_when_every_fine_tile_exists(fake_gsi, tmp_path):
    fake_gsi["avail"] = {"dem1a_png": lambda z, x, y: z == 17}
    paths = sources.fetch_gsi(LAT, LON, 150, str(tmp_path), 1.0, log=quiet)
    assert len(paths) == 1 and {c[1] for c in fake_gsi["calls"]} == {17}


def test_chain_goes_from_zoom_15_to_zoom_14(fake_gsi, tmp_path):
    fake_gsi["avail"] = {"dem5a_png": lambda z, x, y: z == 15 and x % 2 == 0,        # zoom 15 incomplet
                         "dem_png": lambda z, x, y: z == 14}                          # zoom 14 partout
    paths = sources.fetch_gsi(LAT, LON, 900, str(tmp_path), 1.0, log=quiet)
    sizes = [pixel_size(p) for p in paths]
    assert len(paths) == 2 and sizes[0] > sizes[1] > 0                                # zoom 14 puis zoom 15
    assert {c[1] for c in fake_gsi["calls"]} >= {14, 15, 17}


def test_no_data_anywhere_is_still_no_coverage(fake_gsi, tmp_path):
    with pytest.raises(sources.NoCoverage):
        sources.fetch_gsi(LAT, LON, 150, str(tmp_path), 1.0, log=quiet)


def test_zoom_14_has_no_fallback_since_every_layer_is_already_tried(fake_gsi, tmp_path):
    fake_gsi["avail"] = {}
    with pytest.raises(sources.NoCoverage):
        sources.fetch_gsi(LAT, LON, 3000, str(tmp_path), 9.0, log=quiet)
    assert {c[1] for c in fake_gsi["calls"]} == {14}


# --- mélange de résolutions dans warp_window --------------------------------------------------------

def test_source_pixel_size_max_of_a_list(tmp_path):
    a, b = tmp_path / "a.tif", tmp_path / "b.tif"
    write_synthetic_tif(a, np.full((20, 20), 1.0), 600000.0, 5100000.0, 30.0, 32631)
    write_synthetic_tif(b, np.full((20, 20), 1.0), 600600.0, 5100000.0, 5.0, 32631)
    wkt = sources._srs(32631).ExportToWkt()
    assert raster.source_pixel_size_m([str(a), str(b)], wkt, 600300.0, 5099800.0, reduce="max") == pytest.approx(30.0, rel=1e-3)
    assert raster.source_pixel_size_m([str(a), str(b)], wkt, 600300.0, 5099800.0, reduce="min") == pytest.approx(5.0, rel=1e-3)
    with pytest.raises(ValueError):
        raster.source_pixel_size_m(str(a), wkt, 0, 0, reduce="mean")


def test_auto_uses_bilinear_when_any_source_of_a_mixed_list_is_enlarged(tmp_path, capsys):
    rng = np.random.default_rng(3)
    from scipy import ndimage
    coarse, fine = tmp_path / "coarse.tif", tmp_path / "fine.tif"
    write_synthetic_tif(coarse, 500 + ndimage.gaussian_filter(rng.normal(size=(60, 60)), 2) * 60, 600000.0, 5100000.0, 30.0, 32631)
    write_synthetic_tif(fine, 500 + ndimage.gaussian_filter(rng.normal(size=(80, 80)), 2) * 5, 600000.0, 5100000.0 - 1800.0 + 400.0, 5.0, 32631)
    wkt = sources._srs(32631).ExportToWkt()
    cx, cy = 600900.0, 5099100.0
    out = raster.warp_window([str(coarse), str(fine)], wkt, cx, cy, 800.0, 5.0, progress_cb=lambda *a: 1, fill=False)[0]
    assert "ré-échantillonnage bilinear" in capsys.readouterr().out
    region = out[: int(0.25 * out.shape[0])]                                           # zone couverte par la source grossière seule
    gy, gx = np.gradient(np.where(region < -1000, np.nan, region))
    assert np.nanmean((gx == 0) & (gy == 0)) < 0.01


# --- de bout en bout par la ligne de commande -----------------------------------------------------

def test_cli_gsi_source_uses_the_fallback_instead_of_failing(fake_gsi, monkeypatch, tmp_path, capsys):
    import check_stl
    import circle_dem_to_stl
    import mesh_export
    fake_gsi["avail"] = {"dem5a_png": lambda z, x, y: z == 15}                       # aucune tuile DEM1A au zoom 17
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--source", "gsi", "--lat", str(LAT), "--lon", str(LON), "--radius", "150",
                                      "--diameter-mm", "60", "--vexag", "1", "--base-mm", "2",
                                      "--print-spacing-mm", "0.2", "--output", str(tmp_path / "out.stl"),
                                      "--cache-dir", str(tmp_path / "cache")])
    circle_dem_to_stl.main()
    out = capsys.readouterr().out
    assert "repli de résolution au zoom 15" in out and "ré-échantillonnage bilinear" in out
    tris = mesh_export.read_triangles(str(tmp_path / "out.stl"))
    res = {n: ok for n, ok, _ in check_stl.check(tris, 60.0, 2.0, tol_mm=1.5)}
    assert all(res.values()), res
    # La rampe est définie en mètres Web Mercator, étirés de 1/cos(latitude) par rapport aux mètres au sol :
    # pente au sol = hypot(0,3 ; 0,5) / cos(lat). Disque de 300 m à l'échelle 0,2 mm/m : dénivelé attendu en mm.
    expected = float(np.hypot(0.3, 0.5)) / np.cos(np.radians(LAT)) * 300.0 * 0.2
    assert tris[:, :, 2].max() == pytest.approx(expected, rel=0.02)
