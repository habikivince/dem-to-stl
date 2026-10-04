# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests du complément Copernicus (--fill-with-copernicus) : comblement des trous du disque
avec une source secondaire grossière, recalée en altitude. Aucun réseau."""
import sys
from pathlib import Path

import numpy as np
import pytest
from osgeo import osr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import circle_dem_to_stl  # noqa: E402
import download  # noqa: E402
import raster  # noqa: E402
import sources  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402

quiet = lambda *a, **k: None  # noqa: E731
X0, Y0 = 2600000.0, 1204000.0       # coin NO de la zone de test (EPSG:2056)


def terrain(x, y):
    return 1000.0 + 0.05 * (x - X0) + 0.02 * (y - 1200000.0)


def write_plane(path, x_min, x_max, pixel, offset=0.0, epsg=2056):
    n_x, n_y = int((x_max - x_min) / pixel), int(4000 / pixel)
    xs = x_min + (np.arange(n_x) + 0.5) * pixel
    ys = Y0 - (np.arange(n_y) + 0.5) * pixel
    xx, yy = np.meshgrid(xs, ys)
    write_synthetic_tif(path, (terrain(xx, yy) + offset).astype(np.float32), x_min, Y0, pixel, epsg)


def run_warp(inp, radius=1500.0, pixel=8.0):
    wkt, cx, cy = raster.resolve_target_crs_and_center(inp, None, None, X0 + 2000.0, 1202000.0)
    elev, xs, ys, nodata = raster.warp_window(inp, wkt, cx, cy, radius, pixel, progress_cb=lambda *a: 1)
    inside = (xs[None, :] - cx) ** 2 + (ys[:, None] - cy) ** 2 <= radius ** 2
    return elev, xs, ys, nodata, inside, (cx, cy)


def layered(tmp_path, main_path, fill_path):
    calls = []
    inputs = sources.LayeredInputs([str(main_path)])

    def loader():
        calls.append(1)
        return str(fill_path)
    inputs.fill_loader = loader
    return inputs, calls


def test_gap_is_filled_and_offset_removed(tmp_path):
    write_plane(tmp_path / "main.tif", X0, X0 + 2000.0, 2.0)                 # moitié ouest du disque
    write_plane(tmp_path / "fill.tif", X0, X0 + 4000.0, 30.0, offset=2.0)    # Copernicus-like, +2 m
    inputs, calls = layered(tmp_path, tmp_path / "main.tif", tmp_path / "fill.tif")
    elev, xs, ys, nodata, inside, _ = run_warp(inputs)
    assert calls == [1]
    assert not np.any((elev == nodata) & inside)                             # disque plein
    expected = terrain(xs[None, :], ys[:, None])
    assert np.abs(elev - expected)[inside].max() < 0.5                       # +2 m retiré sur les cellules ajoutées


def test_no_secondary_download_when_disc_is_complete(tmp_path):
    write_plane(tmp_path / "main.tif", X0, X0 + 4000.0, 2.0)
    write_plane(tmp_path / "fill.tif", X0, X0 + 4000.0, 30.0)
    inputs, calls = layered(tmp_path, tmp_path / "main.tif", tmp_path / "fill.tif")
    run_warp(inputs)
    assert calls == []


def test_secondary_failure_leaves_disc_partial_without_crash(tmp_path):
    write_plane(tmp_path / "main.tif", X0, X0 + 2000.0, 2.0)
    inputs = sources.LayeredInputs([str(tmp_path / "main.tif")])

    def broken():
        raise RuntimeError("réseau coupé")
    inputs.fill_loader = broken
    elev, xs, ys, nodata, inside, _ = run_warp(inputs)   # ne doit pas lever
    assert np.any((elev == nodata) & inside)


def test_incompatible_altitude_reference_aborts_fill(tmp_path):
    write_plane(tmp_path / "main.tif", X0, X0 + 2000.0, 2.0)
    write_plane(tmp_path / "fill.tif", X0, X0 + 4000.0, 30.0, offset=120.0)
    inputs, calls = layered(tmp_path, tmp_path / "main.tif", tmp_path / "fill.tif")
    elev, xs, ys, nodata, inside, _ = run_warp(inputs)
    assert calls == [1]
    assert np.any((elev == nodata) & inside)          # pas de comblement avec un décalage de 120 m


def test_fill_works_with_geographic_secondary_and_main_without_nodata_tag(tmp_path):
    """Secondaire en EPSG:4326 (cas réel de Copernicus) ; principal sans nodata déclaré."""
    ct = osr.CoordinateTransformation(sources._srs(2056), sources._srs(4326))
    lon0, lat0, _ = ct.TransformPoint(X0 - 500.0, Y0 + 500.0)
    lon1, lat1, _ = ct.TransformPoint(X0 + 4500.0, Y0 - 4500.0)
    step = 0.0003
    n_x, n_y = int((lon1 - lon0) / step), int((lat0 - lat1) / step)
    lons = lon0 + (np.arange(n_x) + 0.5) * step
    lats = lat0 - (np.arange(n_y) + 0.5) * step
    back = osr.CoordinateTransformation(sources._srs(4326), sources._srs(2056))
    # plan linéaire en x,y : on évalue le terrain aux coins et on interpole (suffisant pour le test)
    xs_ = np.array([back.TransformPoint(lo, lats[0])[0] for lo in lons])
    ys_ = np.array([back.TransformPoint(lons[0], la)[1] for la in lats])
    xx, yy = np.meshgrid(xs_, ys_)
    write_synthetic_tif(tmp_path / "cop.tif", (terrain(xx, yy) + 1.0).astype(np.float32), lon0, lat0, step, 4326)
    write_plane(tmp_path / "main.tif", X0, X0 + 2000.0, 2.0)
    import osgeo.gdal as gdal
    ds = gdal.Open(str(tmp_path / "main.tif"), gdal.GA_Update)
    ds.GetRasterBand(1).DeleteNoDataValue()
    ds = None
    inputs, calls = layered(tmp_path, tmp_path / "main.tif", tmp_path / "cop.tif")
    elev, xs, ys, nodata, inside, _ = run_warp(inputs)
    assert nodata == raster.FILL_NODATA and calls == [1]
    assert not np.any((elev == nodata) & inside)
    assert np.abs(elev - terrain(xs[None, :], ys[:, None]))[inside].max() < 2.0


def test_cli_fill_with_input_and_without_latlon(tmp_path, monkeypatch):
    lat, lon = 46.0, 7.6
    x, y = sources._to_epsg(2056, lat, lon)
    main = tmp_path / "main.tif"
    write_synthetic_tif(main, np.full((400, 800), 1500.0), x - 800.0, y + 400.0, 2.0, 2056)   # moitié ouest
    cop = tmp_path / "cop.tif"
    write_synthetic_tif(cop, np.full((200, 200), 1502.0), lon - 0.02, lat + 0.02, 0.0002, 4326)
    called = []
    monkeypatch.setattr(download, "download_copernicus", lambda *a, **k: called.append(1) or [str(cop)])
    out = tmp_path / "o.stl"
    argv = ["x", "--input", str(main), "--lat", str(lat), "--lon", str(lon), "--radius", "500",
            "--diameter-mm", "50", "--vexag", "1", "--base-mm", "2", "--print-spacing-mm", "0.5",
            "--fill-with-copernicus", "--cache-dir", str(tmp_path / "c"), "--output", str(out)]
    monkeypatch.setattr(sys, "argv", argv)
    circle_dem_to_stl.main()
    assert called and out.stat().st_size > 84

    monkeypatch.setattr(sys, "argv", ["x", "--input", str(main), "--cx", str(x), "--cy", str(y), "--radius", "500",
                                      "--fill-with-copernicus", "--output", str(tmp_path / "p.stl")])
    with pytest.raises(SystemExit):
        circle_dem_to_stl.main()
