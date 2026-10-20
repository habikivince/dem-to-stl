# SPDX-License-Identifier: GPL-3.0-or-later
"""Bout en bout : fichier source -> circle_dem_to_stl.main() -> STL réellement écrit puis relu et mesuré.

Chaque cas ci-dessous produit un vrai fichier STL, dont on vérifie la géométrie (étanchéité, orientation des faces,
dimensions), le relief vertical (comparé au gradient connu du terrain synthétique : une unité ou une échelle fausse,
un décalage oublié ou une bande mal choisie le changent) et l'orientation (le sommet doit être au nord-est)."""
import sys
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fmt_helpers import (X0, Y0, assert_stl_ok, cli, hgt_function, plane_array, scientific_dir, srs_wkt,  # noqa: E402
                         utm_gradient, write_asc, write_hgt, write_raster)

GT = (X0, 10.0, 0.0, Y0 + 1000.0, 0.0, -10.0)
CENTRE = dict(cx=X0 + 500.0, cy=Y0 + 500.0)
GRADIENT = (0.05, 0.02)                                   # pente du plan de référence (voir fmt_helpers)


@pytest.fixture(scope="module")
def sci(tmp_path_factory):
    return scientific_dir(tmp_path_factory)


def go(monkeypatch, tmp_path, files, extra=(), **kw):
    kw.setdefault("center", CENTRE)
    out = cli(monkeypatch, tmp_path, files, extra, **kw)
    assert out.exists()
    return out


# --- GeoTIFF ------------------------------------------------------------------------------------------------------

def test_geotiff_float32(monkeypatch, tmp_path):
    f = write_raster(tmp_path / "a.tif", plane_array((100, 100), GT), GT, nodata=-9999.0)
    assert_stl_ok(go(monkeypatch, tmp_path, [f]), 300.0, GRADIENT)


def test_geotiff_int16_with_scale_offset_and_feet(monkeypatch, tmp_path):
    z_ft = plane_array((100, 100), GT, np.float64) / 0.3048
    raw = np.rint((z_ft - 3000.0) / 0.1).astype(np.int16)
    f = write_raster(tmp_path / "p.tif", raw, GT, nodata=-32768.0, scale=0.1, offset=3000.0, unit="ft")
    assert_stl_ok(go(monkeypatch, tmp_path, [f]), 300.0, GRADIENT)           # relief en mètres, pas en pieds


def test_geotiff_unit_option_corrects_a_file_that_does_not_declare_its_unit(monkeypatch, tmp_path):
    f = write_raster(tmp_path / "ft.tif", (plane_array((100, 100), GT, np.float64) / 0.3048).astype(np.float32), GT)
    assert_stl_ok(go(monkeypatch, tmp_path, [f], ["--vertical-unit", "ft"]), 300.0, GRADIENT)


def test_geotiff_multiband_band_two(monkeypatch, tmp_path):
    z = plane_array((100, 100), GT)
    f = write_raster(tmp_path / "mb.tif", np.full_like(z, 777.0), GT, nodata=-9999.0, extra_bands=(z, z + 50))
    assert_stl_ok(go(monkeypatch, tmp_path, [f], ["--band", "2"]), 300.0, GRADIENT)   # la bande 1 est plate


def test_geotiff_tiled_deflate_float64_bigtiff(monkeypatch, tmp_path):
    f = write_raster(tmp_path / "b.tif", plane_array((100, 100), GT, np.float64), GT, nodata=-9999.0,
                     options=["BIGTIFF=YES", "COMPRESS=DEFLATE", "TILED=YES", "BLOCKXSIZE=32", "BLOCKYSIZE=32"])
    assert_stl_ok(go(monkeypatch, tmp_path, [f]), 300.0, GRADIENT)


def test_geotiff_with_a_hole_gives_a_watertight_model(monkeypatch, tmp_path):
    z = plane_array((100, 100), GT)
    z[48:52, 48:52] = -9999.0
    f = write_raster(tmp_path / "h.tif", z, GT, nodata=-9999.0)
    assert_stl_ok(go(monkeypatch, tmp_path, [f]), 300.0, GRADIENT)


# --- ASC ----------------------------------------------------------------------------------------------------------

def test_asc_with_prj(monkeypatch, tmp_path):
    f = write_asc(tmp_path / "a.asc", plane_array((100, 100), GT, np.float64), X0, Y0, 10.0, prj=srs_wkt(2056))
    assert_stl_ok(go(monkeypatch, tmp_path, [f]), 300.0, GRADIENT)


def test_asc_without_prj_needs_source_crs(monkeypatch, tmp_path):
    f = write_asc(tmp_path / "n.asc", plane_array((100, 100), GT, np.float64), X0, Y0, 10.0)
    with pytest.raises(SystemExit) as e:
        go(monkeypatch, tmp_path, [f])
    assert "--source-crs" in str(e.value.code)
    assert_stl_ok(go(monkeypatch, tmp_path, [f], ["--source-crs", "EPSG:2056"], name="ok.stl"), 300.0, GRADIENT)


def test_truncated_asc_stops_with_a_message(monkeypatch, tmp_path):
    f = write_asc(tmp_path / "t.asc", plane_array((100, 100), GT, np.float64), X0, Y0, 10.0, prj=srs_wkt(2056))
    text = Path(f).read_text()
    Path(f).write_text(text[:len(text) // 2])
    with pytest.raises(SystemExit) as e:
        go(monkeypatch, tmp_path, [f])
    assert "tronqué" in str(e.value.code) and not (tmp_path / "out.stl").exists()


# --- SRTM ---------------------------------------------------------------------------------------------------------

def test_srtm_tile_and_two_tiles(monkeypatch, tmp_path):
    a = write_hgt(tmp_path / "N46E007.hgt", 7, 46)
    b = write_hgt(tmp_path / "N46E008.hgt", 8, 46)
    center = dict(lat=46.5, lon=7.5)
    grad = utm_gradient(hgt_function, 32632, *_utm(46.5, 7.5))
    assert_stl_ok(go(monkeypatch, tmp_path, [a], center=center, radius=3000.0), 3000.0, grad)
    center = dict(lat=46.5, lon=8.0)                                       # à cheval sur la couture
    grad = utm_gradient(hgt_function, 32632, *_utm(46.5, 8.0))
    assert_stl_ok(go(monkeypatch, tmp_path, [a, b], center=center, radius=3000.0, name="two.stl"), 3000.0, grad)


def _utm(lat, lon):
    from osgeo import osr
    s, d = osr.SpatialReference(), osr.SpatialReference()
    s.ImportFromEPSG(4326); d.ImportFromEPSG(32632)
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER); d.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    x, y, _ = osr.CoordinateTransformation(s, d).TransformPoint(lon, lat)
    return x, y


# --- GeoPackage ---------------------------------------------------------------------------------------------------

def _gpkg(tmp_path, arr, table, append=False, name="d.gpkg"):
    src = write_raster(tmp_path / f"{table}.tif", arr, GT, nodata=-9999.0)
    opts = ["TILE_FORMAT=TIFF", f"RASTER_TABLE={table}"] + (["APPEND_SUBDATASET=YES"] if append else [])
    out = gdal.Translate(str(tmp_path / name), src, format="GPKG", creationOptions=opts)
    out = None  # noqa: F841
    return str(tmp_path / name)


def test_geopackage_single_table(monkeypatch, tmp_path):
    f = _gpkg(tmp_path, plane_array((100, 100), GT), "mnt")
    assert_stl_ok(go(monkeypatch, tmp_path, [f]), 300.0, GRADIENT)


def test_geopackage_two_tables_need_layer_and_use_the_chosen_one(monkeypatch, tmp_path):
    f = _gpkg(tmp_path, np.full((100, 100), 1000.0, np.float32), "mnt_a")
    f = _gpkg(tmp_path, plane_array((100, 100), GT), "mnt_b", append=True)
    with pytest.raises(SystemExit) as e:
        go(monkeypatch, tmp_path, [f])
    assert "--layer" in str(e.value.code)
    assert_stl_ok(go(monkeypatch, tmp_path, [f], ["--layer", "mnt_b"], name="b.stl"), 300.0, GRADIENT)


# --- NetCDF / HDF5 ------------------------------------------------------------------------------------------------

def _sci_z(lon, lat):
    return 1000.0 + 800.0 * (lon - 7.0) + 1200.0 * (lat - 46.0)


@pytest.mark.parametrize("name, extra", [("simple.nc", []), ("packed.nc", []), ("feet.nc", []), ("time.nc", ["--band", "2"]),
                                         ("netcdf4_as.h5", [])])
def test_netcdf_and_hdf5(monkeypatch, tmp_path, sci, name, extra):
    grad = utm_gradient(_sci_z, 32632, *_utm(46.1, 7.15))
    out = go(monkeypatch, tmp_path, [sci / name], extra, center=dict(lat=46.1, lon=7.15), radius=3000.0)
    assert_stl_ok(out, 3000.0, grad)


def test_netcdf_with_time_dimension_is_refused_without_band(monkeypatch, tmp_path, sci):
    with pytest.raises(SystemExit) as e:
        go(monkeypatch, tmp_path, [sci / "time.nc"], center=dict(lat=46.1, lon=7.15), radius=3000.0)
    assert "--band N" in str(e.value.code) and not (tmp_path / "out.stl").exists()


def test_generic_hdf5_is_reported_as_unsupported(monkeypatch, tmp_path, sci):
    with pytest.raises(SystemExit) as e:
        go(monkeypatch, tmp_path, [sci / "plain.hdf5"], center=dict(lat=46.1, lon=7.15), radius=3000.0)
    assert "HDF5 non géoréférencé" in str(e.value.code)


# --- Compatibilité avec les formes (--sides) et les options existantes ----------------------------------------------

@pytest.mark.parametrize("sides", [4, 6])
def test_regular_polygon_mode_works_with_a_new_format(monkeypatch, tmp_path, sides):
    f = write_asc(tmp_path / "a.asc", plane_array((100, 100), GT, np.float64), X0, Y0, 10.0, prj=srs_wkt(2056))
    out = go(monkeypatch, tmp_path, [f], ["--sides", str(sides), "--circle-mode", "inscribed"])
    assert_stl_ok(out, 300.0, GRADIENT, polygon=sides)


def test_dry_run_with_a_new_format_writes_nothing(monkeypatch, tmp_path):
    import circle_dem_to_stl
    f = write_asc(tmp_path / "a.asc", plane_array((100, 100), GT, np.float64), X0, Y0, 10.0, prj=srs_wkt(2056))
    out = tmp_path / "dry.stl"
    monkeypatch.setattr(sys, "argv", ["x", "--input", f, "--cx", str(X0 + 500), "--cy", str(Y0 + 500), "--radius", "300",
                                      "--diameter-mm", "60", "--vexag", "1", "--base-mm", "2", "--output", str(out),
                                      "--dry-run"])
    circle_dem_to_stl.main()
    assert not out.exists()


def test_unknown_extension_is_still_opened_by_content_and_reported_if_unreadable(monkeypatch, tmp_path):
    (tmp_path / "x.dat").write_bytes(b"nothing useful here" * 10)
    with pytest.raises(SystemExit) as e:
        go(monkeypatch, tmp_path, [tmp_path / "x.dat"])
    assert "illisible" in str(e.value.code)
