# SPDX-License-Identifier: GPL-3.0-or-later
"""Lecteurs de formats sources raster : GeoTIFF, ESRI ASCII Grid, SRTM (.hgt).
(GeoPackage, NetCDF et HDF5 : voir test_input_containers.py.)

Chaque test fabrique de petits fichiers synthétiques et vérifie des VALEURS d'altitude, la taille et l'orientation de
la grille : un fichier qui « s'ouvre » ne prouve rien. Le terrain est un plan incliné (voir fmt_helpers) que le
ré-échantillonnage reproduit exactement ; toute erreur vient donc du lecteur."""
import sys
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal, osr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import formats  # noqa: E402
import raster  # noqa: E402
from fmt_helpers import (X0, Y0, hgt_array, hgt_function, plane, plane_array, prepared, srs_wkt,  # noqa: E402
                         to_geographic, window, write_asc, write_hgt, write_raster)

GT = (X0, 10.0, 0.0, Y0 + 1000.0, 0.0, -10.0)           # 100 x 100 pixels de 10 m, LV95
CENTRE = dict(cx=X0 + 500.0, cy=Y0 + 500.0, radius=300.0, pixel=10.0)   # grille de travail alignée sur les pixels


def truth(xs, ys):
    return plane(xs[None, :], ys[:, None])


def run(names, **kw):
    args = dict(CENTRE)
    args.update(kw)
    return window(names, **args)


# ---------------------------------------------------------------------------------------------------------------
# GeoTIFF
# ---------------------------------------------------------------------------------------------------------------

def test_geotiff_float32_is_read_unchanged_with_correct_orientation(tmp_path):
    path = write_raster(tmp_path / "a.tif", plane_array((100, 100), GT), GT, nodata=-9999.0)
    kind, names = prepared(path)
    assert kind == "raster" and names == [path]                       # aucun VRT : comportement historique
    elev, xs, ys, nodata, _ = run(names)
    assert elev.shape == (60, 60) and nodata == -9999.0
    assert np.abs(elev - truth(xs, ys)).max() < 1e-3
    assert elev[0].mean() > elev[-1].mean()                           # le nord est en haut, il monte vers le nord
    assert elev[:, -1].mean() > elev[:, 0].mean()                     # et vers l'est


@pytest.mark.parametrize("options", [
    ["COMPRESS=DEFLATE", "PREDICTOR=2"], ["COMPRESS=LZW"], ["COMPRESS=PACKBITS"],
    ["COMPRESS=DEFLATE", "TILED=YES", "BLOCKXSIZE=32", "BLOCKYSIZE=32"], ["BIGTIFF=YES"],
    ["COMPRESS=ZSTD"], ["COMPRESS=LERC"],
], ids=lambda o: "+".join(o))
def test_geotiff_compressed_variants_give_identical_values(tmp_path, options):
    codec = next((o.split("=")[1] for o in options if o.startswith("COMPRESS=")), None)
    listing = gdal.GetDriverByName("GTiff").GetMetadataItem("DMD_CREATIONOPTIONLIST") or ""
    if codec in ("ZSTD", "LERC") and codec not in listing:
        pytest.skip(f"cette build de GDAL n'inclut pas le codec {codec}")
    arr = plane_array((100, 100), GT)
    path = write_raster(tmp_path / "c.tif", arr, GT, nodata=-9999.0, options=options)
    elev, xs, ys, _, _ = run(prepared(path)[1])
    assert np.abs(elev - truth(xs, ys)).max() < 1e-3


@pytest.mark.parametrize("dtype", ["int16", "uint16", "int32", "float32", "float64"])
def test_geotiff_numeric_types_keep_their_values(tmp_path, dtype):
    arr = np.rint(plane_array((100, 100), GT, np.float64)).astype(dtype)
    path = write_raster(tmp_path / f"{dtype}.tif", arr, GT)
    elev, xs, ys, _, _ = run(prepared(path)[1])
    expected = arr[50 - 30:50 + 30, 50 - 30:50 + 30].astype(np.float64)       # fenêtre 60 x 60 centrée
    assert elev.shape == expected.shape and np.abs(elev - expected).max() < 1e-6


def test_integer_source_is_resampled_in_floating_point_not_rounded(tmp_path):
    """Damier 0/1 en Int16, moyenné par blocs de 2 x 2 : 0,5 partout. Avec des entiers, GDAL arrondirait à 0 ou 1."""
    i, j = np.indices((100, 100))
    path = write_raster(tmp_path / "check.tif", ((i + j) % 2).astype(np.int16), GT)
    elev, *_ = run(prepared(path)[1], pixel=20.0, resample="average")
    assert elev.shape == (30, 30) and np.allclose(elev, 0.5)


def test_negative_and_decimal_altitudes_are_preserved(tmp_path):
    arr = (plane_array((100, 100), GT, np.float64) - 1400.37).astype(np.float32)     # de -400 à -300, décimales
    path = write_raster(tmp_path / "neg.tif", arr, GT, nodata=-9999.0)
    elev, xs, ys, _, _ = run(prepared(path)[1])
    assert elev.max() < -300 and elev.min() > -401
    assert np.abs(elev - (truth(xs, ys) - 1400.37)).max() < 2e-3


def test_declared_scale_and_offset_are_applied(tmp_path):
    raw = np.rint((plane_array((100, 100), GT, np.float64) - 500.0) / 0.1).astype(np.int16)
    path = write_raster(tmp_path / "packed.tif", raw, GT, scale=0.1, offset=500.0, unit="m")
    kind, names = prepared(path)
    assert names != [path]                                                           # normalisé par un VRT
    elev, xs, ys, _, _ = run(names)
    assert np.abs(elev - truth(xs, ys)).max() < 0.06                                 # précision du quantum 0,1 m
    assert elev.min() > 1000 - 1 and elev.max() < 1100


def test_nodata_never_becomes_zero_and_is_filled_only_when_asked(tmp_path, capsys):
    arr = plane_array((100, 100), GT)
    arr[48:52, 48:52] = -9999.0                                                      # trou de 4 x 4 pixels au centre
    path = write_raster(tmp_path / "hole.tif", arr, GT, nodata=-9999.0)
    names = prepared(path)[1]
    wkt, cx, cy = raster.resolve_target_crs_and_center(names, None, None, CENTRE["cx"], CENTRE["cy"])
    raw, xs, ys, nodata = raster.warp_window(names, wkt, cx, cy, 300.0, 10.0, progress_cb=lambda *a: 1, fill=False)
    assert (raw == nodata).sum() == 16 and not np.any(raw == 0.0)                    # trou conservé tel quel
    assert np.abs(raw[raw != nodata] - truth(xs, ys)[raw != nodata]).max() < 1e-3
    filled, *_ = run(names)                                                          # chemin du programme : comblé...
    assert not np.any(filled == nodata) and np.abs(filled - truth(xs, ys)).max() < 0.5
    assert "16 cellule(s) sans donnée comblée(s)" in capsys.readouterr().out         # ...et le dit


@pytest.mark.parametrize("declared", [None, -9999.0])
def test_area_outside_the_source_footprint_is_nodata_not_zero(tmp_path, declared):
    """La source ne couvre que la moitié ouest de la fenêtre. Sans nodata déclaré, GDAL remplit l'extérieur de
    zéros : une paroi à l'altitude zéro serait fabriquée en silence."""
    gt = (X0, 10.0, 0.0, Y0 + 1000.0, 0.0, -10.0)
    arr = plane_array((100, 50), gt)
    path = write_raster(tmp_path / "half.tif", arr, gt, nodata=declared)
    names = prepared(path)[1]
    wkt, cx, cy = raster.resolve_target_crs_and_center(names, None, None, CENTRE["cx"], CENTRE["cy"])
    elev, xs, ys, nodata = raster.warp_window(names, wkt, cx, cy, 300.0, 10.0, progress_cb=lambda *a: 1, fill=False)
    east = xs > X0 + 500.0
    assert nodata == -9999.0 and np.all(elev[:, east] == nodata) and np.all(elev[:, ~east] != nodata)
    assert not np.any(elev == 0.0)


@pytest.mark.parametrize("declared", [True, False], ids=["nan-declared", "nan-undeclared"])
def test_nan_pixels_are_nodata_never_altitudes(tmp_path, declared):
    arr = plane_array((100, 100), GT)
    arr[50:53, 50:53] = np.nan
    path = write_raster(tmp_path / "nan.tif", arr, GT, nodata=float("nan") if declared else None)
    names = prepared(path)[1]
    wkt, cx, cy = raster.resolve_target_crs_and_center(names, None, None, CENTRE["cx"], CENTRE["cy"])
    raw, xs, ys, nodata = raster.warp_window(names, wkt, cx, cy, 300.0, 10.0, progress_cb=lambda *a: 1, fill=False)
    assert np.isfinite(raw).all() and np.isfinite(nodata)
    assert (raw == nodata).sum() == 9 and np.abs(raw[raw != nodata] - truth(xs, ys)[raw != nodata]).max() < 1e-3


# --- plusieurs bandes ------------------------------------------------------------------------------------------

def _three_bands(tmp_path):
    base = plane_array((100, 100), GT)
    return write_raster(tmp_path / "mb.tif", base, GT, nodata=-9999.0, extra_bands=(base + 100, base + 200))


def test_multiband_file_is_refused_without_band_and_lists_the_bands(tmp_path):
    with pytest.raises(raster.InputError, match=r"3 bandes.*--band N") as e:
        prepared(_three_bands(tmp_path))
    assert "bande 1" in str(e.value) and "bande 3" in str(e.value)


def test_multiband_file_uses_the_band_asked_for(tmp_path):
    elev, xs, ys, _, _ = run(prepared(_three_bands(tmp_path), band=2)[1])
    assert np.abs(elev - (truth(xs, ys) + 100)).max() < 1e-3                         # bande 2, pas la bande 1


def test_band_out_of_range_is_refused(tmp_path):
    with pytest.raises(raster.InputError, match="--band 4.*3 bande"):
        prepared(_three_bands(tmp_path), band=4)


def test_alpha_band_is_not_a_variable(tmp_path):
    arr = plane_array((100, 100), GT)
    path = write_raster(tmp_path / "alpha.tif", arr, GT, nodata=-9999.0, extra_bands=(np.full((100, 100), 255, np.float32),))
    ds = gdal.Open(path, gdal.GA_Update)
    ds.GetRasterBand(2).SetColorInterpretation(gdal.GCI_AlphaBand)
    ds = None
    elev, xs, ys, _, _ = run(prepared(path)[1])
    assert np.abs(elev - truth(xs, ys)).max() < 1e-3


def test_colour_image_is_recognised_as_not_a_dem(tmp_path):
    rgb = np.full((100, 100), 128, np.uint8)
    path = write_raster(tmp_path / "rgb.tif", rgb, GT, extra_bands=(rgb, rgb))
    ds = gdal.Open(path, gdal.GA_Update)
    for i, ci in enumerate((gdal.GCI_RedBand, gdal.GCI_GreenBand, gdal.GCI_BlueBand), 1):
        ds.GetRasterBand(i).SetColorInterpretation(ci)
    ds = None
    with pytest.raises(raster.InputError, match="image couleur"):
        prepared(path)


# --- CRS et unités ---------------------------------------------------------------------------------------------

def test_missing_crs_is_refused_with_the_way_out_and_never_invented(tmp_path):
    path = write_raster(tmp_path / "nocrs.tif", plane_array((100, 100), GT), GT, epsg=None, nodata=-9999.0)
    with pytest.raises(raster.InputError, match="aucun système de coordonnées.*--source-crs"):
        prepared(path)


def test_source_crs_completes_a_file_without_crs(tmp_path):
    path = write_raster(tmp_path / "nocrs.tif", plane_array((100, 100), GT), GT, epsg=None, nodata=-9999.0)
    elev, xs, ys, _, _ = run(prepared(path, source_crs="EPSG:2056")[1])
    assert np.abs(elev - truth(xs, ys)).max() < 1e-3


def test_source_crs_never_replaces_a_declared_crs(tmp_path):
    path = write_raster(tmp_path / "crs.tif", plane_array((100, 100), GT), GT, nodata=-9999.0)
    with pytest.raises(raster.InputError, match="déclare déjà un CRS"):
        prepared(path, source_crs="EPSG:32632")


def test_unknown_source_crs_is_refused(tmp_path):
    path = write_raster(tmp_path / "nocrs.tif", plane_array((100, 100), GT), GT, epsg=None)
    with pytest.raises(raster.InputError, match="CRS non reconnu"):
        prepared(path, source_crs="EPSG:99999999")


def test_missing_georeferencing_is_refused(tmp_path):
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(str(tmp_path / "nogeo.tif"), 10, 10, 1, gdal.GDT_Float32)
    ds.GetRasterBand(1).Fill(5.0)
    ds = None
    with pytest.raises(raster.InputError, match="pas de géoréférencement"):
        prepared(tmp_path / "nogeo.tif")


def test_degrees_declared_for_metre_coordinates_is_refused(tmp_path):
    path = write_raster(tmp_path / "wrong.tif", plane_array((100, 100), GT), GT, epsg=4326)   # 2 600 000 « degrés »
    with pytest.raises(raster.InputError, match="impossible en degrés"):
        prepared(path)


def test_metre_crs_with_degree_like_extent_warns(tmp_path):
    gt = (7.0, 0.001, 0.0, 46.1, 0.0, -0.001)
    log = []
    path = write_raster(tmp_path / "deg.tif", plane_array((100, 100), gt), gt, epsg=2056)
    prepared(path, log=log)
    assert any("allure de degrés" in m for m in log)


def test_undeclared_vertical_unit_is_assumed_metres_with_a_warning(tmp_path):
    log = []
    path = write_raster(tmp_path / "nounit.tif", plane_array((100, 100), GT), GT, nodata=-9999.0)
    names = prepared(path, log=log)[1]
    assert any("non déclarée, supposée en MÈTRES" in m for m in log)
    elev, xs, ys, _, _ = run(names)
    assert np.abs(elev - truth(xs, ys)).max() < 1e-3


@pytest.mark.parametrize("unit, factor", [("ft", 0.3048), ("foot", 0.3048), ("US survey foot", 1200 / 3937),
                                           ("cm", 0.01), ("m", 1.0), ("metre", 1.0)])
def test_declared_vertical_unit_is_converted_to_metres(tmp_path, unit, factor):
    path = write_raster(tmp_path / "u.tif", (plane_array((100, 100), GT, np.float64) / factor).astype(np.float32),
                        GT, nodata=-9999.0, unit=unit)
    elev, xs, ys, _, _ = run(prepared(path)[1])
    assert np.abs(elev - truth(xs, ys)).max() < 0.01


def test_vertical_unit_option_overrides_a_missing_or_wrong_unit(tmp_path):
    log = []
    values = (plane_array((100, 100), GT, np.float64) / 0.3048).astype(np.float32)
    path = write_raster(tmp_path / "ft.tif", values, GT, nodata=-9999.0, unit="m")             # unité déclarée FAUSSE
    elev, xs, ys, _, _ = run(prepared(path, vertical_unit="ft", log=log)[1])
    assert np.abs(elev - truth(xs, ys)).max() < 0.01
    assert any("remplacée par --vertical-unit ft" in m for m in log)


def test_non_length_unit_is_refused(tmp_path):
    path = write_raster(tmp_path / "k.tif", plane_array((100, 100), GT), GT, unit="K")
    with pytest.raises(raster.InputError, match="n'est pas une unité de longueur"):
        prepared(path)


def test_unknown_vertical_unit_option_is_refused(tmp_path):
    path = write_raster(tmp_path / "a.tif", plane_array((100, 100), GT), GT)
    with pytest.raises(raster.InputError, match="--vertical-unit"):
        prepared(path, vertical_unit="fathom")


def test_crs_in_us_feet_is_reprojected_to_metres_and_values_are_right(tmp_path):
    """Source en pieds US (Colorado Nord, EPSG:2231) : le pipeline raisonne en mètres, il passe donc par l'UTM.
    Lire des pieds comme des mètres donnerait un modèle 3,28 fois trop petit ; l'altitude est aussi en pieds US."""
    src = osr.SpatialReference(); src.ImportFromEPSG(4326); src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst = osr.SpatialReference(); dst.ImportFromEPSG(2231); dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    x_ft, y_ft, _ = osr.CoordinateTransformation(src, dst).TransformPoint(-105.0, 40.0)
    gt = (x_ft - 500.0, 10.0, 0.0, y_ft + 500.0, 0.0, -10.0)                    # 100 x 100 pixels de 10 pieds US
    xs_c, ys_c = gt[0] + (np.arange(100) + 0.5) * 10.0, gt[3] - (np.arange(100) + 0.5) * 10.0

    def z_ft(x, y):
        return 5000.0 + 0.2 * (x - gt[0]) + 0.1 * (y - (gt[3] - 1000.0))

    arr = z_ft(xs_c[None, :], ys_c[:, None]).astype(np.float32)
    path = write_raster(tmp_path / "ft.tif", arr, gt, epsg=2231, nodata=-9999.0, unit="US survey foot")
    names = prepared(path)[1]
    elev, xs, ys, _, wkt = window(names, cx=x_ft, cy=y_ft, radius=100.0, pixel=3.0)
    assert "UTM zone 13N" in wkt
    to_src = osr.CoordinateTransformation(osr.SpatialReference(wkt=wkt), dst)
    gx, gy = np.meshgrid(xs, ys)
    pts = np.array(to_src.TransformPoints(np.column_stack([gx.ravel(), gy.ravel()])))
    expected_m = (z_ft(pts[:, 0], pts[:, 1]) * 1200.0 / 3937.0).reshape(gx.shape)
    assert np.abs(elev - expected_m).max() < 0.1
    assert elev.shape == (67, 67)                       # int((2 x 100 m + 1,5 m) / 3 m) cellules : des mètres, pas des pieds


def test_compound_crs_drops_the_vertical_part_without_changing_values(tmp_path):
    log = []
    comp = osr.SpatialReference(); comp.SetFromUserInput("EPSG:2056+5728")
    path = write_raster(tmp_path / "compound.tif", plane_array((100, 100), GT), GT, wkt=comp.ExportToWkt(), nodata=-9999.0)
    names = prepared(path, log=log)[1]
    assert any("AUCUN changement de datum vertical" in m for m in log)
    elev, xs, ys, _, _ = run(names)
    assert np.abs(elev - truth(xs, ys)).max() < 1e-3


def test_corrupted_missing_and_empty_files_give_clear_errors(tmp_path):
    (tmp_path / "junk.tif").write_bytes(b"this is not a tiff" * 20)
    (tmp_path / "empty.tif").write_bytes(b"")
    with pytest.raises(raster.InputError, match="raster illisible"):
        prepared(tmp_path / "junk.tif")
    with pytest.raises(raster.InputError, match="raster illisible"):
        prepared(tmp_path / "empty.tif")
    with pytest.raises(raster.InputError, match="introuvable"):
        prepared(tmp_path / "absent.tif")
    with pytest.raises(raster.InputError, match="n'est pas un fichier"):
        prepared(tmp_path)


def test_layer_option_is_refused_for_formats_without_layers(tmp_path):
    path = write_raster(tmp_path / "a.tif", plane_array((100, 100), GT), GT)
    with pytest.raises(raster.InputError, match="--layer ne s'applique qu'aux"):
        prepared(path, layer="x")


def test_lidar_options_are_refused_for_rasters(tmp_path):
    path = write_raster(tmp_path / "a.tif", plane_array((100, 100), GT), GT)
    with pytest.raises(raster.InputError, match="--lidar-select ne s'applique qu'aux fichiers LAS/LAZ"):
        prepared(path, lidar_select="lowest")


def test_two_adjacent_geotiffs_are_mosaicked(tmp_path):
    g1 = (X0, 10.0, 0.0, Y0 + 1000.0, 0.0, -10.0)
    g2 = (X0 + 500.0, 10.0, 0.0, Y0 + 1000.0, 0.0, -10.0)
    a = write_raster(tmp_path / "w.tif", plane_array((100, 50), g1), g1, nodata=-9999.0)
    b = write_raster(tmp_path / "e.tif", plane_array((100, 50), g2), g2, nodata=-9999.0)
    elev, xs, ys, nodata, _ = run(prepared([a, b])[1])
    assert not np.any(elev == nodata) and np.abs(elev - truth(xs, ys)).max() < 1e-3


def test_raster_and_lidar_cannot_be_mixed(tmp_path):
    a = write_raster(tmp_path / "a.tif", plane_array((100, 100), GT), GT)
    (tmp_path / "x.las").write_bytes(b"LASF" + b"\0" * 300)
    with pytest.raises(raster.InputError, match="mélange"):
        prepared([a, str(tmp_path / "x.las")])


def test_lidar_is_detected_by_signature_not_only_by_extension(tmp_path):
    (tmp_path / "renamed.dat").write_bytes(b"LASF" + b"\0" * 300)
    (tmp_path / "plain.dat").write_bytes(b"II*\0")
    assert formats.detect(str(tmp_path / "renamed.dat")) == "lidar"
    assert formats.detect(str(tmp_path / "plain.dat")) == "raster"
    assert formats.detect("/vsicurl/https://example.org/x.tif") == "raster"


# ---------------------------------------------------------------------------------------------------------------
# ESRI ASCII Grid
# ---------------------------------------------------------------------------------------------------------------

def test_asc_corner_and_centre_registration_agree_and_rows_run_north_to_south(tmp_path):
    arr = plane_array((100, 100), GT)
    a = write_asc(tmp_path / "corner.asc", arr, X0, Y0, 10.0, corner=True, prj=srs_wkt(2056))
    b = write_asc(tmp_path / "centre.asc", arr, X0 + 5.0, Y0 + 5.0, 10.0, corner=False, prj=srs_wkt(2056))
    for path in (a, b):
        elev, xs, ys, _, _ = run(prepared(path)[1])
        assert np.abs(elev - truth(xs, ys)).max() < 2e-3
        assert elev[0].mean() > elev[-1].mean()                                  # première ligne du fichier = nord


def test_asc_nodata_value_and_decimals(tmp_path):
    arr = plane_array((100, 100), GT, np.float64)
    arr[40:44, 40:44] = -9999
    path = write_asc(tmp_path / "h.asc", arr, X0, Y0, 10.0, prj=srs_wkt(2056), fmt="%.4f")
    names = prepared(path)[1]
    wkt, cx, cy = raster.resolve_target_crs_and_center(names, None, None, CENTRE["cx"], CENTRE["cy"])
    raw, xs, ys, nodata = raster.warp_window(names, wkt, cx, cy, 300.0, 10.0, progress_cb=lambda *a: 1, fill=False)
    assert (raw == nodata).sum() == 16 and not np.any(raw == 0.0)
    good = raw != nodata
    assert np.abs(raw[good] - truth(xs, ys)[good]).max() < 1e-3


def test_asc_integer_rows_followed_by_decimal_rows_keep_their_decimals(tmp_path):
    arr = np.full((100, 100), 1000.0)
    arr[60:] += 0.25                                                             # décimales seulement vers la fin
    path = write_asc(tmp_path / "mix.asc", arr, X0, Y0, 10.0, prj=srs_wkt(2056), fmt="%g")
    elev, *_ = run(prepared(path)[1])
    assert set(np.unique(np.round(elev, 3))) == {1000.0, 1000.25}


def test_asc_dx_dy_header_is_understood(tmp_path):
    arr = plane_array((100, 100), GT)
    path = write_asc(tmp_path / "dxdy.asc", arr, X0, Y0, 10.0, prj=srs_wkt(2056), dxdy=(10.0, 10.0))
    elev, xs, ys, _, _ = run(prepared(path)[1])
    assert np.abs(elev - truth(xs, ys)).max() < 2e-3


def test_asc_in_geographic_degrees_is_not_taken_for_metres(tmp_path):
    cell = 0.0002                                                                # ~15 m en longitude
    lon0, lat0 = 7.0, 46.0
    cols = lon0 + (np.arange(100) + 0.5) * cell
    rows = lat0 + 100 * cell - (np.arange(100) + 0.5) * cell
    arr = 1000.0 + 20000.0 * (cols[None, :] - lon0) + 30000.0 * (rows[:, None] - lat0)
    path = write_asc(tmp_path / "deg.asc", arr, lon0, lat0, cell, prj=srs_wkt(4326), fmt="%.4f")
    names = prepared(path)[1]
    elev, xs, ys, _, wkt = window(names, lat=lat0 + 50 * cell, lon=lon0 + 50 * cell, radius=100.0, pixel=5.0)
    assert "UTM zone 32N" in wkt and elev.shape == (40, 40)
    lon, lat = to_geographic(xs, ys, wkt)
    expected = 1000.0 + 20000.0 * (lon - lon0) + 30000.0 * (lat - lat0)
    assert np.abs(elev - expected).max() < 0.05


def test_asc_without_prj_is_refused_then_accepted_with_source_crs(tmp_path):
    path = write_asc(tmp_path / "noprj.asc", plane_array((100, 100), GT), X0, Y0, 10.0)
    with pytest.raises(raster.InputError, match="--source-crs"):
        prepared(path)
    elev, xs, ys, _, _ = run(prepared(path, source_crs="EPSG:2056")[1])
    assert np.abs(elev - truth(xs, ys)).max() < 2e-3


@pytest.mark.parametrize("header, message", [
    (["nrows 3", "xllcorner 0", "yllcorner 0", "cellsize 1"], "champ.*ncols"),
    (["ncols 3", "nrows 3", "xllcorner 0", "cellsize 1"], "yllcorner ou yllcenter"),
    (["ncols 3", "nrows 3", "xllcorner 0", "yllcorner 0"], "cellsize"),
    (["ncols 3", "nrows 3", "xllcorner 0", "yllcorner 0", "cellsize abc"], "cellsize = « abc »"),
    (["ncols 3", "nrows 3", "xllcorner 0", "yllcorner 0", "cellsize -1"], "cellsize doit être strictement positif"),
    (["ncols 0", "nrows 3", "xllcorner 0", "yllcorner 0", "cellsize 1"], "ncols doit être strictement positif"),
    (["ncols 3", "nrows 3", "xllcorner 0", "yllcorner 0", "cellsize 1", "NODATA_value x"], "NODATA_value"),
    (["ncols 3 4", "nrows 3", "xllcorner 0", "yllcorner 0", "cellsize 1"], "mal formé"),
])
def test_asc_malformed_headers_are_refused_with_the_field_named(tmp_path, header, message):
    path = write_asc(tmp_path / "bad.asc", np.ones((3, 3)), 0, 0, 1, header=header)
    with pytest.raises(raster.InputError, match=message):
        prepared(path)


@pytest.mark.parametrize("cut", ["mid_row", "whole_rows", "header_only"])
def test_truncated_asc_is_refused_instead_of_hanging(tmp_path, cut):
    """GDAL ne rend parfois jamais la main sur un ASC dont le corps est plus court que l'en-tête."""
    arr = plane_array((100, 100), GT)
    path = write_asc(tmp_path / "t.asc", arr, X0, Y0, 10.0, prj=srs_wkt(2056))
    text = Path(path).read_text()
    cuts = {"mid_row": len(text) - 25, "whole_rows": text.index("\n", len(text) // 2) + 1,
            "header_only": text.index("\n", text.index("NODATA_value")) + 1}
    Path(path).write_text(text[:cuts[cut]])
    with pytest.raises(raster.InputError, match="tronqué"):
        prepared(path)


def test_asc_with_wrapped_rows_is_still_read(tmp_path):
    """Un rang réparti sur plusieurs lignes est valide : le contrôle compte les valeurs, pas les lignes."""
    arr = plane_array((20, 20), (X0, 10.0, 0.0, Y0 + 200.0, 0.0, -10.0))
    path = tmp_path / "wrap.asc"
    path.write_text("ncols 20\nnrows 20\nxllcorner %r\nyllcorner %r\ncellsize 10\nNODATA_value -9999\n" % (X0, Y0)
                    + "\n".join(" ".join("%.3f" % v for v in arr.ravel()[i:i + 7]) for i in range(0, 400, 7)) + "\n")
    path.with_suffix(".prj").write_text(srs_wkt(2056))
    names = prepared(path)[1]
    ds = gdal.Open(names[0] if not names[0].startswith("/vsimem") else str(path))
    assert np.abs(ds.GetRasterBand(1).ReadAsArray() - arr).max() < 2e-3


# ---------------------------------------------------------------------------------------------------------------
# SRTM (.hgt)
# ---------------------------------------------------------------------------------------------------------------

def test_srtm_tile_values_orientation_and_position(tmp_path):
    path = write_hgt(tmp_path / "N46E007.hgt", 7, 46)
    log = []
    kind, names = prepared(path, log=log)
    assert any("tuile isolée N46E007.hgt" in m for m in log)
    elev, xs, ys, _, wkt = window(names, lat=46.5, lon=7.5, radius=3000.0, pixel=90.0)
    lon, lat = to_geographic(xs, ys, wkt)
    assert np.abs(elev - hgt_function(lon, lat)).max() < 0.7                     # quantum entier de 1 m
    assert elev[0].mean() > elev[-1].mean() and elev[:, -1].mean() > elev[:, 0].mean()


@pytest.mark.parametrize("name, west, north", [("N46E007.hgt", 7, 47), ("S12W077.hgt", -77, -11),
                                                ("n00e000.hgt", 0, 1), ("S01W001.hgt", -1, 0)])
def test_srtm_position_comes_from_the_file_name_with_pixel_is_point_convention(tmp_path, name, west, north):
    lat0 = north - 1
    path = write_hgt(tmp_path / name, west, lat0, func=lambda lon, lat: 0 * lon + 0 * lat + 100)
    ds = gdal.Open(prepared(path)[1][0])
    gt = ds.GetGeoTransform()
    half = 0.5 / 1200.0                                                          # les échantillons sont aux nœuds
    assert gt[0] == pytest.approx(west - half) and gt[3] == pytest.approx(north + half)
    assert (ds.RasterXSize, ds.RasterYSize) == (1201, 1201) and gt[1] == pytest.approx(1 / 1200)


def test_srtm_voids_are_nodata_not_altitudes(tmp_path):
    path = write_hgt(tmp_path / "N46E007.hgt", 7, 46, voids=[(580, 620, 580, 620)])
    names = prepared(path)[1]
    wkt, cx, cy = raster.resolve_target_crs_and_center(names, 46.5, 7.5, None, None)
    raw, xs, ys, nodata = raster.warp_window(names, wkt, cx, cy, 3000.0, 90.0, progress_cb=lambda *a: 1, fill=False)
    assert nodata == -32768.0 and (raw == nodata).sum() > 100
    assert raw[raw != nodata].min() > 0                                          # aucune valeur sentinelle perdue
    filled, *_ = window(names, lat=46.5, lon=7.5, radius=3000.0, pixel=90.0)
    assert filled.min() > 0 and filled.min() > raw[raw != nodata].min() - 1


def test_srtm_name_and_size_are_validated_with_clear_messages(tmp_path):
    good = write_hgt(tmp_path / "N46E007.hgt", 7, 46)
    bad_name = tmp_path / "tile.hgt"
    bad_name.write_bytes(Path(good).read_bytes())
    with pytest.raises(raster.InputError, match="N46E007.hgt"):
        prepared(bad_name)
    truncated = tmp_path / "N47E007.hgt"
    truncated.write_bytes(Path(good).read_bytes()[:-1000])
    with pytest.raises(raster.InputError, match="taille inattendue"):
        prepared(truncated)


def test_srtm_adjacent_tiles_are_assembled_without_a_seam(tmp_path):
    """Référence : UN SEUL GeoTIFF de 2401 colonnes portant exactement les mêmes échantillons, avec la même
    convention « pixel = point » (demi-pixel de décalage). L'assemblage de deux tuiles doit donner la même grille
    à l'arrondi flottant près : tout décalage d'un pixel ou toute colonne dupliquée au bord ferait apparaître un écart."""
    a = write_hgt(tmp_path / "N46E007.hgt", 7, 46)
    b = write_hgt(tmp_path / "N46E008.hgt", 8, 46)
    log = []
    names = prepared([a, b], log=log)[1]
    assert any("assemblage de 2 tuiles" in m for m in log) and not any("manquante" in m for m in log)
    half = 0.5 / 1200.0
    ref = write_raster(tmp_path / "ref.tif", hgt_array(7, 46, ncols=2401).astype(np.float32),
                       (7 - half, 1 / 1200.0, 0.0, 47 + half, 0.0, -1 / 1200.0), epsg=4326, nodata=-32768.0)
    got, xs, ys, nodata, wkt = window(names, lat=46.5, lon=8.0, radius=4000.0, pixel=90.0)
    want, *_ = window([ref], lat=46.5, lon=8.0, radius=4000.0, pixel=90.0)
    lon, lat = to_geographic(xs, ys, wkt)
    assert (lon < 8.0 - 0.01).any() and (lon > 8.0 + 0.01).any()                 # la fenêtre chevauche vraiment la couture
    assert not np.any(got == nodata) and got.shape == want.shape
    assert np.abs(got - want).max() < 1e-2
    assert np.abs(got - hgt_function(lon, lat)).max() < 1.0                      # et reste proche du terrain vrai


def test_srtm_missing_tile_in_the_assembly_is_reported(tmp_path):
    a = write_hgt(tmp_path / "N46E007.hgt", 7, 46)
    b = write_hgt(tmp_path / "N46E009.hgt", 9, 46)
    log = []
    prepared([a, b], log=log)
    assert any("manquante" in m and "N46E008" in m for m in log)


def test_srtm_void_on_a_shared_border_does_not_hide_the_valid_value_of_the_neighbour_tile(tmp_path):
    """La colonne de bordure existe dans les DEUX tuiles : un vide dans l'une ne doit pas effacer la valeur de l'autre."""
    a = write_hgt(tmp_path / "N46E007.hgt", 7, 46, voids=[(500, 520, 1200, 1201)])      # dernière colonne de A
    b = write_hgt(tmp_path / "N46E008.hgt", 8, 46)
    ds = gdal.Open(prepared([a, b])[1][0])
    mosaic = ds.GetRasterBand(1).ReadAsArray()
    assert mosaic.shape == (1201, 2401)
    assert np.all(mosaic[500:520, 1200] == np.rint(hgt_function(8.0, 46.0 + 1.0 - np.arange(500, 520) / 1200.0)))
    assert (mosaic == -32768).sum() == 0                                                  # aucun vide fabriqué


def test_srtm_mosaic_is_one_grid_with_the_expected_extent(tmp_path):
    a = write_hgt(tmp_path / "N46E007.hgt", 7, 46)
    b = write_hgt(tmp_path / "N47E007.hgt", 7, 47)                                        # voisine au NORD
    names = prepared([a, b])[1]
    assert len(names) == 1
    ds = gdal.Open(names[0])
    gt, half = ds.GetGeoTransform(), 0.5 / 1200.0
    assert (ds.RasterXSize, ds.RasterYSize) == (1201, 2401)
    assert gt[0] == pytest.approx(7 - half) and gt[3] == pytest.approx(48 + half)
    assert ds.GetRasterBand(1).GetNoDataValue() == -32768.0


def test_srtm_tiles_of_different_resolutions_are_not_silently_merged(tmp_path):
    a = write_hgt(tmp_path / "N46E007.hgt", 7, 46)
    big = tmp_path / "N46E008.hgt"
    with open(big, "wb") as f:
        f.truncate(3601 * 3601 * 2)                                                       # SRTM1, fichier creux
    log = []
    names = prepared([a, str(big)], log=log)[1]
    assert len(names) == 2 and any("résolutions différentes" in m for m in log)
