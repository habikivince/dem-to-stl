# SPDX-License-Identifier: GPL-3.0-or-later
"""Formats « conteneurs » : GeoPackage raster, NetCDF, HDF5.

Un conteneur peut renfermer plusieurs grilles, plusieurs dates, des variables qui ne sont pas des altitudes : ces
tests vérifient que le programme CHOISIT SEULEMENT quand c'est sans ambiguïté, qu'il refuse sinon en listant les
possibilités, et que les valeurs lues sont les bonnes (paquetage, unités, orientation des axes).

Les NetCDF/HDF5 sont fabriqués par netCDF4 et h5py (écrivains indépendants de GDAL), dans un sous-processus : leurs roues
pip embarquent leur propre HDF5, incompatible avec celui de GDAL dans un même interpréteur. Sans ces deux modules, les
tests concernés sont ignorés avec la raison ; ceux qui n'utilisent que GDAL s'exécutent toujours."""
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal, ogr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import raster  # noqa: E402
from fmt_helpers import (X0, Y0, plane, plane_array, prepared, scientific_dir, srs_wkt,  # noqa: E402
                         to_geographic, window, write_raster)

HERE = Path(__file__).resolve().parent
GT = (X0, 10.0, 0.0, Y0 + 1000.0, 0.0, -10.0)
CENTRE = dict(cx=X0 + 500.0, cy=Y0 + 500.0, radius=300.0, pixel=10.0)


def truth(xs, ys):
    return plane(xs[None, :], ys[:, None])


def run(names, **kw):
    args = dict(CENTRE)
    args.update(kw)
    return window(names, **args)


# ---------------------------------------------------------------------------------------------------------------
# GeoPackage raster
# ---------------------------------------------------------------------------------------------------------------

def make_gpkg(path, arr, table="mnt", tiles="TIFF", append=False, nodata=-9999.0, epsg=2056):
    src = write_raster(Path(path).with_suffix(f".{table}.tif"), arr, GT, epsg=epsg, nodata=nodata)
    opts = [f"TILE_FORMAT={tiles}", f"RASTER_TABLE={table}"] + (["APPEND_SUBDATASET=YES"] if append else [])
    out = gdal.Translate(str(path), src, format="GPKG", creationOptions=opts)
    out = None  # noqa: F841
    return str(path)


def test_gpkg_float32_table_values_and_orientation(tmp_path):
    path = make_gpkg(tmp_path / "dem.gpkg", plane_array((100, 100), GT))
    log = []
    kind, names = prepared(path, log=log)
    assert any("GeoPackage" in m for m in log)
    elev, xs, ys, nodata, _ = run(names)
    assert np.abs(elev - truth(xs, ys)).max() < 1e-3
    assert elev[0].mean() > elev[-1].mean() and elev[:, -1].mean() > elev[:, 0].mean()


def test_gpkg_void_stays_nodata_never_zero(tmp_path):
    arr = plane_array((100, 100), GT)
    arr[48:52, 48:52] = -9999.0
    names = prepared(make_gpkg(tmp_path / "hole.gpkg", arr))[1]
    raw, xs, ys, nodata, _ = run(names, fill=False)
    assert (raw == nodata).sum() == 16 and not np.any(raw == 0.0)


def test_gpkg_integer_tiles_use_the_scale_and_offset_of_the_coverage_ancillary_table(tmp_path):
    """Un MNT GeoPackage entier (extension « tiled gridded coverage ») stocke valeur = brut x échelle + décalage dans la
    table gpkg_2d_gridded_coverage_ancillary ; GDAL les applique à la lecture. Le lecteur ne doit PAS les appliquer
    une deuxième fois. GDAL ne sait pas les écrire : on les pose en SQL (en gardant la convention de GDAL, brut stocké
    = valeur + 32768, et la précision déclarée à 0,1 : avec précision 1, GDAL arrondirait au mètre)."""
    raw = (np.indices((100, 100))[1] * 7 + np.indices((100, 100))[0]).astype(np.int16)
    path = tmp_path / "int16.gpkg"
    ds = gdal.GetDriverByName("GPKG").Create(str(path), 100, 100, 1, gdal.GDT_Int16, ["TILE_FORMAT=PNG", "RASTER_TABLE=mnt16"])
    ds.SetGeoTransform(GT)
    ds.SetProjection(srs_wkt(2056))
    ds.GetRasterBand(1).WriteArray(raw)
    ds = None
    con = sqlite3.connect(path)
    con.execute("UPDATE gpkg_2d_gridded_coverage_ancillary SET scale = 0.1, offset = ?, precision = 0.1",
                (1500.0 - 3276.8,))
    con.commit()
    con.close()
    names = prepared(str(path))[1]
    ds = gdal.Open(names[0])
    got = ds.GetRasterBand(1).ReadAsArray()
    assert np.abs(got - (1500.0 + raw * 0.1)).max() < 1e-3            # une seule application de l'échelle
    elev, xs, ys, _, _ = run(names)
    assert 1500 <= elev.min() and elev.max() < 1500 + 0.1 * (99 * 7 + 99) + 1


def _two_tables(tmp_path):
    a = plane_array((100, 100), GT)
    path = make_gpkg(tmp_path / "two.gpkg", a, table="mnt_a")
    return make_gpkg(path, a + 100.0, table="mnt_b", append=True)


def test_gpkg_with_two_raster_tables_requires_an_explicit_layer(tmp_path):
    with pytest.raises(raster.InputError, match=r"2 tables raster.*--layer NOM") as e:
        prepared(_two_tables(tmp_path))
    assert "mnt_a" in str(e.value) and "mnt_b" in str(e.value)


def test_gpkg_layer_selects_the_table_asked_for(tmp_path):
    path = _two_tables(tmp_path)
    for table, extra in (("mnt_a", 0.0), ("mnt_b", 100.0)):
        elev, xs, ys, _, _ = run(prepared(path, layer=table)[1])
        assert np.abs(elev - (truth(xs, ys) + extra)).max() < 1e-3


def test_gpkg_gdal_subdataset_syntax_is_accepted_as_is(tmp_path):
    path = _two_tables(tmp_path)
    elev, xs, ys, _, _ = run(prepared(f"GPKG:{path}:mnt_b")[1])
    assert np.abs(elev - (truth(xs, ys) + 100.0)).max() < 1e-3


def test_gpkg_unknown_layer_is_refused_with_the_available_tables(tmp_path):
    with pytest.raises(raster.InputError, match="aucune table") as e:
        prepared(_two_tables(tmp_path), layer="nope")
    assert "mnt_a" in str(e.value) and "mnt_b" in str(e.value)


def test_gpkg_vector_layers_next_to_a_single_raster_table_are_ignored(tmp_path):
    path = make_gpkg(tmp_path / "mixed.gpkg", plane_array((100, 100), GT))
    vec = gdal.OpenEx(path, gdal.OF_VECTOR | gdal.OF_UPDATE)
    lyr = vec.CreateLayer("points", None, ogr.wkbPoint)
    vec = None  # noqa: F841
    elev, xs, ys, _, _ = run(prepared(path)[1])
    assert np.abs(elev - truth(xs, ys)).max() < 1e-3


def test_gpkg_with_only_vector_layers_is_not_a_dem(tmp_path):
    path = tmp_path / "vec.gpkg"
    v = ogr.GetDriverByName("GPKG").CreateDataSource(str(path))
    v.CreateLayer("pts", None, ogr.wkbPoint)
    v = None
    with pytest.raises(raster.InputError, match="sans table raster.*pts"):
        prepared(path)


def test_gpkg_image_tiles_are_refused_not_read_as_altitudes(tmp_path):
    src = tmp_path / "rgb.tif"
    ds = gdal.GetDriverByName("GTiff").Create(str(src), 50, 50, 3, gdal.GDT_Byte)
    ds.SetGeoTransform(GT)
    ds.SetProjection(srs_wkt(2056))
    for i, ci in enumerate((gdal.GCI_RedBand, gdal.GCI_GreenBand, gdal.GCI_BlueBand), 1):
        ds.GetRasterBand(i).Fill(40 * i)
        ds.GetRasterBand(i).SetColorInterpretation(ci)
    ds = None
    out = gdal.Translate(str(tmp_path / "ortho.gpkg"), str(src), format="GPKG", creationOptions=["RASTER_TABLE=ortho"])
    out = None  # noqa: F841
    with pytest.raises(raster.InputError, match="tuiles image.*pas des altitudes"):
        prepared(tmp_path / "ortho.gpkg")


def test_corrupted_gpkg_gives_a_clear_error(tmp_path):
    (tmp_path / "bad.gpkg").write_bytes(b"SQLite format 3\0" + b"\xff" * 400)
    (tmp_path / "junk.gpkg").write_bytes(b"not a database" * 40)
    for name in ("bad.gpkg", "junk.gpkg"):
        with pytest.raises(raster.InputError, match="illisible|sans table raster"):
            prepared(tmp_path / name)


# ---------------------------------------------------------------------------------------------------------------
# NetCDF écrit par GDAL (toujours exécuté, même sans netCDF4 ni h5py)
# ---------------------------------------------------------------------------------------------------------------

def test_netcdf_written_by_gdal_round_trips_but_the_variable_must_be_chosen(tmp_path):
    """GDAL nomme sa variable « Band1 », sans standard_name ni unité : rien ne dit que c'est une altitude, donc le
    programme demande --layer au lieu de supposer."""
    gt = (400000.0, 10.0, 0.0, 5100000.0 + 1000.0, 0.0, -10.0)
    src = write_raster(tmp_path / "utm.tif", plane_array((100, 100), gt), gt, epsg=32632, nodata=-9999.0)
    out = gdal.Translate(str(tmp_path / "gdal.nc"), src, format="netCDF")
    out = None  # noqa: F841
    with pytest.raises(raster.InputError, match="ne se déclare pas comme une altitude.*--layer Band1"):
        prepared(tmp_path / "gdal.nc")
    names = prepared(tmp_path / "gdal.nc", layer="Band1")[1]
    got = gdal.Open(names[0])
    arr = got.GetRasterBand(1).ReadAsArray()
    assert arr.shape == (100, 100) and np.abs(arr - plane_array((100, 100), gt)).max() < 1e-3
    assert arr[0].mean() > arr[-1].mean()                                        # nord en haut après lecture


# ---------------------------------------------------------------------------------------------------------------
# NetCDF / HDF5 d'écrivains indépendants
# ---------------------------------------------------------------------------------------------------------------

@pytest.fixture(scope="session")
def sci(tmp_path_factory):
    return scientific_dir(tmp_path_factory)


def sci_truth(lon, lat):
    return 1000.0 + 800.0 * (lon - 7.0) + 1200.0 * (lat - 46.0)


def away_from_voids(elev, nodata, cells=6):
    """Cellules valides à plus de `cells` cellules de travail d'un vide. À côté d'un vide, le ré-échantillonnage
    bilinéaire ne dispose que d'une partie de ses voisins : l'écart y vient de GDAL, pas du lecteur."""
    from scipy import ndimage
    return ~ndimage.binary_dilation(elev == nodata, iterations=cells)


def geo_window(names, fill=True):
    elev, xs, ys, nodata, wkt = window(names, lat=46.1, lon=7.15, radius=3000.0, pixel=100.0, fill=fill)
    lon, lat = to_geographic(xs, ys, wkt)
    return elev, nodata, sci_truth(lon, lat)


def test_netcdf_values_axis_order_and_declared_fill_value(sci):
    """lat croît du sud au nord dans le fichier : GDAL doit retourner la grille ; _FillValue -9999 reste un vide."""
    log = []
    names = prepared(sci / "simple.nc", log=log)[1]
    assert any("WGS84 supposé" in m for m in log)                      # axes en degrés sans datum : dit, pas caché
    elev, nodata, want = geo_window(names, fill=False)
    assert (elev == nodata).sum() > 0 and not np.any(elev == 0.0)
    good = away_from_voids(elev, nodata)
    assert good.sum() > 2000 and np.abs(elev[good] - want[good]).max() < 0.3
    assert elev[0].mean() > elev[-1].mean()                            # le nord est en haut


def test_netcdf_packed_integers_apply_scale_factor_and_add_offset_once(sci):
    names = prepared(sci / "packed.nc")[1]
    elev, nodata, want = geo_window(names, fill=False)
    assert (elev == nodata).sum() > 0 and not np.any(elev == 0.0)      # le trou _FillValue -32768 reste un vide
    good = away_from_voids(elev, nodata)
    assert good.sum() > 2000 and np.abs(elev[good] - want[good]).max() < 0.4     # quantum 0,1 m + ré-échantillonnage
    assert elev[elev != nodata].min() > 1000 - 1 and elev.max() < 1481


def test_netcdf_with_several_variables_picks_the_one_that_declares_itself_an_elevation(sci):
    log = []
    names = prepared(sci / "several.nc", log=log)[1]
    assert any("variable elevation" in m for m in log)
    elev, nodata, want = geo_window(names)
    assert np.abs(elev - want).max() < 0.3                             # ni la température (15) ni l'épaisseur (3)


def test_netcdf_variable_with_a_non_length_unit_is_refused_even_when_asked_for(sci):
    with pytest.raises(raster.InputError, match="« K » n'est pas une unité de longueur"):
        prepared(sci / "several.nc", layer="temperature")


def test_netcdf_explicit_layer_is_honoured_and_unknown_layer_lists_candidates(sci):
    elev, nodata, _ = geo_window(prepared(sci / "several.nc", layer="ice_thickness")[1])
    assert np.allclose(elev, 3.0, atol=1e-3)                           # l'utilisateur a choisi : on obéit
    with pytest.raises(raster.InputError, match="aucune variable de ce nom") as e:
        prepared(sci / "several.nc", layer="nonexistent")
    assert "elevation" in str(e.value) and "ressemble à une altitude" in str(e.value)


def test_netcdf_variables_none_of_which_looks_like_an_elevation_are_ambiguous(sci):
    with pytest.raises(raster.InputError, match=r"2 variables de grille.*--layer NOM") as e:
        prepared(sci / "ambiguous.nc")
    assert "var_a" in str(e.value) and "var_b" in str(e.value)
    elev, nodata, want = geo_window(prepared(sci / "ambiguous.nc", layer="var_b")[1])
    assert np.abs(elev - want).max() < 0.3


def test_netcdf_time_dimension_requires_a_band_and_never_picks_a_slice(sci):
    with pytest.raises(raster.InputError, match=r"3 bandes.*--band N") as e:
        prepared(sci / "time.nc")
    assert "time=0" in str(e.value) and "time=2" in str(e.value)
    for band, extra in ((1, 0.0), (3, 20.0)):
        elev, nodata, want = geo_window(prepared(sci / "time.nc", band=band)[1])
        assert np.abs(elev - (want + extra)).max() < 0.3
    with pytest.raises(raster.InputError, match="--band 9"):
        prepared(sci / "time.nc", band=9)


def test_netcdf_vertical_unit_feet_is_converted(sci):
    log = []
    names = prepared(sci / "feet.nc", log=log)[1]
    assert any("x0.3048" in m for m in log)
    elev, nodata, want = geo_window(names)
    assert np.abs(elev - want).max() < 0.3


def test_netcdf_positive_down_variable_is_turned_into_an_elevation(sci):
    log = []
    names = prepared(sci / "depth.nc", log=log)[1]
    assert any("positive = down" in m for m in log)
    elev, nodata, want = geo_window(names)
    assert np.abs(elev - (want - 1500.0)).max() < 0.3                  # profondeur = -(z - 1500) -> altitude = z - 1500


def test_netcdf_projected_axes_with_grid_mapping_need_no_datum_assumption(sci):
    log = []
    names = prepared(sci / "utm.nc", log=log)[1]
    assert not any("WGS84 supposé" in m for m in log)
    i0, j0 = 40, 35                                                    # centre de la grille, en indices de pixel
    cx, cy = 400000 + 30.0 * i0, 5100000 + 30.0 * j0
    elev, xs, ys, nodata, wkt = window(names, cx=cx, cy=cy, radius=300.0, pixel=10.0)
    assert "UTM zone 32N" in wkt
    want = 1500 + 0.05 * ((xs[None, :] - 400000) / 30.0) + 0.03 * ((ys[:, None] - 5100000) / 30.0)
    assert np.abs(elev - want).max() < 0.02                            # y croissant dans le fichier : bien orienté


def test_netcdf_without_regular_spatial_axes_is_refused(sci):
    with pytest.raises(raster.InputError, match="axes spatiaux réguliers"):
        prepared(sci / "nogeo.nc")


def test_netcdf_source_crs_is_not_allowed_to_replace_the_declared_one(sci):
    with pytest.raises(raster.InputError, match="déclare déjà un CRS"):
        prepared(sci / "utm.nc", source_crs="EPSG:2056")


def test_corrupted_and_truncated_netcdf_hdf5_files_are_reported(sci, tmp_path):
    (tmp_path / "junk.nc").write_bytes(b"not a netcdf" * 50)
    (tmp_path / "junk.h5").write_bytes(b"not an hdf5" * 50)
    data = (sci / "simple.nc").read_bytes()
    (tmp_path / "trunc.nc").write_bytes(data[:len(data) // 3])
    for name, kind in (("junk.nc", "NetCDF"), ("junk.h5", "HDF5"), ("trunc.nc", "NetCDF")):
        with pytest.raises(raster.InputError, match=f"fichier {kind} illisible") as e:
            prepared(tmp_path / name)
        assert "No such file" not in str(e.value)                      # le fichier existe : ne pas laisser croire le contraire


def test_hdf5_netcdf4_model_file_is_read_through_its_cf_structure(sci):
    for name in ("netcdf4_as.h5", "netcdf4_as.hdf5"):
        log = []
        names = prepared(sci / name, log=log)[1]
        assert any("HDF5" in m for m in log)
        elev, nodata, want = geo_window(names, fill=False)
        good = away_from_voids(elev, nodata)
        assert (elev == nodata).sum() > 0 and good.sum() > 2000 and np.abs(elev[good] - want[good]).max() < 0.3


@pytest.mark.parametrize("name, layer", [("plain.hdf5", None), ("groups.h5", None), ("groups.h5", "elevation")])
def test_generic_hdf5_without_georeferencing_is_reported_as_unsupported(sci, name, layer):
    """Un dataset HDF5 brut n'a ni position ni échelle : lui en inventer serait fabriquer un modèle faux. Le message
    liste les datasets trouvés et la marche à suivre."""
    with pytest.raises(raster.InputError, match="HDF5 non géoréférencé") as e:
        prepared(sci / name, layer=layer)
    assert "elevation" in str(e.value) and "gdal_translate" in str(e.value)
