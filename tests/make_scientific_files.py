# SPDX-License-Identifier: GPL-3.0-or-later
"""Crée des NetCDF et HDF5 synthétiques avec les bibliothèques netCDF4 et h5py, ÉCRIVAINS INDÉPENDANTS de GDAL
(un aller-retour GDAL -> GDAL ne prouverait rien sur les conventions : orientation des axes, paquetage, unités).

Exécuté en SOUS-PROCESSUS par les tests (python make_scientific_files.py DOSSIER) : les roues pip de netCDF4 et h5py
embarquent leur propre HDF5, qui entre en conflit avec celui de GDAL si les deux sont chargés dans le même
interpréteur. Sortie : 0 si tout a été créé, 3 si netCDF4 ou h5py est absent (les tests concernés sont alors ignorés).

Terrain commun (lon/lat) : z = 1000 + 800 * (lon - 7.0) + 1200 * (lat - 46.0), lon 7.0..7.3, lat 46.0..46.2."""
import os
import sys

import numpy as np

try:
    import h5py
    import netCDF4
except ImportError:
    sys.exit(3)

D = sys.argv[1]
os.makedirs(D, exist_ok=True)
lat = np.linspace(46.0, 46.2, 51)           # CROISSANT (sud -> nord) : le cas fréquent, GDAL doit retourner la grille
lon = np.linspace(7.0, 7.3, 61)
LON, LAT = np.meshgrid(lon, lat)
z = (1000 + 800 * (LON - 7.0) + 1200 * (LAT - 46.0)).astype(np.float64)


def base(path):
    ds = netCDF4.Dataset(path, "w", format="NETCDF4")
    ds.createDimension("lat", len(lat))
    ds.createDimension("lon", len(lon))
    v = ds.createVariable("lat", "f8", ("lat",)); v[:] = lat; v.units = "degrees_north"; v.standard_name = "latitude"
    v = ds.createVariable("lon", "f8", ("lon",)); v[:] = lon; v.units = "degrees_east"; v.standard_name = "longitude"
    return ds


# 1) variable unique, float32, _FillValue, mètres
ds = base(f"{D}/simple.nc")
v = ds.createVariable("elevation", "f4", ("lat", "lon"), fill_value=-9999.0)
data = z.copy(); data[24:27, 29:32] = -9999.0         # trou de 3 x 3 cellules au centre de la fenêtre des tests
v.set_auto_maskandscale(False); v[:] = data.astype(np.float32); v.units = "m"
v.standard_name = "height_above_reference_ellipsoid"
ds.close()

# 2) entiers 16 bits empaquetés : scale_factor 0.1, add_offset 500, _FillValue -32768
ds = base(f"{D}/packed.nc")
v = ds.createVariable("z", "i2", ("lat", "lon"), fill_value=-32768)
v.scale_factor = 0.1; v.add_offset = 500.0; v.units = "m"
packed = np.rint((z - 500.0) / 0.1); packed[24:27, 29:32] = -32768
v.set_auto_maskandscale(False); v[:] = packed.astype(np.int16)
ds.close()

# 3) plusieurs variables : une altitude, une température (K), une épaisseur de glace (m, pas une altitude)
ds = base(f"{D}/several.nc")
a = ds.createVariable("elevation", "f4", ("lat", "lon")); a[:] = z; a.units = "m"
b = ds.createVariable("temperature", "f4", ("lat", "lon")); b[:] = 15 + 0 * z; b.units = "K"
c = ds.createVariable("ice_thickness", "f4", ("lat", "lon")); c[:] = 3 + 0 * z; c.units = "m"
ds.close()

# 3 bis) plusieurs variables dont AUCUNE ne se déclare comme une altitude
ds = base(f"{D}/ambiguous.nc")
for name in ("var_a", "var_b"):
    a = ds.createVariable(name, "f4", ("lat", "lon")); a[:] = z; a.units = "m"
ds.close()

# 4) dimension temps : elevation(time, lat, lon) = z + 10 * k
ds = base(f"{D}/time.nc")
ds.createDimension("time", 3)
t = ds.createVariable("time", "f8", ("time",)); t[:] = [0, 1, 2]; t.units = "days since 2020-01-01"; t.standard_name = "time"
v = ds.createVariable("elevation", "f4", ("time", "lat", "lon"))
for k in range(3):
    v[k, :, :] = z + 10 * k
v.units = "m"
ds.close()

# 5) unité en pieds
ds = base(f"{D}/feet.nc")
v = ds.createVariable("elevation", "f4", ("lat", "lon")); v[:] = z / 0.3048; v.units = "feet"
ds.close()

# 6) variable orientée vers le bas (profondeur positive) : altitude = -valeur
ds = base(f"{D}/depth.nc")
v = ds.createVariable("depth", "f4", ("lat", "lon")); v[:] = -(z - 1500.0); v.units = "m"; v.positive = "down"
v.standard_name = "sea_floor_depth_below_geoid"
ds.close()

# 7) aucune coordonnée spatiale : dimensions anonymes
ds = netCDF4.Dataset(f"{D}/nogeo.nc", "w")
ds.createDimension("y", 50); ds.createDimension("x", 60)
v = ds.createVariable("elevation", "f4", ("y", "x")); v[:] = z[:50, :60]; v.units = "m"
ds.close()

# 8) projeté UTM 32N : x/y en mètres (y CROISSANT) + grid_mapping CF
x = 400000 + 30.0 * np.arange(80)
y = 5100000 + 30.0 * np.arange(70)
ds = netCDF4.Dataset(f"{D}/utm.nc", "w")
ds.createDimension("y", len(y)); ds.createDimension("x", len(x))
vx = ds.createVariable("x", "f8", ("x",)); vx[:] = x; vx.units = "m"; vx.standard_name = "projection_x_coordinate"
vy = ds.createVariable("y", "f8", ("y",)); vy[:] = y; vy.units = "m"; vy.standard_name = "projection_y_coordinate"
gm = ds.createVariable("crs", "i4")
gm.grid_mapping_name = "transverse_mercator"; gm.longitude_of_central_meridian = 9.0
gm.latitude_of_projection_origin = 0.0; gm.false_easting = 500000.0; gm.false_northing = 0.0
gm.scale_factor_at_central_meridian = 0.9996; gm.semi_major_axis = 6378137.0; gm.inverse_flattening = 298.257223563
gm.spatial_ref = ('PROJCS["WGS 84 / UTM zone 32N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,'
                  '298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
                  'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",9],PARAMETER["scale_factor",0.9996],'
                  'PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1]]')
v = ds.createVariable("dem", "f4", ("y", "x"))
v[:] = 1500 + 0.05 * np.arange(len(x))[None, :] + 0.03 * np.arange(len(y))[:, None]
v.units = "m"; v.grid_mapping = "crs"
ds.close()

# --- HDF5 ---------------------------------------------------------------------------------------------------------
with h5py.File(f"{D}/plain.hdf5", "w") as f:               # HDF5 brut : aucun géoréférencement
    d = f.create_dataset("elevation", data=z.astype(np.float32)); d.attrs["units"] = "m"
with h5py.File(f"{D}/groups.h5", "w") as f:
    f.create_dataset("data/elevation", data=z.astype(np.float32))
    f.create_dataset("data/amplitude", data=(z * 0 + 1).astype(np.float32))
# NetCDF-4 enregistré sous l'extension .h5 : c'est du HDF5 géoréférencé selon le modèle NetCDF-4/CF
import shutil
shutil.copy(f"{D}/simple.nc", f"{D}/netcdf4_as.h5")
shutil.copy(f"{D}/simple.nc", f"{D}/netcdf4_as.hdf5")
print("fichiers créés")
