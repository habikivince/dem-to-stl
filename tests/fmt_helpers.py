# SPDX-License-Identifier: GPL-3.0-or-later
"""Fabrication de petits fichiers sources synthétiques (GDAL seul, aucun téléchargement) et mesures partagées par
test_input_formats.py, test_lidar.py et test_input_formats_cli.py.

Le terrain de référence est un PLAN : z = Z0 + SLOPE_X * (x - X0) + SLOPE_Y * (y - Y0). Il monte vers l'est et vers
le nord, ce qui rend visibles une inversion nord-sud ou est-ouest. Un plan est reproduit exactement par le
ré-échantillonnage moyen ou bilinéaire : toute erreur mesurée vient donc du lecteur, pas de l'interpolation."""
import sys
from pathlib import Path

import numpy as np
from osgeo import gdal, osr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import raster  # noqa: E402

gdal.UseExceptions()
osr.UseExceptions()

X0, Y0 = 2600000.0, 1200000.0        # EPSG:2056 (CH1903+ / LV95), mètres
Z0, SLOPE_X, SLOPE_Y = 1000.0, 0.05, 0.02
quiet = lambda *a, **k: 1  # noqa: E731


def plane(x, y):
    return Z0 + SLOPE_X * (np.asarray(x) - X0) + SLOPE_Y * (np.asarray(y) - Y0)


def srs_wkt(epsg):
    s = osr.SpatialReference()
    s.ImportFromEPSG(epsg)
    return s.ExportToWkt()


def centres(gt, shape):
    """Coordonnées (x, y) des centres de pixels d'une grille géoréférencée nord-en-haut."""
    rows, cols = shape
    return gt[0] + (np.arange(cols) + 0.5) * gt[1], gt[3] + (np.arange(rows) + 0.5) * gt[5]


def plane_array(shape, gt, dtype=np.float32):
    xs, ys = centres(gt, shape)
    return plane(xs[None, :], ys[:, None]).astype(dtype)


def write_raster(path, arr, gt, epsg=2056, fmt="GTiff", options=(), nodata=None, scale=None, offset=None,
                 unit=None, wkt=None, gdal_type=None, extra_bands=()):
    """Écrit `arr` (lignes du nord au sud) ; epsg=None et wkt=None : aucun CRS. Retourne le chemin."""
    gdal_type = gdal_type or gdal.GetDataTypeByName({"float32": "Float32", "float64": "Float64", "int16": "Int16",
                                                     "uint16": "UInt16", "int32": "Int32", "uint8": "Byte"}[arr.dtype.name])
    ds = gdal.GetDriverByName(fmt).Create(str(path), arr.shape[1], arr.shape[0], 1 + len(extra_bands), gdal_type,
                                          list(options))
    ds.SetGeoTransform(gt)
    if wkt or epsg:
        ds.SetProjection(wkt or srs_wkt(epsg))
    for i, a in enumerate((arr, *extra_bands), start=1):
        b = ds.GetRasterBand(i)
        if nodata is not None:
            b.SetNoDataValue(nodata)
        b.WriteArray(a)
    b = ds.GetRasterBand(1)
    if scale is not None:
        b.SetScale(scale)
    if offset is not None:
        b.SetOffset(offset)
    if unit is not None:
        b.SetUnitType(unit)
    ds = None
    return str(path)


def window(inputs, cx=None, cy=None, lat=None, lon=None, radius=60.0, pixel=2.0, resample="auto", fill=True):
    """resolve_target_crs_and_center + warp_window, comme build_mesh : (elev, xs, ys, nodata, dst_wkt).
    fill=False laisse les trous tels que la source les porte (sinon ils sont comblés, comme dans le programme)."""
    wkt, cx2, cy2 = raster.resolve_target_crs_and_center(inputs, lat, lon, cx, cy)
    elev, xs, ys, nodata = raster.warp_window(inputs, wkt, cx2, cy2, radius, pixel, progress_cb=quiet,
                                              resample=resample, fill=fill)
    return elev, xs, ys, nodata, wkt


def to_geographic(xs, ys, wkt):
    """Longitude/latitude (WGS84) des nœuds (xs, ys) d'une grille exprimée dans `wkt`."""
    src = osr.SpatialReference()
    src.ImportFromWkt(wkt)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(4326)
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct = osr.CoordinateTransformation(src, dst)
    gx, gy = np.meshgrid(xs, ys)
    pts = np.array(ct.TransformPoints(np.column_stack([gx.ravel(), gy.ravel()])))
    return pts[:, 0].reshape(gx.shape), pts[:, 1].reshape(gx.shape)


_SCIENTIFIC = {}


def scientific_dir(tmp_path_factory):
    """Dossier des NetCDF/HDF5 fabriqués par make_scientific_files.py (une fois par session). Ignore le test, avec la
    raison, si netCDF4 ou h5py manque."""
    import subprocess

    import pytest
    if "result" not in _SCIENTIFIC:
        d = tmp_path_factory.mktemp("scientific")
        r = subprocess.run([sys.executable, str(Path(__file__).with_name("make_scientific_files.py")), str(d)],
                           capture_output=True, text=True, timeout=120)
        _SCIENTIFIC["result"] = (r.returncode, r.stderr, d)
    code, err, d = _SCIENTIFIC["result"]
    if code == 3:
        pytest.skip("netCDF4 et/ou h5py ne sont pas installés : impossible de fabriquer des NetCDF/HDF5 avec un écrivain "
                    "indépendant de GDAL (pip install netCDF4 h5py)")
    assert code == 0, err
    return d


def prepared(paths, log=None, **opts):
    """formats.prepare avec les options données : (genre, entrées). `log` : liste qui reçoit les messages."""
    import formats
    sink = (lambda msg: log.append(msg)) if log is not None else (lambda msg: None)
    return formats.prepare(paths if isinstance(paths, (list, tuple)) else [paths], formats.SourceOptions(**opts),
                           log=sink)


# --- fichiers ESRI ASCII Grid et SRTM synthétiques --------------------------------------------------------------

def write_asc(path, arr, x, y, cell, corner=True, nodata=-9999, prj=None, header=None, fmt="%.3f", dxdy=None):
    lines = header if header is not None else [
        f"ncols {arr.shape[1]}", f"nrows {arr.shape[0]}",
        f"xll{'corner' if corner else 'center'} {x!r}", f"yll{'corner' if corner else 'center'} {y!r}",
        f"cellsize {cell!r}" if dxdy is None else f"dx {dxdy[0]!r}\ndy {dxdy[1]!r}",
        f"NODATA_value {nodata}"]
    body = "\n".join(" ".join(fmt % v for v in row) for row in arr)
    Path(path).write_text("\n".join(lines) + "\n" + body + "\n")
    if prj:
        Path(path).with_suffix(".prj").write_text(prj)
    return str(path)


def hgt_function(lon, lat, lon0=7.0, lat0=46.0):
    return 500.0 + 800.0 * (lon - lon0) + 1200.0 * (lat - lat0)


def hgt_array(lon0, lat0, func=hgt_function, voids=(), ncols=1201):
    """Échantillons SRTM3 (entiers 16 bits), du nord au sud, placés AUX nœuds : le bord est partagé avec la tuile
    voisine. `ncols` > 1201 sert à fabriquer l'équivalent d'un assemblage de tuiles côte à côte."""
    lon = lon0 + np.arange(ncols) / 1200.0
    lat = lat0 + 1.0 - np.arange(1201) / 1200.0
    data = np.rint(func(lon[None, :] + 0 * lat[:, None], lat[:, None] + 0 * lon[None, :])).astype(">i2")
    for r0, r1, c0, c1 in voids:
        data[r0:r1, c0:c1] = -32768
    return data


def write_hgt(path, lon0, lat0, func=hgt_function, voids=()):
    """Tuile SRTM3 : 1201 x 1201 entiers 16 bits grande-boutiste."""
    Path(path).write_bytes(hgt_array(lon0, lat0, func, voids).tobytes())
    return str(path)


# --- bout en bout : ligne de commande -> STL ---------------------------------------------------------------------

def cli(monkeypatch, tmp_path, files, extra=(), center=None, radius=300.0, diameter=60.0, base=2.0, vexag=1.0,
        spacing=0.5, name="out.stl"):
    """Exécute circle_dem_to_stl.main() comme en ligne de commande (sans interaction) et retourne le chemin du STL."""
    import circle_dem_to_stl
    out = tmp_path / name
    argv = ["x", "--input", *[str(f) for f in files], "--radius", repr(float(radius)), "--diameter-mm", repr(diameter),
            "--vexag", repr(vexag), "--base-mm", repr(base), "--print-spacing-mm", repr(spacing), "--output", str(out),
            "--force", *extra]
    for key, value in (center or {}).items():
        argv += [f"--{key}", repr(float(value))]
    monkeypatch.setattr(sys, "argv", argv)
    circle_dem_to_stl.main()
    return out


def assert_stl_ok(path, radius, gradient, diameter=60.0, base=2.0, vexag=1.0, polygon=0, relief_tol=0.03):
    """Lit le STL écrit et vérifie ce qu'un humain regarderait : fini, étanche, orienté vers l'extérieur, sans triangle
    dégénéré, aux bonnes dimensions, au bon RELIEF (écart d'altitude attendu d'après le gradient connu du terrain,
    donc unités, échelle et décalage vérifiés de bout en bout) et PENCHÉ DU BON CÔTÉ (le sommet est vers le nord-est
    pour un terrain qui monte vers l'est et le nord : détecte toute inversion d'axe).
    polygon : nombre de côtés d'un polygone régulier inscrit dans le cercle de référence (0 : disque)."""
    import check_stl
    tris = check_stl.read_stl(str(path))
    v = tris.reshape(-1, 3)
    assert len(tris) > 1000 and np.isfinite(v).all()
    results = check_stl.check(tris, diameter, base, radius_m=radius, vexag=vexag, relief_m=None)
    failed = [f"{name} : {detail}" for name, ok, detail in results
              if not ok and not (polygon and name.startswith("plus grande dimension"))]
    assert not failed, "; ".join(failed)
    scale = diameter / (2.0 * radius)
    gx, gy = gradient
    expected = 2.0 * radius * float(np.hypot(gx, gy))
    implied = float(v[:, 2].max()) / (scale * vexag)
    if polygon:
        assert 0.0 < implied <= expected + 1.0
        width = v[:, 0].max() - v[:, 0].min()
        height = v[:, 1].max() - v[:, 1].min()
        smallest = diameter * float(np.cos(np.pi / polygon))          # deux fois l'apothème d'un polygone inscrit
        assert smallest - 0.6 <= max(width, height) <= diameter + 0.6
    else:
        assert abs(implied - expected) <= relief_tol * expected + 1.0, \
            f"relief {implied:.1f} m dans le STL, {expected:.1f} m attendu"
        top = v[v[:, 2].argmax()]
        got, want = np.degrees(np.arctan2(top[1], top[0])), np.degrees(np.arctan2(gy, gx))
        assert abs((got - want + 180.0) % 360.0 - 180.0) < 15.0, \
            f"le sommet est à {got:.0f}° du centre, le terrain monte vers {want:.0f}°"
    return tris


def utm_gradient(zfunc, epsg, cx, cy, h=100.0):
    """Gradient (dz/dE, dz/dN) en m/m, dans le CRS projeté epsg, d'un terrain connu en (longitude, latitude)."""
    wkt = srs_wkt(epsg)

    def z(x, y):
        lon, lat = to_geographic(np.array([x]), np.array([y]), wkt)
        return float(zfunc(lon[0, 0], lat[0, 0]))

    return (z(cx + h, cy) - z(cx - h, cy)) / (2 * h), (z(cx, cy + h) - z(cx, cy - h)) / (2 * h)
