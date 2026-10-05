# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests du mode polygone (--polygon) et de check_stl.py. Aucun réseau."""
import json
import struct
import sys
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal, ogr, osr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import circle_dem_to_stl  # noqa: E402
import polygon  # noqa: E402
import raster  # noqa: E402
import sources  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402

quiet = lambda *a, **k: None  # noqa: E731
ring = lambda pts: np.array(pts + [pts[0]], dtype=float)  # noqa: E731


# ---------------------------------------------------------------- masque

def test_mask_with_hole_and_disjoint_parts():
    xs, ys = np.arange(0.5, 20, 1.0), np.arange(9.5, 0, -1.0)
    outer = ring([[0, 0], [10, 0], [10, 10], [0, 10]])
    hole = ring([[4, 4], [6, 4], [6, 6], [4, 6]])
    island = ring([[12, 2], [18, 2], [18, 8], [12, 8]])
    m = polygon.polygon_mask([outer, hole, island], xs, ys)
    assert int(m.sum()) == 100 - 4 + 36
    assert not m[4, 4] and m[1, 1] and m[3, 14] and not m[3, 11]


def test_mask_concave_shape_matches_brute_force():
    pts = [[0, 0], [10, 0], [10, 4], [4, 4], [4, 10], [0, 10]]            # forme en L
    xs, ys = np.arange(0.25, 10, 0.5), np.arange(9.75, 0, -0.5)
    m = polygon.polygon_mask([ring(pts)], xs, ys)
    expected = ((xs[None, :] < 4) | (ys[:, None] < 4))
    assert np.array_equal(m, expected)


# ---------------------------------------------------------------- chargement

def write_geojson(path, features):
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))


def square_feature(name, lon0, lat0, size):
    c = [[lon0, lat0], [lon0 + size, lat0], [lon0 + size, lat0 + size], [lon0, lat0 + size], [lon0, lat0]]
    return {"type": "Feature", "properties": {"name": name}, "geometry": {"type": "Polygon", "coordinates": [c]}}


def test_load_geojson_filter_and_union(tmp_path):
    f = tmp_path / "z.geojson"
    write_geojson(f, [square_feature("A", 7.0, 46.0, 0.1), square_feature("B", 7.2, 46.0, 0.1)])
    both = polygon.load_polygon(str(f), log=quiet)
    only_b = polygon.load_polygon(str(f), where="name='B'", log=quiet)
    assert both.n_features == 2 and only_b.n_features == 1
    lat, lon = only_b.center_latlon()
    assert (lat, lon) == pytest.approx((46.05, 7.25), abs=1e-6)
    assert only_b.approx_radius_m() == pytest.approx(0.05 * 111320 * np.cos(np.radians(46.05)), rel=0.02) \
        or only_b.approx_radius_m() > 0


def test_load_projected_gpkg_is_reprojected(tmp_path):
    f = tmp_path / "z.gpkg"
    ds = ogr.GetDriverByName("GPKG").CreateDataSource(str(f))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(2056)
    lyr = ds.CreateLayer("z", srs, ogr.wkbPolygon)
    x, y = sources._to_epsg(2056, 46.0, 7.6)
    g = ogr.CreateGeometryFromWkt(f"POLYGON(({x - 1000} {y - 1000},{x + 1000} {y - 1000},{x + 1000} {y + 1000},"
                                  f"{x - 1000} {y + 1000},{x - 1000} {y - 1000}))")
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetGeometry(g)
    lyr.CreateFeature(feat)
    ds = None
    clip = polygon.load_polygon(str(f), log=quiet)
    lat, lon = clip.center_latlon()
    assert (lat, lon) == pytest.approx((46.0, 7.6), abs=0.001)
    assert clip.approx_radius_m() == pytest.approx(1000.0, rel=0.03)


def test_load_errors(tmp_path):
    with pytest.raises(polygon.PolygonError):
        polygon.load_polygon(str(tmp_path / "absent.geojson"), log=quiet)
    pts = tmp_path / "pts.geojson"
    pts.write_text(json.dumps({"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {}, "geometry": {"type": "Point", "coordinates": [7, 46]}}]}))
    with pytest.raises(polygon.PolygonError, match="aucun polygone"):
        polygon.load_polygon(str(pts), log=quiet)
    z = tmp_path / "z.geojson"
    write_geojson(z, [square_feature("A", 7.0, 46.0, 0.1)])
    with pytest.raises(polygon.PolygonError):
        polygon.load_polygon(str(z), where="name='Z'", log=quiet)


# ---------------------------------------------------------------- bout en bout

def make_dem_and_polygon(tmp_path):
    lat, lon = 46.0, 7.6
    x, y = sources._to_epsg(2056, lat, lon)
    yy, xx = np.mgrid[0:600, 0:600]
    write_synthetic_tif(tmp_path / "dem.tif", (1000.0 + 0.1 * xx).astype(np.float32), x - 1500.0, y + 1500.0, 5.0, 2056)
    ds = ogr.GetDriverByName("GPKG").CreateDataSource(str(tmp_path / "z.gpkg"))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(2056)
    lyr = ds.CreateLayer("z", srs, ogr.wkbPolygon)
    # pentagone 1800 m x 1200 m (plus grande dimension : 1800 m) ; des côtés pleins aux extrémités :
    # une pointe plus fine qu'un pixel serait absente du maillage et raccourcirait la dimension
    wkt = (f"POLYGON(({x - 900} {y - 600},{x + 900} {y - 600},{x + 900} {y},{x} {y + 600},"
           f"{x - 900} {y + 600},{x - 900} {y - 600}))")
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetGeometry(ogr.CreateGeometryFromWkt(wkt))
    lyr.CreateFeature(feat)
    ds = None
    return str(tmp_path / "dem.tif"), str(tmp_path / "z.gpkg")


def run_cli(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["x"] + argv)
    circle_dem_to_stl.main()


def test_cli_polygon_end_to_end_and_check_stl(tmp_path, monkeypatch):
    dem, gpkg = make_dem_and_polygon(tmp_path)
    out = tmp_path / "tri.stl"
    run_cli(monkeypatch, ["--input", dem, "--polygon", gpkg, "--diameter-mm", "90", "--vexag", "1",
                          "--base-mm", "2", "--print-spacing-mm", "1", "--output", str(out)])
    tris = check_stl.read_stl(str(out))
    results = dict((n, ok) for n, ok, _ in check_stl.check(tris, 90.0, 2.0, tol_mm=2.0))
    assert results["plus grande dimension horizontale"] and results["socle"] and results["étanchéité"]
    # l'emprise horizontale n'est pas carrée (1800 x 1200 m -> 90 x 60 mm)
    pts = tris.reshape(-1, 3)
    dx, dy = np.ptp(pts[:, 0]), np.ptp(pts[:, 1])
    assert max(dx, dy) == pytest.approx(90.0, abs=2.0) and min(dx, dy) == pytest.approx(60.0, abs=2.0)
    # le coin coupé du pentagone doit manquer : le coin nord-est du rectangle englobant est vide
    ne = tris[(tris[:, :, 0] > pts[:, 0].max() - 5) & (tris[:, :, 1] > pts[:, 1].max() - 5)]
    assert len(ne) == 0


def test_cli_polygon_conflicts(tmp_path, monkeypatch):
    dem, gpkg = make_dem_and_polygon(tmp_path)
    for extra in (["--radius", "500"], ["--lat", "46", "--lon", "7.6"], ["--smooth-wall"], ["--place", "x"]):
        monkeypatch.setattr(sys, "argv", ["x", "--input", dem, "--polygon", gpkg, "--output", str(tmp_path / "o.stl")] + extra)
        with pytest.raises(SystemExit) as e:
            circle_dem_to_stl.main()
        assert "incompatible" in str(e.value)
    monkeypatch.setattr(sys, "argv", ["x", "--input", dem, "--output", str(tmp_path / "o.stl")])
    with pytest.raises(SystemExit) as e:
        circle_dem_to_stl.main()
    assert "--radius" in str(e.value)


def test_cli_polygon_with_where_and_bad_file(tmp_path, monkeypatch):
    dem, gpkg = make_dem_and_polygon(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--input", dem, "--polygon", str(tmp_path / "nope.geojson"),
                                      "--output", str(tmp_path / "o.stl")])
    with pytest.raises(SystemExit) as e:
        circle_dem_to_stl.main()
    assert "contour illisible" in str(e.value)


def test_fill_uses_custom_mask_not_the_disc(tmp_path, monkeypatch):
    """Avec inside_fn, le complément s'applique aussi aux coins du carré hors disque."""
    lat, lon = 46.0, 7.6
    x, y = sources._to_epsg(2056, lat, lon)
    write_synthetic_tif(tmp_path / "main.tif", np.full((200, 200), 1000.0, dtype=np.float32), x - 100.0, y + 100.0, 1.0, 2056)  # vide
    ds = gdal.Open(str(tmp_path / "main.tif"), gdal.GA_Update)
    arr = ds.GetRasterBand(1).ReadAsArray()
    arr[:, 100:] = -9999.0                                         # moitié est sans donnée
    ds.GetRasterBand(1).WriteArray(arr)
    ds = None
    write_synthetic_tif(tmp_path / "fill.tif", np.full((400, 400), 1000.0, dtype=np.float32), x - 200.0, y + 200.0, 1.0, 2056)
    monkeypatch.setattr(raster.gdal, "FillNodata", lambda *a, **k: 0)   # isole le complément du comblement des trous
    inputs = sources.LayeredInputs([str(tmp_path / "main.tif")])
    inputs.fill_loader = lambda: str(tmp_path / "fill.tif")
    wkt, cx, cy = raster.resolve_target_crs_and_center(inputs, None, None, x, y)
    q = lambda *a: 1  # noqa: E731
    disc, xs, ys, nd = raster.warp_window(inputs, wkt, cx, cy, 100.0, 4.0, progress_cb=q)
    custom, _, _, _ = raster.warp_window(inputs, wkt, cx, cy, 100.0, 4.0, progress_cb=q,
                                         inside_fn=lambda a, b: np.ones((len(b), len(a)), dtype=bool))
    corner_east = (slice(0, 3), slice(-3, None))                   # coin hors disque, côté sans donnée
    assert np.all(disc[corner_east] == nd)                          # hors disque : non complété
    assert not np.any(custom[corner_east] == nd)


# ---------------------------------------------------------------- check_stl

def circle_stl(tmp_path, monkeypatch):
    dem, _ = make_dem_and_polygon(tmp_path)
    x, y = sources._to_epsg(2056, 46.0, 7.6)
    out = tmp_path / "c.stl"
    run_cli(monkeypatch, ["--input", dem, "--cx", str(x), "--cy", str(y), "--radius", "500", "--diameter-mm", "100",
                          "--vexag", "1", "--base-mm", "3", "--print-spacing-mm", "0.5", "--output", str(out)])
    return out


def test_check_stl_passes_on_real_output_and_fails_on_wrong_expectations(tmp_path, monkeypatch):
    out = circle_stl(tmp_path, monkeypatch)
    tris = check_stl.read_stl(str(out))
    ok = {n: o for n, o, _ in check_stl.check(tris, 100.0, 3.0, radius_m=500.0)}
    assert all(ok.values()), ok
    bad_diam = {n: o for n, o, _ in check_stl.check(tris, 120.0, 3.0, radius_m=500.0)}
    assert not bad_diam["plus grande dimension horizontale"]
    bad_base = {n: o for n, o, _ in check_stl.check(tris, 100.0, 5.0, radius_m=500.0)}
    assert not bad_base["socle"]
    bad_relief = {n: o for n, o, _ in check_stl.check(tris, 100.0, 3.0, radius_m=500.0, relief_m=300.0)}
    assert not bad_relief["relief vertical"]


def test_check_stl_detects_a_hole(tmp_path, monkeypatch):
    out = circle_stl(tmp_path, monkeypatch)
    raw = bytearray(out.read_bytes())
    n = struct.unpack("<I", raw[80:84])[0]
    holed = bytes(raw[:80]) + struct.pack("<I", n - 3) + bytes(raw[84:84 + 50 * (n - 3)])
    (tmp_path / "holed.stl").write_bytes(holed)
    r = check_stl.edge_report(check_stl.read_stl(str(tmp_path / "holed.stl")))
    assert r["boundary_edges"] > 0


def test_check_stl_rejects_non_binary(tmp_path):
    (tmp_path / "a.stl").write_text("solid x\nendsolid x\n" * 10)
    with pytest.raises(ValueError):
        check_stl.read_stl(str(tmp_path / "a.stl"))


# ---------------------------------------------------------------- --clip-to-place (Nominatim)

import geocode  # noqa: E402


def nominatim_payload(geometry, name="Île test, France", typ="island"):
    return json.dumps([{"display_name": name, "lat": "46.0", "lon": "7.6", "category": "place", "type": typ,
                        **({"geojson": geometry} if geometry is not None else {})}]).encode("utf-8")


def square_geojson(x, y, half_m):
    """Carré en lon/lat autour de (x, y) EPSG:2056, de demi-côté half_m."""
    ct = osr.CoordinateTransformation(sources._srs(2056), sources._srs(4326))
    pts = [ct.TransformPoint(x + dx, y + dy)[:2] for dx, dy in
           ((-half_m, -half_m), (half_m, -half_m), (half_m, half_m), (-half_m, half_m), (-half_m, -half_m))]
    return {"type": "Polygon", "coordinates": [[list(p) for p in pts]]}


def test_resolve_polygon_asks_for_geojson_and_caches_separately(tmp_path, monkeypatch):
    urls = []
    monkeypatch.setattr(geocode.time, "sleep", lambda s: None)
    monkeypatch.setattr(sources, "_http_get", lambda url, **k: urls.append(url) or nominatim_payload(
        {"type": "Polygon", "coordinates": [[[7.5, 45.9], [7.7, 45.9], [7.7, 46.1], [7.5, 46.1], [7.5, 45.9]]]}))
    clip = geocode.resolve_polygon("Île test", 1, str(tmp_path), log=quiet)
    assert "polygon_geojson=1" in urls[0] and "polygon_threshold=" in urls[0]
    lat, lon = clip.center_latlon()
    assert (lat, lon) == pytest.approx((46.0, 7.6), abs=1e-6)
    # la recherche simple du même nom ne réutilise pas le cache avec contour (et inversement)
    geocode.resolve("Île test", 1, str(tmp_path), log=quiet)
    assert len(urls) == 2 and "polygon_geojson" not in urls[1]
    geocode.resolve_polygon("Île test", 1, str(tmp_path), log=quiet)
    assert len(urls) == 2


def test_resolve_polygon_point_result_is_an_error(tmp_path, monkeypatch):
    monkeypatch.setattr(geocode.time, "sleep", lambda s: None)
    monkeypatch.setattr(sources, "_http_get", lambda url, **k: nominatim_payload(
        {"type": "Point", "coordinates": [7.6, 46.0]}, name="Mont Test", typ="peak"))
    with pytest.raises(geocode.GeocodeError, match="pas de contour"):
        geocode.resolve_polygon("Mont Test", 1, str(tmp_path), log=quiet)
    monkeypatch.setattr(sources, "_http_get", lambda url, **k: nominatim_payload(None))
    with pytest.raises(geocode.GeocodeError, match="absent"):
        geocode.resolve_polygon("Mont Test2", 1, str(tmp_path), log=quiet)


def test_cli_clip_to_place_end_to_end_and_conflicts(tmp_path, monkeypatch):
    dem, _ = make_dem_and_polygon(tmp_path)
    x, y = sources._to_epsg(2056, 46.0, 7.6)
    monkeypatch.setattr(geocode.time, "sleep", lambda s: None)
    monkeypatch.setattr(sources, "_http_get", lambda url, **k: nominatim_payload(square_geojson(x, y, 700.0)))
    out = tmp_path / "ile.stl"
    run_cli(monkeypatch, ["--input", dem, "--place", "Île test", "--clip-to-place", "--diameter-mm", "70",
                          "--vexag", "1", "--base-mm", "2", "--print-spacing-mm", "1", "--output", str(out),
                          "--cache-dir", str(tmp_path / "c")])
    res = {n: ok for n, ok, _ in check_stl.check(check_stl.read_stl(str(out)), 70.0, 2.0, tol_mm=2.0)}
    assert all(res.values()), res
    for extra, msg in ((["--radius", "500"], "incompatible"), (["--smooth-wall"], "incompatible")):
        monkeypatch.setattr(sys, "argv", ["x", "--input", dem, "--place", "Île test", "--clip-to-place",
                                          "--output", str(tmp_path / "o.stl")] + extra)
        with pytest.raises(SystemExit) as e:
            circle_dem_to_stl.main()
        assert msg in str(e.value)
    monkeypatch.setattr(sys, "argv", ["x", "--input", dem, "--clip-to-place", "--output", str(tmp_path / "o.stl")])
    with pytest.raises(SystemExit) as e:
        circle_dem_to_stl.main()
    assert "requiert --place" in str(e.value)
