# SPDX-License-Identifier: GPL-3.0-or-later
"""Nuages de points LAS/LAZ (lidar.py).

Le nuage synthétique est un RÉSEAU RÉGULIER de points (pas 0,5 m, décalé d'un quart de pas) sur un plan incliné
z = 100 + 0,10 dx + 0,05 dy. Avec des cellules de 2 m alignées sur le réseau, chaque cellule contient 4 x 4 points
disposés symétriquement autour de son centre : la moyenne du sol y vaut EXACTEMENT le plan au centre de la cellule, le
minimum vaut le plan moins 0,1125 m, le maximum le plan plus 0,1125 m. Toute erreur mesurée vient donc de la logique
(classes, doublons, aberrants, comblement), jamais du hasard de l'échantillonnage."""
import struct
import sys
from pathlib import Path

import numpy as np
import pytest
from osgeo import osr

laspy = pytest.importorskip("laspy", reason="laspy n'est pas installé (pip install laspy lazrs) : tests LAS/LAZ ignorés")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import lidar  # noqa: E402
import raster  # noqa: E402
from fmt_helpers import assert_stl_ok, cli, prepared, window  # noqa: E402

X0L, Y0L = 500000.0, 4600000.0                     # WGS 84 / UTM 32N
NODATA = raster.FILL_NODATA
HALF_SPREAD = 0.75 * (0.10 + 0.05)                 # écart du coin du réseau à son centre de cellule : 0,1125 m


def ground_z(x, y):
    return 100.0 + 0.10 * (np.asarray(x) - X0L) + 0.05 * (np.asarray(y) - Y0L)


def lattice(size=200.0, step=0.5, x0=X0L, y0=Y0L):
    i = (np.arange(int(size / step)) + 0.5) * step
    gx, gy = np.meshgrid(x0 + i, y0 + i)
    return gx.ravel(), gy.ravel()


def crs_wkt(code):
    s = osr.SpatialReference()
    s.SetFromUserInput(code)
    return s.ExportToWkt()


def write_las(path, x, y, z, cls, *, version="1.4", fmt=6, wkt="EPSG:32632", geokeys=None, compressed=False,
              withheld=None, scales=(0.001, 0.001, 0.001), evlr=False):
    hdr = laspy.LasHeader(version=version, point_format=fmt)
    hdr.scales = list(scales)
    hdr.offsets = [float(np.floor(x.min())), float(np.floor(y.min())), 0.0]
    las = laspy.LasData(hdr)
    if wkt is not None:
        blob = crs_wkt(wkt).encode() + b"\x00"
        if evlr:
            from laspy.vlrs.vlrlist import VLRList
            las.evlrs = VLRList()
            las.evlrs.append(laspy.VLR(user_id="LASF_Projection", record_id=2112, description="WKT", record_data=blob))
        else:
            las.header.vlrs.append(laspy.VLR(user_id="LASF_Projection", record_id=2112, description="WKT",
                                             record_data=blob))
    if geokeys is not None:
        keys = [(1024, 0, 1, 1)] + [(k, 0, 1, v) for k, v in geokeys.items()]
        blob = struct.pack("<4H", 1, 1, 0, len(keys)) + b"".join(struct.pack("<4H", *k) for k in keys)
        las.header.vlrs.append(laspy.VLR(user_id="LASF_Projection", record_id=34735, description="GeoKeys",
                                         record_data=blob))
    las.x, las.y, las.z = x, y, z
    las.classification = np.asarray(cls, dtype=np.uint8)
    if withheld is not None:
        las.withheld = np.asarray(withheld, dtype=bool)
    las.write(str(path), do_compress=compressed)
    return str(path)


def cloud(*, classified=True, hole=True, veg=True, noise=True, withheld_cell=False, duplicates=False, spike=False):
    """(x, y, z, classe, withheld) du nuage de référence ; les défauts sont ajoutés à la demande."""
    gx, gy = lattice()
    if hole:                                              # « bâtiment » sans point sol : [80, 110] x [80, 110]
        inb = (gx > X0L + 80) & (gx < X0L + 110) & (gy > Y0L + 80) & (gy < Y0L + 110)
        gx, gy = gx[~inb], gy[~inb]
    X, Y, Z = [gx], [gy], [ground_z(gx, gy)]
    C = [np.full(gx.size, 2 if classified else 1)]
    W = [np.zeros(gx.size, bool)]

    def add(x, y, z, c, w=False):
        x, y = (np.atleast_1d(np.asarray(v, dtype=np.float64)) for v in (x, y))
        z = np.broadcast_to(np.asarray(z, dtype=np.float64), x.shape)
        X.append(x); Y.append(y); Z.append(z)
        C.append(np.full(x.size, c)); W.append(np.full(x.size, w, dtype=bool))

    if veg:                                               # végétation +10 m dans [20, 60] x [20, 60]
        inp = (gx > X0L + 20) & (gx < X0L + 60) & (gy > Y0L + 20) & (gy < Y0L + 60)
        add(gx[inp], gy[inp], ground_z(gx[inp], gy[inp]) + 10.0, 5 if classified else 1)
    if noise:                                             # bruit bas (classe 7), 80 m sous le sol
        rng = np.random.default_rng(1)
        nx, ny = X0L + rng.uniform(20, 180, 40), Y0L + rng.uniform(20, 180, 40)
        add(nx, ny, np.full(40, 20.0), 7)
    if withheld_cell:                                     # 5 points « sol » retenus par le fournisseur, +300 m
        add(np.full(5, X0L + 150.6), np.full(5, Y0L + 150.6), ground_z(X0L + 150.6, Y0L + 150.6) + 300.0, 2, True)
    if duplicates:                                        # le coin (+0,75 ; +0,75) de la cellule [60, 62) x [60, 62), x200
        px, py = X0L + 61.75, Y0L + 61.75
        add(np.full(200, px), np.full(200, py), np.full(200, ground_z(px, py)), 2 if classified else 1)
    if spike:                                             # un point « sol » à +400 m
        add(X0L + 131.3, Y0L + 41.3, ground_z(X0L + 131.3, Y0L + 41.3) + 400.0, 2 if classified else 1)
    return (np.concatenate(X), np.concatenate(Y), np.concatenate(Z), np.concatenate(C).astype(np.uint8),
            np.concatenate(W))


def make(path, **kw):
    laz = kw.pop("compressed", False)
    las_kw = {k: kw.pop(k) for k in list(kw) if k in ("version", "fmt", "wkt", "geokeys", "scales", "evlr")}
    x, y, z, c, w = cloud(**kw)
    return write_las(path, x, y, z, c, withheld=w, compressed=laz, **las_kw)


def lidar_window(path, *, cx=X0L + 100, cy=Y0L + 100, radius=90.0, pixel=2.0, fill=True, log=None, **opts):
    kind, inp = prepared(path, log=log, **opts)
    assert kind == "lidar"
    elev, xs, ys, nodata, wkt = window(inp, cx=cx, cy=cy, radius=radius, pixel=pixel, fill=fill)
    return elev, xs, ys, nodata


def truth(xs, ys):
    return ground_z(xs[None, :], ys[:, None])


def in_box(xs, ys, lo, hi):
    return ((xs[None, :] - X0L > lo) & (xs[None, :] - X0L < hi) & (ys[:, None] - Y0L > lo) & (ys[:, None] - Y0L < hi))


HOLE_CELLS = 15 * 15


# ---------------------------------------------------------------------------------------------------------------
# Sélection des points : sol / minimum / maximum
# ---------------------------------------------------------------------------------------------------------------

def test_ground_mode_returns_the_exact_terrain_ignoring_vegetation_and_noise_and_fills_the_hole(tmp_path):
    log = []
    path = make(tmp_path / "c.las")
    elev, xs, ys, nodata = lidar_window(path, log=log)
    assert elev.shape == (90, 90) and nodata == NODATA and np.isfinite(elev).all()
    assert not np.any(elev == NODATA)
    assert np.abs(elev - truth(xs, ys)).max() < 0.01              # moyenne du sol : le plan, trou compris (interpolé)
    assert elev[0].mean() > elev[-1].mean() and elev[:, -1].mean() > elev[:, 0].mean()
    text = "\n".join(log)
    assert "points/m²" in text and "225 cellule(s) sans point comblée(s) par interpolation linear" in text


def test_default_selection_is_ground_even_when_vegetation_dominates_a_cell(tmp_path):
    elev, xs, ys, _ = lidar_window(make(tmp_path / "c.las"))
    patch = in_box(xs, ys, 20, 60)
    assert patch.sum() == 20 * 20
    assert np.abs(elev[patch] - truth(xs, ys)[patch]).max() < 0.01        # pas +10 m de canopée


def test_lowest_and_highest_modes_on_an_unclassified_cloud(tmp_path):
    path = make(tmp_path / "u.las", classified=False)
    log = []
    lo, xs, ys, _ = lidar_window(path, lidar_select="lowest", log=log)
    assert np.abs(lo - (truth(xs, ys) - HALF_SPREAD)).max() < 0.01      # le bruit à 20 m n'est pas le minimum
    assert any("NON classifié" in m and "approximation du sol" in m for m in log)
    hi, xs, ys, _ = lidar_window(path, lidar_select="highest")
    patch = in_box(xs, ys, 20, 60)
    want = truth(xs, ys) + HALF_SPREAD + np.where(patch, 10.0, 0.0)
    assert np.abs(hi - want).max() < 0.01                               # la canopée est bien dans le maximum


def test_unclassified_cloud_with_default_ground_mode_is_refused_with_explicit_alternatives(tmp_path):
    path = make(tmp_path / "u.las", classified=False)
    with pytest.raises(raster.InputError, match="aucun point classé « sol »") as e:
        lidar_window(path)
    msg = str(e.value)
    assert "--lidar-select lowest" in msg and "--lidar-select highest" in msg and "classes présentes" in msg


def test_classified_cloud_with_lowest_mode_says_ground_would_be_better(tmp_path):
    log = []
    lidar_window(make(tmp_path / "c.las"), lidar_select="lowest", log=log)
    assert any("le nuage est classifié" in m and "--lidar-select ground" in m for m in log)


def test_unknown_selection_or_interpolation_is_refused(tmp_path):
    path = make(tmp_path / "c.las")
    with pytest.raises(raster.InputError, match="--lidar-select : valeur inconnue"):
        prepared(path, lidar_select="mean")
    with pytest.raises(raster.InputError, match="--lidar-interp : valeur inconnue"):
        prepared(path, lidar_interp="kriging")


# ---------------------------------------------------------------------------------------------------------------
# Points à écarter : withheld, bruit, doublons, aberrants
# ---------------------------------------------------------------------------------------------------------------

def test_withheld_points_are_excluded(tmp_path):
    """spike_m = 0 : sans cela, le nettoyage des aberrants masquerait l'effet de points mal écartés."""
    path = make(tmp_path / "w.las", withheld_cell=True)
    log = []
    elev, xs, ys, _ = lidar_window(path, lidar_spike_m=0.0, log=log)
    assert np.abs(elev - truth(xs, ys)).max() < 0.01
    assert any("écartés" in m for m in log)


def test_low_noise_class_7_is_excluded_even_in_lowest_mode(tmp_path):
    path = make(tmp_path / "u.las", classified=False)
    elev, xs, ys, _ = lidar_window(path, lidar_select="lowest", lidar_spike_m=0.0)
    assert elev.min() > 90.0                                           # sans l'exclusion : 20 m


def test_exact_duplicate_points_do_not_weigh_on_the_cell_mean(tmp_path):
    log = []
    elev, xs, ys, _ = lidar_window(make(tmp_path / "d.las", duplicates=True), lidar_spike_m=0.0, log=log)
    assert np.abs(elev - truth(xs, ys)).max() < 0.01                  # 200 copies d'un coin déplaceraient la moyenne de ~0,1 m
    assert any("200 doublons exacts" in m for m in log)


def test_an_isolated_aberrant_ground_cell_is_removed_and_interpolated(tmp_path):
    path = make(tmp_path / "s.las", spike=True)
    log = []
    elev, xs, ys, _ = lidar_window(path, log=log)
    assert np.abs(elev - truth(xs, ys)).max() < 0.01
    assert any("1 cellule(s) isolée(s) aberrante(s) écartée(s)" in m for m in log)


def test_aberrant_cell_removal_can_be_switched_off(tmp_path):
    path = make(tmp_path / "s.las", spike=True)
    elev, xs, ys, _ = lidar_window(path, lidar_spike_m=0.0)
    err = elev - truth(xs, ys)
    assert err.max() > 20.0 and np.sum(np.abs(err) > 0.01) == 1       # une seule cellule, +400/17 = +23,5 m


def test_fill_false_returns_the_raw_cells(tmp_path):
    path = make(tmp_path / "s.las", spike=True)
    elev, xs, ys, nodata = lidar_window(path, fill=False)
    assert (elev == nodata).sum() == HOLE_CELLS                       # trou non comblé
    assert (elev - truth(xs, ys))[elev != nodata].max() > 20.0        # et aberrant conservé


# ---------------------------------------------------------------------------------------------------------------
# Trous : méthodes d'interpolation
# ---------------------------------------------------------------------------------------------------------------

def test_hole_without_interpolation_stays_nodata_and_never_zero(tmp_path):
    elev, xs, ys, nodata = lidar_window(make(tmp_path / "c.las"), lidar_interp="none")
    assert (elev == nodata).sum() == HOLE_CELLS and not np.any(elev == 0.0)
    assert np.abs((elev - truth(xs, ys))[elev != nodata]).max() < 0.01


def test_hole_nearest_neighbour_fills_everything_with_a_bounded_error(tmp_path):
    elev, xs, ys, nodata = lidar_window(make(tmp_path / "c.las"), lidar_interp="nearest")
    err = np.abs(elev - truth(xs, ys))
    hole = in_box(xs, ys, 80, 110)
    assert not np.any(elev == nodata) and err[~hole].max() < 0.01
    assert 0.05 < err[hole].max() < 1.6                                # approximation documentée, pas exacte


def test_hole_larger_than_the_gap_limit_keeps_its_centre_nodata(tmp_path):
    log = []
    elev, xs, ys, nodata = lidar_window(make(tmp_path / "c.las"), lidar_max_gap=4.0, log=log)   # 2 cellules de 2 m
    assert elev[47, 42] == nodata                                     # centre du trou (distance > 2 cellules du bord)
    assert elev[47, 35] != nodata and abs(elev[47, 35] - truth(xs, ys)[47, 35]) < 0.01      # bord comblé
    assert any("restent sans donnée" in m for m in log)


def test_linear_fill_never_extrapolates_beyond_the_cloud(tmp_path):
    """Fenêtre qui déborde du nuage de 30 m sur un bord : ces cellules n'ont aucun voisin de l'autre côté."""
    elev, xs, ys, nodata = lidar_window(make(tmp_path / "c.las"), cx=X0L + 130, cy=Y0L + 100, radius=90.0)
    outside = xs[None, :] - X0L > 200
    assert np.all(elev[:, outside[0]] == nodata) and not np.any(elev[:, ~outside[0]] == nodata)


def test_sparse_cloud_warns_that_filled_detail_would_be_invented(tmp_path):
    x, y = lattice(step=4.0)                                          # 1 point par 16 m² pour des cellules de 2 m
    path = write_las(tmp_path / "sparse.las", x, y, ground_z(x, y), np.full(x.size, 2))
    log = []
    lidar_window(path, log=log, lidar_interp="none")
    assert any("moins de la moitié des cellules" in m for m in log)


# ---------------------------------------------------------------------------------------------------------------
# Cellule LiDAR différente du pas de travail
# ---------------------------------------------------------------------------------------------------------------

def test_coarser_lidar_cell_is_resampled_once_and_documented(tmp_path):
    log = []
    elev, xs, ys, nodata = lidar_window(make(tmp_path / "c.las"), lidar_cell=5.0, log=log)
    assert elev.shape == (90, 90)
    assert any("ré-échantillonnée à 2 m (bilinear)" in m for m in log)
    inner = (slice(4, -4), slice(4, -4))
    assert np.abs(elev - truth(xs, ys))[inner].max() < 0.02


def test_finer_lidar_cell_is_averaged_down_and_documented(tmp_path):
    log = []
    elev, xs, ys, nodata = lidar_window(make(tmp_path / "c.las"), lidar_cell=1.0, log=log)
    assert any("ré-échantillonnée à 2 m (average)" in m for m in log)
    assert np.abs(elev - truth(xs, ys)).max() < 0.01


# ---------------------------------------------------------------------------------------------------------------
# Gros nuages, plusieurs fichiers, LAZ, versions
# ---------------------------------------------------------------------------------------------------------------

def test_cloud_is_read_in_chunks_with_identical_result(tmp_path, monkeypatch):
    path = make(tmp_path / "c.las", duplicates=False)
    whole, xs, ys, _ = lidar_window(path)
    calls = []
    original = lidar.LidarInput._consume

    def counting(self, chunk, *a, **k):
        calls.append(len(chunk))
        return original(self, chunk, *a, **k)

    monkeypatch.setattr(lidar, "CHUNK_POINTS", 20_000)
    monkeypatch.setattr(lidar.LidarInput, "_consume", counting)
    chunked, *_ = lidar_window(path)
    assert len(calls) >= 8 and max(calls) <= 20_000                  # 160 000 points / 20 000 par bloc
    assert np.abs(chunked - whole).max() < 1e-9


def test_several_files_are_combined_like_one(tmp_path):
    x, y, z, c, w = cloud()
    left = x < X0L + 100
    a = write_las(tmp_path / "w.las", x[left], y[left], z[left], c[left], withheld=w[left])
    b = write_las(tmp_path / "e.las", x[~left], y[~left], z[~left], c[~left], withheld=w[~left])
    both, xs, ys, _ = _two_files(a, b)
    single, *_ = lidar_window(make(tmp_path / "one.las"))
    assert np.abs(both - single).max() < 1e-9


def _two_files(a, b):
    kind, inp = prepared([a, b])
    elev, xs, ys, nodata, _ = window(inp, cx=X0L + 100, cy=Y0L + 100, radius=90.0, pixel=2.0)
    return elev, xs, ys, nodata


def test_files_that_do_not_touch_the_window_are_skipped_without_reading(tmp_path):
    x, y, z, c, w = cloud()
    near = write_las(tmp_path / "near.las", x, y, z, c, withheld=w)
    far = write_las(tmp_path / "far.las", x + 5000.0, y, z, c, withheld=w)
    both, xs, ys, _ = _two_files(near, far)
    assert np.abs(both - truth(xs, ys)).max() < 0.01


def _laz_available():
    try:
        return bool(laspy.compression.LazBackend.detect_available())
    except Exception:       # noqa: BLE001 - API variable selon la version de laspy
        return False


@pytest.mark.skipif(not _laz_available(), reason="aucun décompresseur LAZ (pip install 'laspy[lazrs]')")
def test_laz_gives_the_same_grid_as_las(tmp_path):
    las, xs, ys, _ = lidar_window(make(tmp_path / "c.las"))
    log = []
    laz, *_ = lidar_window(make(tmp_path / "c.laz", compressed=True), log=log)
    assert any("compressé LAZ" in m for m in log)
    assert np.array_equal(las, laz)


def test_laz_without_decompressor_gives_an_actionable_error(tmp_path, monkeypatch):
    path = make(tmp_path / "c.las")

    def refuse(*a, **k):
        raise laspy.errors.LaspyException("No LazBackend selected, cannot decompress data")

    monkeypatch.setattr(laspy, "open", refuse)
    with pytest.raises(raster.InputError, match=r"compressé \(LAZ\).*pip install 'laspy\[lazrs\]'"):
        prepared(path)


def test_missing_laspy_is_reported_with_the_install_command(monkeypatch, tmp_path):
    path = make(tmp_path / "c.las")
    monkeypatch.setitem(sys.modules, "laspy", None)                    # import laspy -> ImportError
    with pytest.raises(raster.InputError, match="requiert laspy.*pip install laspy"):
        prepared(path)


@pytest.mark.parametrize("version, fmt", [("1.2", 0), ("1.2", 1), ("1.2", 3), ("1.3", 4), ("1.4", 6), ("1.4", 8)])
def test_las_versions_and_point_formats(tmp_path, version, fmt):
    log = []
    elev, xs, ys, _ = lidar_window(make(tmp_path / "v.las", version=version, fmt=fmt), log=log)
    assert any(f"LAS {version} format {fmt}" in m for m in log)
    assert np.abs(elev - truth(xs, ys)).max() < 0.01


# ---------------------------------------------------------------------------------------------------------------
# CRS et unités
# ---------------------------------------------------------------------------------------------------------------

def test_crs_from_geotiff_keys_only(tmp_path):
    log = []
    path = make(tmp_path / "k.las", wkt=None, geokeys={3072: 32632}, version="1.2", fmt=3)
    elev, xs, ys, _ = lidar_window(path, log=log)
    assert any("UTM zone 32N" in m for m in log)
    assert np.abs(elev - truth(xs, ys)).max() < 0.01


def test_crs_wkt_stored_in_an_extended_vlr_is_found(tmp_path):
    log = []
    path = make(tmp_path / "e.las", evlr=True)
    elev, xs, ys, _ = lidar_window(path, log=log)
    assert any("UTM zone 32N" in m for m in log)
    assert np.abs(elev - truth(xs, ys)).max() < 0.01


def test_missing_crs_is_refused_and_never_invented(tmp_path):
    path = make(tmp_path / "n.las", wkt=None)
    with pytest.raises(raster.InputError, match="aucun système de coordonnées.*--source-crs"):
        prepared(path)


def test_source_crs_completes_a_cloud_without_crs_but_never_replaces_one(tmp_path):
    elev, xs, ys, _ = lidar_window(make(tmp_path / "n.las", wkt=None), source_crs="EPSG:32632")
    assert np.abs(elev - truth(xs, ys)).max() < 0.01
    with pytest.raises(raster.InputError, match="déclare déjà un CRS"):
        prepared(make(tmp_path / "c.las"), source_crs="EPSG:32633")
    with pytest.raises(raster.InputError, match="CRS non reconnu"):
        prepared(make(tmp_path / "n2.las", wkt=None), source_crs="EPSG:99999999")


def test_degrees_declared_for_metre_coordinates_are_refused(tmp_path):
    with pytest.raises(raster.InputError, match="impossibles en degrés"):
        prepared(make(tmp_path / "w.las", wkt="EPSG:4326"))


def test_vertical_unit_option_converts_a_cloud_in_feet(tmp_path):
    x, y, z, c, w = cloud()
    path = write_las(tmp_path / "ft.las", x, y, z / 0.3048, c, withheld=w)
    log = []
    elev, xs, ys, _ = lidar_window(path, vertical_unit="ft", log=log)
    assert any("--vertical-unit ft" in m for m in log)
    assert np.abs(elev - truth(xs, ys)).max() < 0.01
    wrong, *_ = lidar_window(path)                                     # sans l'option : mètres (CRS en mètres)
    assert np.abs(wrong - truth(xs, ys)).max() > 100.0


def test_compound_crs_with_feet_vertical_unit_is_honoured_without_a_datum_change(tmp_path):
    x, y, z, c, w = cloud()
    path = write_las(tmp_path / "comp.las", x, y, z / (1200.0 / 3937.0), c, withheld=w, wkt="EPSG:32632+6360")
    log = []
    elev, xs, ys, _ = lidar_window(path, log=log)
    assert any("AUCUN changement de datum vertical" in m for m in log)
    assert np.abs(elev - truth(xs, ys)).max() < 0.01


def test_cloud_in_us_feet_state_plane_is_reprojected_and_converted(tmp_path):
    """Colorado Nord (EPSG:2231), coordonnées ET altitudes en pieds US ; fenêtre en UTM 13N. Lire les pieds comme des
    mètres fausserait l'échelle de 3,28 et les altitudes de ~3 400 m."""
    xc, yc, _ = osr.CoordinateTransformation(*_srs_pair(4326, 2231)).TransformPoint(-105.0, 40.0)
    gx, gy = lattice(size=400.0, step=1.0, x0=xc - 200.0, y0=yc - 200.0)         # 400 ft x 400 ft, pas 1 ft
    z_ft = 5000.0 + 0.10 * (gx - xc) + 0.05 * (gy - yc)
    path = write_las(tmp_path / "co.las", gx, gy, z_ft, np.full(gx.size, 2), wkt="EPSG:2231")
    log = []
    kind, inp = prepared(path, log=log)
    assert any("avertissement : unité verticale non déclarée, supposée identique à l'unité horizontale" in m for m in log)
    elev, xs, ys, nodata, wkt = window(inp, lat=40.0, lon=-105.0, radius=50.0, pixel=2.0)
    assert "UTM zone 13N" in wkt and not np.any(elev == nodata)
    expected = 5000.0 * 1200.0 / 3937.0                                # au centre : 5000 ft US en mètres
    centre = elev[elev.shape[0] // 2 - 2:elev.shape[0] // 2 + 2, elev.shape[1] // 2 - 2:elev.shape[1] // 2 + 2]
    assert abs(centre.mean() - expected) < 0.3
    gradient = np.hypot(np.polyfit(xs, elev.mean(axis=0), 1)[0], np.polyfit(ys, elev.mean(axis=1), 1)[0])
    assert abs(gradient - np.hypot(0.10, 0.05)) < 0.01                 # pente sans dimension : inchangée par les unités


def _srs_pair(a, b):
    s1, s2 = osr.SpatialReference(), osr.SpatialReference()
    s1.ImportFromEPSG(a); s2.ImportFromEPSG(b)
    s1.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER); s2.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return s1, s2


def test_geographic_cloud_is_reprojected_and_the_missing_vertical_unit_is_flagged(tmp_path):
    lon0, lat0 = 7.0, 46.0
    step = np.arange(400) * 0.000004 + 0.000002                        # ≈ 0,3 m en longitude, 0,44 m en latitude
    gx, gy = np.meshgrid(lon0 + step, lat0 + step)
    gx, gy = gx.ravel(), gy.ravel()
    z = 1000.0 + 0.05 * (gx - lon0) / 0.000004 * 0.3 + 0.02 * (gy - lat0) / 0.000004 * 0.44
    path = write_las(tmp_path / "geo.las", gx, gy, z, np.full(gx.size, 2), wkt="EPSG:4326",
                     scales=(1e-8, 1e-8, 0.001))
    log = []
    kind, inp = prepared(path, log=log)
    assert any("CRS géographique sans unité verticale" in m for m in log)
    cx_lon, cx_lat = lon0 + 0.0008, lat0 + 0.0008
    elev, xs, ys, nodata, wkt = window(inp, lat=cx_lat, lon=cx_lon, radius=40.0, pixel=2.0)
    assert "UTM zone 32N" in wkt and not np.any(elev == nodata)
    assert 1000.0 < elev.min() and elev.max() < 1000.0 + 0.05 * 400 * 0.3 + 0.02 * 400 * 0.44 + 1


# ---------------------------------------------------------------------------------------------------------------
# Fichiers inutilisables
# ---------------------------------------------------------------------------------------------------------------

def test_window_outside_the_cloud_is_reported_with_both_extents(tmp_path):
    with pytest.raises(raster.InputError, match="aucun point du nuage dans la zone demandée") as e:
        lidar_window(make(tmp_path / "c.las"), cx=X0L + 90000.0)
    assert "emprise du nuage" in str(e.value)


def test_truncated_las_is_reported(tmp_path):
    path = Path(make(tmp_path / "c.las"))
    data = path.read_bytes()
    path.write_bytes(data[:len(data) // 2])
    with pytest.raises(raster.InputError, match="tronqué|interrompue|illisible|corrompu"):
        lidar_window(str(path))


def test_garbage_with_las_extension_is_reported(tmp_path):
    (tmp_path / "g.las").write_bytes(b"LASF" + bytes(range(256)) * 4)
    with pytest.raises(raster.InputError, match="illisible ou corrompu"):
        prepared(tmp_path / "g.las")
    (tmp_path / "t.las").write_bytes(b"not a point cloud at all")
    with pytest.raises(raster.InputError, match="illisible ou corrompu"):
        prepared(tmp_path / "t.las")


def test_header_with_absurd_sizes_is_rejected_before_laspy_can_allocate(tmp_path):
    """Des octets de bruit lus comme tailles d'enregistrements ont fait demander des gigaoctets à laspy."""
    good = Path(make(tmp_path / "c.las")).read_bytes()
    bad = bytearray(good)
    struct.pack_into("<I", bad, 100, 4_000_000_000)                    # nombre de VLR
    (tmp_path / "vlr.las").write_bytes(bytes(bad))
    with pytest.raises(raster.InputError, match="enregistrements VLR annoncés"):
        prepared(tmp_path / "vlr.las")
    bad = bytearray(good)
    bad[24], bad[25] = 9, 9                                            # version 9.9
    (tmp_path / "ver.las").write_bytes(bytes(bad))
    with pytest.raises(raster.InputError, match="version LAS 9.9 inconnue"):
        prepared(tmp_path / "ver.las")
    (tmp_path / "cut.las").write_bytes(good[:200])                     # coupé dans l'en-tête
    with pytest.raises(raster.InputError, match="en-tête LAS absent ou incomplet"):
        prepared(tmp_path / "cut.las")


def test_empty_cloud_is_reported(tmp_path):
    path = write_las(tmp_path / "empty.las", np.array([X0L]), np.array([Y0L]), np.array([1.0]), np.array([2]))
    las = laspy.read(path)
    las.points = las.points[:0]
    las.write(path)
    with pytest.raises(raster.InputError, match="vide|illisible"):
        prepared(path)


def test_raster_options_that_make_no_sense_for_a_cloud_are_refused(tmp_path):
    path = make(tmp_path / "c.las")
    with pytest.raises(raster.InputError, match="--layer ne s'applique pas aux fichiers LAS/LAZ"):
        prepared(path, layer="x")
    with pytest.raises(raster.InputError, match="--band ne s'applique pas aux fichiers LAS/LAZ"):
        prepared(path, band=2)


def test_clouds_with_different_crs_are_not_silently_combined(tmp_path):
    x, y, z, c, w = cloud()
    a = write_las(tmp_path / "a.las", x, y, z, c, withheld=w)
    b = write_las(tmp_path / "b.las", x, y, z, c, withheld=w, wkt="EPSG:32633")
    with pytest.raises(raster.InputError, match="n'ont pas le même CRS"):
        prepared([a, b])


# ---------------------------------------------------------------------------------------------------------------
# Opérations pures sur la grille
# ---------------------------------------------------------------------------------------------------------------

def test_despike_flags_isolated_spikes_and_pits_but_not_slopes_or_cliffs():
    z = np.full((15, 15), 100.0)
    valid = np.ones_like(z, dtype=bool)
    z[3, 3] += 50.0                                                     # pic isolé
    z[10, 10] -= 50.0                                                   # puits isolé
    flagged = lidar.despike(z, valid, 10.0)
    assert flagged[3, 3] and flagged[10, 10] and flagged.sum() == 2
    ramp = np.tile(np.linspace(0, 400, 15), (15, 1)) + 100.0            # 28 m par cellule, versant raide
    assert not lidar.despike(ramp, valid, 10.0).any()
    cliff = np.where(np.arange(15)[None, :] < 7, 100.0, 140.0) * np.ones((15, 1))
    assert not lidar.despike(cliff, valid, 10.0).any()
    assert not lidar.despike(z, valid, 0.0).any()                       # tol = 0 : désactivé


def test_despike_needs_five_valid_neighbours():
    z = np.full((5, 5), 100.0)
    valid = np.zeros_like(z, dtype=bool)
    valid[2, 2] = True
    valid[1, 1] = valid[1, 2] = True                                    # 2 voisines seulement
    z[2, 2] = 500.0
    assert not lidar.despike(z, valid, 10.0).any()


def test_fill_gaps_linear_reproduces_a_plane_exactly_and_nearest_is_bounded():
    i, j = np.indices((40, 40))
    z = 100.0 + 0.3 * i + 0.1 * j
    valid = np.ones_like(z, dtype=bool)
    valid[15:25, 15:25] = False
    zin = np.where(valid, z, np.nan)
    out, ok = lidar.fill_gaps(zin, valid, "linear", 50)
    assert ok.all() and np.abs(out - z).max() < 1e-9
    out, ok = lidar.fill_gaps(zin, valid, "nearest", 50)
    assert ok.all() and np.abs(out - z).max() <= 0.3 * 5 + 0.1 * 5 + 1e-9


def test_fill_gaps_honours_the_distance_limit_and_never_extrapolates():
    i, j = np.indices((40, 40))
    z = 100.0 + 0.3 * i + 0.1 * j
    valid = np.ones_like(z, dtype=bool)
    valid[10:30, 10:30] = False
    out, ok = lidar.fill_gaps(np.where(valid, z, np.nan), valid, "linear", 3)
    assert ok[10, 10] and not ok[20, 20] and np.isnan(out[20, 20])     # à 3 cellules : rempli ; au centre : non
    far = np.zeros((20, 20), dtype=bool)
    far[:, :10] = True                                                  # moitié valide, moitié vide
    zf = np.where(far, 100.0, np.nan)
    out, ok = lidar.fill_gaps(zf, far, "linear", 50)
    assert not ok[:, 12:].any()                                         # hors de l'enveloppe des cellules valides
    same, ok2 = lidar.fill_gaps(np.full((5, 5), 7.0), np.ones((5, 5), bool), "linear", 50)
    assert ok2.all() and (same == 7.0).all()


# ---------------------------------------------------------------------------------------------------------------
# Bout en bout : nuage -> circle_dem_to_stl.main() -> STL relu et mesuré
# ---------------------------------------------------------------------------------------------------------------

CLI_CENTRE = dict(cx=X0L + 100.0, cy=Y0L + 100.0)
CLI_GRADIENT = (0.10, 0.05)


def go(monkeypatch, tmp_path, files, extra=(), **kw):
    out = cli(monkeypatch, tmp_path, files, extra, center=CLI_CENTRE, radius=80.0, **kw)
    assert out.exists()
    return out


def test_cli_las_ground(monkeypatch, tmp_path):
    assert_stl_ok(go(monkeypatch, tmp_path, [make(tmp_path / "c.las")]), 80.0, CLI_GRADIENT)


@pytest.mark.skipif(not _laz_available(), reason="aucun décompresseur LAZ (pip install 'laspy[lazrs]')")
def test_cli_laz_ground(monkeypatch, tmp_path):
    assert_stl_ok(go(monkeypatch, tmp_path, [make(tmp_path / "c.laz", compressed=True)]), 80.0, CLI_GRADIENT)


def test_cli_unclassified_las_is_refused_then_accepted_with_an_explicit_method(monkeypatch, tmp_path):
    path = make(tmp_path / "u.las", classified=False)
    with pytest.raises(SystemExit) as e:
        go(monkeypatch, tmp_path, [path])
    assert "--lidar-select" in str(e.value.code) and not (tmp_path / "out.stl").exists()
    assert_stl_ok(go(monkeypatch, tmp_path, [path], ["--lidar-select", "lowest"], name="low.stl"), 80.0, CLI_GRADIENT)


def test_cli_las_with_coarser_cell_and_nearest_fill(monkeypatch, tmp_path):
    out = go(monkeypatch, tmp_path, [make(tmp_path / "c.las")], ["--lidar-cell", "4", "--lidar-interp", "nearest"])
    assert_stl_ok(out, 80.0, CLI_GRADIENT)


def test_cli_las_without_crs_needs_source_crs(monkeypatch, tmp_path):
    path = make(tmp_path / "n.las", wkt=None)
    with pytest.raises(SystemExit) as e:
        go(monkeypatch, tmp_path, [path])
    assert "--source-crs" in str(e.value.code)
    assert_stl_ok(go(monkeypatch, tmp_path, [path], ["--source-crs", "EPSG:32632"], name="ok.stl"), 80.0, CLI_GRADIENT)


def test_cli_las_with_the_hole_left_open_still_gives_a_watertight_model(monkeypatch, tmp_path):
    """Sans comblement, le trou devient « sans donnée » dans la grille ; le moteur de maillage le gère (trou exclu)."""
    out = go(monkeypatch, tmp_path, [make(tmp_path / "c.las")], ["--lidar-interp", "none"])
    import check_stl
    tris = check_stl.read_stl(str(out))
    edges = check_stl.edge_report(tris)
    assert edges["nonmanifold_edges"] == 0 and edges["duplicate_directed_edges"] == 0


@pytest.mark.parametrize("sides", [4, 6])
def test_cli_regular_polygon_mode_with_a_cloud(monkeypatch, tmp_path, sides):
    out = go(monkeypatch, tmp_path, [make(tmp_path / "c.las")], ["--sides", str(sides), "--circle-mode", "inscribed"])
    assert_stl_ok(out, 80.0, CLI_GRADIENT, polygon=sides)
