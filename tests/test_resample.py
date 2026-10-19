# SPDX-License-Identifier: GPL-3.0-or-later
"""Choix du rééchantillonnage selon le sens : réduction (average) ou agrandissement (bilinear).

Défaut corrigé : `average` de GDAL, quand le pixel de travail est plus petit que celui de la source, se
comporte comme `nearest` (blocs constants, donc terrasses et pente nulle sur la plupart des cellules)."""
import sys
from pathlib import Path

import numpy as np
import pytest
from osgeo import osr
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import circle_dem_to_stl  # noqa: E402
import raster  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402


def wkt(epsg):
    s = osr.SpatialReference()
    s.ImportFromEPSG(epsg)
    return s.ExportToWkt()


def rough_dem(n, seed=1, amp=60.0):
    return 500.0 + ndimage.gaussian_filter(np.random.default_rng(seed).normal(size=(n, n)), 2.0) * amp


def flat_fraction(a):
    gy, gx = np.gradient(a)
    return float(np.mean((gx == 0) & (gy == 0)))


# --- taille du pixel de la source ------------------------------------------------------------------

def test_source_pixel_size_projected_raster(tmp_path):
    p = tmp_path / "s.tif"
    write_synthetic_tif(p, rough_dem(20), 600000.0, 5100000.0, 30.0, 32631)
    assert raster.source_pixel_size_m(str(p), wkt(32631), 600300.0, 5099700.0) == pytest.approx(30.0, rel=1e-3)


def test_source_pixel_size_geographic_raster_uses_the_coarser_axis(tmp_path):
    p = tmp_path / "g.tif"
    write_synthetic_tif(p, rough_dem(20), 6.0, 46.5, 0.001, 4326)       # 0,001° : ~111 m en y, ~77 m en x à 46°N
    x, y = 6.01, 46.49
    ct = osr.CoordinateTransformation(_srs(4326), _srs(32632))
    cx, cy, _ = ct.TransformPoint(x, y)
    size = raster.source_pixel_size_m(str(p), wkt(32632), cx, cy)
    assert size == pytest.approx(111.2, rel=0.02)


def _srs(epsg):
    s = osr.SpatialReference()
    s.ImportFromEPSG(epsg)
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return s


def test_source_pixel_size_web_mercator_is_ground_size_not_projected_size(tmp_path):
    p = tmp_path / "wm.tif"
    write_synthetic_tif(p, rough_dem(20), 667000.0, 8400000.0, 100.0, 3857)     # ~ latitude 60°N
    ct = osr.CoordinateTransformation(_srs(3857), _srs(32635))
    cx, cy, _ = ct.TransformPoint(667500.0, 8399500.0)
    size = raster.source_pixel_size_m(str(p), wkt(32635), cx, cy)
    assert 45.0 < size < 55.0           # ~ 100 m x cos(60°), pas 100 m


def test_source_pixel_size_takes_the_finest_of_a_list(tmp_path):
    a, b = tmp_path / "a.tif", tmp_path / "b.tif"
    write_synthetic_tif(a, rough_dem(20), 600000.0, 5100000.0, 30.0, 32631)
    write_synthetic_tif(b, rough_dem(20), 600600.0, 5100000.0, 5.0, 32631)
    assert raster.source_pixel_size_m([str(a), str(b)], wkt(32631), 600300.0, 5099800.0) == pytest.approx(5.0, rel=1e-3)


def test_source_pixel_size_unreadable_source_is_none(tmp_path):
    p = tmp_path / "x.tif"
    p.write_bytes(b"pas un raster")
    assert raster.source_pixel_size_m(str(p), wkt(32631), 0.0, 0.0) is None


# --- choix de l'algorithme -------------------------------------------------------------------------

@pytest.mark.parametrize("src,px,expected", [
    (30.0, 5.0, "bilinear"),      # agrandissement x6
    (30.0, 28.0, "bilinear"),     # rapport 1,07 > 1,05
    (30.0, 30.0, "average"),      # même résolution
    (30.0, 29.0, "average"),      # rapport 1,03 : dans la tolérance
    (1.0, 10.0, "average"),       # réduction
    (None, 5.0, "average"),       # source inconnue : comportement historique
])
def test_auto_choice(src, px, expected):
    assert raster.choose_resample_alg("auto", src, px) == expected


@pytest.mark.parametrize("mode", ["average", "bilinear", "cubic", "cubicspline"])
def test_explicit_choice_is_honoured(mode):
    assert raster.choose_resample_alg(mode, 30.0, 5.0) == mode and raster.choose_resample_alg(mode, 1.0, 50.0) == mode


# --- effet sur le terrain ---------------------------------------------------------------------------

def _warp(path, **kw):
    return raster.warp_window(str(path), wkt(32631), 600600.0, 5099400.0, 500.0, 5.0,
                              progress_cb=lambda *a: 1, fill=False, **kw)[0]


@pytest.fixture
def coarse_dem(tmp_path):
    p = tmp_path / "coarse.tif"
    write_synthetic_tif(p, rough_dem(60), 600000.0, 5100000.0, 30.0, 32631)    # pixel de 30 m, grille de travail à 5 m
    return p


def test_auto_removes_the_terraces_that_average_creates_when_the_source_is_enlarged(coarse_dem):
    avg = _warp(coarse_dem, resample="average")
    auto = _warp(coarse_dem)                              # défaut = auto
    assert flat_fraction(avg) > 0.4                       # `average` : blocs constants de 6x6 cellules
    assert flat_fraction(auto) < 0.01                     # auto : surface continue
    assert not np.array_equal(avg, auto)


def test_auto_leaves_downsampling_unchanged(tmp_path):
    p = tmp_path / "fine.tif"
    write_synthetic_tif(p, rough_dem(400), 600000.0, 5100000.0, 2.0, 32631)    # pixel de 2 m, grille à 5 m : réduction
    assert np.array_equal(_warp(p), _warp(p, resample="average"))


def test_bilinear_never_overshoots_source_values(coarse_dem):
    out = _warp(coarse_dem, resample="bilinear")
    src_min, src_max = rough_dem(60).min(), rough_dem(60).max()
    assert out.min() >= src_min - 1e-3 and out.max() <= src_max + 1e-3


def test_unknown_resample_value_is_an_input_error(coarse_dem):
    with pytest.raises(raster.InputError, match="--resample"):
        _warp(coarse_dem, resample="lanczos")


def test_build_mesh_logs_the_choice_only_when_not_average(coarse_dem, tmp_path, capsys):
    kw = dict(lat=None, lon=None, cx=600600.0, cy=5099400.0, radius=500.0, diameter_mm=40.0, vexag=1.0, base_mm=2.0,
              print_spacing_mm=0.8)                      # pixel de travail = 1000 m / 40 mm x 0,8 mm = 20 m < 30 m (source)
    circle_dem_to_stl.build_mesh(str(coarse_dem), **kw)
    assert "ré-échantillonnage bilinear" in capsys.readouterr().out
    fine = tmp_path / "fine.tif"
    write_synthetic_tif(fine, rough_dem(400), 600000.0, 5100000.0, 2.0, 32631)
    circle_dem_to_stl.build_mesh(str(fine), **kw)
    assert "ré-échantillonnage bilinear" not in capsys.readouterr().out


def test_cli_exposes_the_option_with_auto_as_default():
    assert "auto" in raster.RESAMPLE_CHOICES and raster.RESAMPLE_CHOICES[0] == "auto"
    parser_help = open(Path(__file__).resolve().parent.parent / "circle_dem_to_stl.py", encoding="utf-8").read()
    assert '"--resample"' in parser_help and 'default="auto"' in parser_help
