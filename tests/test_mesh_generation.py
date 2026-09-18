"""
Tests d'intégration bout-en-bout : DEM synthétique -> build_mesh() -> STL.

Ces deux propriétés (étanchéité + normales cohérentes) sont celles
vérifiées par stl-audit tout au long du développement ; les régler ici
en tests automatisés évite de repasser par un audit externe à chaque
modification du script.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from circle_dem_to_stl import build_mesh, write_stl_binary  # noqa: E402

from conftest import (  # noqa: E402
    make_cone_dem,
    read_stl_binary,
    assert_watertight,
    signed_volume,
)


def test_cone_watertight_and_normals_consistent(tmp_path, tmp_stl):
    tif_path = tmp_path / "cone.tif"
    cx, cy = make_cone_dem(tif_path)

    tris = build_mesh(
        str(tif_path), cx, cy, radius=140, diameter_mm=56, vexag=3.0,
        base_mm=3.0, print_spacing_mm=0.5,
        sea_level=0.0, min_elevation=-50.0,
    )
    write_stl_binary(tris, str(tmp_stl))

    read_back = read_stl_binary(tmp_stl)
    assert_watertight(read_back)
    assert signed_volume(read_back) > 0, "volume signé négatif : normales incohérentes"


def test_cone_peak_height_matches_scale(tmp_path, tmp_stl):
    """Le point le plus haut du maillage doit correspondre à (altitude
    réelle - altitude min retenue) * échelle * vexag, à la résolution de
    ré-échantillonnage près."""
    tif_path = tmp_path / "cone.tif"
    cx, cy = make_cone_dem(tif_path, peak=80.0, slope=0.6, island_radius=100)

    diameter_mm, radius_m, vexag = 56, 140, 3.0
    scale = (diameter_mm / 2.0) / radius_m

    tris = build_mesh(
        str(tif_path), cx, cy, radius=radius_m, diameter_mm=diameter_mm, vexag=vexag,
        base_mm=3.0, print_spacing_mm=0.5,
        sea_level=0.0, min_elevation=-50.0,
    )
    write_stl_binary(tris, str(tmp_stl))
    read_back = read_stl_binary(tmp_stl)

    z_max = read_back[:, :, 2].max()
    expected = 80.0 * scale * vexag  # min_elev = 0 (niveau de la mer)
    assert abs(z_max - expected) < 1.0, f"Z max = {z_max}, attendu ~{expected}"


def test_base_thickness_applied(tmp_path, tmp_stl):
    tif_path = tmp_path / "cone.tif"
    cx, cy = make_cone_dem(tif_path)

    tris = build_mesh(
        str(tif_path), cx, cy, radius=140, diameter_mm=56, vexag=1.0,
        base_mm=7.5, print_spacing_mm=1.0,
        sea_level=0.0, min_elevation=-50.0,
    )
    write_stl_binary(tris, str(tmp_stl))
    read_back = read_stl_binary(tmp_stl)
    assert read_back[:, :, 2].min() == -7.5
