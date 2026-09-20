"""
Tests du mode île (--sea-level) : les deux cas d'artefacts rencontrés sur
de vraies données (La Réunion) et corrigés au fil du développement.

1. Un artefact isolé (balise, point mal classé) à distance de l'île :
   séparable par composante connexe seule.
2. Du bruit de surface d'eau directement raccordé au littoral : la
   connectivité seule échoue, il faut le critère --land-threshold.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from circle_dem_to_stl import build_mesh
from stl_io import write_stl_binary  # noqa: E402

from conftest import (  # noqa: E402
    make_dem_with_isolated_artifact,
    make_dem_with_coastal_water_noise,
    read_stl_binary,
    assert_watertight,
    signed_volume,
)


def test_isolated_artifact_flattened_to_sea_level(tmp_path, tmp_stl):
    tif_path = tmp_path / "island_artifact.tif"
    cx, cy = make_dem_with_isolated_artifact(tif_path, artifact_dist=140.0)

    tris = build_mesh(
        str(tif_path), None, None, cx, cy, radius=145, diameter_mm=100, vexag=3.0,
        base_mm=3.0, print_spacing_mm=1.0,
        sea_level=0.0, min_elevation=-50.0, land_threshold=1.0,
    )
    write_stl_binary(tris, str(tmp_stl))
    read_back = read_stl_binary(tmp_stl)

    assert_watertight(read_back)
    assert signed_volume(read_back) > 0

    # La zone de l'artefact (100-140m réels) doit être aplatie à Z=0,
    # pas au niveau (20m) de l'artefact.
    top = read_back[read_back[:, :, 2] > -2.9]
    scale = (100 / 2.0) / 145
    r = np.sqrt(top[:, 0] ** 2 + top[:, 1] ** 2)
    band = (r > 105 * scale) & (r < 138 * scale)
    z_band = top[band, 2]
    assert z_band.size > 0
    assert np.allclose(z_band, 0.0), "l'artefact isolé n'a pas été aplati au niveau de la mer"


def test_coastal_water_noise_flattened(tmp_path, tmp_stl):
    tif_path = tmp_path / "island_water_noise.tif"
    cx, cy = make_dem_with_coastal_water_noise(tif_path, island_radius=100, noise_out_to=140.0)

    tris = build_mesh(
        str(tif_path), None, None, cx, cy, radius=145, diameter_mm=100, vexag=3.0,
        base_mm=3.0, print_spacing_mm=1.0,
        sea_level=0.0, min_elevation=-50.0, land_threshold=1.0,
    )
    write_stl_binary(tris, str(tmp_stl))
    read_back = read_stl_binary(tmp_stl)

    assert_watertight(read_back)
    assert signed_volume(read_back) > 0

    # Le bruit (100-140m réels), pourtant raccordé au littoral, doit être
    # aplati : sans --land-threshold, la connectivité seule le laisserait
    # passer puisqu'il touche la vraie terre.
    scale = (100 / 2.0) / 145
    top = read_back[read_back[:, :, 2] > -2.9]
    r = np.sqrt(top[:, 0] ** 2 + top[:, 1] ** 2)
    band = (r > 105 * scale) & (r < 138 * scale)
    z_band = top[band, 2]
    assert z_band.size > 0
    assert np.allclose(z_band, 0.0, atol=1e-4), "le bruit de surface d'eau n'a pas été aplati"

    # Le vrai sommet de l'île (loin du littoral) doit lui rester intact.
    assert read_back[:, :, 2].max() > 20.0 * scale * 3.0 * 0.8  # marge sous la valeur théorique exacte
