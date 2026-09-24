# SPDX-License-Identifier: GPL-3.0-or-later
"""
Tests de download.py limités à ce qui est vérifiable sans accès réseau :
construction des noms/URLs de tuiles, sélection des tuiles par cercle,
et réutilisation du cache. Le téléchargement HTTP réel n'est pas testé
ici (nécessite un accès réseau à amazonaws.com).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from download import tile_name, tile_url, tile_key_for_point, tiles_for_circle, download_copernicus  # noqa: E402


def test_tile_key_matches_known_real_sites():
    # Cas réels rencontrés au fil du projet, valeurs attendues vérifiées manuellement.
    cases = [
        (35.3606, 138.7274, (35, 138)),      # Mont Fuji
        (-3.1616, 35.5878, (-4, 35)),        # Ngorongoro
        (36.975454, -110.096102, (36, -111)),  # Monument Valley
    ]
    for lat, lon, expected in cases:
        assert tile_key_for_point(lat, lon) == expected


def test_tile_name_format():
    assert tile_name(35, 138) == "Copernicus_DSM_COG_10_N35_00_E138_00_DEM"
    assert tile_name(-4, 35) == "Copernicus_DSM_COG_10_S04_00_E035_00_DEM"
    assert tile_name(36, -111) == "Copernicus_DSM_COG_10_N36_00_W111_00_DEM"


def test_tile_url_uses_correct_bucket_per_product():
    assert "copernicus-dem-30m" in tile_url(35, 138, product_m=30)
    assert "copernicus-dem-90m" in tile_url(35, 138, product_m=90)


def test_tiles_for_circle_includes_center_tile():
    tiles = tiles_for_circle(35.3606, 138.7274, 9300)
    assert (35, 138) in tiles


def test_tiles_for_circle_scales_with_radius():
    small = tiles_for_circle(35.3606, 138.7274, 1000)
    large = tiles_for_circle(35.3606, 138.7274, 50000)
    assert len(large) >= len(small)


def test_cached_tile_is_reused_without_network(tmp_path):
    from conftest import write_synthetic_tif
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    name = tile_name(35, 138)
    write_synthetic_tif(cache_dir / f"{name}.tif", __import__("numpy").full((10, 10), 100.0),
                         138.0, 36.0, 0.1, 4326)
    paths = download_copernicus(35.3606, 138.7274, 9300, str(cache_dir))
    assert len(paths) == 1
    assert name in paths[0]
