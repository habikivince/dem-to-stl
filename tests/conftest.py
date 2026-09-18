"""
Utilitaires partagés par les tests : génération de DEM synthétiques,
lecture de STL binaire, vérifications d'étanchéité et de cohérence des
normales.

Ces scénarios reproduisent les cas réels rencontrés (et corrigés) au fil
du développement de l'outil : cône simple, artefact isolé en mer, bruit
de surface d'eau raccordé au littoral, motif en damier créant un
pincement diagonal.
"""
import struct
from collections import Counter

import numpy as np
import pytest
from osgeo import gdal, osr

gdal.UseExceptions()


def write_synthetic_tif(path, elev, origin_x, origin_y, pixel_size, epsg, nodata=-9999.0):
    size_y, size_x = elev.shape
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(str(path), size_x, size_y, 1, gdal.GDT_Float32)
    ds.SetGeoTransform([origin_x, pixel_size, 0, origin_y, 0, -pixel_size])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(elev.astype(np.float32))
    ds = None


def make_cone_dem(path, size=300, cx0=300000.0, cy0=7000000.0, epsg=32740,
                   island_radius=100, peak=80.0, slope=0.6):
    """Île conique simple, mer (nodata) au-delà de island_radius, avec un
    anneau d'artefacts d'interpolation très négatifs juste après le
    littoral (comme observé sur de vraies tuiles LiDAR)."""
    xs = np.arange(size) + cx0 - size / 2
    ys = np.arange(size) + cy0 - size / 2
    X, Y = np.meshgrid(xs, ys[::-1])
    dist = np.sqrt((X - cx0) ** 2 + (Y - cy0) ** 2)
    elev = np.where(dist <= island_radius, peak - dist * slope, -9999.0)
    rim = (dist > island_radius) & (dist <= island_radius + 15)
    elev[rim] = -3000.0 - dist[rim] * 5
    write_synthetic_tif(path, elev, cx0 - size / 2, cy0 + size / 2, 1.0, epsg)
    return cx0, cy0


def make_dem_with_isolated_artifact(path, size=300, cx0=300000.0, cy0=7000000.0,
                                     epsg=32740, island_radius=100, artifact_dist=140.0,
                                     artifact_width=2.0, artifact_elev=20.0):
    """Île conique + une balise/artefact isolé à artifact_dist du centre,
    à une altitude plausible mais sans lien topologique avec l'île."""
    xs = np.arange(size) + cx0 - size / 2
    ys = np.arange(size) + cy0 - size / 2
    X, Y = np.meshgrid(xs, ys[::-1])
    dist = np.sqrt((X - cx0) ** 2 + (Y - cy0) ** 2)
    elev = np.where(dist <= island_radius, 80 - dist * 0.6, -9999.0)
    rim = (dist > island_radius) & (dist <= island_radius + 15)
    elev[rim] = -3000.0 - dist[rim] * 5
    artifact_zone = (dist > artifact_dist - artifact_width) & (dist < artifact_dist + artifact_width)
    elev[artifact_zone] = artifact_elev
    write_synthetic_tif(path, elev, cx0 - size / 2, cy0 + size / 2, 1.0, epsg)
    return cx0, cy0


def make_dem_with_coastal_water_noise(path, size=300, cx0=300000.0, cy0=7000000.0,
                                       epsg=32740, island_radius=100, noise_out_to=140.0,
                                       seed=42):
    """Île à pente linéaire descendant jusqu'à 0 m pile au littoral, avec du
    bruit de surface d'eau (+-0.3m) directement raccordé à la côte au-delà —
    le cas qui a motivé le critère --land-threshold (la connectivité seule
    ne peut pas séparer ce bruit de la vraie terre)."""
    rng = np.random.default_rng(seed)
    xs = np.arange(size) + cx0 - size / 2
    ys = np.arange(size) + cy0 - size / 2
    X, Y = np.meshgrid(xs, ys[::-1])
    dist = np.sqrt((X - cx0) ** 2 + (Y - cy0) ** 2)
    elev = np.where(dist <= island_radius, 80 - dist * 0.8, np.nan)
    noise_zone = (dist > island_radius) & (dist <= noise_out_to)
    elev = np.where(noise_zone, rng.uniform(-0.3, 0.3, size=elev.shape), elev)
    elev = np.where(np.isnan(elev), -9999.0, elev)
    write_synthetic_tif(path, elev, cx0 - size / 2, cy0 + size / 2, 1.0, epsg)
    return cx0, cy0


def read_stl_binary(path):
    with open(path, "rb") as f:
        f.read(80)
        n = struct.unpack("<I", f.read(4))[0]
        tris = []
        for _ in range(n):
            struct.unpack("<3f", f.read(12))  # normale déclarée, ignorée : on la recalcule
            v = [struct.unpack("<3f", f.read(12)) for _ in range(3)]
            f.read(2)
            tris.append(v)
    return np.array(tris)


def edge_histogram(tris, decimals=4):
    """Compte, pour chaque arête (sommets arrondis pour tolérer le bruit
    flottant), combien de triangles la partagent. Un maillage fermé et
    correct a EXACTEMENT 2 pour chaque arête."""
    def key(p):
        return (round(p[0], decimals), round(p[1], decimals), round(p[2], decimals))

    counts = Counter()
    for tri in tris:
        pts = [key(p) for p in tri]
        for a, b in [(0, 1), (1, 2), (2, 0)]:
            edge = tuple(sorted([pts[a], pts[b]]))
            counts[edge] += 1
    return counts


def assert_watertight(tris):
    counts = edge_histogram(tris)
    bad = {e: c for e, c in counts.items() if c != 2}
    assert not bad, f"{len(bad)} arête(s) non partagée(s) par exactement 2 triangles : {list(bad.items())[:5]}"


def signed_volume(tris):
    """Volume signé par le théorème de la divergence. Positif si et
    seulement si les normales sont cohérentes et pointent vers
    l'extérieur sur un maillage fermé."""
    return sum(np.dot(v0, np.cross(v1, v2)) for v0, v1, v2 in tris) / 6.0


@pytest.fixture
def tmp_stl(tmp_path):
    return tmp_path / "out.stl"
