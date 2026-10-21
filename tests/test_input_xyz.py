# SPDX-License-Identifier: GPL-3.0-or-later
"""XYZ texte : une liste « x y z » qui forme une GRILLE RÉGULIÈRE est lue par le pilote GDAL « ASCII Gridded XYZ »
comme n'importe quel raster (c'est le format de livraison de plusieurs MNT 1 m, p. ex. DGM1 de Saxe). Une liste de
points dispersés n'est pas une grille : elle est refusée.

Comme dans test_input_formats.py, le terrain est un plan incliné (voir fmt_helpers) : on vérifie des VALEURS et
l'orientation, pas seulement que le fichier « s'ouvre »."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import raster  # noqa: E402
from fmt_helpers import X0, Y0, centres, plane, prepared, window  # noqa: E402

GT = (X0, 10.0, 0.0, Y0 + 1000.0, 0.0, -10.0)           # 100 x 100 pixels de 10 m, LV95 (EPSG:2056)
CENTRE = dict(cx=X0 + 500.0, cy=Y0 + 500.0, radius=300.0, pixel=10.0)   # grille de travail alignée sur les pixels


def write_xyz(path, order="south_first", skip_row=None, sep=" "):
    """Écrit le plan de test au centre des pixels. order : « south_first » (comme la plupart des livraisons de MNT
    en XYZ) ou « north_first ». skip_row : indice (compté du nord) d'une ligne entière à omettre."""
    xs, ys = centres(GT, (100, 100))                       # ys : du nord au sud
    rows = list(range(100))
    if order == "south_first":
        rows.reverse()
    lines = []
    for i in rows:
        if i == skip_row:
            continue
        lines += [f"{x:.2f}{sep}{ys[i]:.2f}{sep}{plane(x, ys[i]):.3f}" for x in xs]
    Path(path).write_text("\n".join(lines) + "\n")
    return str(path)


def truth(xs, ys):
    return plane(xs[None, :], ys[:, None])


def run(names, **kw):
    args = dict(CENTRE)
    args.update(kw)
    return window(names, **args)


@pytest.mark.parametrize("order", ["south_first", "north_first"])
def test_gridded_xyz_is_read_like_a_raster_with_correct_orientation(tmp_path, order):
    path = write_xyz(tmp_path / f"{order}.xyz", order)
    elev, xs, ys, nodata, _ = run(prepared(path, source_crs="EPSG:2056")[1])
    assert elev.shape == (60, 60)
    assert np.abs(elev - truth(xs, ys)).max() < 2e-3
    assert elev[0].mean() > elev[-1].mean()                           # le plan monte vers le nord (haut du tableau)
    assert elev[:, -1].mean() > elev[:, 0].mean()                     # et vers l'est


def test_gridded_xyz_accepts_comma_separated_columns(tmp_path):
    path = write_xyz(tmp_path / "virgules.xyz", sep=",")
    try:
        names = prepared(path, source_crs="EPSG:2056")[1]
    except raster.InputError:
        pytest.skip("cette version de GDAL ne lit pas les XYZ séparés par des virgules")
    elev, xs, ys, _, _ = run(names)
    assert np.abs(elev - truth(xs, ys)).max() < 2e-3


def test_gridded_xyz_without_crs_is_refused_with_the_way_out(tmp_path):
    path = write_xyz(tmp_path / "nocrs.xyz")
    with pytest.raises(raster.InputError, match="--source-crs"):
        prepared(path)


def test_missing_row_in_a_gridded_xyz_is_nodata_never_zero(tmp_path):
    path = write_xyz(tmp_path / "trou.xyz", skip_row=50)
    elev, xs, ys, nodata, _ = run(prepared(path, source_crs="EPSG:2056")[1], fill=False)
    bad = ~np.isfinite(elev) | (elev == nodata)
    assert bad.any(), "la ligne omise doit rester sans donnée quand le comblement est désactivé"
    assert bad.sum(axis=1).max() == elev.shape[1] and bad.sum(axis=1).min() == 0    # une rangée entière, pas davantage
    ok = ~bad
    assert np.abs(elev[ok] - truth(xs, ys)[ok]).max() < 2e-3          # les autres cellules sont intactes
    assert not np.any(elev[ok] == 0.0)


def test_isolated_missing_point_is_nodata_and_a_real_zero_stays_a_value(tmp_path):
    """Une liste « x y z » en grille avec UN point absent et UN point d'altitude réelle 0,000 (côte) : le point absent
    est sans donnée (comblé ensuite, comme tout trou), le zéro réel reste une valeur. GDAL choisit lui-même une
    valeur NoData absente des données (0 s'il n'y a aucun zéro, sinon -32768) ; ce test le garantit pour dem-to-stl."""
    xs, ys = centres(GT, (100, 100))
    lines, absent, zero = [], (40, 40), (60, 60)           # indices (ligne depuis le nord, colonne)
    for i in range(99, -1, -1):
        for j, x in enumerate(xs):
            if (i, j) == absent:
                continue
            z = 0.0 if (i, j) == zero else float(plane(x, ys[i]))
            lines.append(f"{x:.2f} {ys[i]:.2f} {z:.3f}")
    path = tmp_path / "zero_et_absent.xyz"
    path.write_text("\n".join(lines) + "\n")

    def value_at(elev, wx, wy, i, j):
        r = int(np.argmin(np.abs(wy - ys[i])))
        c = int(np.argmin(np.abs(wx - xs[j])))
        return elev[r, c]

    elev, wx, wy, nodata, _ = run(prepared(str(path), source_crs="EPSG:2056")[1], fill=False)
    assert value_at(elev, wx, wy, *absent) == nodata
    assert value_at(elev, wx, wy, *zero) == 0.0                       # le zéro réel n'est pas pris pour un NoData
    assert value_at(elev, wx, wy, 60, 61) == pytest.approx(float(plane(xs[61], ys[60])), abs=2e-3)
    filled, wx, wy, nodata, _ = run(prepared(str(path), source_crs="EPSG:2056")[1], fill=True)
    assert value_at(filled, wx, wy, *absent) == pytest.approx(float(plane(xs[40], ys[40])), abs=0.5)   # trou comblé
    assert value_at(filled, wx, wy, *zero) == 0.0


def test_scattered_xyz_points_are_refused(tmp_path):
    rng = np.random.default_rng(3)
    x = X0 + rng.uniform(0, 1000, 3000)
    y = Y0 + rng.uniform(0, 1000, 3000)
    path = tmp_path / "disperse.xyz"
    path.write_text("\n".join(f"{a:.2f} {b:.2f} {plane(a, b):.3f}" for a, b in zip(x, y)) + "\n")
    with pytest.raises(raster.InputError, match="Ungridded|illisible|grille"):
        names = prepared(str(path), source_crs="EPSG:2056")[1]
        run(names)


def test_survey_style_xyz_with_header_and_extra_columns_is_refused(tmp_path):
    """Mise en page des sondages bathymétriques NOAA (« survey_id,lat,long,depth,quality_code,active ») : une ligne
    d'en-tête, une colonne texte et des colonnes en plus, des points dispersés. Ce n'est pas une grille XYZ : refusé."""
    rng = np.random.default_rng(5)
    lat = 47.60 + rng.uniform(0, 0.05, 500)
    lon = -122.28 + rng.uniform(0, 0.06, 500)
    lines = ["survey_id,lat,long,depth,quality_code,active"]
    lines += [f"H11292,{a:.8f},{b:.8f},{rng.uniform(1, 60):.2f},1,1" for a, b in zip(lat, lon)]
    path = tmp_path / "H11292.xyz"
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(raster.InputError, match="illisible|not recognized|Ungridded"):
        names = prepared(str(path), source_crs="EPSG:4326")[1]
        run(names)
