# SPDX-License-Identifier: GPL-3.0-or-later
"""Option --base-elevation : altitude imposée au dessus du socle.

Cas visé : couper un massif à une altitude, ou garder le même zéro sur plusieurs pièces. Tout ce qui est plus bas
est aplati à ce niveau ; tout ce qui est plus haut est conservé, y compris un sommet isolé du sommet principal
(le mode île `--sea-level`, lui, ne garde que la plus grande composante de terre et supprime ce sommet).

Les altitudes du STL sont confrontées au calcul : z = (altitude - référence) x échelle x exagération."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import vexag_curve  # noqa: E402
import circle_dem_to_stl  # noqa: E402
from circle_dem_to_stl import build_mesh  # noqa: E402
from conftest import assert_watertight, signed_volume, write_synthetic_tif  # noqa: E402
from test_report_cli import cone, run  # noqa: E402,F401

PIXEL = 10.0
CX, CY = 502000.0, 5002000.0           # centre du disque (m, EPSG:32632)
RADIUS, DIAMETER, BASE = 1900.0, 100.0, 5.0
SCALE = (DIAMETER / 2.0) / RADIUS      # mm par mètre (0,0263)
SMALL_DX, BIG_DX = 1000.0, -1000.0     # position des sommets par rapport au centre (m)


def make_two_summits(path):
    """Fond à 800 m, grand sommet de 1 800 m à l'ouest du centre, petit sommet de 1 300 m à l'est, séparés par un
    fond à 800 m : à 1 000 m, deux zones de terre sans lien entre elles."""
    n = 400
    x0, y0 = 500000.0, 5000000.0
    xs = x0 + (np.arange(n) + 0.5) * PIXEL
    ys = y0 + (np.arange(n)[::-1] + 0.5) * PIXEL
    X, Y = np.meshgrid(xs, ys)

    def bump(cx, cy, h, s):
        return h * np.exp(-(((X - cx) ** 2 + (Y - cy) ** 2) / (2 * s * s)))

    z = 800.0 + bump(CX + BIG_DX, CY, 1000.0, 450.0) + bump(CX + SMALL_DX, CY, 500.0, 300.0)
    write_synthetic_tif(path, z, x0, y0 + n * PIXEL, PIXEL, 32632)
    return z


@pytest.fixture
def two_summits(tmp_path):
    p = tmp_path / "deux_sommets.tif"
    z = make_two_summits(p)
    return str(p), z


def mesh(path, **kw):
    args = dict(radius=RADIUS, diameter_mm=DIAMETER, vexag=1.0, base_mm=BASE, print_spacing_mm=1.0)
    args.update(kw)
    return build_mesh(path, None, None, CX, CY, **args)


def summit_height(tris, dx_m):
    """Hauteur maximale des sommets de maillage situés à moins de 60 m (au sol) de l'axe du sommet."""
    v = tris.reshape(-1, 3)
    near = (np.abs(v[:, 0] - dx_m * SCALE) < 60 * SCALE) & (np.abs(v[:, 1]) < 60 * SCALE)
    return float(v[near, 2].max())


def test_both_summits_are_kept_when_cutting_at_an_altitude(two_summits):
    path, z = two_summits
    tris = mesh(path, base_elevation=1000.0)
    assert summit_height(tris, SMALL_DX) == pytest.approx((1300.0 - 1000.0) * SCALE, abs=0.1)
    assert summit_height(tris, BIG_DX) == pytest.approx((1800.0 - 1000.0) * SCALE, abs=0.1)
    assert_watertight(tris)
    assert signed_volume(tris) > 0


def test_island_mode_drops_the_secondary_summit_which_is_why_the_option_exists(two_summits):
    """Constat sur le mode île, qui n'est pas fait pour couper un massif : la composante secondaire est aplatie."""
    path, _ = two_summits
    tris = mesh(path, sea_level=1000.0, land_threshold=0.5)
    assert summit_height(tris, SMALL_DX) < 0.05                    # supprimé
    assert summit_height(tris, BIG_DX) == pytest.approx((1800.0 - 1000.0) * SCALE, abs=0.1)


def test_everything_below_the_reference_is_flat_at_the_top_of_the_base(two_summits):
    path, _ = two_summits
    tris = mesh(path, base_elevation=1000.0)
    v = tris.reshape(-1, 3)
    top_side = v[v[:, 2] > -1e-6, 2]                               # dessus et arête haute de la paroi
    assert top_side.min() == pytest.approx(0.0, abs=1e-5)
    assert np.mean(np.isclose(top_side, 0.0, atol=1e-5)) > 0.6     # le fond à 800 m est la majeure partie du disque


def test_reference_below_the_lowest_point_thickens_the_base(two_summits):
    path, _ = two_summits
    tris = mesh(path, base_elevation=700.0)
    v = tris.reshape(-1, 3)
    top_side = v[v[:, 2] > -1e-6, 2]
    assert top_side.min() == pytest.approx(100.0 * SCALE, abs=0.02)        # 800 m - 700 m = 100 m de plus de socle
    assert summit_height(tris, BIG_DX) == pytest.approx((1800.0 - 700.0) * SCALE, abs=0.1)
    assert v[:, 2].min() == pytest.approx(-BASE, abs=1e-5)                 # le dessous reste à -base_mm
    assert_watertight(tris)


def test_reference_equal_to_the_lowest_point_changes_nothing(two_summits, capsys):
    path, _ = two_summits
    default = mesh(path)
    capsys.readouterr()
    same = mesh(path, base_elevation=800.0)                         # le fond du terrain de test
    assert "sans effet" in capsys.readouterr().out
    assert default.shape == same.shape and np.allclose(default, same, atol=1e-4)           # 0,1 µm : le fond rééchantillonné vaut 800,0002 m


def test_the_reference_is_a_real_altitude_even_with_an_exaggeration_curve(two_summits):
    """--base-elevation est en mètres réels : la courbe d'exagération s'applique ensuite, au sommet comme au socle."""
    path, _ = two_summits
    curve = vexag_curve.parse("800:800,1000:1000,1800:1400")
    tris = mesh(path, base_elevation=1000.0, vexag_curve=curve)
    assert summit_height(tris, BIG_DX) == pytest.approx((1400.0 - 1000.0) * SCALE, abs=0.1)     # 1800 m -> 1400 m équivalents
    assert summit_height(tris, SMALL_DX) > 0.5                                                  # le petit sommet reste
    assert_watertight(tris)


def test_reference_above_the_highest_point_is_refused(two_summits):
    path, _ = two_summits
    with pytest.raises(SystemExit) as e:
        mesh(path, base_elevation=2000.0)
    assert "--base-elevation" in str(e.value.code) and "aucun relief" in str(e.value.code)


def test_it_cannot_be_combined_with_island_mode(cone, monkeypatch, capsys, tmp_path):
    with pytest.raises(SystemExit) as e:
        run(monkeypatch, capsys, tmp_path, "--base-elevation", "40")        # run() passe déjà --sea-level 0
    assert "s'excluent" in str(e.value.code) and list(tmp_path.iterdir()) == []


def test_cli_end_to_end_height_filename_and_report(cone, monkeypatch, capsys, tmp_path):
    argv = ["circle_dem_to_stl.py", "--input", "synthetique.tif", "--cx", "1000", "--cy", "2000", "--radius", "140",
            "--diameter-mm", "56", "--base-mm", "3", "--vexag", "3", "--print-spacing-mm", "0.5",
            "--base-elevation", "40", "--report", "--output-dir", str(tmp_path)]
    monkeypatch.setattr(sys, "argv", argv)
    circle_dem_to_stl.main()
    out = capsys.readouterr()
    log = out.out + out.err
    assert "--base-elevation 40 m" in log and "aplatis à ce niveau" in log
    stl = next(tmp_path.glob("*.stl"))
    assert "_ref40m" in stl.name
    tris = check_stl.read_stl(str(stl))
    # île de 80 m, échelle 0,2 mm/m, exagération 3 : (80 - 40) x 0,2 x 3 = 24 mm au-dessus du socle
    assert tris[:, :, 2].max() == pytest.approx(24.0, abs=1e-3)
    rep = json.loads(next(tmp_path.glob("*.json")).read_text(encoding="utf-8"))
    assert rep["model"]["base_elevation_m"] == 40.0
    assert rep["relief"]["relief_mm"] == pytest.approx(24.0, abs=1e-3)
