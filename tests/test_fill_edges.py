# SPDX-License-Identifier: GPL-3.0-or-later
"""Comblement des cellules sans donnée : seuls les TROUS sont interpolés. Hors de l'emprise du fichier et autour
d'une grande zone vide (mer laissée sans donnée), gdal.FillNodata prolongerait le relief sur 50 pixels : ce relief
n'existe pas dans la source, il n'est plus créé (sauf --extrapolate-edges). Aucun réseau."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import circle_dem_to_stl  # noqa: E402
import raster  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402

X0, Y0 = 2600000.0, 1202000.0       # coin NO des fichiers de test (EPSG:2056)
NODATA = -9999.0


# --- raster.closed_gaps ------------------------------------------------------------------------------

def grid(h=200, w=200):
    return np.zeros((h, w), dtype=bool)


def test_closed_gaps_fills_a_hole_surrounded_by_data():
    m = grid()
    m[50:60, 50:60] = True
    assert raster.closed_gaps(m).sum() == 100


def test_closed_gaps_ignores_a_hole_wider_than_the_search_distance():
    m = grid(300, 300)
    m[50:250, 50:250] = True              # le centre est à 100 pixels d'une donnée
    assert raster.closed_gaps(m).sum() == 0


def test_closed_gaps_ignores_the_side_of_the_grid_beyond_the_data():
    side, corner = grid(), grid()
    side[:, 150:] = True                  # fichier arrêté à 150 pixels : un seul tronçon sur le pourtour
    corner[:40, 160:] = True
    assert raster.closed_gaps(side).sum() == 0
    assert raster.closed_gaps(corner).sum() == 0


def test_closed_gaps_fills_a_narrow_gap_that_crosses_the_grid():
    river, seam = grid(), grid()
    river[:, 90:110] = True               # rivière, joint entre deux dalles : traverse de part en part
    yy, xx = np.mgrid[0:200, 0:200]
    seam[np.abs((xx - yy) - 100) <= 4] = True      # diagonale, du haut au côté droit
    assert raster.closed_gaps(river).sum() == river.sum()
    assert raster.closed_gaps(seam).sum() == seam.sum()


def test_closed_gaps_degenerate_inputs():
    assert raster.closed_gaps(np.ones((20, 20), dtype=bool)).sum() == 0
    assert raster.closed_gaps(grid()).sum() == 0
    assert raster.closed_gaps(np.array([[True, False, True]])).sum() == 0


# --- raster.warp_window ------------------------------------------------------------------------------

def plane(h, w):
    yy, xx = np.mgrid[0:h, 0:w]
    return 1000.0 + 0.05 * xx + 0.03 * yy


def warp(path, cx, cy, radius, pixel=5.0, **kw):
    wkt, cx, cy = raster.resolve_target_crs_and_center(str(path), None, None, cx, cy)
    return raster.warp_window(str(path), wkt, cx, cy, radius, pixel, progress_cb=lambda *a: 1, **kw)


@pytest.fixture
def east_edge(tmp_path):
    """Fichier de 200 x 200 pixels (1 km) ; la fenêtre de 1 km déborde de 500 m à l'est, soit 100 pixels."""
    write_synthetic_tif(tmp_path / "e.tif", plane(200, 200), X0, Y0, 5.0, 2056)
    return tmp_path / "e.tif", X0 + 1000.0, Y0 - 500.0


def test_nothing_is_invented_beyond_the_file_footprint(east_edge):
    path, cx, cy = east_edge
    stats = {}
    elev, xs, ys, nodata = warp(path, cx, cy, 500.0, stats=stats)
    assert elev.shape == (200, 200)
    assert (elev[:, :100] != nodata).all() and (elev[:, 100:] == nodata).all()
    # ...mais l'ancien comblement aurait prolongé le relief : la zone est signalée dans stats
    edge = stats["edge_cells"]
    assert edge[:, 100:140].all() and not edge[:, :100].any() and not edge[:, 160:].any()


def test_extrapolate_edges_restores_the_old_fill_and_reports_it(east_edge):
    path, cx, cy = east_edge
    stats, strict_stats = {}, {}
    elev, _, _, nodata = warp(path, cx, cy, 500.0, extrapolate_edges=True, stats=stats)
    warp(path, cx, cy, 500.0, stats=strict_stats)
    assert (elev[:, 100:140] != nodata).all() and (elev[:, 160:] == nodata).all()
    assert np.array_equal(stats["edge_cells"], strict_stats["edge_cells"])      # même zone, remplie ou non
    assert stats["edge_cells"].sum() == (elev[:, 100:] != nodata).sum()
    # le relief prolongé reste proche de l'altitude du bord (interpolation pondérée par l'inverse de la distance)
    assert abs(elev[100, 110] - elev[100, 99]) < 5.0


def test_a_small_hole_inside_the_data_is_still_filled(tmp_path):
    z = plane(200, 200)
    z[90:96, 90:96] = NODATA
    write_synthetic_tif(tmp_path / "h.tif", z, X0, Y0, 5.0, 2056)
    stats = {}
    elev, _, _, nodata = warp(tmp_path / "h.tif", X0 + 500.0, Y0 - 500.0, 500.0, stats=stats)
    assert (elev != nodata).all()
    assert abs(elev[92, 92] - plane(200, 200)[92, 92]) < 1.0      # raccordé aux données voisines
    assert not stats["edge_cells"].any()


def test_a_narrow_gap_crossing_the_window_is_still_filled(tmp_path):
    z = plane(200, 200)
    z[:, 90:97] = NODATA                  # rivière sans donnée de 35 m, d'un bord de la fenêtre à l'autre
    write_synthetic_tif(tmp_path / "r.tif", z, X0, Y0, 5.0, 2056)
    stats = {}
    elev, _, _, nodata = warp(tmp_path / "r.tif", X0 + 500.0, Y0 - 500.0, 500.0, stats=stats)
    assert (elev != nodata).all()
    assert not stats["edge_cells"].any()


def test_a_large_empty_zone_inside_the_file_is_left_out(tmp_path):
    z = plane(400, 400)
    z[135:265, 135:265] = NODATA          # « mer » de 130 pixels de côté : son centre est à 65 pixels d'une donnée
    write_synthetic_tif(tmp_path / "s.tif", z, X0, Y0, 5.0, 2056)
    stats = {}
    elev, _, _, nodata = warp(tmp_path / "s.tif", X0 + 1000.0, Y0 - 1000.0, 1000.0, stats=stats)
    assert (elev[135:265, 135:265] == nodata).all()
    assert stats["edge_cells"].sum() > 0              # l'ancien comblement en aurait rempli l'anneau extérieur
    legacy, _, _, _ = warp(tmp_path / "s.tif", X0 + 1000.0, Y0 - 1000.0, 1000.0, extrapolate_edges=True)
    assert (legacy[135:265, 135:265] != nodata).sum() == stats["edge_cells"].sum()


def test_fill_false_is_untouched(east_edge):
    path, cx, cy = east_edge
    stats = {}
    elev, _, _, nodata = warp(path, cx, cy, 500.0, fill=False, extrapolate_edges=True, stats=stats)
    assert (elev[:, 100:] == nodata).all() and "edge_cells" not in stats


# --- ligne de commande -------------------------------------------------------------------------------

COMMON = ["--cx", "2601000", "--cy", "1201000", "--radius", "300", "--diameter-mm", "60", "--vexag", "1",
          "--base-mm", "2", "--print-spacing-mm", "0.5"]


def cli(monkeypatch, tmp_path, capsys, dem, *extra, name="m.stl"):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--input", str(dem), *COMMON, *extra, "--output", name])
    circle_dem_to_stl.main()
    out = capsys.readouterr()
    return tmp_path / "output" / name, out.out + out.err


@pytest.fixture
def cliff_coast(tmp_path):
    """Falaise : terre à 30 m à l'ouest du centre du disque, mer (sans donnée) à l'est. Le disque coupe la côte."""
    yy, xx = np.mgrid[0:400, 0:400]
    z = 30.0 + 0.02 * xx + 0.0 * yy
    z[:, 200:] = NODATA
    write_synthetic_tif(tmp_path / "coast.tif", z, 2600000.0, Y0, 5.0, 2056)
    return tmp_path / "coast.tif"


def test_cli_sea_level_coast_has_no_extrapolated_plateau(monkeypatch, tmp_path, capsys, cliff_coast):
    stl, out = cli(monkeypatch, tmp_path, capsys, cliff_coast, "--sea-level", "0", "--report")
    stl2, out2 = cli(monkeypatch, tmp_path, capsys, cliff_coast, "--sea-level", "0", "--extrapolate-edges",
                     "--report", name="legacy.stl")
    strict, legacy = check_stl.read_stl(str(stl)), check_stl.read_stl(str(stl2))
    # sans l'option, la mer reste à plat jusqu'à la côte ; avec, un plateau de 30 m s'avance dans l'eau
    assert check_stl.signed_volume(legacy) > 1.1 * check_stl.signed_volume(strict)
    assert "laissées sans donnée plutôt que prolongées" in out and "EXTRAPOLÉE" not in out
    assert "EXTRAPOLÉE(S) au-delà des données" in out2 and "ATTENTION" in out2
    rep = json.loads((tmp_path / "output" / "m.json").read_text(encoding="utf-8"))
    rep2 = json.loads((tmp_path / "output" / "legacy.json").read_text(encoding="utf-8"))
    nd, nd2 = rep["source"]["no_data"], rep2["source"]["no_data"]
    assert nd["extrapolate_edges"] is False and nd["cells_extrapolated"] == 0 and nd["cells_left_empty_at_edge"] > 0
    assert nd2["extrapolate_edges"] is True and nd2["cells_left_empty_at_edge"] == 0
    assert nd2["cells_extrapolated"] == nd["cells_left_empty_at_edge"]
    assert rep["options"]["extrapolate_edges"] is False and rep2["options"]["extrapolate_edges"] is True


def test_cli_without_edge_gap_logs_nothing_about_extrapolation(monkeypatch, tmp_path, capsys):
    yy, xx = np.mgrid[0:400, 0:400]
    write_synthetic_tif(tmp_path / "full.tif", 1000.0 + 0.08 * xx + 6.0 * np.sin(yy / 17.0), 2600000.0, Y0, 5.0, 2056)
    _, out = cli(monkeypatch, tmp_path, capsys, tmp_path / "full.tif", "--report")
    assert "EXTRAPOL" not in out and "plutôt que prolongées" not in out
    rep = json.loads((tmp_path / "output" / "m.json").read_text(encoding="utf-8"))
    assert rep["source"]["no_data"]["cells_without_data"] == 0
    assert rep["source"]["no_data"]["cells_left_empty_at_edge"] == 0


@pytest.fixture
def short_dem(tmp_path):
    """Données jusqu'à x = 2601150 seulement : la forme (rayon 300 m) déborde de 150 m, soit 30 pixels."""
    yy, xx = np.mgrid[0:400, 0:400]
    z = 1000.0 + 0.08 * xx + 6.0 * np.sin(yy / 17.0)
    z[:, 230:] = NODATA
    write_synthetic_tif(tmp_path / "short.tif", z, 2600000.0, Y0, 5.0, 2056)
    return tmp_path / "short.tif"


def test_regular_polygon_refused_by_default_and_message_names_the_option(monkeypatch, tmp_path, capsys, short_dem):
    with pytest.raises(SystemExit) as e:
        cli(monkeypatch, tmp_path, capsys, short_dem, "--sides", "6", "--circle-mode", "inscribed")
    msg = str(e.value)
    assert "n'ont pas de donnée valide" in msg and "--sea-level" in msg and "mm imprimés de son bord" in msg
    assert "--extrapolate-edges" in msg


def test_regular_polygon_accepted_with_extrapolate_edges(monkeypatch, tmp_path, capsys, short_dem):
    stl, out = cli(monkeypatch, tmp_path, capsys, short_dem, "--sides", "6", "--circle-mode", "inscribed",
                   "--extrapolate-edges")
    tris = check_stl.read_stl(str(stl))
    assert check_stl.signed_volume(tris) > 0
    assert "EXTRAPOLÉE(S) au-delà des données" in out
