# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 1 de l'audit : écriture STL vectorisée, validation des paramètres, garde-fou de taille,
contrôle précoce de la sortie, erreurs d'entrée claires."""
import struct
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import circle_dem_to_stl  # noqa: E402
import sources  # noqa: E402
import stl_io  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402

quiet = lambda *a, **k: None  # noqa: E731


# ---------------------------------------------------------------- écriture STL

def reference_write(tris, path):
    """Ancienne implémentation (boucle triangle par triangle) conservée comme référence."""
    with open(path, "wb") as f:
        f.write(b"\x00" * 80)
        f.write(struct.pack("<I", len(tris)))
        for t in tris:
            v0, v1, v2 = np.asarray(t[0], float), np.asarray(t[1], float), np.asarray(t[2], float)
            n = np.cross(v1 - v0, v2 - v0)
            norm = np.linalg.norm(n)
            n = n / norm if norm > 0 else np.zeros(3)
            f.write(struct.pack("<3f", *n))
            for v in (v0, v1, v2):
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))


def random_mesh(n, seed=0):
    rng = np.random.default_rng(seed)
    tris = rng.uniform(-100, 100, size=(n, 3, 3))
    tris[3] = tris[3][0]                                   # triangle dégénéré (3 sommets identiques)
    tris[4, 2] = tris[4, 1]                                # triangle plat (2 sommets identiques)
    return tris


def test_vectorized_writer_matches_reference(tmp_path):
    tris = random_mesh(2500)
    reference_write(tris, str(tmp_path / "ref.stl"))
    stl_io.write_stl_binary(tris, str(tmp_path / "new.stl"))
    a, b = (tmp_path / "ref.stl").read_bytes(), (tmp_path / "new.stl").read_bytes()
    assert len(a) == len(b) == 84 + 50 * len(tris)
    ra = np.frombuffer(a[84:], dtype=stl_io.STL_DTYPE)
    rb = np.frombuffer(b[84:], dtype=stl_io.STL_DTYPE)
    assert np.array_equal(ra["v"], rb["v"])                 # sommets strictement identiques
    assert np.abs(ra["n"] - rb["n"]).max() < 1e-5           # normales identiques à l'arrondi float32 près
    assert np.all(rb["n"][3] == 0) and np.all(rb["n"][4] == 0)   # dégénérés : normale nulle
    assert a[:84] == b[:84]


def test_writer_chunking_progress_and_atomicity(tmp_path, monkeypatch):
    monkeypatch.setattr(stl_io, "CHUNK_TRIANGLES", 700)
    calls = []
    tris = random_mesh(2000)
    stl_io.write_stl_binary(tris, str(tmp_path / "c.stl"), progress=lambda i, n, label: calls.append((i, n)))
    assert calls[-1] == (1999, 2000) and len(calls) == 3
    assert check_stl.read_stl(str(tmp_path / "c.stl")).shape == (2000, 3, 3)
    assert not (tmp_path / "c.stl.part").exists()
    # échec en cours d'écriture : ni .part ni fichier tronqué ; un STL existant n'est pas écrasé
    (tmp_path / "keep.stl").write_bytes(b"ancien contenu valide")

    def boom(i, n, label):
        raise RuntimeError("interruption")
    with pytest.raises(RuntimeError):
        stl_io.write_stl_binary(tris, str(tmp_path / "keep.stl"), progress=boom)
    assert (tmp_path / "keep.stl").read_bytes() == b"ancien contenu valide"
    assert not (tmp_path / "keep.stl.part").exists()


def test_empty_mesh_is_valid_file(tmp_path):
    stl_io.write_stl_binary(np.zeros((0, 3, 3)), str(tmp_path / "e.stl"))
    assert (tmp_path / "e.stl").stat().st_size == 84


# ---------------------------------------------------------------- CLI : validations

@pytest.fixture
def dem(tmp_path):
    x, y = 2600000.0, 1201500.0
    yy, xx = np.mgrid[0:300, 0:300]
    write_synthetic_tif(tmp_path / "dem.tif", (1000 + 0.2 * xx + 0.1 * yy).astype(np.float32), x, y, 5.0, 2056)
    return str(tmp_path / "dem.tif")


def run(monkeypatch, dem, tmp_path, **over):
    a = {"input": dem, "cx": "2600750", "cy": "1200750", "radius": "600", "diameter-mm": "40", "vexag": "1",
         "base-mm": "2", "print-spacing-mm": "1", "output": str(tmp_path / "o.stl")}
    a.update({k.replace("_", "-"): v for k, v in over.items()})
    argv = ["x"]
    for k, v in a.items():
        if v is not None:
            argv += [f"--{k}", str(v)]
    monkeypatch.setattr(sys, "argv", argv)
    circle_dem_to_stl.main()


@pytest.mark.parametrize("over,msg", [
    ({"radius": "0"}, "--radius"), ({"radius": "-5"}, "--radius"), ({"radius": "nan"}, "--radius"),
    ({"diameter_mm": "0"}, "--diameter-mm"), ({"diameter_mm": "-40"}, "--diameter-mm"),
    ({"base_mm": "0"}, "--base-mm"), ({"base_mm": "-2"}, "--base-mm"),
    ({"vexag": "-1"}, "--vexag"), ({"vexag": "inf"}, "--vexag"),
    ({"print_spacing_mm": "0"}, "--print-spacing-mm"), ({"print_spacing_mm": "-1"}, "--print-spacing-mm"),
    ({"print_spacing_mm": "500"}, "trop grand"),
])
def test_invalid_parameters_are_clean_errors(monkeypatch, dem, tmp_path, over, msg):
    with pytest.raises(SystemExit) as e:
        run(monkeypatch, dem, tmp_path, **over)
    assert str(e.value).startswith("Erreur :") and msg in str(e.value)
    assert not (tmp_path / "o.stl").exists()


def test_latitude_longitude_and_place_pick_limits(monkeypatch, dem, tmp_path):
    for over, msg in (({"cx": None, "cy": None, "lat": "91", "lon": "7.6"}, "latitude"),
                      ({"cx": None, "cy": None, "lat": "46", "lon": "200"}, "longitude"),
                      ({"place_pick": "0"}, "--place-pick")):
        with pytest.raises(SystemExit) as e:
            run(monkeypatch, dem, tmp_path, **over)
        assert msg in str(e.value)


def test_vexag_zero_only_warns(monkeypatch, dem, tmp_path, capsys):
    run(monkeypatch, dem, tmp_path, vexag="0")
    assert "disque plat" in capsys.readouterr().out
    assert (tmp_path / "o.stl").exists()


# ---------------------------------------------------------------- sortie et entrée

def test_missing_output_folder_is_created_early(monkeypatch, dem, tmp_path):
    out = tmp_path / "nouveau dossier é" / "sous" / "x.stl"
    run(monkeypatch, dem, tmp_path, output=str(out))
    assert out.exists() and not Path(str(out) + ".part").exists()


def test_unwritable_output_fails_before_any_computation(monkeypatch, dem, tmp_path):
    (tmp_path / "fichier").write_text("x")
    called = []
    monkeypatch.setattr(circle_dem_to_stl, "build_mesh", lambda *a, **k: called.append(1))
    with pytest.raises(SystemExit) as e:
        run(monkeypatch, dem, tmp_path, output=str(tmp_path / "fichier" / "x.stl"))
    assert "impossible d'écrire" in str(e.value) and not called
    with pytest.raises(SystemExit) as e:
        run(monkeypatch, dem, tmp_path, output=str(tmp_path))
    assert "est un dossier" in str(e.value)


@pytest.mark.parametrize("kind", ["texte", "vide", "dossier", "absent"])
def test_unreadable_inputs_are_clean_errors(monkeypatch, dem, tmp_path, kind):
    bad = tmp_path / "bad.tif"
    if kind == "texte":
        bad.write_text("ceci n'est pas un tif")
    elif kind == "vide":
        bad.write_bytes(b"")
    elif kind == "dossier":
        bad = tmp_path / "un dossier"
        bad.mkdir()
    with pytest.raises(SystemExit) as e:
        run(monkeypatch, str(bad), tmp_path)
    assert str(e.value).startswith("Erreur : raster illisible")
    assert not (tmp_path / "o.stl").exists()


def test_input_path_with_spaces_and_accents_works(monkeypatch, tmp_path):
    d = tmp_path / "données cervin é"
    d.mkdir()
    yy, xx = np.mgrid[0:300, 0:300]
    write_synthetic_tif(d / "mnt é.tif", (1000 + 0.2 * xx).astype(np.float32), 2600000.0, 1201500.0, 5.0, 2056)
    run(monkeypatch, str(d / "mnt é.tif"), tmp_path, output=str(d / "sortie é.stl"))
    assert (d / "sortie é.stl").stat().st_size > 84


def test_source_failure_is_a_clean_error(monkeypatch, tmp_path):
    def down(*a, **k):
        raise sources.SourceUnavailable("serveur en panne")
    monkeypatch.setitem(sources.SOURCES, "swisstopo", sources.Source("swisstopo", "x", [(5, 45, 11, 48)], down))
    with pytest.raises(SystemExit) as e:
        run(monkeypatch, None, tmp_path, input=None, cx=None, cy=None, lat="46", lon="7.6", source="swisstopo")
    assert "source de données inutilisable" in str(e.value) and "serveur en panne" in str(e.value)


# ---------------------------------------------------------------- garde-fou de taille

def test_estimate_matches_measured_case():
    est = circle_dem_to_stl.estimate_resources(250, 0.2)             # cas mesuré : 1250², 4 908 732 triangles
    assert est["triangles"] == pytest.approx(4_908_732, rel=0.01)
    assert est["ram_bytes"] / 1e9 == pytest.approx(1.05 * 1.0, abs=0.45)   # mesuré : 1,05 Go de pic
    assert circle_dem_to_stl.estimate_resources(250, 0.2, polygon_mode=True)["triangles"] > est["triangles"]


def test_guard_blocks_oversized_grid_unless_forced(capsys):
    with pytest.raises(SystemExit) as e:
        circle_dem_to_stl.check_resources(250, 0.01, False, False, log=quiet, ram_available=lambda: 8_000_000_000)
    assert "--force" in str(e.value)
    circle_dem_to_stl.check_resources(250, 0.01, False, True, log=quiet, ram_available=lambda: 8_000_000_000)
    circle_dem_to_stl.check_resources(250, 0.2, False, False, log=quiet, ram_available=lambda: 8_000_000_000)
    circle_dem_to_stl.check_resources(250, 0.01, False, False, log=quiet, ram_available=lambda: None)  # RAM inconnue : pas de blocage


def test_available_ram_is_readable_on_linux():
    ram = circle_dem_to_stl.available_ram_bytes()
    assert ram is None or ram > 0
