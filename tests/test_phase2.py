# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 2 : dossier de sortie, noms automatiques, formats 3MF et OBJ, cohérence des trois formats,
et non-régression de l'option de téléchargement automatique."""
import os
import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import circle_dem_to_stl  # noqa: E402
import download  # noqa: E402
import mesh_export  # noqa: E402
import sources  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402

gdal.UseExceptions()


@pytest.fixture
def dem(tmp_path):
    yy, xx = np.mgrid[0:300, 0:300]
    write_synthetic_tif(tmp_path / "dem.tif", (1000 + 0.2 * xx + 5 * np.sin(yy / 20.0)).astype(np.float32),
                        2600000.0, 1201500.0, 5.0, 2056)
    return str(tmp_path / "dem.tif")


def run(monkeypatch, tmp_path, dem, *extra, base=True):
    argv = ["x", "--input", dem, "--cx", "2600750", "--cy", "1200750", "--radius", "600", "--diameter-mm", "40",
            "--vexag", "1", "--base-mm", "2", "--print-spacing-mm", "1"] + list(extra)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", argv)
    circle_dem_to_stl.main()


# ---------------------------------------------------------------- dossier de sortie et noms

def test_bare_filename_goes_to_output_folder_and_path_is_respected(monkeypatch, tmp_path, dem):
    run(monkeypatch, tmp_path, dem, "--output", "a.stl")
    assert (tmp_path / "output" / "a.stl").exists() and not (tmp_path / "a.stl").exists()
    run(monkeypatch, tmp_path, dem, "--output", str(tmp_path / "ailleurs" / "b.stl"))
    assert (tmp_path / "ailleurs" / "b.stl").exists()
    run(monkeypatch, tmp_path, dem, "--output", "c.stl", "--output-dir", "mes sorties é")
    assert (tmp_path / "mes sorties é" / "c.stl").exists()


def test_legacy_name_without_known_extension_is_kept_verbatim(monkeypatch, tmp_path, dem):
    run(monkeypatch, tmp_path, dem, "--output", "relief.bin")
    f = tmp_path / "output" / "relief.bin"
    assert f.exists() and check_stl.read_stl(str(f)).shape[1:] == (3, 3)
    assert not (tmp_path / "output" / "relief.bin.stl").exists()


def test_automatic_name_never_overwrites(monkeypatch, tmp_path, dem):
    run(monkeypatch, tmp_path, dem)
    run(monkeypatch, tmp_path, dem)
    names = sorted(p.name for p in (tmp_path / "output").iterdir())
    assert names == ["x2600750_y1200750_40mm.stl", "x2600750_y1200750_40mm_2.stl"]


def test_explicit_output_still_overwrites(monkeypatch, tmp_path, dem):
    run(monkeypatch, tmp_path, dem, "--output", "same.stl")
    first = (tmp_path / "output" / "same.stl").stat().st_mtime_ns
    run(monkeypatch, tmp_path, dem, "--output", "same.stl")
    assert (tmp_path / "output" / "same.stl").stat().st_mtime_ns >= first
    assert [p.name for p in (tmp_path / "output").iterdir()] == ["same.stl"]


def test_slug_and_default_stem():
    class A:
        place, clip_to_place, polygon, lat, lon, cx, cy = "Mont Fuji, Japon", False, None, None, None, 1.0, 2.0
        diameter_mm, vexag = 250.0, 1.0
    assert circle_dem_to_stl.slugify("Île de la Réunion !") == "ile-de-la-reunion"
    assert circle_dem_to_stl.default_stem(A) == "mont-fuji-japon_250mm"
    A.vexag, A.clip_to_place = 2.0, True
    assert circle_dem_to_stl.default_stem(A) == "mont-fuji-japon-contour_250mm_x2"
    A.place, A.lat, A.lon = None, 45.9763, -7.6586
    assert circle_dem_to_stl.default_stem(A) == "lat45.976_lonm7.659_250mm_x2"


# ---------------------------------------------------------------- formats

def test_format_selection_rules(monkeypatch, tmp_path, dem):
    run(monkeypatch, tmp_path, dem, "--output", "t.3mf")                      # l'extension choisit le format
    assert sorted(os.listdir(tmp_path / "output")) == ["t.3mf"]
    run(monkeypatch, tmp_path, dem, "--output", "u", "--format", "stl,obj")
    assert {"u.stl", "u.obj"} <= set(os.listdir(tmp_path / "output"))
    run(monkeypatch, tmp_path, dem, "--output", "v.stl", "--format", "all")    # --format l'emporte, même nom de base
    assert {"v.stl", "v.3mf", "v.obj"} <= set(os.listdir(tmp_path / "output"))
    with pytest.raises(SystemExit) as e:
        run(monkeypatch, tmp_path, dem, "--output", "w", "--format", "ply")
    assert "format inconnu" in str(e.value)


@pytest.fixture
def three_formats(monkeypatch, tmp_path, dem):
    run(monkeypatch, tmp_path, dem, "--output", "m", "--format", "all", "--smooth-wall")
    d = tmp_path / "output"
    return d / "m.stl", d / "m.3mf", d / "m.obj"


def canonical(tris):
    """Multi-ensemble de triangles (sommets arrondis float32), indépendant de l'ordre des triangles."""
    t = np.asarray(tris, dtype=np.float32).reshape(len(tris), 9)
    return t[np.lexsort(t.T[::-1])]


def test_three_formats_have_identical_geometry_and_units(three_formats):
    stl, mf, obj = three_formats
    ts, tm, to = (mesh_export.read_triangles(str(p)) for p in (stl, mf, obj))
    assert len(ts) == len(tm) == len(to) > 1000
    # mêmes sommets, mêmes triangles, même orientation : tableaux identiques triangle par triangle
    assert np.array_equal(ts.astype(np.float32), tm.astype(np.float32))
    assert np.array_equal(ts.astype(np.float32), to.astype(np.float32))
    # mêmes dimensions (mm) et même volume signé positif (normales sortantes)
    for t in (tm, to):
        assert np.allclose(t.reshape(-1, 3).min(axis=0), ts.reshape(-1, 3).min(axis=0))
        assert np.allclose(t.reshape(-1, 3).max(axis=0), ts.reshape(-1, 3).max(axis=0))
    v, f = mesh_export.index_mesh(ts)
    assert mesh_export.signed_volume(v, f) > 0
    assert 38 < np.ptp(ts.reshape(-1, 3)[:, 0]) <= 40.01      # --diameter-mm 40 : unité = millimètre


def test_3mf_package_structure_unit_and_metadata(three_formats):
    _, mf, _ = three_formats
    with zipfile.ZipFile(mf) as z:
        names = z.namelist()
        assert names[0] == "[Content_Types].xml" and "_rels/.rels" in names and "3D/3dmodel.model" in names
        assert z.testzip() is None
    verts, faces, unit, meta = mesh_export.read_3mf(str(mf))
    assert unit == "millimeter"
    assert meta["Title"] == "m" and meta["Application"] == "dem-to-stl" and "40 mm" in meta["Description"]
    assert meta["CreationDate"].endswith("Z") and set(meta) <= {"Title", "Description", "Copyright", "CreationDate",
                                                                 "Application"}
    assert faces.max() == len(verts) - 1 and faces.min() == 0           # indices valides, sommets tous utilisés
    assert len(np.unique(faces)) == len(verts)


def test_3mf_validated_by_official_lib3mf(three_formats):
    lib3mf = pytest.importorskip("lib3mf")
    _, mf, _ = three_formats
    w = lib3mf.get_wrapper()
    model = w.CreateModel()
    model.QueryReader("3mf").ReadFromFile(str(mf))
    assert model.GetUnit() == lib3mf.ModelUnit.MilliMeter
    it = model.GetMeshObjects()
    assert it.MoveNext()
    obj = it.GetCurrentMeshObject()
    verts, faces, _, _ = mesh_export.read_3mf(str(mf))
    assert obj.GetVertexCount() == len(verts) and obj.GetTriangleCount() == len(faces)
    assert obj.IsManifoldAndOriented()


def test_obj_has_only_what_is_relevant(three_formats):
    _, _, obj = three_formats
    text = obj.read_text(encoding="utf-8")
    assert "unité : millimètre" in text
    kinds = {line.split()[0] for line in text.splitlines() if line and not line.startswith("#")}
    assert kinds == {"o", "v", "f"}                                 # ni normales, ni textures, ni matériaux
    assert not obj.with_suffix(".mtl").exists()


def test_check_stl_accepts_all_three_formats(three_formats):
    for p in three_formats:
        tris = mesh_export.read_triangles(str(p))
        res = {n: ok for n, ok, _ in check_stl.check(tris, 40.0, 2.0, tol_mm=1.0)}
        assert all(res.values()), (p.name, res)


def test_index_mesh_merges_vertices_and_keeps_orientation():
    tris = np.array([[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[1, 0, 0], [1, 1, 0], [0, 1, 0]]], dtype=np.float64)
    v, f = mesh_export.index_mesh(tris)
    assert len(v) == 4 and f.shape == (2, 3)
    assert np.array_equal(v[f].astype(np.float64), tris)             # mêmes triangles, même ordre, même orientation
    v0, f0 = mesh_export.index_mesh(np.zeros((0, 3, 3)))
    assert len(v0) == 0 and f0.shape == (0, 3)


def test_attribution_is_written_for_the_data_source(monkeypatch, tmp_path, dem):
    ch = np.full((300, 300), 1500.0, dtype=np.float32)
    x, y = sources._to_epsg(2056, 46.0, 7.6)
    write_synthetic_tif(tmp_path / "ch.tif", ch, x - 750.0, y + 750.0, 5.0, 2056)
    monkeypatch.setitem(sources.SOURCES, "swisstopo", sources.Source("swisstopo", "x", [(5, 45, 11, 48)],
                                                                    lambda *a, **k: [str(tmp_path / "ch.tif")]))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--source", "swisstopo", "--lat", "46", "--lon", "7.6", "--radius", "600",
                                      "--diameter-mm", "40", "--vexag", "1", "--base-mm", "2",
                                      "--print-spacing-mm", "1", "--format", "3mf,obj", "--cache-dir", "c"])
    circle_dem_to_stl.main()
    f = next((tmp_path / "output").glob("*.3mf"))
    assert f.name == "lat46.000_lon7.600_40mm.3mf"
    assert mesh_export.read_3mf(str(f))[3]["Copyright"] == "© swisstopo"
    assert "© swisstopo" in next((tmp_path / "output").glob("*.obj")).read_text(encoding="utf-8")


# ---------------------------------------------------------------- l'option de téléchargement est conservée

def test_download_copernicus_option_still_works(monkeypatch, tmp_path):
    lat, lon = 46.0, 7.6
    n = 400
    lons = lon - 0.05
    yy, xx = np.mgrid[0:n, 0:n]
    write_synthetic_tif(tmp_path / "cop.tif", (1000 + 0.5 * xx + 0.3 * yy).astype(np.float32), lons, lat + 0.05,
                        0.1 / n, 4326)
    calls = []
    monkeypatch.setattr(download, "download_copernicus",
                        lambda *a, **k: calls.append(k.get("product_m")) or [str(tmp_path / "cop.tif")])
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--download-copernicus", "--lat", str(lat), "--lon", str(lon),
                                      "--radius", "1500", "--diameter-mm", "40", "--vexag", "1", "--base-mm", "2",
                                      "--print-spacing-mm", "1", "--output", "cop.stl", "--cache-dir", "cache"])
    circle_dem_to_stl.main()
    assert calls == [30] and (tmp_path / "output" / "cop.stl").stat().st_size > 84
    assert check_stl.read_stl(str(tmp_path / "output" / "cop.stl")).shape[0] > 100


def test_source_copernicus_and_auto_still_download(monkeypatch, tmp_path):
    lat, lon = 0.0, -30.0                       # plein Atlantique : aucune source nationale -> Copernicus
    n = 400
    yy, xx = np.mgrid[0:n, 0:n]
    write_synthetic_tif(tmp_path / "cop.tif", (200 + 0.5 * xx).astype(np.float32), lon - 0.05, lat + 0.05, 0.1 / n, 4326)
    calls = []
    monkeypatch.setattr(download, "download_copernicus", lambda *a, **k: calls.append(1) or [str(tmp_path / "cop.tif")])
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--source", "auto", "--lat", str(lat), "--lon", str(lon), "--radius", "1500",
                                      "--diameter-mm", "40", "--vexag", "1", "--base-mm", "2", "--print-spacing-mm", "1",
                                      "--cache-dir", "cache"])
    circle_dem_to_stl.main()
    assert calls == [1] and len(list((tmp_path / "output").glob("*.stl"))) == 1


# ---------------------------------------------------------------- orientation (--smooth-wall) et nombres

@pytest.mark.parametrize("extra", [[], ["--smooth-wall"]])
@pytest.mark.parametrize("spacing", ["1", "0.5"])
def test_every_mesh_variant_is_watertight_and_consistently_oriented(monkeypatch, tmp_path, dem, extra, spacing):
    """Régression : --smooth-wall produisait des triangles retournés (arêtes mal orientées, volume
    négatif) selon le sens de parcours du contour ; trouvé par le contrôle indépendant de check_stl.py."""
    run(monkeypatch, tmp_path, dem, "--output", "w.stl", "--print-spacing-mm", spacing, *extra)
    tris = check_stl.read_stl(str(tmp_path / "output" / "w.stl"))
    rep = check_stl.edge_report(tris)
    assert rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0 and rep["duplicate_directed_edges"] == 0
    v, f = mesh_export.index_mesh(tris)
    assert mesh_export.signed_volume(v, f) > 0


def test_numbers_roundtrip_exactly_including_tiny_values(tmp_path):
    v = np.array([[0.0, -0.0, 1e-9], [123.456789, -2.9999998, 0.91053063], [1.5e-5, 250.0, -3.3333333]], dtype=np.float32)
    f = np.array([[0, 1, 2]], dtype=np.int64)
    mesh_export.write_3mf(str(tmp_path / "t.3mf"), v, f)
    mesh_export.write_obj(str(tmp_path / "t.obj"), v, f)
    for name, reader in (("t.3mf", lambda p: mesh_export.read_3mf(p)[0]), ("t.obj", lambda p: mesh_export.read_obj(p)[0])):
        back = reader(str(tmp_path / name))
        assert np.array_equal(back, v), name
    assert "e-" not in (tmp_path / "t.obj").read_text(encoding="utf-8")
    assert b"e-" not in zipfile.ZipFile(tmp_path / "t.3mf").read("3D/3dmodel.model")
