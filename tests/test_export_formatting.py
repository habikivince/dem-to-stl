# SPDX-License-Identifier: GPL-3.0-or-later
"""Écriture des sommets du 3MF/OBJ : la voie sans exposant formate chaque valeur distincte une seule fois,
mais le texte produit doit rester identique, caractère pour caractère, à l'ancien formatage valeur par valeur."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import mesh_export  # noqa: E402

VERTEX_FMT = '<vertex x="%.9g" y="%.9g" z="%.9g"/>'
OBJ_FMT = "v %.9g %.9g %.9g\n"


def reference_vertex_lines(vertices, fmt, chunk=200_000):
    """Ancien code de mesh_export._vertex_lines (formatage valeur par valeur), repris tel quel."""
    for s in range(0, len(vertices), chunk):
        block = vertices[s:s + chunk].astype(np.float64).tolist()
        text = "".join(fmt % (x, y, z) for x, y, z in block)
        if "e" in text:
            text = "".join(fmt.replace("%.9g", "%s") % tuple(mesh_export._positional(c) for c in row) for row in block)
        yield text


def _arrays():
    rng = np.random.default_rng(0)
    grid = np.stack(np.meshgrid(np.arange(40) * 0.1 - 2.0, np.arange(40) * 0.1 - 2.0, [0.0])).reshape(3, -1).T
    special = np.array([[0.0, -0.0, 1e-7], [1e-30, -1e-30, 123456.789], [3.4e38, -3.4e38, 0.1],
                        [0.1, 0.2, 0.3], [-0.0, 0.0, 5.0], [1e-5, 99999.99, -2.5e-9]] * 20)
    return {"aléatoire": rng.normal(0, 50, (4000, 3)), "grille": grid, "valeurs particulières": special,
            "un seul sommet": np.array([[1.5, -2.25, 0.0]]), "aucun": np.empty((0, 3))}


@pytest.mark.parametrize("fmt", [VERTEX_FMT, OBJ_FMT])
@pytest.mark.parametrize("name", list(_arrays()))
def test_vertex_text_is_identical_to_the_value_by_value_formatting(name, fmt):
    v = _arrays()[name]
    assert "".join(mesh_export._vertex_lines(v, fmt)) == "".join(reference_vertex_lines(v, fmt))


def test_positional_lines_keep_negative_zero_distinct_from_zero():
    text = mesh_export._positional_lines(np.array([[0.0, -0.0, 1.0], [-0.0, 0.0, 2.0]]), VERTEX_FMT)
    assert 'x="0" y="-0" z="1"' in text and 'x="-0" y="0" z="2"' in text


def test_positional_lines_never_use_exponent_notation():
    v = np.array([[1e-7, -3e-9, 2.5e-30], [1e20, 1e-5, 0.5]])
    text = mesh_export._positional_lines(v, VERTEX_FMT)
    assert "e" not in text.replace("vertex", "")


def test_chunk_boundary_does_not_change_the_text(monkeypatch):
    v = _arrays()["aléatoire"]
    whole = "".join(mesh_export._vertex_lines(v, VERTEX_FMT))
    monkeypatch.setattr(mesh_export, "CHUNK", 777)
    assert "".join(mesh_export._vertex_lines(v, VERTEX_FMT)) == whole
