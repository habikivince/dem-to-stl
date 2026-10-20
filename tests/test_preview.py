# SPDX-License-Identifier: GPL-3.0-or-later
"""preview.py : aperçu ombré PNG (--preview). Sans GDAL : numpy, scipy et la bibliothèque standard.

Le PNG est relu par un décodeur indépendant (signature, CRC de chaque bloc, zlib), pas par la fonction qui l'écrit."""
import struct
import sys
import zlib
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import preview  # noqa: E402


def decode_png(data):
    """Décodeur minimal du seul format écrit ici : 8 bits, gris + opacité, filtre 0. Vérifie les CRC."""
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    pos, chunks = 8, {}
    while pos < len(data):
        (n,) = struct.unpack(">I", data[pos:pos + 4])
        kind, body = data[pos + 4:pos + 8], data[pos + 8:pos + 8 + n]
        (crc,) = struct.unpack(">I", data[pos + 8 + n:pos + 12 + n])
        assert crc == zlib.crc32(kind + body) & 0xFFFFFFFF, f"CRC invalide : {kind!r}"
        chunks[kind] = body
        pos += 12 + n
    w, h, depth, ctype, comp, filt, inter = struct.unpack(">IIBBBBB", chunks[b"IHDR"])
    assert (depth, ctype, comp, filt, inter) == (8, 4, 0, 0, 0) and b"IEND" in chunks
    raw = np.frombuffer(zlib.decompress(chunks[b"IDAT"]), dtype=np.uint8).reshape(h, 1 + 2 * w)
    assert (raw[:, 0] == 0).all()
    return raw[:, 1:].reshape(h, w, 2)


def ramp(h=30, w=30, slope=1.0, toward="se"):
    """Plan incliné : monte vers le sud-est (face au nord-ouest, donc éclairé) ou vers le nord-ouest (dans l'ombre)."""
    rows, cols = np.mgrid[:h, :w].astype(float)
    z = slope * (rows + cols)                      # augmente vers le sud (lignes) et l'est (colonnes)
    return z if toward == "se" else -z


def test_hillshade_flat_ground_is_sine_of_light_altitude():
    sh = preview.hillshade(np.zeros((10, 10)), 1.0)
    assert np.allclose(sh, np.sin(np.radians(preview.LIGHT_ALTITUDE_DEG)))


def test_slope_facing_the_light_is_brighter_than_the_one_facing_away():
    lit = preview.hillshade(ramp(toward="se"), 1.0)[10:20, 10:20].mean()       # monte vers le SE : regarde le NO
    dark = preview.hillshade(ramp(toward="nw"), 1.0)[10:20, 10:20].mean()
    flat = np.sin(np.radians(preview.LIGHT_ALTITUDE_DEG))
    assert dark < flat < lit


def test_north_is_up_and_east_is_right():
    """Colline en pente régulière vers le nord (les lignes vont du nord au sud) : face sud, donc moins éclairée que
    la même pente tournée vers le nord (lumière du nord-ouest) ; idem est/ouest."""
    rows, cols = np.mgrid[:30, :30].astype(float)
    rises_north = -rows                                   # z diminue quand l'indice de ligne (vers le sud) augmente
    rises_south = rows
    assert (preview.hillshade(rises_north, 1.0)[10:20, 10:20].mean()
            < preview.hillshade(rises_south, 1.0)[10:20, 10:20].mean())
    rises_east, rises_west = cols, -cols
    assert (preview.hillshade(rises_east, 1.0)[10:20, 10:20].mean()
            > preview.hillshade(rises_west, 1.0)[10:20, 10:20].mean())


def test_vertical_exaggeration_deepens_the_shading_contrast():
    z = ramp(slope=0.05)
    flat = np.sin(np.radians(preview.LIGHT_ALTITUDE_DEG))
    gentle = abs(preview.hillshade(z, 1.0)[10:20, 10:20].mean() - flat)
    steep = abs(preview.hillshade(z * 8, 1.0)[10:20, 10:20].mean() - flat)
    assert steep > 4 * gentle


def test_render_masks_the_outside_and_keeps_the_inside_opaque():
    n = 40
    elev = ramp(n, n, slope=0.5)
    elev[:, 0] = -9999.0                                                   # colonne sans donnée hors de la pièce
    quad = np.zeros((n - 1, n - 1), dtype=bool)
    quad[5:30, 8:35] = True
    img = preview.render(elev, quad, 1.0)
    assert img.shape == (n, n, 2) and img.dtype == np.uint8
    nodes = np.zeros((n, n), dtype=bool)
    nodes[5:31, 8:36] = True
    assert (img[..., 1][nodes] == 255).all() and (img[..., 1][~nodes] == 0).all()
    assert img[..., 0][nodes].std() == 0 or img[..., 0][nodes].min() > 0       # pente régulière : ombrage uniforme


def test_render_downsamples_large_grids_by_whole_blocks():
    n = 130
    img = preview.render(ramp(n, n, slope=0.2), np.ones((n - 1, n - 1), dtype=bool), 1.0, max_side=40)
    assert max(img.shape[:2]) <= 40 and img.shape[0] == n // 4 and (img[..., 1] == 255).all()


def test_render_refuses_an_empty_piece():
    with pytest.raises(ValueError):
        preview.render(np.zeros((5, 5)), np.zeros((4, 4), dtype=bool), 1.0)


def test_png_roundtrip_is_exact_and_readable_by_an_independent_decoder(tmp_path):
    rng = np.random.default_rng(0)
    img = rng.integers(0, 256, size=(17, 23, 2), dtype=np.uint8)
    assert (decode_png(preview.png_bytes(img)) == img).all()
    path = tmp_path / "p.png"
    elev = ramp(40, 40, slope=0.3)
    w, h = preview.write_preview(path, elev, np.ones((39, 39), dtype=bool), 2.0, vexag=2.0)
    back = decode_png(path.read_bytes())
    assert back.shape == (h, w, 2) == (40, 40, 2)
    assert not list(tmp_path.glob("*.part"))
    try:
        from PIL import Image                      # contrôle croisé avec une vraie bibliothèque, si elle est là
    except ImportError:
        return
    with Image.open(path) as im:
        assert im.mode == "LA" and im.size == (40, 40)
        assert (np.asarray(im) == back).all()


def test_preview_path():
    paths = {"stl": "/tmp/x/cervin_150mm.stl"}
    assert preview.preview_path(paths, "auto") == "/tmp/x/cervin_150mm.png"
    assert preview.preview_path(paths, "/tmp/autre.png") == "/tmp/autre.png"
