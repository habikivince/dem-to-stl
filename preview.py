# SPDX-License-Identifier: GPL-3.0-or-later
"""
preview.py
----------
Aperçu ombré (hillshade) PNG de la zone réellement découpée (--preview), écrit AVANT le maillage : on voit tout de
suite un mauvais centre, des trous ou un relief aplati, sans attendre les minutes de calcul du maillage.

Pur numpy / scipy, sans GDAL ni bibliothèque d'images : le PNG est écrit à la main (zlib + struct).

* Le nord est en haut (les lignes de la grille vont du nord au sud).
* L'ombrage porte sur le relief IMPRIMÉ : les altitudes sont multipliées par l'exagération verticale, comme dans
  le modèle. Un relief plat à l'échelle vraie reste plat sur l'aperçu tant qu'on ne l'exagère pas.
* Hors de la pièce, l'image est transparente. Les nœuds de la pièce sont ceux de report.touched_nodes.
"""
import os
import struct
import zlib

import numpy as np

MAX_SIDE_PX = 1600
LIGHT_AZIMUTH_DEG = 315.0             # lumière venant du nord-ouest, convention cartographique
LIGHT_ALTITUDE_DEG = 45.0


def hillshade(z, pixel_m, azimuth_deg=LIGHT_AZIMUTH_DEG, altitude_deg=LIGHT_ALTITUDE_DEG):
    """Ombrage en [0, 1] d'une grille d'altitudes z (lignes du nord au sud, colonnes d'ouest en est, mêmes unités
    horizontales et verticales que pixel_m). Plat : sin(altitude_deg). z doit être entièrement fini."""
    dz_south, dz_east = np.gradient(z, pixel_m)
    dz_north = -dz_south
    norm = np.sqrt(dz_east ** 2 + dz_north ** 2 + 1.0)
    az, alt = np.radians(azimuth_deg), np.radians(altitude_deg)
    lx, ly, lz = np.sin(az) * np.cos(alt), np.cos(az) * np.cos(alt), np.sin(alt)
    return np.clip((-dz_east * lx - dz_north * ly + lz) / norm, 0.0, 1.0)


def render(elev, quad_ok, pixel_m, vexag=1.0, max_side=MAX_SIDE_PX):
    """Image (H', W', 2) uint8 (gris, opacité) de la pièce. elev (H, W) en mètres, quad_ok (H-1, W-1) : cellules
    incluses. Au-delà de max_side pixels, moyenne par blocs entiers."""
    from scipy import ndimage

    from report import touched_nodes
    node = touched_nodes(quad_ok)
    valid = node & np.isfinite(elev)
    if not valid.any():
        raise ValueError("aucun point valide à représenter")
    # les nœuds hors de la pièce prennent l'altitude du nœud valide le plus proche : pas de falaise factice au bord
    idx = ndimage.distance_transform_edt(~valid, return_distances=False, return_indices=True)
    z = np.asarray(elev, dtype=np.float64)[tuple(idx)] * vexag
    gray = 255.0 * hillshade(z, pixel_m)
    alpha = valid.astype(np.float64) * 255.0
    f = int(np.ceil(max(gray.shape) / max_side))
    if f > 1:
        h, w = (gray.shape[0] // f) * f, (gray.shape[1] // f) * f
        w_blocks = valid[:h, :w].reshape(h // f, f, w // f, f).sum(axis=(1, 3))
        g_sum = (gray[:h, :w] * valid[:h, :w]).reshape(h // f, f, w // f, f).sum(axis=(1, 3))
        gray = np.divide(g_sum, w_blocks, out=np.zeros_like(g_sum), where=w_blocks > 0)
        alpha = 255.0 * (w_blocks > 0)
    out = np.empty(gray.shape + (2,), dtype=np.uint8)
    out[..., 0] = np.clip(np.rint(gray), 0, 255)
    out[..., 1] = np.clip(np.rint(alpha), 0, 255)
    return out


def _chunk(kind, data):
    body = kind + data
    return struct.pack(">I", len(data)) + body + struct.pack(">I", zlib.crc32(body) & 0xFFFFFFFF)


def png_bytes(image):
    """PNG (niveaux de gris + opacité, 8 bits) d'un tableau (H, W, 2) uint8."""
    image = np.ascontiguousarray(image, dtype=np.uint8)
    if image.ndim != 3 or image.shape[2] != 2:
        raise ValueError("image attendue de forme (H, W, 2)")
    h, w = image.shape[:2]
    raw = np.zeros((h, 1 + 2 * w), dtype=np.uint8)           # octet de filtre 0 en tête de chaque ligne
    raw[:, 1:] = image.reshape(h, 2 * w)
    return (b"\x89PNG\r\n\x1a\n" + _chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 4, 0, 0, 0))
            + _chunk(b"IDAT", zlib.compress(raw.tobytes(), 6)) + _chunk(b"IEND", b""))


def preview_path(paths, requested):
    """Chemin de l'aperçu. requested : « auto » (à côté du premier fichier de sortie, même nom en .png) ou un chemin."""
    if requested == "auto":
        return os.path.splitext(next(iter(paths.values())))[0] + ".png"
    return requested


def write_preview(path, elev, quad_ok, pixel_m, vexag=1.0, max_side=MAX_SIDE_PX):
    """Écrit l'aperçu sous un nom temporaire puis le renomme. Retourne (largeur, hauteur) en pixels."""
    image = render(elev, quad_ok, pixel_m, vexag, max_side)
    tmp = str(path) + ".part"
    try:
        with open(tmp, "wb") as f:
            f.write(png_bytes(image))
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    return image.shape[1], image.shape[0]
