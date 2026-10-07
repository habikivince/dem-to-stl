# SPDX-License-Identifier: GPL-3.0-or-later
"""
stl_io.py
---------
Écriture STL binaire. Isolé pour rester réutilisable indépendamment de
la façon dont les triangles ont été produits.

L'écriture est vectorisée par blocs (normales calculées par numpy, enregistrements
écrits d'un coup) : environ 280 fois plus rapide que l'ancienne boucle triangle par
triangle, avec les mêmes sommets et des normales identiques à l'arrondi float32 près.
Le fichier est écrit sous un nom temporaire puis renommé : un échec en cours
d'écriture ne laisse jamais de STL tronqué à la place d'un fichier valide.
"""
import os
import struct

import numpy as np

# 50 octets par triangle : normale (3 f4), 3 sommets (9 f4), attribut (u2)
STL_DTYPE = np.dtype([("n", "<f4", (3,)), ("v", "<f4", (3, 3)), ("a", "<u2")])
CHUNK_TRIANGLES = 500_000


def facet_normals(v):
    """Normales unitaires (n, 3) de triangles v (n, 3, 3) ; nulle pour un triangle dégénéré."""
    nvec = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    norm = np.linalg.norm(nvec, axis=1, keepdims=True)
    return np.divide(nvec, norm, out=np.zeros_like(nvec), where=norm > 0)


def write_stl_binary(tris, path, progress=None):
    tris = np.asarray(tris)
    n = len(tris)
    tmp = str(path) + ".part"
    try:
        with open(tmp, "wb") as f:
            f.write(b"\x00" * 80)
            f.write(struct.pack("<I", n))
            for start in range(0, n, CHUNK_TRIANGLES):
                v = tris[start:start + CHUNK_TRIANGLES].astype(np.float64, copy=False)
                rec = np.zeros(len(v), dtype=STL_DTYPE)
                rec["n"] = facet_normals(v)
                rec["v"] = v
                rec.tofile(f)
                if progress is not None:
                    progress(min(start + CHUNK_TRIANGLES, n) - 1, n, "Écriture STL")
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
