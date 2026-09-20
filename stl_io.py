# SPDX-License-Identifier: GPL-3.0-or-later
"""
stl_io.py
---------
Écriture STL binaire. Isolé pour rester réutilisable indépendamment de
la façon dont les triangles ont été produits.
"""
import struct

import numpy as np


def write_stl_binary(tris, path, progress=None):
    n = len(tris)
    with open(path, "wb") as f:
        f.write(b"\x00" * 80)
        f.write(struct.pack("<I", n))
        for idx, t in enumerate(tris):
            if progress is not None:
                progress(idx, n, "Écriture STL")
            v0, v1, v2 = t[0], t[1], t[2]
            nvec = np.cross(v1 - v0, v2 - v0)
            norm = np.linalg.norm(nvec)
            if norm > 0:
                nvec = nvec / norm
            f.write(struct.pack("<3f", *nvec))
            for v in (v0, v1, v2):
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))
