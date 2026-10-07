# SPDX-License-Identifier: GPL-3.0-or-later
"""
mesh_export.py
--------------
Exports 3MF et OBJ, en complément du STL (stl_io.py), à partir du MÊME maillage.

Principe : les triangles sont d'abord arrondis en float32 (comme dans le STL), puis les
sommets identiques sont fusionnés (maillage indexé, ce qu'exigent 3MF et OBJ). Les
coordonnées sont écrites avec 9 chiffres significatifs, suffisants pour relire exactement un
float32 : relues, elles redonnent exactement les valeurs du STL.
Unités : millimètre partout (explicite dans le 3MF ; l'OBJ n'a pas d'unité, indiquée en
commentaire). Orientation : ordre anti-horaire vu de l'extérieur (normales sortantes), le
même que le STL.

3MF (spécification Core de 3MF Consortium) : paquet ZIP/OPC avec [Content_Types].xml,
_rels/.rels et 3D/3dmodel.model ; unité `millimeter` ; métadonnées Core (Title, Description,
Copyright, CreationDate, Application...). Aucune dépendance en dehors de la bibliothèque
standard et numpy. La bibliothèque officielle `lib3mf` n'est utilisée que par les tests,
comme validateur indépendant, si elle est installée.
"""
import datetime
import os
import re
import xml.etree.ElementTree as ET
import zipfile
from xml.sax.saxutils import escape, quoteattr

import numpy as np

CORE_NS = "http://schemas.microsoft.com/3dmanufacturing/core/2015/02"
FORMATS = ("stl", "3mf", "obj")
CHUNK = 200_000

CONTENT_TYPES = ('<?xml version="1.0" encoding="UTF-8"?>\n'
                 '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                 '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                 '<Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>'
                 '</Types>')
RELS = ('<?xml version="1.0" encoding="UTF-8"?>\n'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Target="/3D/3dmodel.model" Id="rel0" '
        'Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>')


def index_mesh(tris):
    """(sommets float32 (V, 3), faces int64 (F, 3)) : fusionne les sommets identiques après
    arrondi float32, sans changer l'ordre ni l'orientation des triangles."""
    flat = np.ascontiguousarray(np.asarray(tris).reshape(-1, 3), dtype=np.float32)
    n = len(flat)
    if n == 0:
        return flat, np.zeros((0, 3), dtype=np.int64)
    order = np.lexsort((flat[:, 2], flat[:, 1], flat[:, 0]))
    pts = flat[order]
    new = np.ones(n, dtype=bool)
    new[1:] = np.any(pts[1:] != pts[:-1], axis=1)
    ids = np.empty(n, dtype=np.int64)
    ids[order] = np.cumsum(new) - 1
    return pts[new], ids.reshape(-1, 3)


def signed_volume(vertices, faces):
    """Volume signé (mm³) : positif si les normales sont sortantes."""
    v = vertices.astype(np.float64)
    a, b, c = v[faces[:, 0]], v[faces[:, 1]], v[faces[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def _positional(x):
    """Écriture sans exposant (certains lecteurs OBJ/3MF s'en accommodent mal) qui relit le même float32."""
    return np.format_float_positional(np.float32(x), unique=True, trim="-")


def _vertex_lines(vertices, fmt):
    """Coordonnées en %.9g : 9 chiffres significatifs suffisent pour relire exactement un float32,
    donc la même géométrie que le STL. Les rares valeurs notées avec exposant (|x| < 1e-4)
    sont réécrites sans exposant."""
    for s in range(0, len(vertices), CHUNK):
        block = vertices[s:s + CHUNK].astype(np.float64).tolist()
        text = "".join(fmt % (x, y, z) for x, y, z in block)
        if "e" in text:
            text = "".join(fmt.replace("%.9g", "%s") % tuple(_positional(c) for c in row) for row in block)
        yield text


def _face_lines(faces, fmt, offset):
    for s in range(0, len(faces), CHUNK):
        yield "".join(fmt % (a + offset, b + offset, c + offset) for a, b, c in faces[s:s + CHUNK].tolist())


def write_3mf(path, vertices, faces, title="dem-to-stl", description="", copyright_text="", application="dem-to-stl",
              created=None):
    created = created or datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    meta = [("Title", title), ("Description", description), ("Copyright", copyright_text),
            ("CreationDate", created), ("Application", application)]
    head = ('<?xml version="1.0" encoding="UTF-8"?>\n'
            f'<model unit="millimeter" xml:lang="en-US" xmlns="{CORE_NS}">')
    head += "".join(f"<metadata name={quoteattr(k)}>{escape(v)}</metadata>" for k, v in meta if v)
    head += f"<resources><object id=\"1\" type=\"model\" name={quoteattr(title)}><mesh><vertices>"
    tmp = str(path) + ".part"
    try:
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED, compresslevel=4) as z:
            z.writestr("[Content_Types].xml", CONTENT_TYPES)
            z.writestr("_rels/.rels", RELS)
            with z.open("3D/3dmodel.model", "w", force_zip64=True) as m:
                m.write(head.encode("utf-8"))
                for block in _vertex_lines(vertices, '<vertex x="%.9g" y="%.9g" z="%.9g"/>'):
                    m.write(block.encode("ascii"))
                m.write(b"</vertices><triangles>")
                for block in _face_lines(faces, '<triangle v1="%d" v2="%d" v3="%d"/>', 0):
                    m.write(block.encode("ascii"))
                m.write(b'</triangles></mesh></object></resources><build><item objectid="1"/></build></model>')
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def write_obj(path, vertices, faces, title="dem-to-stl", comment=""):
    tmp = str(path) + ".part"
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as f:
            f.write(f"# dem-to-stl : {title}\n")
            f.write("# unité : millimètre (le format OBJ n'a pas d'unité ; axe Z vers le haut)\n")
            if comment:
                f.write(f"# {comment}\n")
            f.write(f"# {len(vertices)} sommets, {len(faces)} triangles\n")
            f.write("o " + re.sub(r"\s+", "_", title) + "\n")
            for block in _vertex_lines(vertices, "v %.9g %.9g %.9g\n"):
                f.write(block)
            for block in _face_lines(faces, "f %d %d %d\n", 1):
                f.write(block)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


# ------------------------------------------------------------------ relecture (tests, check_stl.py)

def read_obj(path):
    verts, faces = [], []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("v "):
                verts.append([float(t) for t in line.split()[1:4]])
            elif line.startswith("f "):
                faces.append([int(t.split("/")[0]) - 1 for t in line.split()[1:]])
    return np.array(verts, dtype=np.float32), np.array(faces, dtype=np.int64)


def read_3mf(path, flush=200_000):
    """(sommets, faces, unité, métadonnées) d'un 3MF, relu en flux (analyse XML incrémentale) :
    la mémoire reste raisonnable même pour des millions de triangles."""
    vchunks, fchunks, vbuf, fbuf, meta = [], [], [], [], {}
    unit, parents = None, {}
    with zipfile.ZipFile(path) as z, z.open("3D/3dmodel.model") as fh:
        for event, el in ET.iterparse(fh, events=("start", "end")):
            tag = el.tag.rsplit("}", 1)[-1]
            if event == "start":
                if tag == "model":
                    unit = el.get("unit")
                elif tag in ("vertices", "triangles"):
                    parents[tag] = el
                continue
            if tag == "vertex":
                vbuf.append((float(el.get("x")), float(el.get("y")), float(el.get("z"))))
                if len(vbuf) >= flush:
                    vchunks.append(np.array(vbuf, dtype=np.float32)); vbuf.clear(); parents["vertices"].clear()
            elif tag == "triangle":
                fbuf.append((int(el.get("v1")), int(el.get("v2")), int(el.get("v3"))))
                if len(fbuf) >= flush:
                    fchunks.append(np.array(fbuf, dtype=np.int64)); fbuf.clear(); parents["triangles"].clear()
            elif tag == "metadata" and el.get("name"):
                meta[el.get("name")] = el.text or ""
    if vbuf:
        vchunks.append(np.array(vbuf, dtype=np.float32))
    if fbuf:
        fchunks.append(np.array(fbuf, dtype=np.int64))
    verts = np.concatenate(vchunks) if vchunks else np.zeros((0, 3), dtype=np.float32)
    faces = np.concatenate(fchunks) if fchunks else np.zeros((0, 3), dtype=np.int64)
    return verts, faces, unit, meta


def read_triangles(path):
    """Triangles (n, 3, 3) float64 d'un STL binaire, d'un OBJ ou d'un 3MF (d'après l'extension)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".stl":
        import check_stl
        return check_stl.read_stl(path)
    if ext == ".obj":
        v, f = read_obj(path)
    elif ext == ".3mf":
        v, f, _, _ = read_3mf(path)
    else:
        raise ValueError(f"extension non gérée : {ext} (attendu .stl, .obj ou .3mf)")
    return v.astype(np.float64)[f]
