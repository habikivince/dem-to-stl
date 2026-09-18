#!/usr/bin/env python3
"""
island_dem_to_stl.py
---------------------
Convertit un DEM (GeoTIFF ou VRT, CRS projeté en mètres) couvrant une île
entière en un STL solide et étanche, découpé selon le contour réel des
données valides (= le littoral) plutôt qu'un cercle arbitraire — utile
quand la donnée en mer est absente ou aberrante (artefacts d'interpolation
bien au-delà du nodata déclaré).

Usage :
    python island_dem_to_stl.py --input reunion_mosaic.vrt --output reunion.stl \
        --size-mm 250 --vexag 3.0 --min-elevation -50
"""
import argparse
import struct
import sys

import numpy as np
from osgeo import gdal

gdal.UseExceptions()


def log(msg):
    print(f"[ISLAND2STL] {msg}", flush=True)


def erode3x3(m):
    p = np.pad(m, 1, constant_values=False)
    out = m.copy()
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            out &= p[1 + di:1 + di + m.shape[0], 1 + dj:1 + dj + m.shape[1]]
    return out


def dilate3x3(m):
    p = np.pad(m, 1, constant_values=False)
    out = np.zeros_like(m)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            out |= p[1 + di:1 + di + m.shape[0], 1 + dj:1 + dj + m.shape[1]]
    return out


def fix_diagonal_pinches(quad_ok):
    for _ in range(4):
        diag1 = quad_ok[:-1, :-1] & quad_ok[1:, 1:] & ~quad_ok[:-1, 1:] & ~quad_ok[1:, :-1]
        diag2 = quad_ok[:-1, 1:] & quad_ok[1:, :-1] & ~quad_ok[:-1, :-1] & ~quad_ok[1:, 1:]
        if not (diag1.any() or diag2.any()):
            break
        fix = diag1 | diag2
        quad_ok[:-1, 1:] |= fix
        quad_ok[1:, :-1] |= fix
        quad_ok[:-1, :-1] |= diag2
        quad_ok[1:, 1:] |= diag2
    return quad_ok


def build_mesh(input_path, size_mm, vexag, base_mm, print_spacing_mm, min_elevation):
    src = gdal.Open(input_path)
    gt0 = src.GetGeoTransform()
    width_m = abs(gt0[1]) * src.RasterXSize
    height_m = abs(gt0[5]) * src.RasterYSize
    xmin = gt0[0]
    ymax = gt0[3]
    xmax = xmin + gt0[1] * src.RasterXSize
    ymin = ymax + gt0[5] * src.RasterYSize
    nodata = src.GetRasterBand(1).GetNoDataValue()
    log(f"emprise source : {width_m/1000:.2f} x {height_m/1000:.2f} km")

    scale = size_mm / max(width_m, height_m)
    pixel_size_m = print_spacing_mm / scale
    log(f"échelle : {scale:.8f} mm/m | résolution de ré-échantillonnage : {pixel_size_m:.2f} m/pixel")

    warped = gdal.Warp(
        "", input_path, format="MEM",
        outputBounds=(xmin, ymin, xmax, ymax),
        xRes=pixel_size_m, yRes=pixel_size_m,
        resampleAlg="average",
        dstNodata=nodata,
    )
    band = warped.GetRasterBand(1)
    gdal.FillNodata(band, None, maxSearchDist=50, smoothingIterations=0)
    elev = band.ReadAsArray().astype(np.float64)
    wgt = warped.GetGeoTransform()
    H, W = elev.shape
    log(f"grille ré-échantillonnée : {W} x {H} points")

    xs = wgt[0] + (np.arange(W) + 0.5) * wgt[1]
    ys = wgt[3] + (np.arange(H) + 0.5) * wgt[5]

    valid = (elev != nodata) if nodata is not None else np.ones_like(elev, dtype=bool)
    mask = valid & (elev > min_elevation)
    n_excluded_artifacts = int((valid & (elev <= min_elevation)).sum())
    if n_excluded_artifacts:
        log(f"{n_excluded_artifacts} pixels exclus comme artefacts implausibles (<= {min_elevation} m, hors nodata déclaré)")
    if not mask.any():
        sys.exit("Erreur : aucun pixel valide trouvé (seuil --min-elevation trop haut ?).")

    min_elev = elev[mask].min()
    z_top_mm = (elev - min_elev) * scale * vexag

    Xmm = xs * scale
    Ymm = ys * scale

    quad_ok = mask[:-1, :-1] & mask[1:, :-1] & mask[:-1, 1:] & mask[1:, 1:]
    quad_ok = dilate3x3(erode3x3(quad_ok))
    quad_ok = fix_diagonal_pinches(quad_ok)
    if not quad_ok.any():
        sys.exit("Erreur : masque vide après nettoyage morphologique.")
    log(f"{int(quad_ok.sum())} cellules incluses (silhouette terre)")

    ii, jj = np.nonzero(quad_ok)

    def corner(di, dj):
        i, j = ii + di, jj + dj
        return np.stack([Xmm[j], Ymm[i], z_top_mm[i, j]], axis=1)

    c00, c10, c01, c11 = corner(0, 0), corner(1, 0), corner(0, 1), corner(1, 1)

    top_tris = np.concatenate([
        np.stack([c00, c10, c01], axis=1),
        np.stack([c10, c11, c01], axis=1),
    ])

    zb = -base_mm
    c00b, c10b, c01b, c11b = c00.copy(), c10.copy(), c01.copy(), c11.copy()
    for c in (c00b, c10b, c01b, c11b):
        c[:, 2] = zb

    bot_tris = np.concatenate([
        np.stack([c00b, c01b, c10b], axis=1),
        np.stack([c10b, c01b, c11b], axis=1),
    ])

    wall_tris = []

    def add_wall(p1t, p2t, outward):
        p1b = (p1t[0], p1t[1], zb)
        p2b = (p2t[0], p2t[1], zb)
        if outward:
            wall_tris.append((p1t, p2t, p1b))
            wall_tris.append((p2t, p2b, p1b))
        else:
            wall_tris.append((p2t, p1t, p1b))
            wall_tris.append((p2b, p2t, p1b))

    H1, W1 = quad_ok.shape
    for i in range(H1):
        for j in range(W1 + 1):
            left = quad_ok[i, j - 1] if j > 0 else False
            right = quad_ok[i, j] if j < W1 else False
            if left != right:
                p1 = (Xmm[j], Ymm[i], z_top_mm[i, j])
                p2 = (Xmm[j], Ymm[i + 1], z_top_mm[i + 1, j])
                add_wall(p1, p2, outward=bool(left))
    for j in range(W1):
        for i in range(H1 + 1):
            up = quad_ok[i - 1, j] if i > 0 else False
            down = quad_ok[i, j] if i < H1 else False
            if up != down:
                p1 = (Xmm[j], Ymm[i], z_top_mm[i, j])
                p2 = (Xmm[j + 1], Ymm[i], z_top_mm[i, j + 1])
                add_wall(p1, p2, outward=bool(down))

    wall_tris = np.array(wall_tris, dtype=np.float64)
    all_tris = np.concatenate([top_tris, bot_tris, wall_tris])
    log(f"{len(top_tris)} triangles (dessus) + {len(bot_tris)} (dessous) + {len(wall_tris)} (paroi) "
        f"= {len(all_tris)} triangles au total")
    return all_tris


def write_stl_binary(tris, path):
    n = len(tris)
    with open(path, "wb") as f:
        f.write(b"\x00" * 80)
        f.write(struct.pack("<I", n))
        for t in tris:
            v0, v1, v2 = t[0], t[1], t[2]
            nvec = np.cross(v1 - v0, v2 - v0)
            norm = np.linalg.norm(nvec)
            if norm > 0:
                nvec = nvec / norm
            f.write(struct.pack("<3f", *nvec))
            for v in (v0, v1, v2):
                f.write(struct.pack("<3f", *v))
            f.write(struct.pack("<H", 0))


def main():
    p = argparse.ArgumentParser(description="Découpe par silhouette (littoral) + export STL solide.")
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--size-mm", type=float, default=250.0, help="Taille cible de la plus grande dimension, en mm")
    p.add_argument("--vexag", type=float, default=1.0)
    p.add_argument("--base-mm", type=float, default=3.0)
    p.add_argument("--print-spacing-mm", type=float, default=0.2)
    p.add_argument("--min-elevation", type=float, default=-50.0,
                    help="Seuil sous lequel une valeur (hors nodata déclaré) est traitée comme un artefact, en m")
    args = p.parse_args()

    tris = build_mesh(args.input, args.size_mm, args.vexag, args.base_mm,
                       args.print_spacing_mm, args.min_elevation)
    write_stl_binary(tris, args.output)
    log(f"Terminé : {args.output}")


if __name__ == "__main__":
    main()
