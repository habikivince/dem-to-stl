#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""
compare_dems.py
---------------
Compare deux DEM (fichiers, VRT ou listes de fichiers) sur la même zone et à la
résolution de travail du script principal, pour vérifier qu'une source automatique
ne perd rien par rapport à des données téléchargées à la main.

Les deux rasters passent par raster.warp_window (même fenêtre, même CRS cible,
même ré-échantillonnage 'average'), mais SANS le comblement des trous (FillNodata)
du pipeline : sinon la bordure d'un raster qui couvre moins de surface est
extrapolée sur jusqu'à 50 pixels et fausse la comparaison. Seules les cellules
réellement mesurées dans les deux rasters sont comparées. La résolution de travail est celle affichée par circle_dem_to_stl.py
(« résolution de ré-échantillonnage »).

Exemple :
  python compare_dems.py manuel_05m.vrt swisstopo_2m.vrt --lat 45.9763 --lon 7.6586 \
         --radius 3000 --pixel-size 8
"""
import argparse
import sys

import numpy as np

import raster


def compare_dems(a, b, lat, lon, cx, cy, radius, pixel_size):
    """Statistiques de B - A sur les cellules valides des deux rasters, dans le disque."""
    wkt, cx, cy = raster.resolve_target_crs_and_center(a, lat, lon, cx, cy)
    quiet = lambda *args: 1  # noqa: E731
    ea, xs, ys, na = raster.warp_window(a, wkt, cx, cy, radius, pixel_size, progress_cb=quiet, fill=False)
    eb, _, _, nb = raster.warp_window(b, wkt, cx, cy, radius, pixel_size, progress_cb=quiet, fill=False)
    inside = (xs[None, :] - cx) ** 2 + (ys[:, None] - cy) ** 2 <= radius ** 2

    def valid(e, nd):
        ok = np.isfinite(e)
        return ok & (e != nd) if nd is not None else ok

    va, vb = valid(ea, na) & inside, valid(eb, nb) & inside
    both = va & vb
    if not both.any():
        raise ValueError("aucune cellule valide commune aux deux rasters dans le disque")
    d = eb[both] - ea[both]
    return {
        "cells_disc": int(inside.sum()), "cells_a": int(va.sum()), "cells_b": int(vb.sum()),
        "cells_both": int(both.sum()),
        "cells_only_a": int((va & ~vb).sum()), "cells_only_b": int((vb & ~va).sum()),
        "relief_a_m": float(ea[va].max() - ea[va].min()),
        "mean": float(d.mean()), "std": float(d.std()), "rms": float(np.sqrt(np.mean(d ** 2))),
        "mean_abs": float(np.abs(d).mean()), "p95_abs": float(np.percentile(np.abs(d), 95)),
        "max_abs": float(np.abs(d).max()),
        "share_gt_1m": float(np.mean(np.abs(d) > 1.0)), "share_gt_5m": float(np.mean(np.abs(d) > 5.0)),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("a", help="DEM de référence (ex. téléchargé à la main, haute résolution)")
    p.add_argument("b", help="DEM à comparer (ex. fichier mis en cache par --source)")
    p.add_argument("--lat", type=float)
    p.add_argument("--lon", type=float)
    p.add_argument("--cx", type=float)
    p.add_argument("--cy", type=float)
    p.add_argument("--radius", type=float, required=True)
    p.add_argument("--pixel-size", type=float, required=True, help="m/pixel de comparaison")
    args = p.parse_args()
    if (args.lat is None) == (args.cx is None) or (args.lat is None) != (args.lon is None) \
            or (args.cx is None) != (args.cy is None):
        sys.exit("Erreur : fournis --lat/--lon OU --cx/--cy (dans le CRS de A).")
    try:
        s = compare_dems(args.a, args.b, args.lat, args.lon, args.cx, args.cy, args.radius, args.pixel_size)
    except ValueError as e:
        sys.exit(f"Erreur : {e}")
    print(f"Cellules dans le disque : {s['cells_disc']} ; valides A : {s['cells_a']} ; "
          f"valides B : {s['cells_b']} ; communes : {s['cells_both']} "
          f"(seulement A : {s['cells_only_a']}, seulement B : {s['cells_only_b']})")
    print(f"Relief de A sur la zone : {s['relief_a_m']:.1f} m")
    print(f"Écart B - A : moyenne {s['mean']:+.3f} m (décalage systématique éventuel, ex. référence "
          f"d'altitude), écart-type {s['std']:.3f} m, RMS {s['rms']:.3f} m")
    print(f"|écart| : moyen {s['mean_abs']:.3f} m, 95e centile {s['p95_abs']:.3f} m, max {s['max_abs']:.3f} m")
    print(f"Part des cellules > 1 m : {100 * s['share_gt_1m']:.2f} % ; > 5 m : {100 * s['share_gt_5m']:.2f} %")


if __name__ == "__main__":
    main()
