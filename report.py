# SPDX-License-Identifier: GPL-3.0-or-later
"""
report.py
---------
Hauteur du relief imprimé et rapport JSON d'accompagnement (--report). Sans GDAL : seulement numpy et la
bibliothèque standard, donc testable sans fichier ni lecteur de raster.

* relief_stats / describe_relief : altitude minimale et maximale du modèle, hauteur de relief en mm au-dessus du
  socle, avec et sans exagération verticale (ce que le journal affiche avant le maillage).
* suggest_vexag : exagération INDICATIVE. Jamais appliquée : elle n'est qu'affichée. L'usage (fraise, épaisseur
  de planche, rendu voulu) décide.
* build_report / write_report : fichier JSON écrit à côté du modèle, avec la source, sa mention de licence, l'échelle,
  l'exagération et la ligne de commande, de quoi refaire ou citer la pièce.
"""
import datetime
import json
import math
import os
import sys

import numpy as np

SCHEMA_VERSION = 1

# Repère indicatif : relief imprimé d'environ 22,5 % du diamètre (45 mm pour 200 mm), exagération bornée à 1..8.
TARGET_RELIEF_RATIO = 0.225
MIN_SUGGESTED_VEXAG = 1.0
MAX_SUGGESTED_VEXAG = 8.0
HINT_MIN_GAP = 0.25                   # écart relatif minimal entre l'exagération demandée et le repère pour l'afficher


def touched_nodes(quad_ok):
    """Masque (H, W) des nœuds de grille qui sont un coin d'au moins une cellule de quad_ok (H-1, W-1) : ce sont les
    seuls dont l'altitude entre dans le modèle."""
    cp = np.pad(quad_ok, 1)
    return cp[:-1, :-1] | cp[:-1, 1:] | cp[1:, :-1] | cp[1:, 1:]


def relief_stats(min_elev_m, max_elev_m, scale_mm_per_m, vexag, diameter_mm, equivalent_relief_m=None):
    """Hauteur du relief au-dessus du socle. min_elev_m : altitude du plan z = 0 du modèle (le minimum, ou le niveau de
    la mer en mode île). equivalent_relief_m : dénivelé après une courbe d'exagération (--vexag-curve), qui remplace
    alors le dénivelé réel pour les hauteurs imprimées. Retourne un dict de float (JSON-compatible)."""
    relief_m = max(0.0, float(max_elev_m) - float(min_elev_m))
    base_m = relief_m if equivalent_relief_m is None else max(0.0, float(equivalent_relief_m))
    true_mm = base_m * scale_mm_per_m
    return {"min_elevation_m": float(min_elev_m), "max_elevation_m": float(max_elev_m), "relief_m": relief_m,
            "relief_equivalent_m": None if equivalent_relief_m is None else base_m,
            "relief_true_mm": true_mm, "relief_mm": true_mm * vexag, "vexag": float(vexag),
            "diameter_percent": 100.0 * true_mm * vexag / diameter_mm if diameter_mm > 0 else None,
            "suggested_vexag": suggest_vexag(true_mm, diameter_mm)}


def suggest_vexag(relief_true_mm, diameter_mm):
    """Exagération qui donnerait un relief d'environ TARGET_RELIEF_RATIO x diamètre, bornée à 1..8 ; None si le relief
    est nul. Indicatif seulement."""
    if not (relief_true_mm > 0 and diameter_mm > 0):
        return None
    v = TARGET_RELIEF_RATIO * diameter_mm / relief_true_mm
    return float(min(max(v, MIN_SUGGESTED_VEXAG), MAX_SUGGESTED_VEXAG))


def describe_relief(stats):
    """Lignes de journal (liste de str) : une ligne de fait, et un repère si l'exagération s'en écarte nettement."""
    if stats["relief_m"] <= 0:
        return [f"relief : nul (altitude {stats['min_elevation_m']:.1f} m partout), disque plat"]
    pct = stats["diameter_percent"]
    lines = [f"relief : altitudes de {stats['min_elevation_m']:.1f} à {stats['max_elevation_m']:.1f} m "
             f"(dénivelé {stats['relief_m']:.0f} m), soit {stats['relief_mm']:.1f} mm au-dessus du socle "
             f"(exagération x{stats['vexag']:g} ; {stats['relief_true_mm']:.1f} mm "
             + ("à l'échelle vraie" if stats["relief_equivalent_m"] is None
                else f"avec la courbe, sans l'exagération (dénivelé équivalent {stats['relief_equivalent_m']:.0f} m)")
             + (f" ; {pct:.0f} % du diamètre)" if pct is not None else ")")]
    s = stats["suggested_vexag"]
    if s is not None and stats["vexag"] > 0 and abs(s - stats["vexag"]) / stats["vexag"] >= HINT_MIN_GAP:
        lines.append(f"  repère indicatif (relief de {100 * TARGET_RELIEF_RATIO:.1f} % du diamètre, exagération "
                     f"bornée entre {MIN_SUGGESTED_VEXAG:g} et {MAX_SUGGESTED_VEXAG:g}) : --vexag {s:.1f} donnerait "
                     f"environ {s * stats['relief_true_mm']:.1f} mm. À ajuster selon l'usage.")
    return lines


def report_path(paths, requested):
    """Chemin du rapport. requested : « auto » (à côté du premier fichier de sortie, même nom en .json) ou un chemin."""
    if requested == "auto":
        return os.path.splitext(next(iter(paths.values())))[0] + ".json"
    return requested


def _clean(value):
    """Rend récursivement une valeur JSON-compatible : types numpy convertis, non-finis remplacés par None."""
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return _clean(value.tolist())
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def build_report(*, outputs, model, relief, source, mesh_info, location, options, argv=None, created=None):
    """Assemble le dictionnaire du rapport. Aucun calcul : chaque section vient de l'appelant.

    outputs : {format: chemin} ; model : forme, diamètre, socle, exagération, échelle ; relief : relief_stats() ;
    source : {used, input, fill_with_copernicus, attribution} ; mesh_info : nombres de triangles ; location :
    lat/lon/cx/cy/radius_m connus ; options : options de maillage ; argv : ligne de commande."""
    created = created or datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0)
    return _clean({
        "schema": SCHEMA_VERSION,
        "generator": "dem-to-stl (circle_dem_to_stl.py)",
        "created_utc": created.isoformat().replace("+00:00", "Z"),
        "outputs": {k: os.path.abspath(v) for k, v in outputs.items()},
        "model": model,
        "relief": relief,
        "source": source,
        "mesh": mesh_info,
        "location": location,
        "options": options,
        "argv": list(sys.argv[1:] if argv is None else argv),
    })


def write_report(path, data):
    """Écrit le JSON sous un nom temporaire puis le renomme : jamais de rapport tronqué. Lève OSError/ValueError."""
    text = json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    tmp = str(path) + ".part"
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            f.write(text)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
