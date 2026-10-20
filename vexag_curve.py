# SPDX-License-Identifier: GPL-3.0-or-later
"""
vexag_curve.py
--------------
Courbe d'exagération non linéaire (--vexag-curve) : au lieu d'un facteur unique, une correspondance « altitude réelle
-> altitude équivalente » appliquée AVANT l'échelle et l'exagération habituelles (--vexag continue de multiplier).

    --vexag-curve "0:0,500:500,2000:1100,4000:1800"

signifie : jusqu'à 500 m, relief à l'échelle vraie ; de 500 à 2000 m, 600 m équivalents pour 1500 m réels (facteur
0,4) ; de 2000 à 4000 m, 700 pour 2000 (0,35). Utile pour comprimer des sommets très hauts, ou dilater un relief bas, sur
un format limité.

Interpolation PCHIP (Fritsch-Carlson) : monotone et sans dépassement entre deux points, donc ni creux ni pic
artificiel ; la hauteur relative des points reste ordonnée (le plus haut reste le plus haut). Hors des points donnés,
prolongement linéaire de la pente du premier ou du dernier tronçon (courbe continue). Une partie plate (deux altitudes
équivalentes égales) écrête le relief.

Ce n'est PAS une réplique fidèle : les hauteurs relatives entre altitudes sont modifiées. Option facultative, signalée
dans le journal.

Pur numpy / scipy, sans GDAL.
"""
import math

import numpy as np
from scipy.interpolate import PchipInterpolator


class CurveError(ValueError):
    """Courbe invalide (message destiné à l'utilisateur)."""


class AltitudeCurve:
    def __init__(self, points):
        pts = [(float(a), float(b)) for a, b in points]
        if len(pts) < 2:
            raise CurveError("la courbe demande au moins deux points « altitude:altitude_équivalente ».")
        if not all(math.isfinite(v) for p in pts for v in p):
            raise CurveError("la courbe ne doit contenir que des nombres finis.")
        for (a0, b0), (a1, b1) in zip(pts, pts[1:]):
            if not a1 > a0:
                raise CurveError(f"les altitudes doivent être strictement croissantes ({a0:g} puis {a1:g}).")
            if b1 < b0:
                raise CurveError(f"les altitudes équivalentes ne peuvent pas décroître ({b0:g} puis {b1:g}) : "
                                 f"un relief inversé passerait sous le socle.")
        if pts[-1][1] == pts[0][1]:
            raise CurveError("toutes les altitudes équivalentes sont égales : le relief serait nul.")
        self.points = pts
        self._x = np.array([p[0] for p in pts])
        self._y = np.array([p[1] for p in pts])
        self._pchip = PchipInterpolator(self._x, self._y, extrapolate=False)
        self._slope_lo = (self._y[1] - self._y[0]) / (self._x[1] - self._x[0])
        self._slope_hi = (self._y[-1] - self._y[-2]) / (self._x[-1] - self._x[-2])

    def __call__(self, altitude):
        """Altitude équivalente (m), scalaire ou tableau. NaN reste NaN."""
        a = np.asarray(altitude, dtype=np.float64)
        inside = np.nan_to_num(a, nan=self._x[0])                         # évalué partout, NaN réinjecté ensuite
        out = np.where(inside < self._x[0], self._y[0] + (inside - self._x[0]) * self._slope_lo,
                       np.where(inside > self._x[-1], self._y[-1] + (inside - self._x[-1]) * self._slope_hi,
                                self._pchip(np.clip(inside, self._x[0], self._x[-1]))))
        out = np.where(np.isnan(a), np.nan, out)
        return float(out) if out.ndim == 0 else out

    def segments(self, vexag=1.0):
        """[(altitude début, altitude fin, facteur moyen)] : facteur = tronçon équivalent / tronçon réel x vexag."""
        return [(a0, a1, (b1 - b0) / (a1 - a0) * vexag) for (a0, b0), (a1, b1) in zip(self.points, self.points[1:])]

    def describe(self, vexag=1.0):
        """Lignes de journal : avertissement de non-fidélité et exagération moyenne par tronçon."""
        seg = " | ".join(f"{a0:g} à {a1:g} m : x{f:.2f}" for a0, a1, f in self.segments(vexag))
        return ["ATTENTION : --vexag-curve déforme le relief de façon non uniforme : les hauteurs relatives entre "
                "altitudes ne sont plus à l'échelle vraie, la pièce n'est pas une réplique fidèle.",
                f"  courbe (PCHIP, exagération moyenne par tronçon, --vexag compris) : {seg}"]

    def as_list(self):
        return [[a, b] for a, b in self.points]


def parse(text):
    """« alt:équivalente,alt:équivalente,... » (séparateurs : virgule, point-virgule ou espace ; point décimal)."""
    items = [t for t in text.replace(";", ",").replace(" ", ",").split(",") if t.strip()]
    pts = []
    for item in items:
        a, sep, b = item.partition(":")
        try:
            if not sep:
                raise ValueError
            pts.append((float(a), float(b)))
        except ValueError:
            raise CurveError(f"point illisible : {item!r} (attendu « altitude:altitude_équivalente », en mètres, "
                             f"avec un point décimal).") from None
    return AltitudeCurve(pts)
