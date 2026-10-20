# SPDX-License-Identifier: GPL-3.0-or-later
"""
regular_polygon.py
------------------
Polygones réguliers convexes associés à un cercle de référence, avec arrondi exact des sommets.
Géométrie pure (numpy seul, aucun accès disque ni GDAL) : testable sans le moindre fichier, comme mesh.py.

Quatre longueurs, qu'il ne faut jamais confondre (toutes dans la même unité, quelle qu'elle soit) :

  ref_radius     rayon du cercle de RÉFÉRENCE, tel que donné par l'utilisateur ;
  apothem        apothème a : distance du centre au milieu d'un côté (rayon du cercle inscrit) ;
  circumradius   rayon circonscrit Rc : distance du centre à un sommet vif (rayon du cercle circonscrit) ;
  corner_radius  rayon d'arrondi r d'un sommet (rayon du congé, sans rapport avec le cercle de référence).

Mode « inscribed »     : le polygone est inscrit dans le cercle de référence, ses sommets vifs sont SUR le
                         cercle : Rc = ref_radius, a = Rc cos(pi/n).
Mode « circumscribed » : le polygone est circonscrit au cercle de référence, ses côtés lui sont TANGENTS :
                         a = ref_radius, Rc = a / cos(pi/n).

Le cercle de référence décrit toujours le polygone à angles vifs sous-jacent, arrondi ou non : l'arrondi
est une modification du polygone, il ne change pas la signification du rayon de référence.

Arrondi : chaque sommet est remplacé par un congé circulaire tangent aux deux côtés voisins, de rayon r.
Le centre du congé est à r / cos(pi/n) du sommet, sur la bissectrice ; les points de tangence sont à
r tan(pi/n) du sommet ; l'arc balaie 2 pi / n. Le contour arrondi est le polygone d'apothème (a - r) élargi
d'un disque de rayon r : il reste convexe, et r est limité par l'apothème (r = a donne un cercle de rayon a).

Orientation par défaut : le polygone est « debout » (un sommet vers le nord, +y, quand n est impair ; un côté
vers le nord quand n est pair), puis tourné de ``rotation_deg`` degrés dans le sens antihoraire.
"""
import math
from dataclasses import dataclass

import numpy as np

MIN_SIDES = 3
MAX_SIDES = 360
CIRCLE_MODES = ("inscribed", "circumscribed")
MAX_ARC_SEGMENTS = 20000          # garde-fou : tolérance absurdement fine
REL_TOL = 1e-9                    # tolérance relative sur les comparaisons de longueurs
BOUNDARY_TOL = 1e-9               # tolérance relative de contains() : frontière incluse de façon déterministe

_SIDE_NAMES = {3: "triangle", 4: "carré", 5: "pentagone", 6: "hexagone", 7: "heptagone", 8: "octogone",
               9: "ennéagone", 10: "décagone", 12: "dodécagone"}


class ShapeError(ValueError):
    """Paramètres de forme invalides."""


def _finite(v):
    return isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool) and math.isfinite(v)


@dataclass(frozen=True)
class RegularPolygon:
    n: int
    circle_mode: str
    ref_radius: float
    corner_radius: float = 0.0
    rotation_deg: float = 0.0

    def __post_init__(self):
        if isinstance(self.n, bool) or not isinstance(self.n, (int, np.integer)):
            raise ShapeError(f"le nombre de côtés doit être un entier (reçu : {self.n!r}).")
        if not MIN_SIDES <= self.n <= MAX_SIDES:
            raise ShapeError(f"le nombre de côtés doit être compris entre {MIN_SIDES} et {MAX_SIDES} "
                             f"(reçu : {self.n}).")
        if self.circle_mode not in CIRCLE_MODES:
            raise ShapeError(f"mode de cercle inconnu : {self.circle_mode!r} (attendu : {', '.join(CIRCLE_MODES)}).")
        if not _finite(self.ref_radius) or self.ref_radius <= 0:
            raise ShapeError(f"le rayon du cercle de référence doit être un nombre fini strictement positif "
                             f"(reçu : {self.ref_radius!r}).")
        if not _finite(self.corner_radius) or self.corner_radius < 0:
            raise ShapeError(f"le rayon d'arrondi doit être un nombre fini positif ou nul "
                             f"(reçu : {self.corner_radius!r}).")
        if not _finite(self.rotation_deg):
            raise ShapeError(f"la rotation doit être un nombre fini (reçu : {self.rotation_deg!r}).")
        if self.corner_radius > self.apothem * (1.0 + REL_TOL):
            raise ShapeError(f"le rayon d'arrondi ({self.corner_radius:g}) dépasse l'apothème du polygone "
                             f"({self.apothem:g}) : les congés de deux sommets voisins se chevaucheraient.")

    # ------------------------------------------------------------------ grandeurs dérivées
    @property
    def half_angle(self):
        """pi / n : demi-angle au centre d'un côté."""
        return math.pi / self.n

    @property
    def circumradius(self):
        """Rayon circonscrit Rc : centre -> sommet du polygone à angles vifs."""
        if self.circle_mode == "inscribed":
            return float(self.ref_radius)
        return float(self.ref_radius) / math.cos(self.half_angle)

    @property
    def apothem(self):
        """Apothème a : centre -> milieu d'un côté (inchangé par l'arrondi)."""
        if self.circle_mode == "circumscribed":
            return float(self.ref_radius)
        return float(self.ref_radius) * math.cos(self.half_angle)

    @property
    def side(self):
        """Longueur d'un côté du polygone à angles vifs."""
        return 2.0 * self.apothem * math.tan(self.half_angle)

    @property
    def max_corner_radius(self):
        """Plus grand rayon d'arrondi admis : l'apothème (le contour est alors un cercle de rayon a)."""
        return self.apothem

    @property
    def phase(self):
        """Angle (radians, depuis +x, antihoraire) du sommet numéro 0."""
        base = math.pi / 2.0 if self.n % 2 else math.pi / 2.0 - self.half_angle
        return base + math.radians(self.rotation_deg)

    @property
    def fillet_center_distance(self):
        """Distance du centre du polygone au centre d'un congé : (a - r) / cos(pi/n)."""
        return max(self.apothem - self.corner_radius, 0.0) / math.cos(self.half_angle)

    @property
    def extreme_radius(self):
        """Distance maximale du centre au contour réel (arrondi compris)."""
        return self.fillet_center_distance + self.corner_radius

    @property
    def flat_length(self):
        """Longueur rectiligne restante d'un côté après arrondi des deux extrémités."""
        return 2.0 * max(self.apothem - self.corner_radius, 0.0) * math.tan(self.half_angle)

    def area(self):
        a, r, h = self.apothem, self.corner_radius, self.half_angle
        return self.n * a * a * math.tan(h) - self.n * r * r * (math.tan(h) - h)

    def perimeter(self):
        a, r, h = self.apothem, self.corner_radius, self.half_angle
        return 2.0 * self.n * max(a - r, 0.0) * math.tan(h) + 2.0 * math.pi * r

    def bbox(self):
        """(xmin, xmax, ymin, ymax) exacts du contour (fonction d'appui du polygone intérieur + r)."""
        ang = self.phase + 2.0 * math.pi * np.arange(self.n) / self.n
        c, r = self.fillet_center_distance, self.corner_radius
        cx, cy = c * np.cos(ang), c * np.sin(ang)
        return (float(cx.min() - r), float(cx.max() + r), float(cy.min() - r), float(cy.max() + r))

    def half_extent(self):
        """Demi-côté du plus petit carré centré sur le centre du polygone et qui le contient."""
        return max(abs(v) for v in self.bbox())

    def scaled(self, k):
        """Même forme, toutes les longueurs multipliées par k (changement d'unité)."""
        return RegularPolygon(self.n, self.circle_mode, self.ref_radius * k, self.corner_radius * k,
                              self.rotation_deg)

    def name(self):
        return _SIDE_NAMES.get(self.n, f"polygone à {self.n} côtés")

    # ------------------------------------------------------------------ discrétisation des arcs
    def arc_segments(self, tol=None, segments=None):
        """Nombre de segments par arc d'arrondi (0 sans arrondi). ``segments`` l'impose ; sinon ``tol`` est la
        flèche maximale admise (écart entre une corde et son arc, même unité que les rayons) : les sommets sont
        sur l'arc exact, les cordes sont à l'intérieur, l'erreur est d'un seul côté et bornée par ``tol``."""
        r = self.corner_radius
        if r <= 0:
            return 0
        if segments is not None:
            if isinstance(segments, bool) or not isinstance(segments, (int, np.integer)) or segments < 1:
                raise ShapeError(f"le nombre de segments par arc doit être un entier >= 1 (reçu : {segments!r}).")
            k = int(segments)
        else:
            if not _finite(tol) or tol <= 0:
                raise ShapeError(f"la tolérance de flèche doit être un nombre fini strictement positif (reçu : {tol!r}).")
            if tol >= r:
                k = 1
            else:
                # angle au centre d'une corde de flèche tol : r (1 - cos(theta/2)) = tol, soit
                # theta = 4 asin(sqrt(tol / 2r)) (forme stable pour tol très petit, où acos(1 - tol/r) vaudrait 0)
                theta = 4.0 * math.asin(math.sqrt(tol / (2.0 * r)))
                k = max(1, math.ceil(2.0 * self.half_angle / theta))
        if k > MAX_ARC_SEGMENTS:
            raise ShapeError(f"{k} segments par arc seraient nécessaires (maximum {MAX_ARC_SEGMENTS}) : "
                             f"la tolérance de flèche est trop fine pour ce rayon d'arrondi.")
        return k

    def outline_points(self, tol=None, segments=None):
        """Sommets exacts du contour, dans le sens antihoraire, tableau (N, 2). Sans arrondi : les n sommets
        vifs. Avec arrondi : les arcs discrétisés, tous sur le vrai cercle de congé, extrémités de tangence
        comprises ; deux points consécutifs confondus (r = a, côtés réduits à un point) sont fusionnés."""
        n = self.n
        theta = self.phase + 2.0 * math.pi * np.arange(n) / n
        if self.corner_radius <= 0:
            return np.stack([self.circumradius * np.cos(theta), self.circumradius * np.sin(theta)], axis=1)
        k = self.arc_segments(tol, segments)
        r, c, h = self.corner_radius, self.fillet_center_distance, self.half_angle
        gamma = theta[:, None] - h + (2.0 * h) * np.arange(k + 1)[None, :] / k          # (n, k + 1)
        pts = np.stack([c * np.cos(theta)[:, None] + r * np.cos(gamma),
                        c * np.sin(theta)[:, None] + r * np.sin(gamma)], axis=-1).reshape(-1, 2)
        eps = 1e-12 * max(self.apothem, 1e-300)
        keep = np.ones(len(pts), dtype=bool)
        keep[1:] = np.linalg.norm(pts[1:] - pts[:-1], axis=1) > eps
        if len(pts) > 1 and np.linalg.norm(pts[0] - pts[-1]) <= eps:
            keep[0] = False
        return pts[keep]

    def discretize(self, tol=None, segments=None):
        return DiscretizedOutline(self, self.outline_points(tol, segments), self.arc_segments(tol, segments))

    # ------------------------------------------------------------------ appartenance et projection radiale
    def _local(self, x, y):
        """Coordonnées (qx, qy) du point dans le repère du secteur de sommet le plus proche (axe qx = direction
        du sommet), et angle local beta dans [-pi/n, pi/n]."""
        x, y = np.broadcast_arrays(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))
        psi = np.arctan2(y, x)
        step = 2.0 * self.half_angle
        beta = psi - self.phase
        beta = beta - np.rint(beta / step) * step
        rho = np.hypot(x, y)
        return rho * np.cos(beta), rho * np.sin(beta), beta

    def contains(self, x, y):
        """Masque booléen (broadcast de x et y) : points dans le contour, frontière comprise. Test exact, sans
        discrétisation des arcs. Dans le secteur d'un sommet, seuls comptent les deux côtés adjacents et, dans le
        cône issu du centre du congé, le disque d'arrondi."""
        qx, qy, _ = self._local(x, y)
        h, a, r = self.half_angle, self.apothem, self.corner_radius
        ch, sh = math.cos(h), math.sin(h)
        a_tol = a * (1.0 + BOUNDARY_TOL)          # un point sur la frontière, au bruit flottant près, est dedans
        inside = (qx * ch + qy * sh <= a_tol) & (qx * ch - qy * sh <= a_tol)
        if r > 0:
            dx, dy = qx - self.fillet_center_distance, qy
            in_cone = dx * math.tan(h) >= np.abs(dy)
            inside &= ~in_cone | (dx * dx + dy * dy <= (r * (1.0 + BOUNDARY_TOL)) ** 2)
        return inside

    def radial_distance(self, psi):
        """Distance du centre au contour dans la direction d'angle ``psi`` (radians). Le contour étant convexe et
        contenant le centre, elle est unique : c'est ce qui permet la projection radiale."""
        psi = np.asarray(psi, dtype=np.float64)
        step = 2.0 * self.half_angle
        beta = psi - self.phase
        beta = beta - np.rint(beta / step) * step
        h, a, r = self.half_angle, self.apothem, self.corner_radius
        rho_edge = a / np.cos(h - np.abs(beta))
        if r <= 0:
            return rho_edge
        c = self.fillet_center_distance
        sb, cb = np.sin(beta), np.cos(beta)
        disc = r * r - (c * sb) ** 2
        t_disk = c * cb + np.sqrt(np.maximum(disc, 0.0))             # sortie du rayon hors du disque d'arrondi
        dx, dy = t_disk * cb - c, t_disk * sb
        on_arc = (disc >= 0) & (dx * math.tan(h) >= np.abs(dy))      # le point de sortie est bien sur l'arc
        return np.where(on_arc, t_disk, rho_edge)

    def project(self, x, y):
        """Projette les points (x, y) sur le point du contour le PLUS PROCHE (projection perpendiculaire).

        Pourquoi pas radialement : le long d'un côté oblique, un escalier de cellules avance dans le sens du côté
        à chaque marche, mais pas dans le sens angulaire vu du centre : des marches « vers le centre » font reculer
        l'angle, deux nœuds distincts se retrouvent sur le même rayon et donc au même point, et la paroi repasse
        deux fois sur la même arête. Le pied de la perpendiculaire ne dépend que de la coordonnée le long du côté,
        qui, elle, ne recule jamais. (Pour un cercle, rayon et perpendiculaire coïncident.)

        Dans un secteur de sommet (repère local, axe qx = direction du sommet), le point le plus proche d'un point
        intérieur est : le pied sur la droite du côté adjacent le plus proche ; ou, dans le cône issu du centre du
        congé, le point du cercle d'arrondi aligné avec le centre du congé. Un point hors de la forme (rare : bord
        d'une cellule comblée par le nettoyage du masque) est projeté radialement."""
        x, y = np.broadcast_arrays(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))
        qx, qy, beta = self._local(x, y)
        h, a, r = self.half_angle, self.apothem, self.corner_radius
        ch, sh = math.cos(h), math.sin(h)
        d = a - (qx * ch + np.abs(qy) * sh)                      # distance à la droite du côté adjacent le plus proche
        fx = qx + d * ch
        fy = qy + d * np.where(qy >= 0, sh, -sh)
        if r > 0:
            c = self.fillet_center_distance
            dx, dy = qx - c, qy
            dist = np.hypot(dx, dy)
            safe = np.where(dist > 0, dist, 1.0)
            ux = np.where(dist > 0, dx / safe, 1.0)              # au centre exact du congé : direction du sommet
            uy = np.where(dist > 0, dy / safe, 0.0)
            in_cone = dx * math.tan(h) >= np.abs(dy)
            fx = np.where(in_cone, c + r * ux, fx)
            fy = np.where(in_cone, r * uy, fy)
        ang = np.arctan2(y, x) - beta                            # angle de l'axe local (direction du sommet)
        gx = fx * np.cos(ang) - fy * np.sin(ang)
        gy = fx * np.sin(ang) + fy * np.cos(ang)
        outside = ~self.contains(x, y)
        if outside.any():
            psi = np.arctan2(y, x)
            rho = self.radial_distance(psi)
            gx = np.where(outside, rho * np.cos(psi), gx)
            gy = np.where(outside, rho * np.sin(psi), gy)
        return gx, gy

    # ------------------------------------------------------------------ description
    def describe(self, unit="", digits=3):
        mode = ("inscrit : sommets sur le cercle de référence" if self.circle_mode == "inscribed"
                else "circonscrit : côtés tangents au cercle de référence")
        f = lambda v: f"{v:.{digits}f} {unit}".rstrip()  # noqa: E731
        s = (f"{self.name()} régulier {mode} ; rayon de référence {f(self.ref_radius)}, "
             f"apothème {f(self.apothem)}, rayon circonscrit {f(self.circumradius)}, côté {f(self.side)}")
        if self.corner_radius > 0:
            s += (f" ; sommets arrondis au rayon {f(self.corner_radius)} (rayon extrême du contour "
                  f"{f(self.extreme_radius)}, partie rectiligne {f(self.flat_length)} par côté)")
        return s


@dataclass(frozen=True)
class DiscretizedOutline:
    """Contour d'un RegularPolygon avec ses sommets exacts : ce que mesh.build_outline_rim consomme."""
    shape: RegularPolygon
    points: np.ndarray            # (N, 2), antihoraire
    segments_per_arc: int

    def project(self, x, y):
        return self.shape.project(x, y)
