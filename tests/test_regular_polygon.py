# SPDX-License-Identifier: GPL-3.0-or-later
"""Polygones réguliers inscrits ou circonscrits à un cercle de référence, à angles vifs ou arrondis (--sides).

Trois niveaux, tous fondés sur des propriétés géométriques calculées à part du code testé :
  1. regular_polygon.py seul : relations rayon / apothème / côté, tangence, aire et périmètre en forme
     fermée comparés à la formule du lacet sur un contour finement discrétisé, appartenance, projection sur
     le point le plus proche (force brute), flèche des arcs, paramètres invalides ;
  2. mesh.build_mesh_from_arrays(outline=...) sur des grilles synthétiques : étanchéité, orientation,
     triangles non dégénérés, volume d'un relief plat = socle x aire exacte (à l'erreur de flèche près),
     sommets de la paroi exactement sur le contour ;
  3. la ligne de commande : validations, --dry-run, STL de bout en bout, rayon de téléchargement."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_stl  # noqa: E402
import circle_dem_to_stl  # noqa: E402
import download  # noqa: E402
import mesh  # noqa: E402
import regular_polygon  # noqa: E402
from conftest import write_synthetic_tif  # noqa: E402
from regular_polygon import RegularPolygon, ShapeError  # noqa: E402

gdal.UseExceptions()

MODES = ("inscribed", "circumscribed")


# ---------------------------------------------------------------- outils indépendants du code testé

def shoelace(P):
    x, y = P[:, 0], P[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def polyline_length(P):
    return float(np.linalg.norm(np.roll(P, -1, axis=0) - P, axis=1).sum())


def dist_to_polyline(Q, P):
    """Distance de chaque point de Q (m, 2) à la ligne fermée P (k, 2)."""
    A, B = P, np.roll(P, -1, axis=0)
    AB = B - A
    best = np.full(len(Q), np.inf)
    for a, ab in zip(A, AB):
        t = np.clip(((Q - a) @ ab) / (ab @ ab), 0.0, 1.0)
        best = np.minimum(best, np.linalg.norm(Q - (a + t[:, None] * ab), axis=1))
    return best


def inside_convex(Q, P):
    """Appartenance au polygone convexe antihoraire P (k, 2) : produit vectoriel >= 0 pour chaque arête."""
    ok = np.ones(len(Q), dtype=bool)
    for a, b in zip(P, np.roll(P, -1, axis=0)):
        ok &= (b[0] - a[0]) * (Q[:, 1] - a[1]) - (b[1] - a[1]) * (Q[:, 0] - a[0]) >= 0.0
    return ok


def make(n, mode, frac=0.0, rot=0.0, ref=10.0):
    """Polygone dont le rayon d'arrondi vaut `frac` fois l'apothème."""
    a = RegularPolygon(n, mode, ref, 0.0, rot).apothem
    return RegularPolygon(n, mode, ref, frac * a, rot)


def sharp_vertices(s):
    ang = s.phase + 2.0 * math.pi * np.arange(s.n) / s.n
    return np.stack([s.circumradius * np.cos(ang), s.circumradius * np.sin(ang)], axis=1)


def fine(s, tol=1e-6):
    return s.outline_points(tol=tol)


SIDES = (3, 4, 5, 6, 7, 8, 12, 50)
FRACS = (0.0, 0.05, 0.5, 0.999, 1.0)


# ================================================================ 1. géométrie pure

@pytest.mark.parametrize("n", range(3, 41))
def test_radius_apothem_circumradius_side_relations(n):
    R, h = 10.0, math.pi / n
    ins, cir = RegularPolygon(n, "inscribed", R), RegularPolygon(n, "circumscribed", R)
    assert ins.circumradius == R and ins.apothem == pytest.approx(R * math.cos(h), rel=1e-14)
    assert ins.side == pytest.approx(2 * R * math.sin(h), rel=1e-12)
    assert cir.apothem == R and cir.circumradius == pytest.approx(R / math.cos(h), rel=1e-14)
    assert cir.side == pytest.approx(2 * R * math.tan(h), rel=1e-12)
    for s in (ins, cir):                       # côté = distance réelle entre deux sommets voisins
        V = sharp_vertices(s)
        assert np.linalg.norm(V[1] - V[0]) == pytest.approx(s.side, rel=1e-12)
        assert s.area() == pytest.approx(shoelace(V), rel=1e-12)
        assert s.perimeter() == pytest.approx(n * s.side, rel=1e-12)


@pytest.mark.parametrize("n", SIDES)
def test_inscribed_vertices_lie_on_the_reference_circle(n):
    s = RegularPolygon(n, "inscribed", 10.0)
    V = sharp_vertices(s)
    assert np.allclose(np.hypot(V[:, 0], V[:, 1]), 10.0, rtol=0, atol=1e-12)
    # le polygone est dans le cercle : tout point du contour est à distance <= R, avec égalité aux seuls sommets
    P = fine(s)
    assert np.hypot(P[:, 0], P[:, 1]).max() == pytest.approx(10.0, abs=1e-12)
    mids = 0.5 * (V + np.roll(V, -1, axis=0))
    assert np.allclose(np.hypot(mids[:, 0], mids[:, 1]), 10.0 * math.cos(math.pi / n), atol=1e-12)


@pytest.mark.parametrize("n", SIDES)
def test_circumscribed_sides_are_tangent_to_the_reference_circle(n):
    s = RegularPolygon(n, "circumscribed", 10.0)
    V = sharp_vertices(s)
    mids = 0.5 * (V + np.roll(V, -1, axis=0))
    assert np.allclose(np.hypot(mids[:, 0], mids[:, 1]), 10.0, atol=1e-12)         # le milieu d'un côté touche le cercle
    # chaque côté est tangent : la distance du centre à la droite du côté vaut R, et au plus près
    for v0, v1 in zip(V, np.roll(V, -1, axis=0)):
        e = v1 - v0
        assert abs(e[0] * v0[1] - e[1] * v0[0]) / np.linalg.norm(e) == pytest.approx(10.0, abs=1e-12)
    # le cercle est dans le polygone et ne touche le contour qu'en n points
    t = np.linspace(0, 2 * math.pi, 3600, endpoint=False)
    C = 10.0 * np.stack([np.cos(t), np.sin(t)], axis=1)
    assert s.contains(C[:, 0], C[:, 1]).all()
    # et chaque côté ne touche le cercle qu'en un seul point, son milieu (distance minimale du centre au côté = R)
    u = np.linspace(0.0, 1.0, 2001)
    for v0, v1 in zip(V, np.roll(V, -1, axis=0)):
        pts = v0[None, :] + u[:, None] * (v1 - v0)[None, :]
        rho = np.hypot(pts[:, 0], pts[:, 1])
        assert rho.min() == pytest.approx(10.0, abs=1e-12) and int(np.argmin(rho)) == 1000
        assert (rho > 10.0 + 1e-9).sum() >= 1998


@pytest.mark.parametrize("n", (3, 4, 5, 6, 8, 12))
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("frac", (0.05, 0.3, 0.7, 0.999))
def test_rounded_area_and_perimeter_match_closed_forms_and_shoelace(n, mode, frac):
    s = make(n, mode, frac)
    a, r, h = s.apothem, s.corner_radius, math.pi / n
    area = n * a * a * math.tan(h) - n * r * r * (math.tan(h) - h)
    perimeter = 2 * n * (a - r) * math.tan(h) + 2 * math.pi * r
    assert s.area() == pytest.approx(area, rel=1e-13) and s.perimeter() == pytest.approx(perimeter, rel=1e-13)
    P = fine(s, 1e-7)                                  # polygone inscrit dans les arcs : borne par défaut
    assert shoelace(P) == pytest.approx(area, rel=1e-6) and shoelace(P) <= area * (1 + 1e-12)
    assert polyline_length(P) == pytest.approx(perimeter, rel=1e-6)


@pytest.mark.parametrize("n", SIDES)
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("frac", FRACS)
def test_outline_is_convex_counterclockwise_and_every_point_is_exactly_on_the_boundary(n, mode, frac):
    s = make(n, mode, frac, rot=23.0)
    P = s.outline_points(tol=1e-4)
    e = np.roll(P, -1, axis=0) - P
    cross = e[:, 0] * np.roll(e, -1, axis=0)[:, 1] - e[:, 1] * np.roll(e, -1, axis=0)[:, 0]
    assert (cross >= -1e-12).all(), "contour non convexe ou horaire"
    assert shoelace(P) > 0
    rho = np.hypot(P[:, 0], P[:, 1])
    assert np.allclose(rho, s.radial_distance(np.arctan2(P[:, 1], P[:, 0])), rtol=1e-12, atol=1e-12)
    d = np.linalg.norm(np.roll(P, -1, axis=0) - P, axis=1)
    assert d.min() > 0, "deux sommets consécutifs confondus"


@pytest.mark.parametrize("n", (3, 4, 5, 6, 7, 8))
@pytest.mark.parametrize("mode", MODES)
def test_rounded_corners_are_tangent_arcs_of_radius_r(n, mode):
    """Les extrémités d'arc sont sur les côtés du polygone d'origine (tangence), à la distance r tan(pi/n) du
    sommet, et tous les points d'un arc sont à r du centre du congé."""
    s = make(n, mode, 0.4)
    r, h = s.corner_radius, math.pi / n
    V = sharp_vertices(s)
    k = s.arc_segments(tol=1e-3)
    P = s.outline_points(tol=1e-3).reshape(n, k + 1, 2)
    centers = (s.fillet_center_distance) * np.stack(
        [np.cos(s.phase + 2 * math.pi * np.arange(n) / n), np.sin(s.phase + 2 * math.pi * np.arange(n) / n)], axis=1)
    assert np.allclose(np.linalg.norm(P - centers[:, None, :], axis=2), r, atol=1e-12)
    for i in range(n):
        for end, w in ((P[i, 0], V[i - 1]), (P[i, -1], V[(i + 1) % n])):
            assert np.linalg.norm(end - V[i]) == pytest.approx(r * math.tan(h), rel=1e-10)
            # l'extrémité est sur le côté d'origine [V_i, w]
            e = w - V[i]
            t = np.dot(end - V[i], e) / np.dot(e, e)
            assert 0 <= t <= 1 and np.linalg.norm(V[i] + t * e - end) < 1e-10
    assert s.extreme_radius == pytest.approx((s.apothem - r) / math.cos(h) + r, rel=1e-13)
    assert np.hypot(*fine(s, 1e-7).T).max() == pytest.approx(s.extreme_radius, abs=2e-7)


@pytest.mark.parametrize("n", (3, 4, 6, 9))
@pytest.mark.parametrize("mode", MODES)
def test_corner_radius_equal_to_apothem_gives_a_circle(n, mode):
    s = make(n, mode, 1.0)
    P = fine(s, 1e-6)
    a = s.apothem
    assert np.allclose(np.hypot(P[:, 0], P[:, 1]), a, atol=1e-9)
    assert shoelace(P) == pytest.approx(math.pi * a * a, rel=1e-6)
    assert s.area() == pytest.approx(math.pi * a * a, rel=1e-12)
    assert len(P) == len(np.unique(np.round(P, 12), axis=0)), "points confondus quand les côtés se réduisent à un point"
    t = np.linspace(0, 2 * math.pi, 997)
    Q = np.stack([np.cos(t), np.sin(t)], axis=1)
    assert s.contains(0.999 * a * Q[:, 0], 0.999 * a * Q[:, 1]).all()
    assert not s.contains(1.001 * a * Q[:, 0], 1.001 * a * Q[:, 1]).any()


@pytest.mark.parametrize("n", (3, 5, 6, 8, 50))
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("frac", (0.0, 0.3, 1.0))
@pytest.mark.parametrize("rot", (0.0, 17.0, -40.0))
def test_contains_matches_a_polygon_test_on_the_fine_outline(n, mode, frac, rot):
    s = make(n, mode, frac, rot)
    tol = 1e-4
    P = fine(s, tol)
    rng = np.random.default_rng(n)
    E = 1.2 * s.half_extent()
    Q = rng.uniform(-E, E, size=(3000, 2))
    ref = inside_convex(Q, P)
    got = s.contains(Q[:, 0], Q[:, 1])
    far = dist_to_polyline(Q, P) > 2 * tol               # hors de la bande où la discrétisation (flèche <= tol) compte
    assert far.sum() > 2000 and (ref[far] == got[far]).all()
    assert got.sum() > 100 and (~got).sum() > 100


@pytest.mark.parametrize("n", (3, 4, 6, 7))
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("frac", (0.0, 0.5))
def test_contains_includes_the_boundary_and_excludes_just_outside(n, mode, frac):
    s = make(n, mode, frac, rot=11.0)
    P = s.outline_points(tol=1e-4)
    assert s.contains(P[:, 0], P[:, 1]).all()                    # sommets exacts et points de tangence : dedans
    mids = 0.5 * (P + np.roll(P, -1, axis=0))                    # milieux de segments : à l'intérieur (convexe)
    assert s.contains(mids[:, 0], mids[:, 1]).all()
    assert s.contains(0.0, 0.0)
    rho = np.hypot(P[:, 0], P[:, 1])
    out = P * (1 + 1e-6)                  # contour convexe contenant le centre : tout agrandissement radial sort
    assert not s.contains(out[:, 0], out[:, 1]).any()
    far = P * 1.01
    assert not s.contains(far[:, 0], far[:, 1]).any()
    assert rho.max() <= s.extreme_radius + 1e-9


@pytest.mark.parametrize("n", (3, 4, 5, 6, 8, 12))
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("frac", (0.0, 0.2, 0.9))
def test_project_returns_the_nearest_boundary_point(n, mode, frac):
    s = make(n, mode, frac, rot=7.0)
    tol = 1e-5
    P = fine(s, tol)
    rng = np.random.default_rng(100 + n)
    E = s.half_extent()
    Q = rng.uniform(-E, E, size=(4000, 2))
    Q = Q[s.contains(Q[:, 0], Q[:, 1])]
    px, py = s.project(Q[:, 0], Q[:, 1])
    # 1. le point projeté est exactement sur le contour
    assert np.allclose(np.hypot(px, py), s.radial_distance(np.arctan2(py, px)), rtol=1e-11, atol=1e-11)
    # 2. c'est le plus proche : la distance au point projeté égale la distance à la polyligne fine, à la flèche
    #    près (la polyligne est inscrite dans les arcs, donc à au plus `tol` du vrai contour)
    d_true = np.hypot(px - Q[:, 0], py - Q[:, 1])
    d_poly = dist_to_polyline(Q, P)
    assert len(Q) > 500
    assert (d_true <= d_poly + tol + 1e-9).all()
    assert (d_poly <= d_true + tol + 1e-9).all()


@pytest.mark.parametrize("n", (3, 6, 11))
def test_project_keeps_the_position_along_a_side_for_staircase_like_points(n):
    """Raison d'être de la projection perpendiculaire : des nœuds en escalier le long d'un côté oblique (qui
    avancent le long du côté tout en s'enfonçant ou en remontant vers le centre) doivent se projeter dans l'ordre,
    à l'abscisse exacte qu'ils ont le long du côté. La projection radiale, elle, peut les faire reculer ou
    les confondre."""
    s = RegularPolygon(n, "circumscribed", 10.0, 0.0, 13.0)
    V = sharp_vertices(s)
    v0, v1 = V[0], V[1]
    length = np.linalg.norm(v1 - v0)
    e = (v1 - v0) / length
    nrm = np.array([e[1], -e[0]])                                   # normale sortante (polygone antihoraire)
    t = np.linspace(0.3, 0.7, 200) * length
    depth = 0.5 * (1 + np.sign(np.sin(np.arange(200) * 0.9))) * 1.0 + 0.05    # marches : 0,05 puis 1,05 du bord
    Q = v0[None, :] + t[:, None] * e[None, :] - depth[:, None] * nrm[None, :]
    px, py = s.project(Q[:, 0], Q[:, 1])
    along = (px - v0[0]) * e[0] + (py - v0[1]) * e[1]
    assert np.allclose(along, t, atol=1e-10) and (np.diff(along) > 0).all()
    assert np.allclose((px - v0[0]) * nrm[0] + (py - v0[1]) * nrm[1], 0.0, atol=1e-10)     # et sur le côté


@pytest.mark.parametrize("r", (0.5, 2.0, 7.0))
@pytest.mark.parametrize("tol", (0.05, 0.01, 0.002, 1e-4))
def test_arc_segments_respect_the_sagitta_tolerance_and_are_minimal(r, tol):
    n = 6
    s = RegularPolygon(n, "circumscribed", 50.0, r, 0.0)
    k = s.arc_segments(tol=tol)
    h = math.pi / n

    def sagitta(kk):
        return r * (1 - math.cos(h / kk))

    assert sagitta(k) <= tol * (1 + 1e-12)
    assert k == 1 or sagitta(k - 1) > tol
    # vérification géométrique directe : écart maximal entre un arc et sa corde
    P = s.outline_points(tol=tol).reshape(n, k + 1, 2)
    c = s.fillet_center_distance
    center0 = c * np.array([math.cos(s.phase), math.sin(s.phase)])
    a0, a1 = P[0, 0] - center0, P[0, 1] - center0
    assert np.linalg.norm(a0) == pytest.approx(r, rel=1e-12) and np.linalg.norm(a1) == pytest.approx(r, rel=1e-12)
    assert r - np.linalg.norm(0.5 * (a0 + a1)) <= tol * (1 + 1e-9)       # écart corde-arc mesuré, pas calculé
    assert s.outline_points(tol=tol).shape[0] == n * (k + 1)


def test_arc_segments_special_cases():
    s = RegularPolygon(6, "inscribed", 10.0, 1.0)
    assert RegularPolygon(6, "inscribed", 10.0).arc_segments(tol=0.01) == 0
    assert s.arc_segments(tol=5.0) == 1                      # tolérance supérieure au rayon : une corde
    assert s.arc_segments(segments=7, tol=1e-9) == 7         # le nombre imposé l'emporte
    assert len(s.outline_points(segments=3)) == 6 * 4
    for bad in (0, -2, 2.5, True):
        with pytest.raises(ShapeError):
            s.arc_segments(segments=bad)
    for bad in (0, -1.0, float("nan"), float("inf"), None):
        with pytest.raises(ShapeError):
            s.arc_segments(tol=bad)
    with pytest.raises(ShapeError, match="maximum"):
        RegularPolygon(6, "inscribed", 10.0, 4.0).arc_segments(tol=1e-300)


@pytest.mark.parametrize("kw", [
    dict(n=2), dict(n=361), dict(n=3.5), dict(n=True), dict(n=0), dict(n=-6), dict(n="6"),
    dict(circle_mode="tangent"), dict(circle_mode=None),
    dict(ref_radius=0), dict(ref_radius=-1.0), dict(ref_radius=float("nan")), dict(ref_radius=float("inf")),
    dict(corner_radius=-0.1), dict(corner_radius=float("nan")), dict(corner_radius=float("inf")),
    dict(corner_radius=20.0), dict(rotation_deg=float("nan")), dict(rotation_deg=float("inf")),
])
def test_invalid_parameters_are_rejected(kw):
    base = dict(n=6, circle_mode="inscribed", ref_radius=10.0)
    base.update(kw)
    with pytest.raises(ShapeError):
        RegularPolygon(**base)


def test_valid_limits_are_accepted():
    for n in (regular_polygon.MIN_SIDES, regular_polygon.MAX_SIDES):
        for mode in MODES:
            s = RegularPolygon(n, mode, 10.0, 0.0)
            assert s.area() > 0
    s = RegularPolygon(6, "inscribed", 10.0)
    RegularPolygon(6, "inscribed", 10.0, s.apothem)               # maximum exact
    RegularPolygon(6, "inscribed", 10.0, s.apothem * (1 + 1e-12))  # bruit flottant toléré


@pytest.mark.parametrize("n", (3, 4, 5, 6, 7, 8, 9, 12))
def test_default_orientation_is_upright_and_rotation_is_counterclockwise(n):
    s = RegularPolygon(n, "circumscribed", 10.0)
    V = sharp_vertices(s)
    if n % 2:                                       # impair : un sommet exactement au nord
        top = V[np.argmax(V[:, 1])]
        assert top[0] == pytest.approx(0.0, abs=1e-12) and top[1] == pytest.approx(s.circumradius)
    else:                                           # pair : un côté horizontal au nord (deux sommets à la même hauteur)
        ys = np.sort(V[:, 1])[::-1]
        assert ys[0] == pytest.approx(ys[1], abs=1e-12) and ys[0] == pytest.approx(s.apothem)
        assert abs(V[np.argmax(V[:, 1])][0]) <= s.side / 2 + 1e-12
    # rotation de 360/n : mêmes sommets ; de 90 degrés : le polygone tourne dans le sens antihoraire
    R = RegularPolygon(n, "circumscribed", 10.0, 0.0, 360.0 / n)
    assert np.allclose(np.sort(np.round(sharp_vertices(R), 9), axis=0), np.sort(np.round(V, 9), axis=0), atol=1e-8)
    q = RegularPolygon(n, "circumscribed", 10.0, 0.0, 90.0)
    c, sn = 0.0, 1.0
    rot = np.stack([c * V[:, 0] - sn * V[:, 1], sn * V[:, 0] + c * V[:, 1]], axis=1)
    assert np.allclose(sharp_vertices(q), rot, atol=1e-12)


@pytest.mark.parametrize("n", SIDES)
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("frac", (0.0, 0.4, 1.0))
@pytest.mark.parametrize("rot", (0.0, 29.0))
def test_bbox_is_exact(n, mode, frac, rot):
    s = make(n, mode, frac, rot)
    tol = 1e-7
    P = fine(s, tol)
    x0, x1, y0, y1 = s.bbox()
    assert (x0, y0) == pytest.approx((P[:, 0].min(), P[:, 1].min()), abs=2 * tol)
    assert (x1, y1) == pytest.approx((P[:, 0].max(), P[:, 1].max()), abs=2 * tol)
    assert s.half_extent() == max(abs(x0), abs(x1), abs(y0), abs(y1))


def test_scaled_changes_unit_only():
    s = make(7, "circumscribed", 0.3, 5.0, ref=2.0)
    k = 3.5
    t = s.scaled(k)
    assert t.area() == pytest.approx(s.area() * k * k, rel=1e-13)
    assert t.perimeter() == pytest.approx(s.perimeter() * k, rel=1e-13)
    assert t.extreme_radius == pytest.approx(s.extreme_radius * k, rel=1e-13)
    assert np.allclose(np.array(t.bbox()), k * np.array(s.bbox()), rtol=1e-13)
    assert t.arc_segments(tol=0.01 * k) == s.arc_segments(tol=0.01)


def test_describe_mentions_the_main_figures():
    txt = RegularPolygon(6, "circumscribed", 30.0, 2.0).describe("mm", 2)
    for part in ("hexagone", "circonscrit", "30.00 mm", "apothème", "arrondis", "2.00 mm"):
        assert part in txt


# ================================================================ 2. maillage avec contour exact

SCALE, BASE = 0.1, 2.0            # mm imprimés par unité de grille ; épaisseur du socle


def shape_mesh(n, mode, frac, rot, *, ref=30.0, vexag=1.5, relief=True, seed=0, tol=0.005, **kw):
    """Maillage d'un polygone sur une grille synthétique au pas de 1 (unité de terrain) ; renvoie aussi la forme
    en mm et le contour discrétisé."""
    spec = make(n, mode, frac, rot, ref=ref)
    half = int(math.ceil(spec.half_extent())) + 3
    xs = np.arange(-half, half + 1, dtype=float)
    ys = xs[::-1].copy()
    X, Y = np.meshgrid(xs, ys)
    clean, quad_ok = mesh.clean_mask(spec.contains(X, Y))
    if relief:
        elev = 500.0 + ndimage.gaussian_filter(np.random.default_rng(seed).normal(size=X.shape), 3) * 400
    else:
        elev = np.full(X.shape, 500.0)
    outline = spec.scaled(SCALE).discretize(tol=tol)
    tris, nt, nb, nw = mesh.build_mesh_from_arrays(elev, quad_ok, xs, ys, SCALE, vexag, BASE, elev[clean].min(),
                                                   outline=outline, **kw)
    return tris, spec.scaled(SCALE), outline


def assert_valid_solid(tris):
    rep = check_stl.edge_report(tris)
    assert rep["boundary_edges"] == 0 and rep["nonmanifold_edges"] == 0 and rep["duplicate_directed_edges"] == 0, rep
    assert check_stl.signed_volume(tris) > 0, "faces tournées vers l'intérieur"
    assert check_stl.degenerate_count(tris) == 0, "triangles dégénérés"
    assert tris[:, :, 2].min() == pytest.approx(-BASE)


def wall_vertices(tris):
    nrm = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    norm = np.linalg.norm(nrm, axis=1)
    vertical = np.abs(nrm[:, 2]) <= 1e-9 * np.maximum(norm, 1e-300)
    return np.unique(tris[vertical].reshape(-1, 3), axis=0)


@pytest.mark.parametrize("n", (3, 4, 6, 7))
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("frac", (0.0, 0.5, 1.0))
@pytest.mark.parametrize("rot", (0.0, 37.0))
def test_mesh_is_closed_oriented_and_has_no_degenerate_triangle(n, mode, frac, rot):
    tris, _, _ = shape_mesh(n, mode, frac, rot, seed=n)
    assert_valid_solid(tris)


@pytest.mark.parametrize("kw", [
    {"bottom": "merged"}, {"bottom": "merged-rows"}, {"diagonal": "adaptive"},
    {"top": "adaptive", "top_tol_mm": 0.05},
])
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("n,frac,rot", [(5, 0.0, 0.0), (6, 0.5, 37.0), (3, 0.3, 11.0), (8, 1.0, 0.0)])
def test_every_triangulation_mode_stays_a_valid_solid_with_exact_outline(kw, mode, n, frac, rot):
    tris, _, _ = shape_mesh(n, mode, frac, rot, seed=3, **kw)
    assert_valid_solid(tris)


@pytest.mark.parametrize("n", (3, 4, 5, 6, 8, 12, 40))
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("frac", (0.0, 0.2, 0.9))
@pytest.mark.parametrize("rot", (0.0, 13.0))
def test_flat_relief_volume_is_base_times_exact_area_within_the_sagitta_bound(n, mode, frac, rot):
    """Relief plat : volume = socle x aire de l'empreinte. L'empreinte du maillage est le polygone inscrit dans
    les arcs exacts ; elle ne peut manquer que les segments de cercle (au plus 2 pi r eps) plus la tolérance de
    fusion des points (1e-4 mm le long du périmètre)."""
    eps = 0.005
    tris, shp, _ = shape_mesh(n, mode, frac, rot, relief=False, tol=eps)
    assert_valid_solid(tris)
    vol = check_stl.signed_volume(tris)
    area_mesh = vol / BASE
    miss = shp.area() - area_mesh
    slack = shp.perimeter() * 2e-4
    assert -slack <= miss <= 2 * math.pi * shp.corner_radius * eps + slack, (miss, slack)


@pytest.mark.parametrize("n", (3, 4, 5, 6, 9))
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("frac", (0.0, 0.35))
def test_wall_is_vertical_and_lies_exactly_on_the_outline(n, mode, frac):
    tris, shp, outline = shape_mesh(n, mode, frac, 21.0, seed=5)
    W = wall_vertices(tris)
    assert len(W) > 0
    rho = np.hypot(W[:, 0], W[:, 1])
    exact = shp.radial_distance(np.arctan2(W[:, 1], W[:, 0]))
    assert np.allclose(rho, exact, rtol=0, atol=2e-6)                    # sommets de paroi exactement sur le contour
    # chaque sommet exact du contour (angles vifs, extrémités de tangence, points d'arc) est sur la paroi, à la
    # tolérance de fusion de mesh.build_outline_rim (1e-4 mm) près : soit c'est un sommet de la paroi, soit la
    # corde entre deux sommets voisins, eux-mêmes sur l'arc, passe à moins de 1e-4 mm de lui
    top = W[W[:, 2] > -BASE + 1e-9]
    chain = top[np.argsort(np.arctan2(top[:, 1], top[:, 0])), :2]          # contour étoilé : tri par angle polaire
    assert dist_to_polyline(outline.points, chain).max() <= 1e-4 + 1e-9
    # parois verticales : tout sommet haut a son sommet bas à la même position horizontale
    low = W[np.isclose(W[:, 2], -BASE)]
    high = W[~np.isclose(W[:, 2], -BASE)]
    assert len(low) > 0 and len(high) > 0
    for p in high[:200]:
        assert np.isclose(np.hypot(low[:, 0] - p[0], low[:, 1] - p[1]), 0.0, atol=1e-6).any()


@pytest.mark.parametrize("n", (3, 4, 6))
@pytest.mark.parametrize("mode", MODES)
def test_sharp_corners_are_present_in_the_wall(n, mode):
    tris, shp, _ = shape_mesh(n, mode, 0.0, 9.0, seed=2)
    W = wall_vertices(tris)
    V = sharp_vertices(shp)
    d = np.linalg.norm(V[:, None, :] - W[None, :, :2], axis=2).min(axis=1)
    assert (d <= 1e-4 + 1e-6).all(), d


def test_shape_with_a_hole_is_refused_with_a_clear_message():
    spec = make(6, "inscribed", 0.0, 0.0, ref=30.0)
    half = int(math.ceil(spec.half_extent())) + 3
    xs = np.arange(-half, half + 1, dtype=float)
    ys = xs[::-1].copy()
    X, Y = np.meshgrid(xs, ys)
    m = spec.contains(X, Y) & ~((np.abs(X - 3) < 3) & (np.abs(Y) < 3))
    clean, quad_ok = mesh.clean_mask(m)
    elev = np.full(X.shape, 10.0)
    with pytest.raises(SystemExit, match="un seul contour"):
        mesh.build_mesh_from_arrays(elev, quad_ok, xs, ys, SCALE, 1.0, BASE, 10.0,
                                    outline=spec.scaled(SCALE).discretize(tol=0.005))


# ================================================================ 3. ligne de commande

@pytest.fixture
def dem(tmp_path):
    """MNT synthétique 2 km x 2 km (pixel 5 m, CRS suisse), relief doux."""
    yy, xx = np.mgrid[0:400, 0:400]
    elev = 1000.0 + 0.08 * xx + 6.0 * np.sin(yy / 17.0) + 4.0 * np.cos(xx / 23.0)
    write_synthetic_tif(tmp_path / "dem.tif", elev, 2600000.0, 1202000.0, 5.0, 2056)
    return str(tmp_path / "dem.tif")


COMMON = ["--cx", "2601000", "--cy", "1201000", "--radius", "300", "--diameter-mm", "60", "--vexag", "1",
          "--base-mm", "2", "--print-spacing-mm", "0.5"]


def cli(monkeypatch, tmp_path, dem, *extra, output="p.stl"):
    argv = ["x", "--input", dem, *COMMON, *extra]
    if output:
        argv += ["--output", output]
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", argv)
    circle_dem_to_stl.main()
    return tmp_path / "output" / output if output else None


def fails(monkeypatch, tmp_path, dem, *extra):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--input", dem, *COMMON, "--dry-run", *extra])
    with pytest.raises(SystemExit) as e:
        circle_dem_to_stl.main()
    return str(e.value)


@pytest.mark.parametrize("n,mode,corner,rot", [
    (6, "circumscribed", 0, 0), (6, "circumscribed", 2, 0), (6, "inscribed", 2, 15),
    (3, "inscribed", 0, 0), (3, "circumscribed", 3, 0), (4, "inscribed", 1, 45), (5, "circumscribed", 0, 0),
    (12, "inscribed", 0, 0),
])
def test_cli_produces_a_watertight_oriented_stl_with_the_exact_printed_size(monkeypatch, tmp_path, dem, n, mode,
                                                                          corner, rot):
    args = ["--sides", str(n), "--circle-mode", mode, "--corner-radius", str(corner), "--rotation", str(rot)]
    out = cli(monkeypatch, tmp_path, dem, *args)
    tris = check_stl.read_stl(str(out))
    assert_valid_solid(tris)
    # dimensions imprimées : celles du polygone exact en mm (échelle 0,1 mm/m), au bruit float32 près
    ref_mm = 30.0
    shp = RegularPolygon(n, mode, ref_mm, min(corner, RegularPolygon(n, mode, ref_mm).apothem), rot)
    x0, x1, y0, y1 = shp.bbox()
    v = tris.reshape(-1, 3)
    # angles vifs : sommets exacts, donc dimensions exactes. Arrondis : les sommets du maillage sont sur l'arc mais
    # le point extrême de l'arc peut tomber entre deux sommets : le maillage est dedans, d'au plus la flèche (0,005 mm)
    slack = 0.005 if corner else 2e-4
    for ax, lo, hi in ((0, x0, x1), (1, y0, y1)):
        assert lo - 2e-4 <= v[:, ax].min() <= lo + slack
        assert hi - slack <= v[:, ax].max() <= hi + 2e-4
    assert v[:, 2].min() == pytest.approx(-2.0)                     # plancher du socle
    # empreinte : section horizontale du solide = aire exacte, à la flèche et à la tolérance de fusion près
    assert check_stl.signed_volume(tris) > shp.area() * 2.0 * (1 - 1e-3)


def test_cli_dry_run_prints_the_real_printed_dimensions(monkeypatch, tmp_path, dem, capsys):
    cli(monkeypatch, tmp_path, dem, "--sides", "6", "--circle-mode", "circumscribed", "--corner-radius", "2",
        "--dry-run", output=None)
    out = capsys.readouterr().out
    shp = RegularPolygon(6, "circumscribed", 30.0, 2.0)
    x0, x1, y0, y1 = shp.bbox()
    assert f"pièce : {x1 - x0:.2f} x {y1 - y0:.2f} mm" in out
    assert "68.66 x 60.00 mm" in out
    assert "hexagone" in out and "circonscrit" in out and "segment(s) par arc" in out
    assert "6cotes-circonscrit_r2mm" in out
    assert not (tmp_path / "output").exists() or not list((tmp_path / "output").glob("*.stl"))


def test_cli_inscribed_piece_fits_in_the_reference_circle_and_circumscribed_exceeds_it(monkeypatch, tmp_path, dem):
    ins_t = check_stl.read_stl(str(cli(monkeypatch, tmp_path, dem, "--sides", "7", "--circle-mode", "inscribed",
                                       output="i.stl")))
    cir_t = check_stl.read_stl(str(cli(monkeypatch, tmp_path, dem, "--sides", "7", "--circle-mode", "circumscribed",
                                       output="c.stl")))
    ins, cir = wall_vertices(ins_t), wall_vertices(cir_t)                 # sommets de la paroi = le contour imprimé
    assert np.hypot(ins[:, 0], ins[:, 1]).max() == pytest.approx(30.0, abs=2e-4)     # sommets sur le cercle de Ø 60
    assert np.hypot(ins[:, 0], ins[:, 1]).min() == pytest.approx(30.0 * math.cos(math.pi / 7), abs=2e-4)
    # circonscrit : côtés tangents au cercle de Ø 60 (plus petite distance centre-paroi = 30), sommets au-delà
    assert np.hypot(cir[:, 0], cir[:, 1]).min() >= 30.0 - 2e-4
    assert np.hypot(cir[:, 0], cir[:, 1]).max() == pytest.approx(30.0 / math.cos(math.pi / 7), abs=2e-4)


def test_cli_corner_radius_equal_to_apothem_matches_the_smooth_disc(monkeypatch, tmp_path, dem):
    """r = apothème donne un cercle de rayon a : avec un hexagone circonscrit, a = rayon de référence, donc le
    même disque que le chemin historique --smooth-wall. Mêmes dimensions ; volumes très proches."""
    poly = check_stl.read_stl(str(cli(monkeypatch, tmp_path, dem, "--sides", "6", "--circle-mode", "circumscribed",
                                      "--corner-radius", "30", output="a.stl")))
    disc = check_stl.read_stl(str(cli(monkeypatch, tmp_path, dem, "--smooth-wall", output="b.stl")))
    for t in (poly, disc):
        assert_valid_solid(t)
    pv, dv = poly.reshape(-1, 3), disc.reshape(-1, 3)
    for ax in (0, 1):
        assert pv[:, ax].min() == pytest.approx(dv[:, ax].min(), abs=2e-3)
        assert pv[:, ax].max() == pytest.approx(dv[:, ax].max(), abs=2e-3)
    assert check_stl.signed_volume(poly) == pytest.approx(check_stl.signed_volume(disc), rel=5e-3)


def test_cli_relief_follows_the_dem_inside_the_shape(monkeypatch, tmp_path, dem):
    """Altitude du dessus (mm) = (z - zmin) x échelle (0,1 mm/m, exagération 1) : l'écart avec le relief analytique
    du MNT synthétique est constant (= -zmin x échelle) pour tous les sommets intérieurs de la pièce."""
    tris = check_stl.read_stl(str(cli(monkeypatch, tmp_path, dem, "--sides", "8", "--circle-mode", "inscribed",
                                      "--corner-radius", "1")))
    v = np.unique(tris.reshape(-1, 3), axis=0)
    shp = RegularPolygon(8, "inscribed", 30.0, 1.0)
    inner = v[(np.hypot(v[:, 0], v[:, 1]) < 0.9 * shp.apothem) & (v[:, 2] > -2.0 + 1e-6)]
    assert len(inner) > 1000
    X, Y = 2601000.0 + inner[:, 0] / 0.1, 1201000.0 + inner[:, 1] / 0.1          # mm imprimés -> mètres de terrain
    xx, yy = (X - 2600000.0) / 5.0 - 0.5, (1202000.0 - Y) / 5.0 - 0.5            # indices de pixel (centres)
    analytic = (1000.0 + 0.08 * xx + 6.0 * np.sin(yy / 17.0) + 4.0 * np.cos(xx / 23.0)) * 0.1
    diff = inner[:, 2] - analytic
    assert diff.max() - diff.min() < 0.01                                          # mm : relief et pente respectés
    assert inner[:, 2].max() - inner[:, 2].min() > 1.0                             # et le relief est bien là


@pytest.mark.parametrize("extra,msg", [
    (["--circle-mode", "inscribed"], "requiert --sides"),
    (["--corner-radius", "1"], "requiert --sides"),
    (["--rotation", "10"], "requiert --sides"),
    (["--corner-tol-mm", "0.01"], "requiert --sides"),
    (["--corner-segments", "4"], "requiert --sides"),
    (["--sides", "6"], "--circle-mode"),
    (["--sides", "2", "--circle-mode", "inscribed"], "entre 3 et 360"),
    (["--sides", "361", "--circle-mode", "inscribed"], "entre 3 et 360"),
    (["--sides", "6", "--circle-mode", "inscribed", "--smooth-wall"], "incompatible avec --smooth-wall"),
    (["--sides", "6", "--circle-mode", "inscribed", "--polygon", "x.geojson"], "incompatible avec --polygon"),
    (["--sides", "6", "--circle-mode", "inscribed", "--place", "Paris", "--clip-to-place"], "incompatible"),
    (["--sides", "6", "--circle-mode", "inscribed", "--corner-radius", "-1"], "positif ou nul"),
    (["--sides", "6", "--circle-mode", "inscribed", "--corner-radius", "nan"], "positif ou nul"),
    (["--sides", "6", "--circle-mode", "inscribed", "--corner-radius", "26"], "dépasse l'apothème"),
    (["--sides", "6", "--circle-mode", "inscribed", "--rotation", "inf"], "--rotation"),
    (["--sides", "6", "--circle-mode", "inscribed", "--corner-radius", "1", "--corner-tol-mm", "0"], "--corner-tol-mm"),
    (["--sides", "6", "--circle-mode", "inscribed", "--corner-radius", "1", "--corner-segments", "0"],
     "--corner-segments"),
    (["--sides", "6", "--circle-mode", "inscribed", "--corner-radius", "1", "--corner-tol-mm", "1e-12"],
     "--corner-tol-mm"),
])
def test_cli_rejects_invalid_shape_options(monkeypatch, tmp_path, dem, extra, msg):
    assert msg in fails(monkeypatch, tmp_path, dem, *extra)


@pytest.fixture
def half_dem(tmp_path):
    """MNT sans donnée à droite de x = centre - 250 m. Le ré-échantillonnage comble les trous jusqu'à 50 pixels
    (250 m) des données valides : il reste donc sans donnée toute la moitié droite de la forme (x > centre)."""
    yy, xx = np.mgrid[0:400, 0:400]
    elev = 1000.0 + 0.08 * xx + 6.0 * np.sin(yy / 17.0)
    elev[:, 150:] = -9999.0
    write_synthetic_tif(tmp_path / "half.tif", elev, 2600000.0, 1202000.0, 5.0, 2056)
    return str(tmp_path / "half.tif")


@pytest.mark.parametrize("extra", [["--sides", "6", "--circle-mode", "inscribed"],
                                   ["--sides", "6", "--circle-mode", "circumscribed", "--corner-radius", "2"],
                                   ["--sides", "3", "--circle-mode", "circumscribed"]])
def test_cli_refuses_a_dem_that_does_not_cover_the_shape(monkeypatch, tmp_path, half_dem, extra):
    """Sans cette garde, le bord exact projetterait les nœuds du masque (arrêté au centre) sur un contour situé à
    plusieurs dizaines de pas : maillage fermé en apparence, mais étiré et faux."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--input", half_dem, *COMMON, *extra, "--output", "h.stl"])
    with pytest.raises(SystemExit) as e:
        circle_dem_to_stl.main()
    msg = str(e.value)
    assert "n'ont pas de donnée valide" in msg and "--sea-level" in msg and "mm imprimés de son bord" in msg
    assert not list((tmp_path / "output").glob("*.stl")) if (tmp_path / "output").exists() else True


def test_cli_sea_level_lets_a_partial_dem_fill_the_whole_shape(monkeypatch, tmp_path, half_dem):
    out = cli(monkeypatch, tmp_path, half_dem, "--sides", "6", "--circle-mode", "circumscribed", "--corner-radius",
              "2", "--sea-level", "900")
    tris = check_stl.read_stl(str(out))
    assert_valid_solid(tris)
    shp = RegularPolygon(6, "circumscribed", 30.0, 2.0)
    assert check_stl.signed_volume(tris) > shp.area() * 2.0 * 0.999                # toute la forme est remplie


@pytest.mark.parametrize("mode,frac", [("inscribed", 0.0), ("circumscribed", 0.2)])
def test_shape_coverage_guard_tolerates_missing_cells_near_the_edge_only(mode, frac):
    shp = make(5, mode, frac, 10.0, ref=100.0)
    xs = np.arange(-130.0, 130.1, 5.0)
    ys = xs[::-1].copy()
    gx, gy = np.meshgrid(xs, ys)
    inside = shp.contains(gx, gy)
    px, py = shp.project(gx[inside], gy[inside])
    depth = np.zeros(gx.shape)
    depth[inside] = np.hypot(px - gx[inside], py - gy[inside])
    # cellules manquantes jusqu'à 12 m du bord (moins de 3 pas de 5 m) : accepté
    circle_dem_to_stl.check_shape_coverage(shp, inside, inside & (depth > 12.0), xs, ys, 15.0, 0.1)
    # une cellule manquante à plus de 15 m du bord : refusé, avec sa profondeur en mm imprimés
    kept = inside & (depth > 12.0)
    ii, jj = np.unravel_index(np.argmax(depth), depth.shape)
    kept[ii, jj] = False
    with pytest.raises(SystemExit, match=r"1 point\(s\) de la grille.*mm imprimés de son bord"):
        circle_dem_to_stl.check_shape_coverage(shp, inside, kept, xs, ys, 15.0, 0.1)


def test_cli_corner_radius_at_the_exact_apothem_is_accepted(monkeypatch, tmp_path, dem, capsys):
    apothem_mm = 30.0 * math.cos(math.pi / 6)
    cli(monkeypatch, tmp_path, dem, "--sides", "6", "--circle-mode", "inscribed", "--corner-radius",
        repr(apothem_mm), "--dry-run", output=None)
    assert "pièce : 51.96 x 51.96 mm" in capsys.readouterr().out


def test_cli_warns_when_the_corner_is_below_the_print_spacing(monkeypatch, tmp_path, dem, capsys):
    cli(monkeypatch, tmp_path, dem, "--sides", "6", "--circle-mode", "inscribed", "--corner-radius", "0.2",
        "--dry-run", output=None)
    assert "inférieur au pas d'impression" in capsys.readouterr().out


def test_cli_segments_per_arc_can_be_forced(monkeypatch, tmp_path, dem, capsys):
    cli(monkeypatch, tmp_path, dem, "--sides", "6", "--circle-mode", "inscribed", "--corner-radius", "3",
        "--corner-segments", "5", "--dry-run", output=None)
    out = capsys.readouterr().out
    assert "5 segment(s) par arc" in out and "imposé par --corner-segments" in out
    tris = check_stl.read_stl(str(cli(monkeypatch, tmp_path, dem, "--sides", "6", "--circle-mode", "inscribed",
                                      "--corner-radius", "3", "--corner-segments", "5")))
    assert_valid_solid(tris)


def test_cli_shape_exports_the_same_solid_in_stl_obj_and_3mf(monkeypatch, tmp_path, dem):
    pytest.importorskip("lib3mf")
    import mesh_export
    cli(monkeypatch, tmp_path, dem, "--sides", "5", "--circle-mode", "inscribed", "--corner-radius", "1.5",
        "--format", "all", output="f")
    vols = []
    for ext in ("stl", "obj", "3mf"):
        t = mesh_export.read_triangles(str(tmp_path / "output" / f"f.{ext}"))
        assert_valid_solid(t)
        vols.append(check_stl.signed_volume(t))
    assert vols[1] == pytest.approx(vols[0], rel=1e-5) and vols[2] == pytest.approx(vols[0], rel=1e-5)


def test_check_stl_verifies_the_exact_footprint_of_a_shape(monkeypatch, tmp_path, dem, capsys):
    out = cli(monkeypatch, tmp_path, dem, "--sides", "6", "--circle-mode", "circumscribed", "--corner-radius", "2")
    tris = check_stl.read_stl(str(out))
    shp = RegularPolygon(6, "circumscribed", 30.0, 2.0)
    res = check_stl.check(tris, 60.0, 2.0, radius_m=300.0, shape_mm=shp, tol_mm=0.01)
    assert all(ok for _, ok, _ in res), res
    assert any("dimensions de la pièce" in name for name, _, _ in res)
    wrong = check_stl.check(tris, 60.0, 2.0, shape_mm=RegularPolygon(6, "inscribed", 30.0, 2.0), tol_mm=0.01)
    assert not wrong[0][1]                                        # un autre polygone est refusé
    assert not check_stl.check(tris, 60.0, 2.0, tol_mm=0.01)[0][1]      # et sans forme, on attend le diamètre : 68,66 != 60
    # la ligne de commande de check_stl.py
    monkeypatch.setattr(sys, "argv", ["check_stl.py", str(out), "--diameter-mm", "60", "--base-mm", "2", "--radius",
                                      "300", "--sides", "6", "--circle-mode", "circumscribed", "--corner-radius", "2"])
    with pytest.raises(SystemExit) as e:
        check_stl.main()
    assert e.value.code == 0 and "[OK ] dimensions de la pièce" in capsys.readouterr().out
    monkeypatch.setattr(sys, "argv", ["check_stl.py", str(out), "--diameter-mm", "60", "--base-mm", "2", "--sides", "6"])
    with pytest.raises(SystemExit, match="--circle-mode"):
        check_stl.main()


def test_default_stem_and_estimate_with_a_shape():
    class A:
        place, clip_to_place, polygon, lat, lon, cx, cy = None, False, None, None, None, 1.0, 2.0
        diameter_mm, vexag = 60.0, 1.0
        sides, circle_mode, corner_radius = 6, "circumscribed", 2.0
    assert circle_dem_to_stl.default_stem(A) == "x1_y2_60mm_6cotes-circonscrit_r2mm"
    A.circle_mode, A.corner_radius = "inscribed", None
    assert circle_dem_to_stl.default_stem(A) == "x1_y2_60mm_6cotes-inscrit"
    A.sides = None                                                  # sans forme : nom historique inchangé
    assert circle_dem_to_stl.default_stem(A) == "x1_y2_60mm"
    disc = circle_dem_to_stl.estimate_resources(60.0, 0.5)
    shp = circle_dem_to_stl.estimate_resources(60.0, 0.5, area_mm2=math.pi * 30.0 ** 2, size_mm=60.0)
    assert shp["triangles"] == pytest.approx(disc["triangles"], rel=1e-12)       # un disque : même estimation
    sq = circle_dem_to_stl.estimate_resources(60.0, 0.5, area_mm2=60.0 * 60.0, size_mm=60.0)
    assert sq["triangles"] == pytest.approx(4 * 60.0 * 60.0 / 0.25) and sq["grid"] == pytest.approx(120.0)


@pytest.mark.parametrize("n,mode,frac,rot", [(4, "circumscribed", 0.0, 0.0), (6, "circumscribed", 0.0, 0.0),
                                             (4, "inscribed", 0.0, 0.0), (5, "circumscribed", 0.3, 20.0)])
def test_cli_downloads_data_for_the_disc_that_encloses_the_whole_piece(monkeypatch, tmp_path, n, mode, frac, rot):
    """Un polygone circonscrit dépasse le cercle de référence : les données demandées à la source doivent
    couvrir toute la pièce (rayon extrême du contour), pas seulement le cercle de référence."""
    lat, lon = 46.0, 7.6
    nn = 400
    yy, xx = np.mgrid[0:nn, 0:nn]
    write_synthetic_tif(tmp_path / "cop.tif", (1000 + 0.5 * xx + 0.3 * yy).astype(np.float32), lon - 0.05, lat + 0.05,
                        0.1 / nn, 4326)
    radii = []
    monkeypatch.setattr(download, "download_copernicus",
                        lambda la, lo, radius_m, *a, **k: radii.append(radius_m) or [str(tmp_path / "cop.tif")])
    monkeypatch.chdir(tmp_path)
    a = RegularPolygon(n, mode, 1500.0).apothem
    monkeypatch.setattr(sys, "argv", ["x", "--download-copernicus", "--lat", str(lat), "--lon", str(lon),
                                      "--radius", "1500", "--diameter-mm", "40", "--vexag", "1", "--base-mm", "2",
                                      "--print-spacing-mm", "1", "--output", "cop.stl", "--cache-dir", "cache",
                                      "--sides", str(n), "--circle-mode", mode, "--corner-radius",
                                      str(frac * a * 40.0 / 3000.0), "--rotation", str(rot)])
    circle_dem_to_stl.main()
    expected = make(n, mode, frac, rot, ref=1500.0).extreme_radius
    assert radii == [pytest.approx(expected, rel=1e-9)]
    assert check_stl.edge_report(check_stl.read_stl(str(tmp_path / "output" / "cop.stl")))["boundary_edges"] == 0


def test_cli_without_sides_still_requests_the_reference_radius(monkeypatch, tmp_path):
    lat, lon = 46.0, 7.6
    nn = 400
    yy, xx = np.mgrid[0:nn, 0:nn]
    write_synthetic_tif(tmp_path / "cop.tif", (1000 + 0.5 * xx).astype(np.float32), lon - 0.05, lat + 0.05, 0.1 / nn,
                        4326)
    radii = []
    monkeypatch.setattr(download, "download_copernicus",
                        lambda la, lo, radius_m, *a, **k: radii.append(radius_m) or [str(tmp_path / "cop.tif")])
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["x", "--download-copernicus", "--lat", str(lat), "--lon", str(lon),
                                      "--radius", "1500", "--diameter-mm", "40", "--vexag", "1", "--base-mm", "2",
                                      "--print-spacing-mm", "1", "--output", "cop.stl", "--cache-dir", "cache"])
    circle_dem_to_stl.main()
    assert radii == [1500.0]
