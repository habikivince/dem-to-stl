# SPDX-License-Identifier: GPL-3.0-or-later
"""
lidar.py
--------
Nuages de points LiDAR (LAS 1.0-1.4, LAZ) -> grille d'altitudes régulière pour le moteur de maillage commun.

Un nuage n'est PAS un raster : il n'a ni grille, ni nodata, ni résolution native. Il n'est donc pas forcé dans
un lecteur raster. Le nuage est lu par blocs (mémoire constante, quelle que soit sa taille) et rastérisé
directement sur la grille de travail de raster.warp_window (même origine, même pas, même CRS cible), de sorte
qu'il n'y a AUCUN ré-échantillonnage supplémentaire (sauf si --lidar-cell demande une cellule différente du
pas de travail : la grille LiDAR est alors ré-échantillonnée une fois, bilinéaire ou moyenne).

Chaîne, par bloc de points :
  1. emprise : les points hors de la fenêtre demandée sont écartés avant tout calcul ;
  2. CRS : transformation vers le CRS de travail (UTM, ou le CRS projeté du nuage) si nécessaire ;
  3. écarts : points « withheld » (retenus par le fournisseur) et classes de bruit ASPRS 7 et 18 ;
  4. sélection : sol (classe 2), ou minimum / maximum par cellule pour un nuage non classifié ;
  5. doublons exacts (x, y, z identiques) retirés au sein du bloc ;
  6. agrégation par cellule (moyenne du sol, minimum ou maximum).
Puis, sur la grille complète : écart des cellules isolées aberrantes (pics et puits), puis comblement
documenté des cellules sans point (triangulation linéaire ou plus proche voisin, limité en distance ; jamais
d'extrapolation hors de l'emprise du nuage en mode linéaire). Les cellules restantes sont « sans donnée ».

Unités : les coordonnées sont celles du CRS du nuage ; l'altitude est convertie en mètres d'après l'unité
verticale du CRS (clé GeoTIFF 4099 ou CRS composé), sinon supposée identique à l'unité horizontale, avec
avertissement. Aucun changement de datum vertical n'est appliqué.

Dépendances : laspy (obligatoire pour LAS/LAZ), lazrs ou laszip (pour .laz). scipy sert au comblement.
"""
import math
import os
import struct

import numpy as np
from osgeo import gdal, osr
from scipy import ndimage

try:
    from scipy.spatial import QhullError
except ImportError:                     # scipy plus ancien
    from scipy.spatial.qhull import QhullError

import raster
from raster import FILL_NODATA, InputError

gdal.UseExceptions()
osr.UseExceptions()

CHUNK_POINTS = 2_000_000
GROUND_CLASS = 2
NOISE_CLASSES = (7, 18)              # ASPRS : bruit bas (7) et bruit haut (18)
DEFAULT_SPIKE_M = 10.0
DEFAULT_MAX_GAP_CELLS = 50           # comme gdal.FillNodata(maxSearchDist=50) du chemin raster
SELECTIONS = ("ground", "lowest", "highest")
INTERPOLATIONS = ("linear", "nearest", "none")

# codes EPSG d'unités de longueur (clé VerticalUnitsGeoKey)
_EPSG_UNITS = {9001: 1.0, 9002: 0.3048, 9003: 1200.0 / 3937.0}
_CLASS_NAMES = {0: "jamais classé", 1: "non classé", 2: "sol", 3: "végétation basse", 4: "végétation moyenne",
                5: "végétation haute", 6: "bâtiment", 7: "bruit bas", 8: "réservé", 9: "eau", 10: "rail",
                11: "route", 12: "réservé", 13: "fil", 14: "fil", 15: "pylône", 16: "connexion de fil",
                17: "pont", 18: "bruit haut"}


def _laspy():
    try:
        import laspy
    except ImportError:
        raise InputError("la lecture des fichiers LAS/LAZ requiert laspy, absent de cet environnement : "
                         "pip install laspy   (pour les .laz : pip install 'laspy[lazrs]')")
    return laspy


def _n(value):
    """Entier groupé par milliers (espace), pour des messages lisibles."""
    return f"{int(value):,}".replace(",", " ")


def _first_line(exc):
    text = str(exc).strip()
    return text.splitlines()[0] if text else type(exc).__name__


# ---------------------------------------------------------------------------
# En-tête et CRS
# ---------------------------------------------------------------------------

def _evlr_bytes(path, start, count):
    """Enregistrements EVLR (LAS 1.4) lus à la main : le WKT du CRS y est parfois stocké."""
    out = []
    try:
        with open(path, "rb") as f:
            offset = start
            for _ in range(count):
                f.seek(offset)
                head = f.read(60)
                if len(head) < 60:
                    break
                _reserved, user, record_id, length = struct.unpack_from("<H16sHQ", head)
                out.append((user.rstrip(b"\x00").decode("ascii", "ignore"), record_id, f.read(length)))
                offset += 60 + length
    except (OSError, struct.error):
        pass
    return out


def _wkt_text(blob):
    return blob.split(b"\x00", 1)[0].decode("utf-8", "ignore").strip()


def _geokeys(blob):
    """GeoKeyDirectoryTag (34735) -> {clé: valeur} pour les clés dont la valeur est directe."""
    if len(blob) < 8:
        return {}
    _version, _rev, _minor, n = struct.unpack_from("<4H", blob, 0)
    keys = {}
    for i in range(n):
        if 8 + 8 * i + 8 > len(blob):
            break
        key, location, _count, value = struct.unpack_from("<4H", blob, 8 + 8 * i)
        if location == 0:
            keys[key] = value
    return keys


def _crs_from_records(records):
    """(srs horizontal ou None, mètres par unité verticale ou None, description du datum vertical ou None)
    d'après les enregistrements (user_id, record_id, octets). Le WKT (2112) est prioritaire sur les clés
    GeoTIFF (34735). Un CRS défini par l'utilisateur dans les clés GeoTIFF (code 32767) est inutilisable."""
    wkt = next((b for u, r, b in records if u == "LASF_Projection" and r == 2112), None)
    if wkt:
        text = _wkt_text(wkt)
        try:
            srs = osr.SpatialReference()
            srs.ImportFromWkt(text)
            srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            vfactor, vname = None, None
            if srs.IsCompound():
                try:
                    vfactor = float(srs.GetTargetLinearUnits("VERT_CS"))
                except RuntimeError:
                    pass
                vname = srs.GetAttrValue("VERT_CS")
                srs.StripVertical()
            return srs, vfactor, vname
        except RuntimeError:
            pass
    blob = next((b for u, r, b in records if u == "LASF_Projection" and r == 34735), None)
    if blob:
        keys = _geokeys(blob)
        code = keys.get(3072) if keys.get(3072) not in (None, 0, 32767) else keys.get(2048)
        if code not in (None, 0, 32767):
            srs = osr.SpatialReference()
            try:
                srs.ImportFromEPSG(int(code))
            except RuntimeError:
                return None, None, None
            srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            vfactor = _EPSG_UNITS.get(keys.get(4099))
            vname = None
            if keys.get(4096) not in (None, 0, 32767):
                vname = f"EPSG:{keys[4096]}"
                if vfactor is None:
                    vs = osr.SpatialReference()
                    try:
                        vs.ImportFromEPSG(int(keys[4096]))
                        vfactor = float(vs.GetLinearUnits())
                    except RuntimeError:
                        pass
            return srs, vfactor, vname
    return None, None, None


class _FileInfo:
    def __init__(self, path, count, mins, maxs, version, point_format, compressed, srs, vfactor, vname):
        self.path, self.count, self.mins, self.maxs = path, count, mins, maxs
        self.version, self.point_format, self.compressed = version, point_format, compressed
        self.srs, self.vfactor, self.vname = srs, vfactor, vname


def _check_header(path):
    """Contrôle minimal de l'en-tête AVANT de laisser laspy l'analyser : sur un fichier corrompu, laspy peut tenter
    d'allouer des quantités de mémoire absurdes à partir de tailles lues dans des octets de bruit."""
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as f:
            head = f.read(375)
    except OSError as e:
        raise InputError(f"fichier LAS/LAZ illisible : {path} ({e.strerror or e})")
    if len(head) < 227 or head[:4] != b"LASF":
        raise InputError(f"fichier LAS/LAZ illisible ou corrompu : {path} (en-tête LAS absent ou incomplet)")
    major, minor = head[24], head[25]
    header_size, = struct.unpack_from("<H", head, 94)
    point_offset, n_vlr = struct.unpack_from("<II", head, 96)
    problem = None
    if major != 1 or minor > 4:
        problem = f"version LAS {major}.{minor} inconnue"
    elif not 227 <= header_size <= min(size, 375):
        problem = f"taille d'en-tête incohérente ({header_size} octets)"
    elif not header_size <= point_offset <= size:
        problem = f"début des points ({point_offset}) hors du fichier de {size} octets : fichier tronqué ?"
    elif n_vlr * 54 > point_offset - header_size:
        problem = f"{n_vlr} enregistrements VLR annoncés, plus qu'il n'en tient avant les points"
    if problem:
        raise InputError(f"fichier LAS/LAZ illisible ou corrompu : {path} ({problem})")


def _inspect(path, laspy):
    """Ouvre l'en-tête, lit un premier bloc de points (détecte un LAZ sans décompresseur, un fichier corrompu)."""
    _check_header(path)
    try:
        with laspy.open(path) as f:
            h = f.header
            records = [(v.user_id, v.record_id, bytes(v.record_data_bytes())) for v in h.vlrs]
            n_evlr = getattr(h, "number_of_evlrs", 0) or 0
            if n_evlr and not any(u == "LASF_Projection" for u, _r, _b in records):
                records += _evlr_bytes(path, h.start_of_first_evlr, n_evlr)
            info = _FileInfo(path, int(h.point_count), tuple(map(float, h.mins)), tuple(map(float, h.maxs)),
                             f"{h.version.major}.{h.version.minor}", int(h.point_format.id),
                             bool(h.are_points_compressed), *_crs_from_records(records))
            if info.count > 0:
                next(f.chunk_iterator(min(1000, info.count)))
    except InputError:
        raise
    except Exception as e:      # laspy lève des types variés (LaspyException, struct.error, ValueError, OSError)
        text = _first_line(e)
        if "LazBackend" in text or "laz" in text.lower() and "backend" in text.lower():
            raise InputError(f"{os.path.basename(path)} est compressé (LAZ) mais aucun décompresseur n'est "
                             f"installé : pip install 'laspy[lazrs]'  ({text[:80]})")
        raise InputError(f"fichier LAS/LAZ illisible ou corrompu : {path} ({text[:120]})")
    if info.count <= 0:
        raise InputError(f"nuage de points vide : {path}")
    if not all(math.isfinite(v) for v in info.mins + info.maxs):
        raise InputError(f"en-tête LAS invalide (emprise non finie) : {path}")
    return info


class LidarInput:
    """Entrée LAS/LAZ validée (en-têtes seulement : aucun point n'est chargé avant grid_window)."""

    is_lidar = True

    def __init__(self, paths, opts, log=raster.log):
        self.log = log
        self.paths = list(paths)
        self.select = opts.lidar_select or "ground"
        self.interp = opts.lidar_interp or "linear"
        self.cell = opts.lidar_cell
        self.max_gap_m = opts.lidar_max_gap
        self.spike_m = DEFAULT_SPIKE_M if opts.lidar_spike_m is None else float(opts.lidar_spike_m)
        if self.select not in SELECTIONS:
            raise InputError(f"--lidar-select : valeur inconnue {self.select!r} (choix : {', '.join(SELECTIONS)})")
        if self.interp not in INTERPOLATIONS:
            raise InputError(f"--lidar-interp : valeur inconnue {self.interp!r} "
                             f"(choix : {', '.join(INTERPOLATIONS)})")
        laspy = _laspy()
        self.infos = [_inspect(p, laspy) for p in self.paths]

        first = self.infos[0]
        if first.srs is not None and opts.source_crs is not None:
            raise InputError(f"--source-crs : {os.path.basename(first.path)} déclare déjà un CRS ; l'option ne "
                             f"sert qu'aux fichiers qui n'en ont pas.")
        if first.srs is None:
            if opts.source_crs is None:
                raise InputError(f"{os.path.basename(first.path)} ne déclare aucun système de coordonnées "
                                 f"(ni WKT ni clés GeoTIFF exploitables) : degrés ou mètres, datum ? Indique-le "
                                 f"avec --source-crs (ex. EPSG:2056).")
            srs = osr.SpatialReference()
            try:
                srs.SetFromUserInput(opts.source_crs)
            except RuntimeError as e:
                raise InputError(f"--source-crs : CRS non reconnu « {opts.source_crs} » ({_first_line(e)[:100]})")
            srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            first.srs = srs
            log(f"CRS du nuage imposé par --source-crs : {srs.GetName()}")
        self.srs = first.srs
        for info in self.infos[1:]:
            other = info.srs if info.srs is not None else self.srs
            if not other.IsSame(self.srs):
                raise InputError(f"les fichiers LAS/LAZ n'ont pas le même CRS ({os.path.basename(first.path)} : "
                                 f"{self.srs.GetName()} ; {os.path.basename(info.path)} : {other.GetName()}) : "
                                 f"reprojette-les avec PDAL avant de les assembler")
        self.srs_wkt = self.srs.ExportToWkt()
        self._check_extent()
        self.vfactor, vorigin = self._vertical_unit(opts, first)
        total = sum(i.count for i in self.infos)
        mins = [min(i.mins[k] for i in self.infos) for k in range(3)]
        maxs = [max(i.maxs[k] for i in self.infos) for k in range(3)]
        unit = self.srs.GetLinearUnitsName() if self.srs.IsProjected() else "degrés"
        log(f"entrée LAS/LAZ : {len(self.infos)} fichier(s), {_n(total)} points, LAS {first.version} format "
            f"{first.point_format}{' (compressé LAZ)' if first.compressed else ''}, CRS {self.srs.GetName()}, "
            f"emprise X {mins[0]:.2f}..{maxs[0]:.2f} Y {mins[1]:.2f}..{maxs[1]:.2f} ({unit}), "
            f"Z {mins[2] * self.vfactor:.2f}..{maxs[2] * self.vfactor:.2f} m ({vorigin})")
        if first.vname:
            log(f"  datum vertical déclaré : {first.vname} — AUCUN changement de datum vertical appliqué")

    def _check_extent(self):
        xs = [v for i in self.infos for v in (i.mins[0], i.maxs[0])]
        ys = [v for i in self.infos for v in (i.mins[1], i.maxs[1])]
        if self.srs.IsGeographic() and (max(map(abs, xs)) > 360.0 or max(map(abs, ys)) > 90.0):
            raise InputError(f"le CRS du nuage est géographique (degrés) mais ses coordonnées "
                             f"[{min(xs):g}..{max(xs):g}] x [{min(ys):g}..{max(ys):g}] sont impossibles en degrés : "
                             f"coordonnées projetées (mètres) ? Corrige le CRS (--source-crs).")

    def _vertical_unit(self, opts, info):
        from formats import VERTICAL_UNIT_CHOICES
        if opts.vertical_unit is not None:
            return VERTICAL_UNIT_CHOICES[opts.vertical_unit], f"--vertical-unit {opts.vertical_unit}"
        if info.vfactor:
            return float(info.vfactor), "unité verticale déclarée par le CRS"
        if self.srs.IsProjected():
            f = float(self.srs.GetLinearUnits())
            if abs(f - 1.0) < 1e-9:
                return 1.0, "mètres, comme l'unité horizontale"
            self.log(f"avertissement : unité verticale non déclarée, supposée identique à l'unité horizontale "
                     f"({self.srs.GetLinearUnitsName()}, 1 unité = {f:.7g} m) ; --vertical-unit m|ft|us-ft pour "
                     f"la corriger")
            return f, f"supposée {self.srs.GetLinearUnitsName()}"
        self.log("avertissement : CRS géographique sans unité verticale : altitudes supposées en mètres "
                 "(--vertical-unit m|ft|us-ft pour corriger)")
        return 1.0, "supposée en mètres"

    # ------------------------------------------------------------------
    # Rastérisation
    # ------------------------------------------------------------------

    def grid_window(self, dst_srs_wkt, cx, cy, radius, pixel_size_m, progress_cb=None, inside_fn=None,
                    fill=True, resample="auto"):
        """Même contrat que raster.warp_window : (elev, xs, ys, nodata) sur la fenêtre carrée de côté 2*radius
        centrée en (cx, cy), au pas pixel_size_m, dans dst_srs_wkt. fill=False : cellules brutes, sans écart des
        aberrants ni comblement."""
        px = float(pixel_size_m)
        xmin, ymax = cx - radius, cy + radius
        W = max(1, int((2.0 * radius + px / 2.0) / px))      # mêmes dimensions que gdal.Warp(outputBounds, xRes)
        H = W
        cell = float(self.cell) if self.cell else px
        Wc = max(1, int(math.ceil(W * px / cell - 1e-9)))
        Hc = Wc
        dst = osr.SpatialReference()
        dst.ImportFromWkt(dst_srs_wkt)
        dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        stats = _Stats()
        grid, count = self._accumulate(dst, xmin, ymax, cell, Wc, Hc, progress_cb, stats)
        self._report(stats, count, cell, px)

        valid = count > 0
        z = np.where(valid, grid, np.nan)
        if fill:
            z, valid = self._clean(z, valid, cell)
        if abs(cell - px) > 1e-9 * max(cell, px) or (Wc, Hc) != (W, H):
            z, valid = self._to_working_grid(z, valid, dst_srs_wkt, xmin, ymax, cell, px, W, H)
        elev = np.where(valid, z, FILL_NODATA).astype(np.float64)
        if not np.isfinite(elev).all():
            raise InputError("la grille LiDAR contient des valeurs non finies (NaN ou infini) : nuage corrompu ?")
        xs = xmin + (np.arange(W) + 0.5) * px
        ys = ymax - (np.arange(H) + 0.5) * px
        return elev, xs, ys, FILL_NODATA

    def _accumulate(self, dst, xmin, ymax, cell, Wc, Hc, progress_cb, stats):
        laspy = _laspy()
        xmax, ymin = xmin + Wc * cell, ymax - Hc * cell
        total = sum(i.count for i in self.infos)
        select = self.select
        if select == "ground":
            acc_sum, acc_cnt = np.zeros(Hc * Wc), np.zeros(Hc * Wc, dtype=np.int64)
        else:
            acc_sum = np.full(Hc * Wc, np.inf if select == "lowest" else -np.inf)
            acc_cnt = np.zeros(Hc * Wc, dtype=np.int64)
        done = 0
        if progress_cb is None and total >= 20_000_000:
            progress_cb = _LogProgress(self.log)
        for info in self.infos:
            same = dst.IsSame(self.srs)
            if same:
                box = (xmin, ymin, xmax, ymax)
                to_dst = None
            else:
                back = osr.CoordinateTransformation(dst, self.srs)
                to_dst = osr.CoordinateTransformation(self.srs, dst)
                edge = np.linspace(0.0, 1.0, 41)
                ring = np.concatenate([
                    np.column_stack([xmin + (xmax - xmin) * edge, np.full_like(edge, ymin)]),
                    np.column_stack([xmin + (xmax - xmin) * edge, np.full_like(edge, ymax)]),
                    np.column_stack([np.full_like(edge, xmin), ymin + (ymax - ymin) * edge]),
                    np.column_stack([np.full_like(edge, xmax), ymin + (ymax - ymin) * edge])])
                t = np.array(back.TransformPoints(ring))[:, :2]
                pad = 0.01 * max(np.ptp(t[:, 0]), np.ptp(t[:, 1]))
                box = (t[:, 0].min() - pad, t[:, 1].min() - pad, t[:, 0].max() + pad, t[:, 1].max() + pad)
            if info.maxs[0] < box[0] or info.mins[0] > box[2] or info.maxs[1] < box[1] or info.mins[1] > box[3]:
                done += info.count
                stats.files_skipped += 1
                self._progress(progress_cb, done, total)
                continue
            read = 0
            try:
                with laspy.open(info.path) as f:
                    for chunk in f.chunk_iterator(CHUNK_POINTS):
                        read += len(chunk)
                        self._consume(chunk, box, to_dst, xmin, ymax, cell, Wc, Hc, acc_sum, acc_cnt, stats)
                        self._progress(progress_cb, done + read, total)
            except InputError:
                raise
            except Exception as e:
                raise InputError(f"lecture LAS/LAZ interrompue dans {os.path.basename(info.path)} après "
                                 f"{_n(read)} points sur {_n(info.count)} (fichier tronqué ou corrompu : "
                                 f"{_first_line(e)[:100]})")
            if read < info.count:
                raise InputError(f"fichier LAS/LAZ tronqué : {os.path.basename(info.path)} annonce "
                                 f"{_n(info.count)} points, {_n(read)} lus")
            stats.read += read
            done += info.count
        if stats.in_window == 0:
            lo = [min(i.mins[k] for i in self.infos) for k in range(2)]
            hi = [max(i.maxs[k] for i in self.infos) for k in range(2)]
            raise InputError(f"aucun point du nuage dans la zone demandée (fenêtre X {xmin:.1f}..{xmax:.1f}, "
                             f"Y {ymin:.1f}..{ymax:.1f} dans le CRS de travail ; emprise du nuage "
                             f"X {lo[0]:.1f}..{hi[0]:.1f}, Y {lo[1]:.1f}..{hi[1]:.1f} dans son CRS) : "
                             f"centre ou rayon hors du nuage ?")
        if stats.used == 0:
            raise self._no_points_error(stats)
        if select == "ground":
            with np.errstate(invalid="ignore", divide="ignore"):
                grid = np.where(acc_cnt > 0, acc_sum / np.maximum(acc_cnt, 1), np.nan)
        else:
            grid = np.where(acc_cnt > 0, acc_sum, np.nan)
        return grid.reshape(Hc, Wc), acc_cnt.reshape(Hc, Wc)

    def _no_points_error(self, stats):
        present = ", ".join(f"{c} ({_CLASS_NAMES.get(c, 'autre')}) : {_n(n)}" for c, n in stats.top_classes(6))
        if self.select == "ground":
            return InputError(
                f"aucun point classé « sol » (classe 2) dans la zone : le nuage n'est pas classifié (classes "
                f"présentes : {present}). Classe-le d'abord (PDAL filters.csf / filters.smrf, LAStools "
                f"lasground), ou choisis explicitement une méthode : --lidar-select lowest (minimum par cellule, "
                f"approximation du sol, imprécise sous la végétation dense et le bâti) ou --lidar-select highest "
                f"(surface avec végétation et bâtiments).")
        return InputError(f"aucun point exploitable dans la zone (classes présentes : {present} ; les points "
                          f"« withheld » et les classes de bruit 7 et 18 sont toujours écartés)")

    @staticmethod
    def _progress(cb, done, total):
        if cb is not None and total:
            cb(min(1.0, done / total), "", None)

    def _consume(self, chunk, box, to_dst, xmin, ymax, cell, Wc, Hc, acc_sum, acc_cnt, stats):
        x = np.asarray(chunk.x, dtype=np.float64)
        y = np.asarray(chunk.y, dtype=np.float64)
        sel = np.flatnonzero((x >= box[0]) & (x <= box[2]) & (y >= box[1]) & (y <= box[3]))
        if sel.size == 0:
            return
        x, y = x[sel], y[sel]
        z = np.asarray(chunk.z, dtype=np.float64)[sel] * self.vfactor
        cls = np.asarray(chunk.classification)[sel].astype(np.int64)
        try:
            withheld = np.asarray(chunk.withheld, dtype=bool)[sel]
        except (AttributeError, ValueError):
            withheld = np.zeros(sel.size, dtype=bool)
        if to_dst is not None:
            t = np.array(to_dst.TransformPoints(np.column_stack([x, y])))
            x, y = t[:, 0], t[:, 1]
        ix = np.floor((x - xmin) / cell)
        iy = np.floor((ymax - y) / cell)
        inb = (ix >= 0) & (ix < Wc) & (iy >= 0) & (iy < Hc) & np.isfinite(ix) & np.isfinite(iy)
        if not inb.any():
            return
        x, y, z, cls, withheld = x[inb], y[inb], z[inb], cls[inb], withheld[inb]
        flat = (iy[inb].astype(np.int64) * Wc + ix[inb].astype(np.int64))
        stats.in_window += x.size
        stats.add_classes(cls)

        keep = np.isfinite(z) & ~withheld & ~np.isin(cls, NOISE_CLASSES)
        stats.dropped_flags += int((~keep).sum())
        if self.select == "ground":
            keep &= cls == GROUND_CLASS
        if not keep.any():
            return
        x, y, z, flat = x[keep], y[keep], z[keep], flat[keep]
        # doublons exacts (x, y, z) dans le bloc
        order = np.lexsort((z, y, x))
        x, y, z, flat = x[order], y[order], z[order], flat[order]
        dup = np.zeros(x.size, dtype=bool)
        dup[1:] = (x[1:] == x[:-1]) & (y[1:] == y[:-1]) & (z[1:] == z[:-1])
        if dup.any():
            stats.duplicates += int(dup.sum())
            x, y, z, flat = x[~dup], y[~dup], z[~dup], flat[~dup]
        stats.used += z.size

        order = np.argsort(flat, kind="stable")
        flat, z = flat[order], z[order]
        starts = np.flatnonzero(np.r_[True, flat[1:] != flat[:-1]])
        cells = flat[starts]
        counts = np.diff(np.r_[starts, flat.size])
        if self.select == "ground":
            acc_sum[cells] += np.add.reduceat(z, starts)
        elif self.select == "lowest":
            acc_sum[cells] = np.minimum(acc_sum[cells], np.minimum.reduceat(z, starts))
        else:
            acc_sum[cells] = np.maximum(acc_sum[cells], np.maximum.reduceat(z, starts))
        acc_cnt[cells] += counts

    # ------------------------------------------------------------------
    # Bilan, nettoyage, comblement
    # ------------------------------------------------------------------

    def _report(self, stats, count, cell, px):
        log = self.log
        classes = ", ".join(f"{c} ({_CLASS_NAMES.get(c, 'autre')}) : {_n(n)}" for c, n in stats.top_classes(6))
        log(f"LiDAR : {_n(stats.read)} points lus par blocs de {_n(CHUNK_POINTS)}, {_n(stats.in_window)} dans la "
            f"zone (classes : {classes})")
        classified = {c for c, _ in stats.top_classes(256)} - {0, 1}
        if GROUND_CLASS not in classified and self.select != "ground":
            log("avertissement : nuage NON classifié (aucun point de classe 2 dans la zone) : "
                + ("le minimum par cellule n'est qu'une approximation du sol, trop haut sous la végétation dense "
                   "et le bâti si la cellule est petite devant ces objets" if self.select == "lowest" else
                   "la surface comprend végétation et bâtiments"))
        elif GROUND_CLASS in classified and self.select != "ground":
            log(f"remarque : le nuage est classifié ; --lidar-select ground donnerait le sol réel au lieu du "
                f"{'minimum' if self.select == 'lowest' else 'maximum'} par cellule")
        log(f"  points retenus ({self.select}) : {_n(stats.used)} ; écartés : {_n(stats.dropped_flags)} "
            f"(withheld, bruit 7/18 ou altitude non finie) et {_n(stats.duplicates)} doublons exacts (au sein "
            f"d'un bloc)")
        n_cells = count.size
        occupied = int((count > 0).sum())
        density = stats.used / (n_cells * cell * cell)
        log(f"  grille LiDAR {count.shape[1]} x {count.shape[0]}, cellule {cell:.3g} m : {occupied} cellule(s) "
            f"avec point ({100.0 * occupied / n_cells:.1f} %), densité retenue ≈ {density:.2f} points/m²")
        if occupied / n_cells < 0.5:
            log(f"avertissement : moins de la moitié des cellules contient un point : le nuage est trop clairsemé "
                f"pour une cellule de {cell:.3g} m ; le détail comblé serait inventé. Augmente --lidar-cell "
                f"(≈ {math.sqrt(4.0 / max(density, 1e-9)):.2g} m donnerait ~4 points par cellule).")

    def _clean(self, z, valid, cell):
        log = self.log
        if self.spike_m > 0:
            flagged = despike(z, valid, self.spike_m)
            n = int(flagged.sum())
            if n:
                log(f"{n} cellule(s) isolée(s) aberrante(s) écartée(s) (au-delà de {self.spike_m:g} m + l'étendue "
                    f"des 8 voisines ; --lidar-spike-m 0 pour désactiver)")
                valid = valid & ~flagged
                z = np.where(valid, z, np.nan)
        gap_cells = (DEFAULT_MAX_GAP_CELLS if self.max_gap_m is None else self.max_gap_m / cell)
        if self.interp != "none" and gap_cells > 0 and not valid.all() and valid.any():
            before = int(valid.sum())
            z, valid = fill_gaps(z, valid, self.interp, gap_cells)
            n = int(valid.sum()) - before
            left = int((~valid).sum())
            log(f"{n} cellule(s) sans point comblée(s) par interpolation {self.interp} (distance ≤ "
                f"{gap_cells * cell:.3g} m = {gap_cells:.3g} cellules) ; {left} restent sans donnée"
                + (" (pas d'extrapolation hors de l'emprise du nuage)" if self.interp == "linear" else ""))
        elif self.interp == "none":
            log(f"aucun comblement (--lidar-interp none) : {int((~valid).sum())} cellule(s) sans donnée")
        return z, valid

    def _to_working_grid(self, z, valid, dst_srs_wkt, xmin, ymax, cell, px, W, H):
        """Cellule LiDAR différente du pas de travail : un seul ré-échantillonnage GDAL, documenté."""
        Hc, Wc = z.shape
        src = gdal.GetDriverByName("MEM").Create("", Wc, Hc, 1, gdal.GDT_Float64)
        src.SetGeoTransform((xmin, cell, 0.0, ymax, 0.0, -cell))
        src.SetProjection(dst_srs_wkt)
        band = src.GetRasterBand(1)
        band.SetNoDataValue(FILL_NODATA)
        band.WriteArray(np.where(valid, z, FILL_NODATA))
        alg = "bilinear" if cell > px else "average"
        self.log(f"grille LiDAR de {cell:.3g} m ré-échantillonnée à {px:.3g} m ({alg}) : "
                 + ("la cellule choisie est plus grossière que le pas de travail" if cell > px else
                    "la cellule choisie est plus fine que le pas de travail"))
        out = gdal.Warp("", src, format="MEM", outputBounds=(xmin, ymax - H * px, xmin + W * px, ymax),
                        xRes=px, yRes=px, resampleAlg=alg, dstNodata=FILL_NODATA, outputType=gdal.GDT_Float64)
        arr = out.GetRasterBand(1).ReadAsArray()
        ok = np.isfinite(arr) & (arr != FILL_NODATA)
        return np.where(ok, arr, np.nan), ok


class _LogProgress:
    """Avancement par dixièmes dans le journal, pour les gros nuages (lecture de plusieurs minutes)."""

    def __init__(self, log):
        self.log, self.last = log, 0

    def __call__(self, fraction, _msg, _data):
        step = int(fraction * 10)
        if step > self.last:
            self.last = step
            self.log(f"LiDAR : {step * 10} % des points lus")
        return 1


class _Stats:
    def __init__(self):
        self.read = 0
        self.in_window = 0
        self.used = 0
        self.dropped_flags = 0
        self.duplicates = 0
        self.files_skipped = 0
        self._classes = np.zeros(256, dtype=np.int64)

    def add_classes(self, cls):
        self._classes += np.bincount(np.clip(cls, 0, 255), minlength=256)

    def top_classes(self, n):
        order = np.argsort(self._classes)[::-1]
        return [(int(c), int(self._classes[c])) for c in order[:n] if self._classes[c] > 0]


# ---------------------------------------------------------------------------
# Opérations sur la grille (pures, testées isolément)
# ---------------------------------------------------------------------------

_NEIGHBOURS = np.ones((3, 3), dtype=bool)
_NEIGHBOURS[1, 1] = False


def despike(z, valid, tol):
    """Masque des cellules isolées aberrantes : au-dessus (ou au-dessous) de TOUTES leurs voisines valides de
    plus de `tol` mètres en plus de l'étendue (max - min) de ces voisines. Exige au moins 5 voisines valides ;
    un relief rugueux élargit donc la marge, un pic ou un puits isolé (point de bruit) est écarté, un versant
    ou une falaise ne l'est pas. tol <= 0 : rien n'est écarté."""
    if tol <= 0:
        return np.zeros(z.shape, dtype=bool)
    with np.errstate(invalid="ignore"):
        hi = ndimage.maximum_filter(np.where(valid, z, -np.inf), footprint=_NEIGHBOURS, mode="constant",
                                    cval=-np.inf)
        lo = ndimage.minimum_filter(np.where(valid, z, np.inf), footprint=_NEIGHBOURS, mode="constant",
                                    cval=np.inf)
        n = ndimage.convolve(valid.astype(np.int8), _NEIGHBOURS.astype(np.int8), mode="constant", cval=0)
        spread = np.where(n >= 5, hi - lo, np.inf)
        return valid & ((z > hi + spread + tol) | (z < lo - spread - tol))


def fill_gaps(z, valid, method, gap_cells):
    """Comble les cellules invalides situées à moins de `gap_cells` cellules d'une cellule valide. `nearest` :
    valeur de la cellule valide la plus proche. `linear` : interpolation linéaire sur la triangulation de
    Delaunay des cellules valides voisines des trous (jamais hors de leur enveloppe convexe : pas
    d'extrapolation). Retourne (z, valid) ; z reste NaN là où rien n'est comblé."""
    out, ok = z.copy(), valid.copy()
    dist, (ri, ci) = ndimage.distance_transform_edt(~valid, return_indices=True)
    target = ~valid & (dist <= gap_cells)
    if not target.any():
        return out, ok
    if method == "nearest":
        out[target] = z[ri[target], ci[target]]
        ok[target] = True
        return out, ok
    from scipy.interpolate import LinearNDInterpolator
    ring = valid & (ndimage.distance_transform_edt(valid) <= 2.0)
    pts = np.column_stack(np.nonzero(ring))
    if len(pts) < 3:
        return out, ok
    try:
        interp = LinearNDInterpolator(pts.astype(np.float64), z[ring])
        tpts = np.column_stack(np.nonzero(target)).astype(np.float64)
        vals = interp(tpts)
    except (QhullError, ValueError):
        return out, ok
    good = np.isfinite(vals)
    rows, cols = np.nonzero(target)
    out[rows[good], cols[good]] = vals[good]
    ok[rows[good], cols[good]] = True
    return out, ok
