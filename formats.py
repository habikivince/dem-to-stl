# SPDX-License-Identifier: GPL-3.0-or-later
"""
formats.py
----------
Détection du format source et normalisation des rasters d'altitude avant le pipeline commun
(raster.warp_window -> mesh -> export). Chaîne :

    détection (extension, signature) -> lecteur GDAL adapté -> normalisation -> raster.warp_window ...

La normalisation ne touche jamais au moteur de maillage : elle rend, pour chaque fichier, un nom lisible
par GDAL qui désigne UNE bande d'altitudes en MÈTRES, avec son géoréférencement et son CRS connus. C'est
le fichier d'origine lui-même quand il n'y a rien à corriger (cas du GeoTIFF mono-bande habituel : aucun
changement de comportement), sinon un VRT en mémoire (/vsimem/) qui applique, sans copier les données :
  - le choix de la bande (jamais deviné) ;
  - le facteur d'échelle / décalage déclaré (entiers « empaquetés » des NetCDF, GeoPackage 16 bits...) ;
  - la conversion d'unité verticale vers le mètre (déclarée par le fichier ou --vertical-unit) ;
  - l'inversion de signe d'une variable `positive = down` (profondeurs) ;
  - le CRS fourni par --source-crs quand le fichier n'en déclare aucun ;
  - le retrait de la composante verticale d'un CRS composé (aucun changement de datum vertical).

Les fichiers LAS/LAZ (nuages de points) ne sont PAS des rasters : ils passent par lidar.py.

Formats : GeoTIFF (.tif/.tiff), ESRI ASCII Grid (.asc), SRTM (.hgt), GeoPackage raster (.gpkg),
NetCDF (.nc), HDF5 (.h5/.hdf5, uniquement s'il est géoréférencé : NetCDF-4/CF ou sous-jeu de données
reconnu par GDAL), VRT et tout autre raster que GDAL sait lire (comportement historique), LAS/LAZ.
"""
import math
import os
import re
import uuid
from dataclasses import dataclass

from osgeo import gdal, osr

import raster
from raster import InputError

gdal.UseExceptions()
osr.UseExceptions()

RASTER_EXTENSIONS = (".tif", ".tiff", ".asc", ".hgt", ".gpkg", ".nc", ".h5", ".hdf5", ".vrt")
LIDAR_EXTENSIONS = (".las", ".laz")
SUPPORTED_EXTENSIONS = RASTER_EXTENSIONS + LIDAR_EXTENSIONS

# mètres par unité, d'après le nom normalisé (minuscules, sans espace, tiret ni souligné)
_LENGTH_UNITS = {
    "m": 1.0, "metre": 1.0, "metres": 1.0, "meter": 1.0, "meters": 1.0,
    "ft": 0.3048, "foot": 0.3048, "feet": 0.3048, "internationalfoot": 0.3048, "internationalfeet": 0.3048,
    "usft": 1200.0 / 3937.0, "ftus": 1200.0 / 3937.0, "footus": 1200.0 / 3937.0, "ussurveyfoot": 1200.0 / 3937.0,
    "ussurveyfeet": 1200.0 / 3937.0, "usfoot": 1200.0 / 3937.0,
    "km": 1000.0, "kilometre": 1000.0, "kilometer": 1000.0, "kilometres": 1000.0, "kilometers": 1000.0,
    "cm": 0.01, "centimetre": 0.01, "centimeter": 0.01, "centimetres": 0.01, "centimeters": 0.01,
    "mm": 0.001, "millimetre": 0.001, "millimeter": 0.001, "millimetres": 0.001, "millimeters": 0.001,
}
VERTICAL_UNIT_CHOICES = {"m": 1.0, "ft": 0.3048, "us-ft": 1200.0 / 3937.0}

# CF : standard_name d'une variable qui est une altitude (la hauteur de l'eau, la température... n'en sont pas)
_ELEVATION_STANDARD_NAMES = {
    "altitude", "height", "surface_altitude", "height_above_mean_sea_level", "height_above_reference_ellipsoid",
    "height_above_geopotential_datum", "geoid_height_above_reference_ellipsoid", "land_surface_altitude",
    "sea_floor_depth_below_geoid", "sea_floor_depth_below_sea_level", "sea_floor_depth_below_sea_surface",
    "sea_floor_depth_below_reference_ellipsoid", "bedrock_altitude", "surface_altitude_above_geoid",
}
_ELEVATION_NAMES = {
    "elevation", "altitude", "height", "z", "dem", "dtm", "dsm", "topo", "topography", "relief", "orography",
    "orog", "elev", "alt", "hgt", "terrain", "bathymetry", "depth", "surface_altitude",
}
_LATITUDE_UNITS = {"degreesnorth", "degreenorth", "degreesn", "degreen", "degrees_north", "degree_north"}
_LONGITUDE_UNITS = {"degreeseast", "degreeeast", "degreese", "degreee", "degrees_east", "degree_east"}


@dataclass
class SourceOptions:
    """Options de lecture propres aux formats. Toute valeur à None = non demandée par l'utilisateur."""
    layer: str = None            # table GeoPackage, variable NetCDF, dataset HDF5
    band: int = None             # bande / tranche à utiliser (fichier multibande, dimension temps...)
    source_crs: str = None       # CRS à attribuer à un fichier qui n'en déclare aucun
    vertical_unit: str = None    # m | ft | us-ft : remplace l'unité verticale déclarée ou absente
    lidar_select: str = None     # ground | lowest | highest
    lidar_cell: float = None     # taille de cellule de la grille LiDAR, en mètres (défaut : pixel de travail)
    lidar_interp: str = None     # linear | nearest | none
    lidar_max_gap: float = None  # distance maximale de comblement, en mètres (défaut : 50 cellules)
    lidar_spike_m: float = None  # tolérance d'isolement des pics/puits, en mètres (0 = désactivé)

    def raster_options(self):
        return [n for n, v in (("--layer", self.layer), ("--band", self.band), ("--source-crs", self.source_crs),
                               ("--vertical-unit", self.vertical_unit)) if v is not None]

    def lidar_options(self):
        return [n for n, v in (("--lidar-select", self.lidar_select), ("--lidar-cell", self.lidar_cell),
                               ("--lidar-interp", self.lidar_interp), ("--lidar-max-gap", self.lidar_max_gap),
                               ("--lidar-spike-m", self.lidar_spike_m)) if v is not None]


# ---------------------------------------------------------------------------
# Détection
# ---------------------------------------------------------------------------

_VIRTUAL_PREFIXES = ("/vsi", "NETCDF:", "HDF5:", "GPKG:", "http://", "https://")


def _is_virtual(path):
    return path.startswith(_VIRTUAL_PREFIXES)


def detect(path):
    """'lidar' (nuage de points LAS/LAZ) ou 'raster' (tout ce que GDAL lit). Détection par l'extension puis par
    la signature du fichier (les 4 premiers octets d'un LAS/LAZ valent « LASF »), pour qu'un nuage renommé ne
    soit pas confondu avec un raster."""
    if _is_virtual(path):
        return "raster"
    is_las = os.path.splitext(path)[1].lower() in LIDAR_EXTENSIONS
    what = "fichier LAS/LAZ" if is_las else "raster"           # « raster illisible » : formulation historique du programme
    if not os.path.exists(path):
        raise InputError(f"{what} illisible : {path} (fichier introuvable)")
    if not os.path.isfile(path):
        raise InputError(f"{what} illisible : {path} (ce n'est pas un fichier)")
    if is_las:
        return "lidar"
    try:
        with open(path, "rb") as f:
            if f.read(4) == b"LASF":
                return "lidar"
    except OSError as e:
        raise InputError(f"{what} illisible : {path} ({e.strerror or e})")
    return "raster"


def prepare(paths, opts=None, log=raster.log):
    """Point d'entrée unique. Retourne (genre, entrées) :
      - ('raster', [noms GDAL]) : un nom par fichier, désignant une bande d'altitudes en mètres ;
      - ('lidar', lidar.LidarInput) : nuage de points, rastérisé plus tard sur la grille de travail.
    Lève InputError, avec un message exploitable, pour tout fichier inutilisable : ambigu (plusieurs bandes,
    variables, couches), sans CRS, tronqué, ou qui n'est pas un MNT."""
    opts = opts or SourceOptions()
    paths = [str(p) for p in (paths if isinstance(paths, (list, tuple)) else [paths])]
    if not paths:
        raise InputError("aucun fichier d'entrée")
    kinds = {detect(p) for p in paths}
    if len(kinds) > 1:
        raise InputError("--input mélange des nuages de points (LAS/LAZ) et des rasters : un seul genre à la fois")
    if kinds == {"lidar"}:
        bad = [n for n in opts.raster_options() if n in ("--layer", "--band")]   # --source-crs et --vertical-unit valent aussi pour LAS/LAZ
        if bad:
            raise InputError(f"{', '.join(bad)} ne s'applique pas aux fichiers LAS/LAZ (voir les options --lidar-*)")
        import lidar
        return "lidar", lidar.LidarInput(paths, opts, log=log)
    bad = opts.lidar_options()
    if bad:
        raise InputError(f"{', '.join(bad)} ne s'applique qu'aux fichiers LAS/LAZ")
    if opts.vertical_unit is not None and opts.vertical_unit not in VERTICAL_UNIT_CHOICES:
        raise InputError(f"--vertical-unit : valeur inconnue {opts.vertical_unit!r} "
                         f"(choix : {', '.join(VERTICAL_UNIT_CHOICES)})")
    names = [_prepare_raster(p, opts, log) for p in paths]
    _report_srtm_assembly(paths, log)
    return "raster", _mosaic_srtm(paths, names, log)


# ---------------------------------------------------------------------------
# Ouverture par format
# ---------------------------------------------------------------------------

def _first_line(exc):
    text = str(exc).strip()
    return text.splitlines()[0] if text else "format non reconnu"


def _gdal_open(path, **kw):
    try:
        return gdal.OpenEx(path, gdal.OF_RASTER, **kw) if kw else gdal.Open(path)
    except RuntimeError as e:
        raise InputError(f"raster illisible : {path} ({_first_line(e)})")


_HGT_NAME = re.compile(r"^([NS])(\d{2})([EW])(\d{3})\.hgt$", re.IGNORECASE)


def _check_hgt(path):
    """Une tuile SRTM n'a pas d'en-tête : sa position vient de son NOM (N46E007.hgt = coin sud-ouest) et sa
    taille donne la résolution. GDAL échoue sans explication sur un nom ou une taille inattendus."""
    name = os.path.basename(path)
    if not _HGT_NAME.match(name):
        raise InputError(f"tuile SRTM {name} : le nom doit être de la forme N46E007.hgt (coin sud-ouest de la "
                         f"tuile, 1° x 1°) car un fichier .hgt ne contient aucun géoréférencement")
    size = os.path.getsize(path)
    side = math.isqrt(size // 2)
    if size % 2 or side * side * 2 != size or side not in (1201, 3601):
        raise InputError(f"tuile SRTM {name} : taille inattendue ({size} octets) ; seules les tuiles SRTM3 "
                         f"(1201 x 1201 échantillons, 2 884 802 octets) et SRTM1 (3601 x 3601, 25 934 402 octets) "
                         f"sont reconnues. Fichier tronqué ou autre produit ?")


def _check_asc_header(path):
    """En-tête ESRI ASCII Grid vérifié AVANT l'ouverture par GDAL, dont le refus (« format non reconnu ») ne dit
    pas quel champ manque. Retourne les champs lus (clés en minuscules)."""
    keys = {}
    header_lines = 0
    try:
        with open(path, "r", errors="replace") as f:
            for _ in range(16):
                line = f.readline()
                if not line:
                    break
                parts = line.split()
                if not parts:
                    header_lines += 1
                    continue
                if not re.match(r"^[A-Za-z_]+$", parts[0]):
                    break
                if len(parts) != 2:
                    raise InputError(f"en-tête ASC mal formé dans {path} : ligne « {line.strip()[:60]} »")
                keys[parts[0].lower()] = parts[1]
                header_lines += 1
    except OSError as e:
        raise InputError(f"fichier illisible : {path} ({e.strerror or e})")
    missing = [k for k in ("ncols", "nrows") if k not in keys]
    if "xllcorner" not in keys and "xllcenter" not in keys:
        missing.append("xllcorner ou xllcenter")
    if "yllcorner" not in keys and "yllcenter" not in keys:
        missing.append("yllcorner ou yllcenter")
    if "cellsize" not in keys and not ("dx" in keys and "dy" in keys):
        missing.append("cellsize")
    if missing:
        raise InputError(f"en-tête ASC incomplet dans {path} : champ(s) manquant(s) : {', '.join(missing)}")

    def num(key, kind=float):
        try:
            value = kind(keys[key])
        except ValueError:
            raise InputError(f"en-tête ASC invalide dans {path} : {key} = « {keys[key]} » n'est pas un nombre")
        if kind is float and not math.isfinite(value):
            raise InputError(f"en-tête ASC invalide dans {path} : {key} n'est pas fini")
        return value

    for key in ("ncols", "nrows"):
        if num(key, int) <= 0:
            raise InputError(f"en-tête ASC invalide dans {path} : {key} doit être strictement positif")
    for key in ("cellsize", "dx", "dy"):
        if key in keys and num(key) <= 0:
            raise InputError(f"en-tête ASC invalide dans {path} : {key} doit être strictement positif")
    for key in ("xllcorner", "xllcenter", "yllcorner", "yllcenter"):
        if key in keys:
            num(key)
    if "nodata_value" in keys:
        try:
            float(keys["nodata_value"])
        except ValueError:
            raise InputError(f"en-tête ASC invalide dans {path} : NODATA_value = « {keys['nodata_value']} »")
    _check_asc_body(path, int(keys["ncols"]), int(keys["nrows"]), header_lines)
    return keys


def _check_asc_body(path, ncols, nrows, header_lines):
    """Le corps contient-il ncols x nrows valeurs ? Contrôle fait SANS GDAL : sur un fichier tronqué, la lecture
    GDAL d'un ASC peut ne jamais rendre la main. Voie rapide : nombre de lignes = nrows et dernière ligne à
    ncols valeurs (fichier normal, un rang par ligne) ; sinon, comptage complet des valeurs."""
    try:
        size = os.path.getsize(path)
        newlines = 0
        ends_with_newline = True
        with open(path, "rb") as f:
            while True:
                chunk = f.read(8 << 20)
                if not chunk:
                    break
                newlines += chunk.count(b"\n")
                ends_with_newline = chunk.endswith(b"\n")
            tail_len = min(size, 4 << 20)
            f.seek(size - tail_len)
            tail = f.read(tail_len).rstrip().splitlines()
        data_lines = newlines - header_lines + (0 if ends_with_newline else 1)
        if tail and data_lines == nrows and len(tail[-1].split()) == ncols:
            return
        values = 0
        carry = b""
        skipped = 0
        with open(path, "rb") as f:
            while skipped < header_lines:
                f.readline()
                skipped += 1
            while True:
                chunk = f.read(8 << 20)
                if not chunk:
                    break
                chunk = carry + chunk
                cut = max(chunk.rfind(b" "), chunk.rfind(b"\n"), chunk.rfind(b"\t"), chunk.rfind(b"\r"))
                carry, chunk = chunk[cut + 1:], chunk[:cut + 1]
                values += len(chunk.split())
            values += len(carry.split())
    except OSError as e:
        raise InputError(f"fichier illisible : {path} ({e.strerror or e})")
    if values < ncols * nrows:
        raise InputError(f"fichier ASC tronqué : {path} annonce {ncols} x {nrows} = {ncols * nrows} valeurs "
                         f"(en-tête) mais n'en contient que {values}")


def _check_asc_data(path, ds, keys):
    """Lecture de la dernière ligne de données : GDAL n'indexe les lignes qu'à la demande, et un fichier tronqué
    ne ferait sinon échouer qu'en plein ré-échantillonnage."""
    if (ds.RasterXSize, ds.RasterYSize) != (int(keys["ncols"]), int(keys["nrows"])):
        raise InputError(f"fichier ASC incohérent : en-tête {keys['ncols']} x {keys['nrows']}, "
                         f"GDAL lit {ds.RasterXSize} x {ds.RasterYSize}")
    try:
        ds.GetRasterBand(1).ReadRaster(0, ds.RasterYSize - 1, ds.RasterXSize, 1)
    except RuntimeError as e:
        raise InputError(f"fichier ASC tronqué ou mal formé : {path} (lignes de données manquantes ou courtes "
                         f"par rapport à l'en-tête {keys['ncols']} x {keys['nrows']} : {_first_line(e)[:120]})")


def _entry_name(sub_name):
    """Nom lisible d'un sous-jeu de données GDAL : NETCDF:"f.nc":elevation -> elevation, GPKG:f.gpkg:t -> t."""
    return sub_name.rsplit(":", 1)[-1].strip('"')


def _match_entry(entries, wanted):
    """entries : [(nom_gdal, description, ...)]. Retourne l'entrée désignée par --layer : nom complet, ou dernier
    composant du chemin (« elevation » pour /data/elevation) s'il est unique ; None sinon."""
    w = wanted.strip()
    exact = [e for e in entries if _entry_name(e[0]) in (w, "/" + w.lstrip("/"), w.lstrip("/"))]
    if len(exact) == 1:
        return exact[0]
    tail = [e for e in entries if _entry_name(e[0]).rsplit("/", 1)[-1] == w]
    return tail[0] if len(tail) == 1 else None


def _kind_of(layer_meta):
    """Métadonnées CF d'une variable -> (standard_name, long_name, units, positive)."""
    return tuple((layer_meta.get(k) or "").strip() for k in ("standard_name", "long_name", "units", "positive"))


def _variable_metadata(ds):
    """Métadonnées CF de la variable désignée par la bande 1 d'un jeu NetCDF/HDF5 (dict vide si absentes)."""
    band = ds.GetRasterBand(1)
    bmd = band.GetMetadata()
    var = bmd.get("NETCDF_VARNAME", "")
    md = {k: v for k, v in bmd.items()}
    for key in ("standard_name", "long_name", "units", "positive"):
        v = ds.GetMetadataItem(f"{var}#{key}")
        if v and key not in md:
            md[key] = v
    return var, md


def _looks_like_elevation(name, md):
    standard, long_name, _units, _ = _kind_of(md)
    if standard.lower() in _ELEVATION_STANDARD_NAMES:
        return True
    if name.rsplit("/", 1)[-1].lower() in _ELEVATION_NAMES:
        return True
    return any(w in long_name.lower() for w in ("elevation", "altitude", "topograph", "bathymetr", "height above"))


def _describe_entries(entries):
    lines = []
    for name, desc, probable in entries:
        lines.append(f"    {_entry_name(name)}  {desc.split('] ', 1)[-1] if '] ' in desc else desc}"
                     + ("   <- ressemble à une altitude" if probable else ""))
    return "\n".join(lines)


def _require_georeferenced(kind, path, ds, names):
    """HDF5 : aucune convention de géoréférencement n'est universelle. Un jeu sans transformation affine reconnue
    est refusé AVEC la liste de ses datasets, plutôt que de deviner une emprise. (Un NetCDF-4 enregistré en .h5,
    avec ses axes lat/lon ou x/y CF, est géoréférencé et passe.)"""
    if ds.GetGeoTransform(can_return_null=True) is not None:
        return
    if kind == "HDF5":
        raise InputError(f"HDF5 non géoréférencé : {path} (datasets : {', '.join(names)}). Seuls sont pris en "
                         f"charge les fichiers HDF5 qui suivent le modèle NetCDF-4/CF, avec des axes spatiaux "
                         f"réguliers (lat/lon ou x/y) et un CRS ; une grille brute sans position ne permet ni "
                         f"l'échelle ni le placement. Ajoute le géoréférencement (gdal_translate -a_srs -a_ullr, "
                         f"ou conversion en NetCDF/GeoTIFF) puis relance.")
    raise InputError(f"NetCDF sans axes spatiaux réguliers : {path} (variable(s) : {', '.join(names)}). Il faut des "
                     f"coordonnées lat/lon ou x/y à pas constant ; une grille irrégulière ou sans coordonnées "
                     f"n'est pas prise en charge.")


def _open_scientific(path, opts):
    """NetCDF (.nc) et HDF5 (.h5/.hdf5) par le pilote netCDF de GDAL, forcé par le préfixe NETCDF: (un NetCDF-4
    renommé .h5 est sinon pris pour un HDF5 générique, sans géoréférencement). Retourne (ds, nom GDAL, variable).
    Une variable n'est jamais prise pour une altitude sur sa seule présence : sans --layer, il faut qu'elle se
    déclare comme telle (standard_name CF ou nom usuel), et une seule candidate."""
    ext = os.path.splitext(path)[1].lower()
    kind = "NetCDF" if ext == ".nc" else "HDF5"
    gdal_name = f'NETCDF:"{path}"'
    try:
        ds = gdal.Open(gdal_name)
    except RuntimeError as e:
        # Pour un fichier présent mais invalide, GDAL répond « No such file or directory » : trompeur, on le tait.
        reason = _first_line(e)
        detail = "" if os.path.exists(path) and "No such file" in reason else f" ({reason})"
        raise InputError(f"fichier {kind} illisible : {path}{detail}. "
                         f"Corrompu, tronqué, ou bibliothèque HDF5/netCDF de GDAL sans ce format ?")
    subs = ds.GetSubDatasets()
    if ds.RasterCount == 0:
        if not subs:
            raise InputError(f"{kind} sans variable de grille 2D exploitable : {path}")
        entries = []
        for name, desc in subs:
            sub = gdal.Open(name)
            var, md = _variable_metadata(sub) if sub.RasterCount else (_entry_name(name), {})
            entries.append((name, desc, _looks_like_elevation(_entry_name(name), md)))
        if opts.layer is not None:
            chosen = _match_entry(entries, opts.layer)
            if chosen is None:
                raise InputError(f"--layer {opts.layer} : aucune variable de ce nom dans {path}. Disponibles :\n"
                                 + _describe_entries(entries))
        else:
            likely = [e for e in entries if e[2]]
            if len(likely) != 1:
                why = ("aucune ne se déclare comme une altitude" if not likely
                       else f"{len(likely)} se déclarent comme une altitude")
                raise InputError(f"{kind} {os.path.basename(path)} contient {len(entries)} variables de grille et "
                                 f"{why} : choisis-en une avec --layer NOM. Disponibles :\n"
                                 + _describe_entries(entries))
            chosen = likely[0]
        ds = gdal.Open(chosen[0])
        _require_georeferenced(kind, path, ds, [_entry_name(n) for n, _ in subs])
        return ds, chosen[0], _entry_name(chosen[0])
    var, md = _variable_metadata(ds)
    if opts.layer is not None:
        if opts.layer.strip().lstrip("/").rsplit("/", 1)[-1] != var.lstrip("/").rsplit("/", 1)[-1]:
            raise InputError(f"--layer {opts.layer} : le fichier ne contient que la variable « {var} »")
    elif not _looks_like_elevation(var, md):
        standard, long_name, units, _ = _kind_of(md)
        raise InputError(f"la variable « {var} » de {os.path.basename(path)} ne se déclare pas comme une altitude "
                         f"(standard_name « {standard or '-'} », long_name « {long_name or '-'} », "
                         f"unités « {units or '-'} »). Si c'en est une, ajoute --layer {var}.")
    _require_georeferenced(kind, path, ds, [var])
    return ds, gdal_name, var


def _open_gpkg(path, opts):
    """GeoPackage RASTER uniquement. Une table de tuiles couleur (fond de carte) n'est pas un MNT : seules les
    tables d'altitude (Int16/UInt16/Int32/Float32/Float64, extension « gridded coverage » ou PNG 16 bits)
    sont acceptées. Plusieurs tables raster : --layer obligatoire."""
    try:
        ds = gdal.OpenEx(path, gdal.OF_RASTER)
    except RuntimeError:
        try:
            vec = gdal.OpenEx(path, gdal.OF_VECTOR)
        except RuntimeError as e:
            raise InputError(f"GeoPackage illisible : {path} ({_first_line(e)})")
        names = [vec.GetLayer(i).GetName() for i in range(vec.GetLayerCount())]
        raise InputError(f"GeoPackage sans table raster : {path} ne contient que des couches vectorielles "
                         f"({', '.join(names) or 'aucune'}) ; un MNT doit être une table raster d'altitudes")
    subs = ds.GetSubDatasets()
    if ds.RasterCount == 0:
        if not subs:
            raise InputError(f"GeoPackage sans table raster : {path}")
        entries = [(n, d, False) for n, d in subs]
        if opts.layer is None:
            raise InputError(f"GeoPackage {os.path.basename(path)} : {len(entries)} tables raster, choisis-en une "
                             f"avec --layer NOM. Disponibles :\n" + _describe_entries(entries))
        chosen = _match_entry(entries, opts.layer)
        if chosen is None:
            raise InputError(f"--layer {opts.layer} : aucune table raster de ce nom dans {path}. Disponibles :\n"
                             + _describe_entries(entries))
        ds, name, layer = gdal.Open(chosen[0]), chosen[0], _entry_name(chosen[0])
    elif opts.layer is not None:
        try:
            ds = gdal.OpenEx(path, gdal.OF_RASTER, open_options=[f"TABLE={opts.layer}"])
        except RuntimeError as e:
            raise InputError(f"--layer {opts.layer} : table raster inconnue dans {path} ({_first_line(e)[:100]})")
        name, layer = f"GPKG:{path}:{opts.layer}", opts.layer
    else:
        name, layer = path, None
    bands = [ds.GetRasterBand(i + 1) for i in range(ds.RasterCount)]
    colors = {gdal.GetColorInterpretationName(b.GetColorInterpretation()) for b in bands}
    if ds.RasterCount > 1 or bands[0].DataType == gdal.GDT_Byte or colors & {"Red", "Green", "Blue", "Palette"}:
        raise InputError(f"la table raster {layer or os.path.basename(path)} de {path} contient des tuiles image "
                         f"({ds.RasterCount} bande(s) {gdal.GetDataTypeName(bands[0].DataType)}), pas des altitudes : "
                         f"un MNT GeoPackage est une grille Int16/UInt16/Float32 à une bande")
    return ds, name, layer


def _open(path, opts):
    """Retourne (ds, nom GDAL, couche, format) pour un fichier raster, en appliquant les règles du format."""
    ext = os.path.splitext(path)[1].lower()
    virtual = _is_virtual(path)
    if ext == ".gpkg" and not virtual:
        ds, name, layer = _open_gpkg(path, opts)
        return ds, name, layer, "GeoPackage"
    if ext in (".nc", ".h5", ".hdf5") and not virtual:
        ds, name, layer = _open_scientific(path, opts)
        return ds, name, layer, "NetCDF" if ext == ".nc" else "HDF5"
    if opts.layer is not None:
        raise InputError("--layer ne s'applique qu'aux fichiers GeoPackage (.gpkg), NetCDF (.nc) et HDF5 "
                         "(.h5, .hdf5)")
    if ext == ".hgt" and not virtual:
        _check_hgt(path)
    asc_keys = _check_asc_header(path) if ext == ".asc" and not virtual else None
    ds = _gdal_open(path)
    if ds.RasterCount < 1:
        if ds.GetSubDatasets():
            raise InputError(f"{path} est un conteneur de sous-jeux de données ; donne le nom d'un sous-jeu "
                             f"(ex. {ds.GetSubDatasets()[0][0]})")
        raise InputError(f"raster sans bande : {path}")
    if asc_keys is not None:
        _check_asc_data(path, ds, asc_keys)
    fmt = {".tif": "GeoTIFF", ".tiff": "GeoTIFF", ".asc": "ESRI ASCII Grid", ".hgt": "SRTM"}.get(ext)
    return ds, path, None, fmt or ds.GetDriver().LongName


# ---------------------------------------------------------------------------
# Bande, géoréférencement, unités
# ---------------------------------------------------------------------------

def _band_label(ds, i):
    b = ds.GetRasterBand(i)
    parts = [gdal.GetDataTypeName(b.DataType), gdal.GetColorInterpretationName(b.GetColorInterpretation())]
    if b.GetDescription():
        parts.insert(0, b.GetDescription())
    parts += [f"{k[len('NETCDF_DIM_'):]}={v}" for k, v in b.GetMetadata().items()
              if k.startswith("NETCDF_DIM_") and k != "NETCDF_DIM_EXTRA"]
    if b.GetUnitType():
        parts.append(f"unité {b.GetUnitType()}")
    return f"    bande {i} : " + ", ".join(parts)


def _select_band(ds, opts, path):
    """Jamais de choix arbitraire : une seule bande utile -> elle ; sinon --band est obligatoire. Une bande alpha
    (transparence) n'est pas une variable."""
    n = ds.RasterCount
    if opts.band is not None:
        if not 1 <= opts.band <= n:
            raise InputError(f"--band {opts.band} : {os.path.basename(path)} a {n} bande(s) (1 à {n})")
        return opts.band
    useful = [i for i in range(1, n + 1)
              if ds.GetRasterBand(i).GetColorInterpretation() != gdal.GCI_AlphaBand]
    if len(useful) == 1:
        return useful[0]
    cis = {gdal.GetColorInterpretationName(ds.GetRasterBand(i).GetColorInterpretation()) for i in useful}
    hint = (" Ressemble à une image couleur (rouge/vert/bleu), pas à un MNT."
            if {"Red", "Green", "Blue"} <= cis else "")
    shown = "\n".join(_band_label(ds, i) for i in useful[:12]) + ("\n    ..." if len(useful) > 12 else "")
    raise InputError(f"{os.path.basename(path)} contient {len(useful)} bandes (variables, dates ou canaux) et le "
                     f"programme ne choisit pas à ta place.{hint} Indique celle à utiliser avec --band N :\n{shown}")


def _unit_factor(text):
    return _LENGTH_UNITS.get(re.sub(r"[\s_\-.]", "", text.lower()))


def _crs_vertical_factor(srs):
    """Mètres par unité de la composante verticale d'un CRS composé ; None si le CRS n'en a pas."""
    if srs is None or not srs.IsCompound():
        return None
    try:
        return float(srs.GetTargetLinearUnits("VERT_CS"))
    except RuntimeError:
        return None


def _vertical_unit(band, srs, opts, label, log):
    """(facteur vers le mètre, origine de l'information). Ordre : --vertical-unit, unité de la bande/variable,
    composante verticale du CRS composé ; sinon mètres SUPPOSÉS, avec avertissement (jamais en silence)."""
    declared = (band.GetUnitType() or "").strip()
    if opts.vertical_unit is not None:
        factor = VERTICAL_UNIT_CHOICES[opts.vertical_unit]
        known = _unit_factor(declared) if declared else None
        if known is not None and abs(known - factor) > 1e-9:
            log(f"{label} : avertissement : unité verticale déclarée « {declared} » remplacée par "
                f"--vertical-unit {opts.vertical_unit}")
        return factor, f"--vertical-unit {opts.vertical_unit}"
    if declared:
        factor = _unit_factor(declared)
        if factor is None:
            raise InputError(f"{label} : l'unité de la grille « {declared} » n'est pas une unité de longueur : "
                             f"ce n'est probablement pas une altitude. Si c'en est une, précise son unité avec "
                             f"--vertical-unit m|ft|us-ft.")
        return factor, f"déclarée « {declared} »"
    factor = _crs_vertical_factor(srs)
    if factor is not None:
        return factor, "composante verticale du CRS"
    log(f"{label} : avertissement : unité verticale non déclarée, supposée en MÈTRES "
        f"(précise-la avec --vertical-unit m|ft|us-ft si ce n'est pas le cas)")
    return 1.0, "non déclarée, supposée en mètres"


def _cf_geographic(ds):
    """Vrai si le fichier déclare des axes CF en degrés nord et est (variables `lat#units`...), sans CRS."""
    units = {(v or "").strip().lower().replace(" ", "") for k, v in ds.GetMetadata().items() if k.endswith("#units")}
    return bool(units & _LATITUDE_UNITS) and bool(units & _LONGITUDE_UNITS)


def _srs_from_user(text):
    srs = osr.SpatialReference()
    try:
        srs.SetFromUserInput(text)
    except RuntimeError as e:
        raise InputError(f"--source-crs : CRS non reconnu « {text} » ({_first_line(e)[:100]}). "
                         f"Exemples : EPSG:4326, EPSG:2056, EPSG:32632")
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return srs


def _check_extent_vs_crs(ds, gt, srs, label, log):
    """Garde-fou « degrés ou mètres » : une emprise incompatible avec le type du CRS est refusée (géographique)
    ou signalée (projeté)."""
    xs = [gt[0], gt[0] + ds.RasterXSize * gt[1] + ds.RasterYSize * gt[2]]
    ys = [gt[3], gt[3] + ds.RasterXSize * gt[4] + ds.RasterYSize * gt[5]]
    if srs.IsGeographic():
        if max(abs(v) for v in xs) > 360.0 + 1e-6 or max(abs(v) for v in ys) > 90.0 + 1e-6:
            raise InputError(f"{label} : le CRS est géographique (degrés) mais l'emprise "
                             f"[{min(xs):g}..{max(xs):g}] x [{min(ys):g}..{max(ys):g}] est impossible en degrés : "
                             f"coordonnées en réalité projetées (mètres) ? Corrige le CRS (--source-crs).")
    elif srs.IsProjected():
        if max(abs(v) for v in xs) <= 360.0 and max(abs(v) for v in ys) <= 90.0 and abs(gt[1]) < 0.5:
            log(f"{label} : avertissement : le CRS est projeté ({srs.GetLinearUnitsName()}) mais l'emprise et le "
                f"pixel ({gt[1]:g}) ont l'allure de degrés : vérifie le CRS (--source-crs)")


def _georef(ds, opts, label, log):
    """Retourne (srs, wkt_à_imposer ou None). Un fichier sans transformation affine ou sans CRS est refusé (ou,
    sans CRS, complété par --source-crs) : le programme ne fabrique pas de géoréférencement."""
    gt = ds.GetGeoTransform(can_return_null=True)
    if gt is None or tuple(gt) == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0):
        raise InputError(f"{label} : pas de géoréférencement (transformation affine absente"
                         + (" ; géolocalisation par tableaux ou points d'appui non prise en charge" if
                            ds.GetGCPCount() or ds.GetMetadata("GEOLOCATION") else "")
                         + "). Sans position, l'échelle du relief est inconnue : géoréférence le fichier "
                           "(gdal_edit.py, gdalwarp) avant de l'utiliser.")
    if gt[1] == 0.0 or gt[5] == 0.0:
        raise InputError(f"{label} : pixel de taille nulle ({gt[1]:g} x {gt[5]:g})")
    wkt = ds.GetProjection()
    override = None
    if wkt:
        if opts.source_crs is not None:
            raise InputError(f"--source-crs : {label} déclare déjà un CRS ; l'option ne sert qu'aux fichiers qui "
                             f"n'en ont pas (un CRS déclaré n'est jamais remplacé en silence).")
        srs = osr.SpatialReference()
        srs.ImportFromWkt(wkt)
    elif opts.source_crs is not None:
        srs = _srs_from_user(opts.source_crs)
        override = srs.ExportToWkt()
        log(f"{label} : CRS imposé par --source-crs : {srs.GetName()}")
    elif _cf_geographic(ds):
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        override = srs.ExportToWkt()
        log(f"{label} : avertissement : axes en degrés (lat/lon CF) mais datum non déclaré : WGS84 supposé "
            f"(écart possible de ~1-2 m avec NAD83/ETRS89, négligeable pour l'impression ; --source-crs pour "
            f"préciser)")
    else:
        raise InputError(f"{label} : aucun système de coordonnées (CRS). Impossible de savoir si les "
                         f"coordonnées sont des degrés ou des mètres, ni dans quel datum. Indique-le avec "
                         f"--source-crs (ex. --source-crs EPSG:2056), ajoute le fichier .prj voisin (ASC), ou "
                         f"géoréférence le fichier.")
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    _check_extent_vs_crs(ds, gt, srs, label, log)
    return srs, override


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

def _float_type(dtype):
    """Type de travail sans perte : Float64 pour les sources à 32 bits entiers ou plus, Float32 sinon."""
    return gdal.GDT_Float64 if dtype in (gdal.GDT_Int32, gdal.GDT_UInt32, gdal.GDT_Float64,
                                         getattr(gdal, "GDT_Int64", -1), getattr(gdal, "GDT_UInt64", -1)) \
        else gdal.GDT_Float32


def _nodata_arg(band):
    nd = band.GetNoDataValue()
    if nd is None:
        return "nan" if band.DataType in (gdal.GDT_Float32, gdal.GDT_Float64) else None
    return "nan" if math.isnan(nd) else repr(float(nd))


def _build_vrt(ds, band_index, scale, offset, out_wkt):
    """VRT en mémoire : une bande en virgule flottante, valeur = brut * scale + offset, nodata conservé."""
    band = ds.GetRasterBand(band_index)
    kw = {"format": "VRT", "bandList": [band_index], "outputType": _float_type(band.DataType)}
    if scale != 1.0 or offset != 0.0:
        kw["scaleParams"] = [[0.0, 1.0, offset, offset + scale]]
    nodata = _nodata_arg(band)
    if nodata is not None:
        kw["noData"] = nodata
    if out_wkt:
        kw["outputSRS"] = out_wkt
    name = f"/vsimem/dem2stl_{uuid.uuid4().hex}.vrt"
    out = gdal.Translate(name, ds, **kw)
    out = None  # noqa: F841 - ferme et écrit le VRT
    return name


def _prepare_raster(path, opts, log):
    label = os.path.basename(path) if not _is_virtual(path) else path
    if not _is_virtual(path):
        path = os.path.abspath(path)   # le VRT en mémoire référence le fichier : chemin indépendant du cwd
    ds, name, layer, fmt = _open(path, opts)
    band_index = _select_band(ds, opts, path)
    band = ds.GetRasterBand(band_index)
    if band.GetColorTable() is not None and opts.band is None:
        raise InputError(f"{label} : image indexée (palette de couleurs), pas une grille d'altitudes")
    srs, crs_override = _georef(ds, opts, label, log)
    factor, unit_origin = _vertical_unit(band, srs, opts, label, log)

    # échelle / décalage déclarés par le format (unscaled = brut * échelle + décalage), puis unité, puis signe
    declared_scale = band.GetScale()
    declared_offset = band.GetOffset()
    scale = 1.0 if declared_scale is None else float(declared_scale)
    offset = 0.0 if declared_offset is None else float(declared_offset)
    if not (math.isfinite(scale) and scale != 0.0 and math.isfinite(offset)):
        raise InputError(f"{label} : facteur d'échelle ({declared_scale}) ou décalage ({declared_offset}) invalide")
    sign = 1.0
    if (ds.GetMetadataItem(f"{band.GetMetadata().get('NETCDF_VARNAME', '')}#positive") or "").lower() == "down" \
            or (band.GetMetadata().get("positive") or "").lower() == "down":
        sign = -1.0
        log(f"{label} : variable orientée vers le bas (positive = down) : signe inversé pour obtenir des altitudes")
    scale, offset = scale * factor * sign, offset * factor * sign

    out_wkt = crs_override
    if srs.IsCompound():
        horizontal = srs.Clone()
        horizontal.StripVertical()
        out_wkt = horizontal.ExportToWkt()
        log(f"{label} : CRS composé ({srs.GetName()}) : composante verticale ignorée, "
            f"AUCUN changement de datum vertical appliqué")

    needs_vrt = (band_index != 1 or ds.RasterCount > 1 or scale != 1.0 or offset != 0.0 or out_wkt is not None)
    result = _build_vrt(_gdal_open(name), band_index, scale, offset, out_wkt) if needs_vrt else name

    bits = [fmt + (f" (table/variable {layer})" if layer else ""),
            f"{ds.RasterXSize} x {ds.RasterYSize}", gdal.GetDataTypeName(band.DataType),
            "CRS " + (srs.GetName() or "?"), f"altitudes en mètres ({unit_origin})" if factor == 1.0
            else f"unité verticale convertie en mètres x{factor:.6g} ({unit_origin})"]
    if ds.RasterCount > 1:
        bits.append(f"bande {band_index}/{ds.RasterCount}")
    if declared_scale not in (None, 1.0) or declared_offset not in (None, 0.0):
        bits.append(f"valeur = brut x {declared_scale if declared_scale is not None else 1:g} "
                    f"+ {declared_offset if declared_offset is not None else 0:g} appliqué")
    log(f"entrée {label} : " + ", ".join(bits))
    if needs_vrt:
        log(f"  lecture normalisée par un VRT en mémoire (aucune copie des données)")
    return result


# ---------------------------------------------------------------------------
# SRTM : tuile isolée ou assemblage
# ---------------------------------------------------------------------------

def _tile_key(path):
    m = _HGT_NAME.match(os.path.basename(path))
    if not m:
        return None
    lat = int(m.group(2)) * (1 if m.group(1).upper() == "N" else -1)
    lon = int(m.group(4)) * (1 if m.group(3).upper() == "E" else -1)
    return lat, lon


def _mosaic_srtm(paths, names, log):
    """Plusieurs tuiles .hgt de même résolution -> UN seul VRT mosaïque en mémoire (aucune copie des données).

    Les tuiles voisines partagent une rangée/colonne d'échantillons : elles se recouvrent donc d'un pixel. Si on les
    donnait séparément à gdal.Warp, chaque cellule de sortie à cheval sur ce recouvrement serait moyennée par
    chaque tuile sur sa seule part de la cellule, et la dernière tuile écraserait l'autre : un écart d'une pente de
    pixel apparaît le long de la couture. Vu comme une seule grille, le recouvrement porte des valeurs identiques et
    le résultat est celui d'un raster unique. Les vides (-32768) d'une tuile ne masquent pas une valeur de l'autre.
    Des résolutions mixtes (SRTM1 + SRTM3) ne sont PAS fusionnées : cela imposerait un ré-échantillonnage silencieux."""
    idx = [i for i, p in enumerate(paths) if os.path.splitext(p)[1].lower() == ".hgt"]
    if len(idx) < 2:
        return names
    sizes = set()
    for i in idx:
        sizes.add(os.path.getsize(paths[i]))
    if len(sizes) > 1:
        log("SRTM : avertissement : tuiles de résolutions différentes (SRTM1 et SRTM3) ; elles ne sont pas fusionnées "
            "(pas de ré-échantillonnage caché), la couture entre elles peut donc différer d'un pixel")
        return names
    vrt = f"/vsimem/dem2stl_{uuid.uuid4().hex}.vrt"
    try:
        out = gdal.BuildVRT(vrt, [names[i] for i in idx])
    except RuntimeError as e:
        raise InputError(f"assemblage des tuiles SRTM impossible : {_first_line(e)}")
    if out is None:
        raise InputError("assemblage des tuiles SRTM impossible : GDAL n'a pas pu construire la mosaïque")
    out = None  # noqa: F841 - ferme et écrit le VRT
    kept = [n for i, n in enumerate(names) if i not in idx]
    return [vrt] + kept


def _report_srtm_assembly(paths, log):
    tiles = [(p, _tile_key(p)) for p in paths if os.path.splitext(p)[1].lower() == ".hgt"]
    if not tiles:
        return
    keys = sorted({k for _, k in tiles if k})
    if len(keys) == 1:
        log(f"SRTM : tuile isolée {os.path.basename(tiles[0][0])} (1° x 1°) ; tout ce qui sort de la tuile est "
            f"sans donnée : fournis les tuiles voisines si le disque la déborde")
        return
    lats = [k[0] for k in keys]
    lons = [k[1] for k in keys]
    absent = [(la, lo) for la in range(min(lats), max(lats) + 1) for lo in range(min(lons), max(lons) + 1)
              if (la, lo) not in keys]
    log(f"SRTM : assemblage de {len(keys)} tuiles ; les échantillons de bordure partagés entre tuiles voisines "
        f"(1201/3601 points pour 1°) se recouvrent exactement, sans décalage ni couture")
    if absent:
        names = ", ".join(f"{'N' if la >= 0 else 'S'}{abs(la):02d}{'E' if lo >= 0 else 'W'}{abs(lo):03d}"
                          for la, lo in absent[:6])
        log(f"SRTM : avertissement : tuile(s) manquante(s) dans le rectangle couvert ({names}"
            f"{'...' if len(absent) > 6 else ''}) : leurs cellules seront sans donnée")
