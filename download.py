# SPDX-License-Identifier: GPL-3.0-or-later
"""
download.py
-----------
Téléchargement intégré du Copernicus DEM GLO-30 (30 m) / GLO-90 (90 m)
depuis le compartiment S3 public géré par Sinergise/AWS Open Data —
accès HTTPS simple, sans authentification, un objet par tuile de 1x1°.

Convention de nommage confirmée sur https://copernicus-dem-30m.s3.amazonaws.com/readme.html :
    Copernicus_DSM_COG_<résolution en arc-secondes>_<NS><lat>_00_<EW><lon>_00_DEM/...DEM.tif
    résolution = 10 pour GLO-30, 30 pour GLO-90 ; lat/lon = coin sud-ouest de la tuile (partie entière).
"""
import math
import os
import shutil
import urllib.request
import urllib.error

BUCKETS = {30: "copernicus-dem-30m", 90: "copernicus-dem-90m"}
RESOLUTIONS_ARCSEC = {30: 10, 90: 30}
HTTP_TIMEOUT = 60


def tile_name(lat_deg, lon_deg, product_m=30):
    """Nom de la tuile (sans extension) couvrant le coin sud-ouest entier
    (lat_deg, lon_deg) — lat_deg/lon_deg doivent déjà être des entiers
    (utiliser tile_key_for_point pour un point quelconque)."""
    ns = f"N{lat_deg:02d}" if lat_deg >= 0 else f"S{abs(lat_deg):02d}"
    ew = f"E{lon_deg:03d}" if lon_deg >= 0 else f"W{abs(lon_deg):03d}"
    res = RESOLUTIONS_ARCSEC[product_m]
    return f"Copernicus_DSM_COG_{res}_{ns}_00_{ew}_00_DEM"


def tile_url(lat_deg, lon_deg, product_m=30):
    name = tile_name(lat_deg, lon_deg, product_m)
    bucket = BUCKETS[product_m]
    return f"https://{bucket}.s3.amazonaws.com/{name}/{name}.tif"


def tile_key_for_point(lat, lon):
    """Tuile entière (coin SW) contenant le point (lat, lon)."""
    return math.floor(lat), math.floor(lon)


def tiles_for_circle(center_lat, center_lon, radius_m, margin_m=2000.0):
    """Détermine l'ensemble des tuiles 1°x1° dont l'emprise recoupe le
    cercle (centre, rayon), avec une marge de sécurité. Renvoie une liste
    de (lat_deg, lon_deg) triée, sans doublon."""
    lat_span = (radius_m + margin_m) / 111_320.0
    lon_span = (radius_m + margin_m) / (111_320.0 * max(math.cos(math.radians(center_lat)), 0.01))

    lat_min = math.floor(center_lat - lat_span)
    lat_max = math.floor(center_lat + lat_span)
    lon_min = math.floor(center_lon - lon_span)
    lon_max = math.floor(center_lon + lon_span)

    return sorted({
        (la, lo)
        for la in range(lat_min, lat_max + 1)
        for lo in range(lon_min, lon_max + 1)
    })


def download_copernicus(center_lat, center_lon, radius_m, cache_dir, product_m=30, log=print):
    """Télécharge (ou réutilise depuis le cache) toutes les tuiles
    Copernicus DEM nécessaires pour couvrir le cercle demandé. Retourne
    la liste des chemins locaux des fichiers .tif obtenus.

    Une tuile qui n'existe pas dans le compartiment (HTTP 404 — cas des
    tuiles océaniques, absentes par convention, cf. readme.html : "ocean
    areas do not have tiles, there one can assume height values equal to
    zero") est ignorée sans faire échouer l'ensemble.
    """
    os.makedirs(cache_dir, exist_ok=True)
    needed = tiles_for_circle(center_lat, center_lon, radius_m)
    log(f"{len(needed)} tuile(s) Copernicus GLO-{product_m} nécessaire(s) pour ce rayon")

    paths = []
    failed = []
    for lat_deg, lon_deg in needed:
        name = tile_name(lat_deg, lon_deg, product_m)
        local_path = os.path.join(cache_dir, f"{name}.tif")
        if os.path.exists(local_path):
            log(f"  {name} : déjà en cache")
            paths.append(local_path)
            continue

        url = tile_url(lat_deg, lon_deg, product_m)
        part = local_path + ".part"
        try:
            log(f"  {name} : téléchargement...")
            req = urllib.request.Request(url, headers={"User-Agent": "dem-to-stl"})
            with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as resp, open(part, "wb") as out:
                shutil.copyfileobj(resp, out)
            os.replace(part, local_path)  # écriture atomique : jamais de .tif partiel en cache
            paths.append(local_path)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                log(f"  {name} : absente (probablement océan) — ignorée")
            else:
                log(f"  {name} : échec HTTP {e.code}")
                failed.append(name)
        except OSError as e:  # URLError, timeout, connexion coupée
            log(f"  {name} : échec réseau ({e})")
            failed.append(name)
        finally:
            if os.path.exists(part):
                os.remove(part)

    if failed:
        # Une tuile manquante pour une autre raison qu'un 404 laisserait un trou
        # silencieux dans le relief : on préfère échouer explicitement.
        raise RuntimeError("Téléchargement Copernicus incomplet, tuile(s) en échec : " + ", ".join(failed))
    if not paths:
        raise RuntimeError("Aucune tuile Copernicus n'a pu être obtenue pour cette zone.")
    return paths
