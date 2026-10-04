# SPDX-License-Identifier: GPL-3.0-or-later
"""
geocode.py
----------
Résolution d'un nom de lieu (« Mont Fuji », « Everest »...) en coordonnées WGS84
via Nominatim (OpenStreetMap). Utilisé par `--place`.

Politique d'usage du serveur public Nominatim (rappelée par les pages officielles
et plusieurs sources tierces) : 1 requête/seconde au maximum, User-Agent qui
identifie l'application (surcharge possible avec DEM2STL_USER_AGENT, idéalement
avec un contact), pas de géocodage en masse, mise en cache des résultats, mention
« © OpenStreetMap contributors » (licence ODbL). Ce module n'émet qu'une requête par
lieu et met le résultat en cache sur disque.

On prend le PREMIER résultat selon l'ordre de pertinence de Nominatim, sans
heuristique supplémentaire (une préférence pour les sommets choisirait de mauvais
lieux pour des noms ambigus). Les autres candidats sont affichés ; --place-pick N
permet de choisir le N-ième.
"""
import hashlib
import json
import os
import time
import urllib.parse

import sources

NOMINATIM_URL = os.environ.get("DEM2STL_GEOCODER_URL", "https://nominatim.openstreetmap.org/search")
MAX_RESULTS = 5
MIN_INTERVAL_S = 1.0
ADMIN_TYPES = {"administrative", "state", "region", "county", "country", "province", "municipality",
               "city", "town", "village", "suburb", "island", "archipelago", "continent"}

_last_request = [0.0]


class GeocodeError(RuntimeError):
    """Lieu introuvable ou service de géocodage indisponible."""


def _cache_path(cache_dir, query):
    key = hashlib.sha1(" ".join(query.split()).casefold().encode("utf-8")).hexdigest()[:16]
    return os.path.join(cache_dir, "geocode", f"{key}.json")


def search(query, cache_dir, log=print):
    """Liste de candidats [{name, lat, lon, category, type}] (ordre Nominatim)."""
    query = query.strip()
    if not query:
        raise GeocodeError("nom de lieu vide")
    path = _cache_path(cache_dir, query)
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    wait = MIN_INTERVAL_S - (time.monotonic() - _last_request[0])
    if wait > 0:
        time.sleep(wait)
    url = NOMINATIM_URL + "?" + urllib.parse.urlencode({
        "q": query, "format": "jsonv2", "limit": MAX_RESULTS, "dedupe": 1, "accept-language": "fr,en"})
    try:
        data = sources._http_get(url)
    except sources.SourceUnavailable as e:
        raise GeocodeError(f"géocodage indisponible ({e})")
    finally:
        _last_request[0] = time.monotonic()
    if data is None:
        raise GeocodeError("géocodage : réponse 404")
    try:
        raw = json.loads(data.decode("utf-8"))
        results = [{"name": r["display_name"], "lat": float(r["lat"]), "lon": float(r["lon"]),
                    "category": r.get("category", ""), "type": r.get("type", "")} for r in raw]
    except (ValueError, KeyError, TypeError) as e:
        raise GeocodeError(f"géocodage : réponse invalide ({e})")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False)
    return results


def resolve(query, pick, cache_dir, log=print):
    """(lat, lon) du pick-ième résultat (1 = premier) ; journalise les candidats."""
    results = search(query, cache_dir, log)
    if not results:
        raise GeocodeError(f"aucun résultat pour « {query} » (essaie un nom plus précis ou --lat/--lon)")
    if not 1 <= pick <= len(results):
        raise GeocodeError(f"--place-pick {pick} hors limites (1 à {len(results)} résultats)")
    chosen = results[pick - 1]
    log(f"Lieu « {query} » -> {chosen['name']} ({chosen['category']}/{chosen['type']}) "
        f"lat={chosen['lat']:.5f} lon={chosen['lon']:.5f} [résultat {pick}/{len(results)}]")
    if chosen["type"] in ADMIN_TYPES:
        log("Lieu : ATTENTION, c'est le centre d'une zone administrative ou d'une île, pas forcément le "
            "point voulu (sommet, site) ; précise le nom ou utilise --lat/--lon.")
    if len(results) > 1:
        for i, r in enumerate(results, 1):
            if i != pick:
                log(f"  autre candidat {i} : {r['name']} ({r['category']}/{r['type']}) "
                    f"lat={r['lat']:.5f} lon={r['lon']:.5f}")
        log("  (choisis un autre candidat avec --place-pick N)")
    log("Géocodage : données © OpenStreetMap contributors (ODbL), via Nominatim.")
    return chosen["lat"], chosen["lon"]
