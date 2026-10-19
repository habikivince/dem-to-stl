# Audit comparatif dem-to-stl — Rapport 4 : lot 3 (sources nationales et données réelles)

Date : 10 octobre 2026. **Aucun fichier de votre projet n'a été modifié.** Toutes les mesures ont été faites sur des copies de travail, avec des scripts écrits pour l'occasion (je peux vous les fournir).

Convention : **[Fait]** constaté dans le code ou mesuré ; **[Déduction]** ; **[Hypothèse — à vérifier par benchmark]** ; **[Opinion]**.

---

## 0. Résumé

Les 13 dépôts de ce lot n'apportent **aucune méthode de maillage** : ce sont des chaînes de récupération de données. Leur intérêt est ailleurs : en les comparant à `sources.py` et `raster.py`, j'ai trouvé **quatre faiblesses réelles de `dem-to-stl`, toutes vérifiées par un test** (pas seulement par lecture) :

1. **`average` se comporte comme `nearest` en suréchantillonnage** (vérifié avec GDAL 3.12.2). Avec un DEM plus grossier que la grille d'impression, la surface est faite de terrasses : de 25 à 77 % des cellules ont une pente exactement nulle. L'erreur de pente est à peu près divisée par deux avec `bilinear` ou `cubic` sur un vrai DEM. **C'est le seul point du lot 3 qui améliore la qualité du STL à partir du même DEM.** Classe B.
2. **Un téléchargement tronqué lève une exception brute non rattrapée**, et **une page d'erreur reçue en HTTP 200 est enregistrée dans le cache puis réutilisée** sans contrôle. Classe B, correction faible.
3. **GSI : à partir d'un pixel de travail ≲ 3 m (zooms 16 et 17), seule la couche DEM1A est interrogée.** Si elle n'a pas la tuile, aucune couche plus grossière n'est essayée. Classe B.
4. **Aucun contrôle de débit ni `Retry-After`** (le dépôt `Geo3Dprint` en a un, pour un fournisseur dont il cite l'usage équitable). Classe C.

Deux idées sont à écarter ou à différer (fusion planaire avec tolérance, manifeste JSON) et une liste de sources candidates est donnée en section 3, avec leur niveau de vérification.

## 1. Dépôts du lot 3 : niveau de lecture et verdict

« Lu » = code pertinent lu en détail ; « partiel » = fonctions clés lues ; « survol » = README, signatures, recherches ciblées. **Ce qui n'est pas lu n'est pas couvert par mes conclusions.**

| Dépôt | Licence (détectée) | Lecture | Apport | Verdict |
|---|---|---|---|---|
| **nrw-terrain-print** (fkhr79) | **PolyForm Noncommercial 1.0.0** | `wcs.py` et `limits.py` lus ; `mesh.py` (913 l.), `boundary.py` (807 l.), `lod2.py` (1 913 l.) non lus | Cache validé (TIFF + SHA-256), limites de réponse ; sources NRW et Basse-Saxe | **Fiche 19, 24** (idées seulement, licence non commerciale) |
| **Geo3Dprint** (aareblau) | **aucune** | `http_utils.py` lu ; `terrain_pipeline.py` (1 940 l.) : liste des fonctions + début du fusionneur planaire | Débit poli, `Retry-After`, fusion planaire à tolérance | **Fiche 21, 22** (idées seulement) |
| **jp_geo_generator** | MIT | `main.py` lu à moitié (120/263 lignes) | Rien (tuiles texte, `np.max`, `pyfqmr`) | Écarter |
| **japan-3d-terrain** | MIT | README + signatures | Fusion automatique DEM5A/5B/5C + DEM10B, découpe par limite administrative | Fiche 20 (inspiration) |
| **dem2stl-relief** | **aucune** | `fetch_data.py` : lignes autour du `Warp` ; `mesh.py` non lu | **Rééchantillonnage selon le sens** | **Fiche 18** (idée seulement) |
| **portland-topo-map** | docs : Unlicense/CC0 ; **code : licence non détectée** | recherches ciblées | Tuiles, découpe à la limite de la ville, sources USGS/OpenTopography | Écarter |
| **geo-printing** (stretchyboy) | **aucune** | README + recherches | OS Terrain 50 (UK) via l'API OS Data Hub | Fiche 24 |
| **DEMtoPrint** (Spanholz) | MIT | README | MNT 1 m de Saxe, XYZ → GeoTIFF | Fiche 24 |
| **TerrainToSTL** (Unity) | **aucune** | recherche TNMAccess | NED 1/3″ seulement : moins riche que votre `fetch_usgs` | Écarter |
| **PeakPrinter** | GPL-3 | survol | Copernicus via CDSE | Écarter (vous avez Copernicus) |
| **elevation-to-stl** (katie-jones) | MIT | survol | SRTM 1″ sur AWS | Écarter |
| **convert-dem-to-stl** | MIT | survol | Conteneur Docker | Écarter |
| **turkeytopgraphyDEM2STL** | MIT | survol | OpenTopography ; `nan_to_num` → 0 | Écarter |

## 2. Fiches comparatives

### Fiche 18 — Rééchantillonnage selon le sens (réduction ou agrandissement) · Classe B

**Projet concurrent.** *dem2stl-relief* ([Fait] commentaire dans `fetch_data.py`) : « downsampling → average ; upsampling (print finer than the DEM) → cubic, avoids blocky DEM pixels ».

**Mon projet.** [Fait] `raster.py` utilise `resampleAlg="average"` pour la fenêtre principale, quel que soit le rapport entre le pixel de travail et le pixel de la source. Seul le complément Copernicus choisit `bilinear` sous 25 m de pixel. [Fait] Aucun traitement du suréchantillonnage dans `raster.py`, `circle_dem_to_stl.py` ni le README.

**Mesures.**

1. **[Fait, mesuré avec GDAL 3.12.2 via `rasterio`]** Source de 128×128 pixels de pas 8, grille d'impression au pas 1 (facteur 8) :

| Algorithme | Cellules à pente nulle | Bloc 8×8 constant |
|---|---|---|
| `average` | **56,5 %** | oui |
| `nearest` | 56,5 % | oui |
| `bilinear` | 0 % | non |
| `cubic` | 0 % | non |

`average` et `nearest` donnent ici **le même résultat**.

2. **[Mesure sur terrain synthétique]** Part de cellules à pente exactement nulle avec `nearest`/`average` : 25 % (facteur 4), 56 % (×8), 77 % (×16). L'erreur de hauteur (RMSE) est presque identique entre les trois méthodes (écarts < 4 %) : les terrasses ne sont **pas** une erreur d'altitude.

3. **[Mesure sur SRTM3 réel lissé]** Erreur de **pente** (RMSE du gradient relatif au gradient vrai) :

| Facteur | `nearest`/`average` | `bilinear` | `cubic` |
|---|---|---|---|
| ×4 | 0,930 | 0,471 | 0,409 |
| ×8 | 1,299 | 0,726 | 0,682 |

Sur un terrain synthétique très fin (fBm), le gain est plus faible (0,935 → 0,863 à ×4), car le détail fin n'est pas récupérable.

**Différence.** [Déduction] Pour une source grossière suréchantillonnée, `average` produit des paliers plats séparés de marches, avec une pente fausse (nulle dans les paliers, concentrée sur les bords). `bilinear` ou `cubic` donnent une surface continue et divisent par environ 2 l'erreur de pente sur ce DEM réel.

**Quand cela arrive chez vous.** [Fait] Les sources nationales « ne récupèrent que la résolution nécessaire » (README), donc elles réduisent plutôt qu'elles n'agrandissent. [Déduction] Le suréchantillonnage apparaît avec : `--input` d'un DEM global grossier (Copernicus 30 m, SRTM, GEBCO, ETOPO) sur un petit rayon, le repli Copernicus de `--source auto`, USGS 1/3″ (~10 m) en repli, GSI à zoom plafonné.

**Gain.** Important pour ces cas, nul pour les autres. [Opinion] C'est la seule amélioration de ce lot qui change le **rendu de la surface** pour un DEM donné.

**Difficulté.** Faible : comparer le pixel de travail à la taille du pixel de la source exprimée en mètres (la source peut être en degrés, ce qui demande la latitude) et choisir l'algorithme. **Risque.** Faible à moyen :
- [Hypothèse — à vérifier] `cubic` peut dépasser localement (« ringing ») près de falaises ; `bilinear` ne dépasse jamais mais laisse des arêtes de pente (surface « facettée »). Choix à trancher sur vos DEM.
- Avec plusieurs sources de pas différents, la règle doit tenir compte de la plus grossière.
- Les résultats changent seulement dans le cas suréchantillonné ; les tests existants sur données synthétiques fines ne devraient pas bouger.

**Verdict : intégrer.** À proposer avec un réglage `--resample auto|average|bilinear|cubic`.

### Fiche 19 — Validation des téléchargements et du cache · Classe B

**Projet concurrent.** *nrw-terrain-print* ([Fait] `wcs.py`) : avant de mettre une tuile en cache, il contrôle le type de contenu, refuse une taille annoncée supérieure à une limite, écrit la réponse en flux avec plafond, vérifie les premiers octets TIFF, puis consigne l'URL, la date et l'empreinte SHA-256 dans un manifeste JSON ; à chaque lecture du cache il appelle `_validate_dgm_cache_hit` ([Déduction] revalidation de l'empreinte ; je n'ai pas lu le corps de la fonction). La clé de cache est dérivée de l'URL.

**Mon projet.** [Fait] `_download` écrit `resp.read()` entier dans un fichier `.part`, puis le renomme, **sans contrôle de contenu** ; si le fichier existe, il est renvoyé tel quel. Le chemin Kartverket définit des octets magiques TIFF (`TIFF_MAGICS`) pour reconnaître la réponse ; le chemin GSI écrit l'image dans le cache **avant** de la décoder.

**Mesures** (j'ai extrait vos fonctions `_http_get` et `_download` **telles quelles** et les ai exécutées contre un serveur local) :

1. **[Fait]** Réponse tronquée (500 octets reçus sur 1 000 annoncés) : exception **`http.client.IncompleteRead`**, qui n'est ni `SourceUnavailable` ni rattrapée par les `except` de `_http_get`. **Une seule tentative**, aucune nouvelle. Cause : `IncompleteRead` dérive de `HTTPException`, pas de `OSError`.
2. **[Fait]** Réponse HTTP 200 dont le corps est un rapport d'erreur XML (87 octets) : `_download` l'écrit dans le cache comme une tuile. Le deuxième appel, serveur arrêté, **renvoie ce fichier sans le revalider**.

**Gain.** Moyen : un téléchargement interrompu ou un service qui répond 200 avec une erreur (cela arrive avec les WCS) bloque ensuite l'utilisateur jusqu'à vidage manuel du cache. **Difficulté** faible : rattraper `HTTPException`, vérifier le type du fichier (octets magiques TIFF/PNG, ou ouverture GDAL) avant le renommage, contrôler la longueur annoncée. **Risque** faible.

**Verdict : intégrer** (contrôle de contenu + nouvelle tentative sur troncature). Le manifeste SHA-256 de *nrw-terrain-print* est **inutile** pour votre usage : le contrôle de contenu suffit.

### Fiche 20 — Repli de résolution par tuile (GSI) · Classe B

**Projet concurrent.** *japan-3d-terrain* ([Fait] README) fusionne automatiquement des mailles DEM5A/5B/5C (5 m) et DEM10B (10 m) selon ce qui existe.

**Mon projet.** [Fait] `_gsi_tile` essaie les couches dans l'ordre `dem1a_png` (zoom ≤ 17), `dem5a_png`, `dem5b_png` (≤ 15), `dem_png` (≤ 14), mais saute toute couche dont le zoom maximum est inférieur au zoom demandé.

**Mesure** (réseau simulé : toutes les tuiles absentes ; la logique de sélection de couches est la vôtre, GDAL remplacé par un stub) :

| Zoom demandé | Couches interrogées |
|---|---|
| 17 | `dem1a_png` seulement |
| 16 | `dem1a_png` seulement |
| 15 | `dem1a_png`, `dem5a_png`, `dem5b_png` |
| 14 | les quatre |

**[Fait]** Le zoom 16 correspond à un pixel de travail d'environ 2 à 3 m, le zoom 17 à ≤ 1,5 m (à 35,4° N). **[Hypothèse — à vérifier]** La couverture DEM1A serait partielle ; ce point n'est pas vérifié. Si c'est exact, une demande de pixel fin hors zone DEM1A retourne « aucune tuile » alors que du 5 m existe, et `--source auto` se rabat sur Copernicus 30 m.

**Gain.** Moyen pour le Japon à haute résolution. **Difficulté** moyenne : demander la tuile à `min(z, zmax)` de chaque couche et la suréchantillonner dans la mosaïque (en lien direct avec la fiche 18). **Risque** faible. **Verdict : intégrer** après la fiche 18, avec test réel sur une zone hors DEM1A.

### Fiche 21 — Débit poli et `Retry-After` · Classe C

**Projet concurrent.** *Geo3Dprint* ([Fait] `http_utils.py`) : au plus 40 requêtes par minute vers `geo.admin.ch` (constante nommée « fair use »), verrou partagé, `Retry-After` interprété (secondes ou date), attente exponentielle jusqu'à 5 tentatives, agent utilisateur explicite.

**Mon projet.** [Fait] 3 tentatives, attente exponentielle plafonnée à 8 s, agent utilisateur configurable (`DEM2STL_USER_AGENT`) ; pas de limitation de débit ni lecture de `Retry-After`.

**Différence.** [Déduction] Pour une grande zone swisstopo (centaines de tuiles) ou des WCS découpés en jusqu'à 400 fenêtres, des requêtes en rafale risquent de déclencher des 429. [Hypothèse] Je n'ai pas observé de blocage réel ; les conditions d'usage équitable du fournisseur ne sont pas vérifiées par moi (je ne cite que la constante du dépôt).

**Gain** faible. **Difficulté** faible. **Verdict : étudier**, seulement si vous rencontrez des 429.

### Fiche 22 — Fusion planaire à tolérance (grille régulière) · Classe D

**Projet concurrent.** *Geo3Dprint* ([Fait]) fusionne des rectangles de cellules « planes » avec une tolérance (par défaut le plus grand de 10 % du pas horizontal et 0,025 % de l'amplitude verticale), détecte des parois raides, écrit un fond en éventail, et ([Déduction] d'après les noms `_write_adaptive_grid_conforming` et `_patch_ring_vertices`) raccorde les patchs sans jonction en T. Conçu pour un modèle de **surface** (toits, bâtiments).

**Mesure** (simulation par quadtree, **optimiste** : 2 triangles par bloc, sans surcoût de raccord ; part des triangles de la grille complète ; SRTM3 lissé, 513×513) :

| Lissage | Tolérance 0,10 % | 0,25 % | 0,50 % |
|---|---|---|---|
| σ=1 px | 86–90 % | 48–56 % | 24–29 % |
| σ=2 px | 60–71 % | 27–34 % | 14–18 % |
| σ=4 px | 32–44 % | 14–19 % | 7–10 % |

**Lecture.** [Déduction] À même tolérance et mêmes données, c'est de l'ordre de 1,3 à 1,5 fois moins efficace que RTIN (rapport 3, fiche 12 : 44–53 % / 20–25 % / 10–13 % pour σ=2), et ce sans compter le surcoût de raccord. [Opinion] L'approche convient aux toits et aux plaines parfaitement plates, pas aux reliefs de montagne. **Verdict : ne pas intégrer** pour votre usage.

### Fiche 23 — Manifeste JSON de résultat · Classe D

**Projet concurrent.** *nrw-terrain-print* écrit, à côté des STL/3MF, un JSON avec entrées, paramètres, contrôles de maillage et nombre d'arêtes ouvertes.

**Mon projet.** [Fait] Aucune sortie JSON dans `circle_dem_to_stl.py` (journal affiché seulement).

**Verdict : étudier seulement si vous voulez archiver vos reliefs avec leurs paramètres.** [Opinion] Gain de traçabilité, aucun effet sur la qualité du STL.

### Fiche 24 — Sources candidates pour `--source`

**Statut de vérification : aucune de ces sources n'a été testée par moi.** Mon environnement n'a pas accès à ces services. Les informations viennent des dépôts du lot (README, notices tierces) et sont **déclaratives**.

| Source | Pays | Résolution | Accès décrit par le dépôt | Licence déclarée par le dépôt | Remarque |
|---|---|---|---|---|---|
| Geobasis NRW DGM1 | Allemagne (Rhénanie du Nord-Westphalie) | 1 m | WCS `https://www.wcs.nrw.de/geobasis/wcs_nw_dgm`, couverture `nw_dgm` | dl-de/zero-2.0 | Même type de service que votre Kartverket/IGN |
| LGLN Basse-Saxe DGM1 | Allemagne | 1 m | STAC + WCS | **CC BY 4.0** (attribution) | Le dépôt traite Basse-Saxe en « première extension » |
| DEM1 Saxe | Allemagne | 1 m | téléchargement de fichiers XYZ par le portail `geodaten.sachsen.de` | non précisée par le dépôt | Pas d'API décrite : téléchargement manuel |
| OS Terrain 50 | Royaume-Uni | 50 m | API OS Data Hub (`api.os.uk/downloads`), archive ASCII d'environ 2 Go | non précisée par le dépôt | Résolution grossière |
| Copernicus GLO-30 sur AWS | Monde | 30 m | tuiles COG `copernicus-dem-30m.s3.amazonaws.com` | non précisée par le dépôt | Vous avez déjà Copernicus |
| AWS Terrain Tiles (terrarium) | Monde | variable (non précisée par le dépôt) | `elevation-tiles-prod.s3.amazonaws.com` | non précisée par le dépôt | Sans clé |

[Opinion] Les seules vraies nouveautés pour vous sont NRW et la Basse-Saxe (1 m, accès WCS/STAC, licences ouvertes selon le dépôt). Leur intérêt dépend des reliefs que vous prévoyez en Allemagne. Le dépôt *nrw-terrain-print* étant sous licence non commerciale, je recommande de lire les **documentations officielles des services** et non son code.

## 3. Constats transversaux

**Nodata chez les autres.**

| Dépôt | Valeur ou traitement |
|---|---|
| dem2stl-relief | −9999 ; **nodata → 0 et valeurs négatives écrêtées** (pas de bathymétrie) |
| geo-printing | −9999 → NaN ; puis `nan_to_num` → 0 |
| turkeytopgraphyDEM2STL | `nan_to_num` → **0** (trous remplis par le niveau de la mer, sans signalement) |
| PeakPrinter | `nan_to_num` → 0, puis écarte les triangles contenant la valeur nodata (comportement exact non vérifié) |
| jp_geo_generator | « e » → −1 |
| DEM.Net (rapport 3) | altitude 0 = invalide |
| **dem-to-stl** | −9999 ou valeur de la source, comblement GDAL, masque, ouverture, repli Copernicus recalé en altitude |

[Opinion] Le traitement de vos trous est nettement le plus rigoureux : les dépôts du lot remplacent silencieusement les vides par 0, ce qui fabrique des fosses ou du faux niveau de la mer.

**Altitudes de référence.** [Déduction, d'après une recherche de mots-clés geoid/datum/EGM] Aucun autre dépôt du lot ne semble traiter la différence de référence verticale entre sources ; votre repli Copernicus la **recale sur la zone commune** avec la source principale. [Opinion] Point fort à conserver.

**Gros DEM.** [Fait] Aucun dépôt du lot ne traite des fichiers plus gros que la mémoire ; la plupart téléchargent des tuiles et les assemblent en mémoire. Votre `gdal.Warp` fenêtré (avec `/vsicurl/`) est plus robuste.

**Dépendances.** [Fait] *nrw-terrain-print* recommande Docker car GDAL et `rasterio` imposent des bibliothèques natives ; les autres utilisent `requests`, `numpy`, `rasterio`. Vous : GDAL seulement. Rien à ajouter.

**Licences.** MIT : jp_geo_generator, japan-3d-terrain, DEMtoPrint, elevation-to-stl, convert-dem-to-stl, turkeytopgraphyDEM2STL. GPL-3 : PeakPrinter. **Non commerciale** : nrw-terrain-print. **Aucune détectée** : Geo3Dprint, dem2stl-relief, geo-printing, TerrainToSTL, portland-topo-map (code). Les fiches 18, 19, 21 s'implémentent **sans reprise de code** (idées simples, tests faits sur votre propre code).

## 4. Mesures effectuées dans ce rapport

1. **`_http_get` et `_download`** : extraction des fonctions par analyse syntaxique de `sources.py`, exécution contre un serveur HTTP local (troncature ; HTTP 200 avec XML d'erreur).
2. **Rééchantillonnage** : `rasterio` 1.5.2 / GDAL 3.12.2 (`average`, `nearest`, `bilinear`, `cubic`), puis simulation numpy/scipy pour l'erreur de pente (fBm synthétique et SRTM3 réel).
3. **`_gsi_tile`** : appel direct avec réseau simulé et GDAL remplacé par un stub (seule la sélection de couches est évaluée).
4. **Fusion planaire** : quadtree (optimiste) sur SRTM3 lissé.

**Non mesuré** : comportement réel des services (aucun accès), couverture DEM1A de la GSI, impact sur vos DEM, suite de tests de `dem-to-stl` (GDAL absent chez moi).

## 5. Parties non lues du lot 3

`mesh.py`, `boundary.py` et `lod2.py` de nrw-terrain-print ; `terrain_pipeline.py` de Geo3Dprint (au-delà du début du fusionneur planaire) ; `mesh.py` et `pipeline.py` de dem2stl-relief ; la moitié de `main.py` de jp_geo_generator ; la chaîne complète de portland-topo-map.
