# Audit comparatif dem-to-stl — Rapport 6 : les 48 dépôts restants

Date : 11 octobre 2026. **Aucun fichier de votre projet n'a été modifié par ce rapport.**

Convention : **[Fait]** constaté dans le code ou mesuré ; **[Déduction]** ; **[Hypothèse — à vérifier par benchmark]** ; **[Opinion]**.

## 0. Résumé

1. **Aucun de ces 48 dépôts n'apporte une méthode de maillage ou de traitement de données supérieure à la vôtre.** Plusieurs **confirment** des constats précédents : une deuxième implémentation indépendante de la couverture gloutonne par rectangles (`ginoclement/terrain`), le fond par `earcut` (`MapMaker`, `Heightmap2STL`, `Generate-Terrain-with-GPX`), le TIN Martini (`MapMaker`), le `INTER_CUBIC` pour agrandir une carte de hauteurs (`Heightmap2STL`), la boucle de pincements jusqu'à stabilité (`ginoclement/terrain`).
2. **Une erreur de ma part est corrigée** : **5 dépôts avaient été classés « sans code exploitable » à tort**, car leurs archives (téléchargées par moi) avaient été tronquées à 60 Mo par mon plafond de taille. Je les ai clonés et relus (section 2). Aucune autre archive n'était invalide (vérifié sur les 93).
3. **Une mesure nouvelle** compare quatre façons de trianguler le fond d'un disque (section 3). Elle a conduit à remplacer ma première proposition (rectangles gloutons) par un quadtree équilibré dans le prototype livré (rapport 7).
4. Deux autres dépôts n'ont pas de code dans l'archive uploadée (README seul) : ils sont **non évaluables**, pas « sans intérêt ».

## 1. Les 48 dépôts

« lu » = code ou documentation du mécanisme concerné lu ; « survol » = README et signaux de code (recherches ciblées) seulement. **Ce qui n'est que survolé n'est pas couvert par mes conclusions** au-delà de la catégorie.

| # | Dépôt | Langages (LOC) | Licence détectée | Lecture | Constat | Verdict |
|---|---|---|---|---|---|---|
| 1 | Wenshuishi0528/topotile-studio | Python 12012, HTML 2927, Shell 9 (14948) | CC-BY (LICENSE) | lu (liste des fonctions de mesh_repair.py et printability.py) | Modèles de villes et de terrains depuis OSM ; réparation avec trimesh, contrôles d'imprimabilité. Licence **CC BY-NC-SA 4.0** (détectée à tort comme « CC-BY » au rapport 1). | E — hors cœur, licence non commerciale |
| 2 | gora-print-router-main | Python 13739 (13739) | AGPL (LICENSE) | survol (README, signaux) | Trace GPX transformée en modèle imprimable (routes, découpe en tuiles). Application complète en 13 700 lignes. | E — hors cœur ; AGPL |
| 3 | agrignon/MapMaker | TS 12412, HTML 13 (12425) | ABSENTE | lu (en-têtes de `solid.ts`, `terrain.ts`) | Application web : terrain par **Martini RTIN** (maillage adaptatif à erreur bornée), parois sur les bords, **fond triangulé par earcut sur le contour**. Confirme que le TIN borné est utilisé pour un solide imprimable (rapport 3, fiche 12). README : « All rights reserved ». | D — idée déjà retenue ; ne rien copier |
| 4 | ginoclement/terrain | JS 5526, HTML 376 (5902) | ABSENTE | lu (`mesh.js` : en-tête, nettoyage des pincements, fusion plate) | Application web (formes libres, cercle, lettres). **Deuxième implémentation indépendante** de la couverture gloutonne par rectangles + éventails (« optimizeFlat », activé par défaut) et de la suppression des cellules en pincement **en boucle jusqu'à stabilité** (confirme la correction apportée). Fond possible en « heightmap inversée » (gravure de socle). | D — confirme les fiches 1 et 5 |
| 5 | terraprint-main | Python 3882, HTML 1082, Shell 406 (5375) | ABSENTE | survol | Deux outils : planificateur de missions (PWA) et générateur de terrain. | E — hors cœur |
| 6 | ModelRift/terrain-to-3d | TS 4413, JS 23, HTML 15 (4451) | ABSENTE | survol (signaux : opérations booléennes) | Éditeur web de terrains imprimables. | E — hors cœur |
| 7 | MattMullin/openmoon | JS 2247, Python 1786, HTML 214 (4247) | ABSENTE | survol | Navigateur de terrain lunaire (Cesium + FastAPI) avec export. | E — hors cœur |
| 8 | formscope/gps-terrain-stl | Python 3104, HTML 891 (3995) | ABSENTE | survol (signaux : Delaunay) | Trace GPX/IGC vers un modèle de terrain. | E — hors cœur |
| 9 | andrewsheerin/topographic-maps | Python 2386, JS 1282, HTML 12 (3680) | ABSENTE | survol | Application « Topo2STL » : zone dessinée, DEM, STL. | E — hors cœur |
| 10 | cjorgensenmd/dem-to-3d-print-blender | Python 3584 (3584) | MIT (LICENSE) | survol | Extension Blender (MIT) ; dépend de Blender. | E — hors cœur |
| 11 | hicatlab/DEM-Cropper | JS 2093, Python 826, HTML 258 (3177) | MIT (LICENSE.md) | survol | Application web de découpe de DEM (MIT). | E — hors cœur |
| 12 | chris17453/terraprint3d | Python 2759 (2759) | README.md: MIT | survol | Application Python de terrains imprimables (MIT d'après le README). | E — hors cœur |
| 13 | Lord0fBytes/gpx-terrain-studio | TS 2332, HTML 13 (2345) | ABSENTE | survol | Trace GPX, terrain encadré étanche (TS). | E — hors cœur |
| 14 | 3d-relief-map-main | JS 2042, HTML 153 (2195) | ABSENTE | survol | Rectangle choisi sur une carte, relief exporté (JS). | E — hors cœur |
| 15 | map-generator-main | Python 2046 (2046) | README.md: licensed under [CC BY-N | survol | URL Google Maps vers STL (licence CC d'après le README). | E — hors cœur |
| 16 | niyoseris/3droute | Rust 834, JS 690, HTML 268 (1792) | MIT (LICENSE) | survol | Visualiseur de trajets 3D (Rust + JS, MIT). | E — hors cœur |
| 17 | jkoether/Terrain_Trails | Python 1668 (1668) | GPL-3 (LICENSE) | survol (signaux : booléens) | Modèles de sentiers de course : traces OSM incrustées dans le terrain. | E — hors cœur ; GPL-3 |
| 18 | PaulZC/SRTM_to_STL | Python 1539 (1539) | MIT (LICENSE) | survol | Chaîne SRTM → CSV → STL (MIT). | E — redondant |
| 19 | ThatcherC/Terrain2STL | JS 536, C 407, HTML 274 (1299) | package.json: "license": "ISC | survol | Service web HGT → STL (JS + C). | E — redondant |
| 20 | JGGroot/ODTT | JS 1042, HTML 38 (1080) | ABSENTE | survol | Application web de terrains imprimables. | E — hors cœur |
| 21 | spelfrog/3d-aachen | Python 1004 (1004) | MIT (LICENSE) | survol | Convertisseur TIFF → STL, quartiers étiquetés (MIT). | E — redondant |
| 22 | alb0/3d-greenland | Notebook 440, Matlab 149 (589) | Apache-2.0 (LICENSE.txt) | survol | Notebook et notes de production d'un modèle du Groenland (Apache-2.0). | E — hors cœur |
| 23 | coelacant1/ElevationDataTo3DModel | Python 580 (580) | GPL-3 (LICENSE) | lu (grep ciblé) | Outils Python (pyvista, numpy-stl) ; décimation **par pas** (`decimate_array`), nodata sous l'eau mis à 0. | E — redondant |
| 24 | terrainprint-main | HTML 578 (578) | ABSENTE | survol | Page HTML unique de 578 lignes. | E — hors cœur |
| 25 | r-barnes/DEMto3D | Python 547 (547) | MIT (LICENSE) | lu (grep ciblé) | CDED → STL : **tuiles à ergots d'assemblage** (paramètres `hsep`, `anchorsize`, `tab_dep`), lissage gaussien, `multiprocessing` par tuile. | D — tuilage modulaire hors cœur |
| 26 | IdeaShakerLab/dem-to-stl | Python 542 (542) | MIT (LICENSE) | survol (même lignée que r-barnes) | Variante du précédent (MIT) avec le même générateur de facettes. | D — idem |
| 27 | billyauhk/dem2stl | C 364, C++ 151 (515) | ABSENTE | survol | Implémentation C/C++ (avertissement dans le README). | E — redondant |
| 28 | lopoc/dem2stl | Python 504 (504) | ABSENTE | lu (grep ciblé) | « GeoSTL » : lecture fenêtrée `rasterio.windows.from_bounds`, STL étanche ; équivalent de votre `gdal.Warp` fenêtré. | E — redondant |
| 29 | GeoffWilliams/dem2stl | Python 469 (469) | LGPL (LICENSE) | survol | Image de DEM (PIL) → STL (LGPL). | E — redondant |
| 30 | Hands3DLab/Generate-Terrain-with-GPX | Python 469 (469) | README.md: MIT | survol (grep : earcut, trimesh) | Trace GPX sur terrain ; **fond par `mapbox_earcut`**. | D — idée earcut mesurée (voir §3) |
| 31 | ScreamingHawk/geotiff-to-mesh | Python 454 (454) | MIT (LICENSE) | survol | Fusion de plusieurs GeoTIFF en un STL (numpy-stl, MIT). | E — redondant |
| 32 | 41pha1/Heightmap2STL | Python 440 (440) | ABSENTE | lu (grep ciblé) | Image de relief → STL : **fond par `mapbox_earcut`** sur le contour, **`cv.resize(..., INTER_CUBIC)`** pour agrandir la carte de hauteurs (même constat que la fiche 18), gradient de Sobel. | D — deux idées déjà retenues |
| 33 | cvr/dem2stl | Python 364 (364) | GPL-3 (LICENSE) | survol | Script GDAL/OGR de l'ère Python 2 (GPL-3). | E — obsolète |
| 34 | oliver-rew/topo_tool | Python 287, Shell 6 (293) | ABSENTE | survol | Convertisseur GeoTIFF → STL « basé en grande partie sur » un autre outil. | E — redondant |
| 35 | BrendanLuke15/DTM-2-STL | Python 203, Matlab 59 (262) | MIT (LICENSE) | survol | MNT planétaires (Astropedia) → STL via PIL, Matlab (MIT). | E — hors cœur |
| 36 | Usuy-Leon/Topo3D-world | Python 242 (242) | MIT (LICENSE) | survol | Reprojection `rasterio` puis STL (MIT). | E — redondant |
| 37 | aparshin/srtm-to-stl | JS 202 (202) | package.json: "license": "ISC | survol | Script Node de 202 lignes. | E — redondant |
| 38 | ewandennis/geotiff2stl | Python 162 (162) | ABSENTE | survol | GeoTIFF flottant → STL (scikit-image). | E — redondant |
| 39 | aligundogdu/TerraStl | TS 145, JS 12 (157) | README.md: licensed under the MIT  | survol | Application web locale (TS). | E — hors cœur |
| 40 | GeoTiff2STL-main | Python 92 (92) | ABSENTE | survol | Script de 92 lignes. | E — redondant |
| 41 | F-Reid-McLain/macon-3d-print | — (0) | ABSENTE | **relu après clonage** (README, liste des scripts) | Relief imprimable du centre de Macon (Géorgie) : terrain USGS, bâtiments OSM et Microsoft, routes gravées, 21 tuiles à l'échelle 1:12 500. | E — hors cœur (OSM, ville) |
| 42 | terrain-relief-konfigurator-main | — (0) | ABSENTE | archive réduite à un README (160 octets) | Générateur web de reliefs modulaires, aperçu et export 3MF (d'après le README seul). | Non évaluable (aucun code dans l'archive) |
| 43 | terrain_map_3dprinting-main | — (0) | ABSENTE | archive réduite à un README et un fichier de 2 Ko | « Carte de terrain imprimée en 3D avec les données SRTM ». | Non évaluable |
| 44 | UnintelligibleMaker/PyTerrainModeler | — (0) | ABSENTE | **relu après clonage** (README, fonctions) | Données MapZen/OpenTopography ou NOAA XYZ ; classes de faces et normales écrites à la main (MIT). | E — redondant |
| 45 | danny-kumpf/gpx-tracks-3d | — (0) | ABSENTE | **relu après clonage** (README, fonctions de `auto_cut.py`) | Grave une trace GPX dans une tuile STL en subdivisant les triangles proches de la trace. | E — hors cœur |
| 46 | furrysalamander/mini-map-maker | — (0) | ABSENTE | **relu après clonage** (README, fonctions) | Nuages de points LiDAR USGS (.laz) → DEM via LAStools → STL (MIT) ; même générateur de facettes que r-barnes. | E — hors cœur |
| 47 | marcinbator/asc-to-3d-model-blender-script | — (0) | ABSENTE | **relu après clonage** | Script Blender de 65 lignes pour fichiers ASC (MIT). | E — redondant |
| 48 | riggsd/ArcGIS-Pro-DEM-3D-Print-Toolbox | — (0) | GPL-3 (LICENSE.txt) | lu (liste des fonctions de `DEM3DPrint.pyt`, 2 522 lignes) | Boîte à outils ArcGIS Pro : DEM → STL/3MF, découpe par polygones, **couleurs par couches** (écrivain 3MF), **parois vectorisées** (`_ns_wall`, `_ew_wall`). Dépend d'ArcPy (propriétaire). | E — dépendance propriétaire |

## 2. Correction de mon triage

[Fait] Mon script de téléchargement limitait chaque archive à 60 Mo. Cinq dépôts dépassent cette taille (de 90 à 245 Mo : STL, images et données dans l'historique), leurs archives étaient tronquées et mon extraction les a silencieusement ignorées, d'où « 0 ligne de code » et la catégorie « sans code exploitable » (rapport 1). Je les ai clonés avec `git clone --depth 1` : `macon-3d-print`, `PyTerrainModeler`, `gpx-tracks-3d`, `mini-map-maker`, `asc-to-3d-model-blender-script`. Ils contiennent de 65 à 1 100 lignes de Python, relues en README et liste de fonctions : verdict **hors cœur ou redondant**, comme les autres. Le rapport 1 a été annoté.
[Déduction] Les 88 autres archives sont valides : j'ai ouvert chacune avec un lecteur ZIP. Le tableau du rapport 1 reste donc correct pour elles, hormis la licence de `topotile-studio` (CC BY-NC-SA 4.0, non « CC-BY »).

## 3. Mesure : quatre fonds pour un disque

Disque de 1 250 cellules de diamètre (1 223 300 cellules), pas 0,2 mm (cas Ø250 mm). Même disque pour les quatre lignes. Le « allongement » est le rapport L²/(2·aire) du pire triangle (équilatéral : 1,15 ; demi-carré : 2).

| Fond | Triangles du fond | Total (dessus + fond) | Allongement médian / max | Durée de construction |
|---|---|---|---|---|
| Grille (actuel) | 2 446 600 | 4 893 200 | 2,0 / 2,0 | — |
| Rectangles gloutons (`my3dmaps`, `ginoclement/terrain`) | 6 452 | 2 453 052 (−49,9 %) | 1 120 / **600 608** | 1,2 s |
| `earcut` sur le contour (`MapMaker`, `Heightmap2STL`, `Hands3DLab`) | **4 990** | 2 451 590 (−49,9 %) | 10 / **30 000 000** | 0,01 s |
| **Quadtree équilibré** (prototype livré) | 28 462 | 2 475 062 (−49,4 %) | **2,0 / 4,0** | 1,1 s |

[Fait] Mesures sur un disque ; le fond `earcut` a aussi été essayé sur un disque de 400 cellules : **solide valide** (0 arête ouverte, volume identique, 0 triangle d'aire exactement nulle, accepté par `admesh` et `prusa-slicer --info`), mais l'aire minimale d'un triangle y est de 9,5×10⁻⁸ mm² (5 triangles sous 10⁻⁶ mm²).
[Opinion] `earcut` ne fait gagner qu'environ 1,4 % de plus que le quadtree sur le total (le dessus domine), au prix de triangles quasi dégénérés. Le quadtree est le meilleur compromis ; `earcut` n'est pas retenu.

## 4. Idées et constats transversaux

- **Fonctions hors cœur** (aucune amélioration de la fidélité) : traces GPX gravées dans le terrain (`gora-print-router`, `formscope/gps-terrain-stl`, `Terrain_Trails`, `gpx-terrain-studio`, `Generate-Terrain-with-GPX`, `gpx-tracks-3d`, `3droute`) ; villes et OSM (`topotile-studio`, `macon-3d-print`, `map-generator`, `Topo2STL`, `terraprint`) ; tuiles à ergots d'assemblage (`DEMto3D` de r-barnes et ses dérivés) ; couleurs par couches en 3MF (boîte à outils ArcGIS) ; terrains planétaires (`openmoon`, `DTM-2-STL`).
- **Nodata** : [Fait] `coelacant1/ElevationDataTo3DModel` met le nodata sous l'eau à 0 ; `dem2stl-relief`, `geo-printing`, `turkeytopgraphyDEM2STL` aussi (rapport 4). Rien ne remplace votre traitement.
- **Gros fichiers** : [Fait] `lopoc/dem2stl` utilise la lecture fenêtrée de `rasterio` (équivalente à votre `gdal.Warp` fenêtré) ; `coelacant1` décime par pas (sans moyenne). Rien de supérieur.
- **Nuages de points** : [Fait] `mini-map-maker` produit un DEM à partir de fichiers `.laz` avec LAStools (logiciel non libre). [Opinion] Une entrée LAZ n'apporterait rien à vos sources actuelles, qui fournissent déjà des modèles numériques de terrain.
- **Parois vectorisées** : [Fait] la boîte à outils ArcGIS (`_ns_wall`, `_ew_wall`) confirme que les parois se construisent sans boucle Python (déjà retenu, rapport 2 fiche 3).
- **Licences** : parmi les 48, seules des licences MIT (PaulZC, spelfrog, r-barnes, IdeaShakerLab, ScreamingHawk, Usuy-Leon, BrendanLuke15, PyTerrainModeler, mini-map-maker, marcinbator, DEM-Cropper, cjorgensenmd, 3droute) permettraient un emprunt de code ; **aucune idée retenue n'en demande** (ce sont des idées déjà ré-implémentées ou des observations). `MapMaker` : « All rights reserved » ; `topotile-studio` : CC BY-NC-SA ; `gora-print-router` : AGPL ; 24 dépôts sans licence détectée (sur 48) : ne rien copier.

## 5. Réponses aux dix questions pour ce lot

1. **Meilleures idées algorithmiques** : aucune nouvelle ; confirmations seulement.
2. **Meilleurs STL** : aucun dépôt ne produit une surface plus juste que la vôtre.
3. **Gros DEM** : aucun mieux que `dem-to-stl`.
4. **Données imparfaites** : `dem-to-stl`.
5. **Architecture** : rien à reprendre.
6. **Meilleure idée** : le quadtree équilibré pour le fond (voir §3) est une conclusion de ce lot, pas une idée d'un dépôt.
7. **Meilleur rapport gain/coût** : inchangé (rapport 5).
8. **Plus difficile mais intéressante** : TIN borné (confirmé par `MapMaker`).
9. **Ne valent pas la complexité** : tuiles à ergots, GPX, OSM, couleurs par couches, LAZ, courbure terrestre.
10. **Remise en cause de l'architecture** : non.

## 6. Parties non lues

La quasi-totalité du code des 33 dépôts « survol » ; `topotile-studio` au-delà de deux fichiers ; `gora-print-router` ; `MapMaker` au-delà de l'en-tête de `solid.ts` et de `terrain.ts` ; `ginoclement/terrain` hors `mesh.js`. Je n'ai exécuté aucun de ces dépôts.