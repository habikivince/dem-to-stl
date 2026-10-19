# Audit comparatif dem-to-stl — Rapport 1 (état intermédiaire)

Date : 9 octobre 2026. **Aucun fichier de votre projet n'a été modifié.** Ce document est un état d'avancement : baseline, inventaire, triage automatisé et premières constatations. L'audit approfondi des dépôts n'est **pas** terminé.

Convention : **[Fait]** constaté dans le code ou mesuré ; **[Déduction]** ; **[Hypothèse — à vérifier par benchmark]** ; **[Opinion]**.

## 1. Baseline `dem-to-stl` (lue dans l'archive fournie, 9 fichiers de code du cœur)

- [Fait] Pipeline : `circle_dem_to_stl.py` (CLI, orchestration) → `raster.py` (GDAL : `gdal.Warp` vers un dataset MEM, fenêtre carrée autour du centre, rééchantillonnage `average` au pas `print_spacing/échelle`, reprojection UTM automatique, `gdal.FillNodata` maxSearchDist=50) → `mesh.py` (NumPy/SciPy) → `stl_io.py` / `mesh_export.py`.
- [Fait] Conséquence structurelle : la taille de la grille de travail dépend de `--diameter-mm / --print-spacing-mm`, **pas de la taille du DEM source**. La mémoire est bornée par le maillage de sortie ; la gestion des gros DEM sources est déléguée à GDAL (VRT, `/vsicurl/`, lecture partielle). Les techniques de tuilage/streaming d'autres dépôts ont donc peu de prise sur le goulot réel.
- [Fait] Maillage : surface supérieure = 2 triangles par cellule avec **diagonale fixe** ; **fond plat = même grille, même densité que le dessus** (`bot_tris` construit à partir des mêmes quads) ; parois verticales construites par **double boucle Python** sur les arêtes de bord (`build_mesh_from_arrays`) ; soupe de triangles float64 (N×3×3), pas de maillage indexé natif.
- [Fait] Robustesse du masque : ouverture morphologique 3×3 (`erode3x3`/`dilate3x3`), comblement des pincements diagonaux (`fix_diagonal_pinches`, **4 itérations max**), composante connexe principale en mode île (`--sea-level`), bord circulaire optionnel (`--smooth-wall`, projection radiale du contour).
- [Fait] Export : STL binaire écrit par blocs vectorisés (fichier `.part` puis renommage), 3MF (unité mm, écrit avec la bibliothèque standard), OBJ ; `index_mesh` fusionne les sommets par `lexsort` après arrondi float32.
- [Fait] Validation : `check_stl.py` relit le fichier indépendamment du code de maillage (dimensions, socle, relief, arêtes ouvertes / non-manifold / mal orientées) ; 127 fonctions de test sur 12 fichiers ; garde-fou RAM (≈220 octets/triangle, constantes mesurées sur un cas de 4,9 M de triangles) et `--dry-run`.
- [Fait] Dépendances : GDAL, NumPy, SciPy (+ pytest). Licence GPL-3.0-or-later. Seule la version Linux est utilisée pour l'estimation de RAM (`/proc/meminfo`).
- [Limite de cette baseline] Je n'ai pas exécuté la suite de tests ni de benchmark : GDAL n'est pas installé dans mon environnement (NumPy et SciPy le sont). Les chiffres de performance cités viennent de votre README, pas de mes mesures.

## 2. Inventaire des dépôts

- **93 dépôts uniques** analysés (91 GitHub + 2 Codeberg), alors que vous en annonciez 88 : la différence est à clarifier (votre liste contenait en plus 3 doublons — `fgebhart/mapa`, `suchanek/satprint`, `gora-print-router` — et `rayshader-master.zip` a été envoyé deux fois, octet pour octet identique).
- 74 dépôts GitHub téléchargés par moi (archive de la branche par défaut, API GitHub inutilisable : quota épuisé) ; 19 archives uploadées par vous (dont les 2 dépôts Codeberg).
- Extraction sans images ni données volumineuses ; les comptes de lignes (LOC) ne concernent que le code source.
- Hors cœur (appli web, GPX, OSM, ville, multi-couleur) : **28**
- Lot 1 — maillage / robustesse : **17**
- Redondant probable (petit script grille→STL) : **16**
- Lot 3 — sources nationales / données : **13**
- Lot 2 — gros DEM / autres écosystèmes / CNC : **12**
- Sans code exploitable détecté : **7**

## 3. Premières constatations (provisoires)

### Fond du socle allégé

[Fait] Terrain2STL remplace le fond par un éventail unique quand la plaque est un rectangle plein, et son commentaire indique que le fond par cellule doublait le nombre de triangles (≈4,3 M contre ≈2,2 M). [Fait] `dem-to-stl` triangule son fond à la même densité que le dessus, donc environ la moitié de ses triangles est un plan plat. [Déduction] Pour un disque simple, le contour du fond est étoilé par rapport au centre : un éventail est possible, avec repli sur la grille complète si le masque a des trous ou plusieurs parties. [Hypothèse — à vérifier par benchmark] Les slicers (PrusaSlicer, Cura) acceptent sans problème de longs triangles très fins ; gain attendu : ~−45 % de triangles, de taille de fichier et de RAM, **sans effet sur la qualité de la surface**. Classe provisoire : B. Difficulté moyenne, risque moyen.

### Parois et indexation vectorisées

[Fait] `gpx3mf` (MIT) et Terrain2STL construisent parois et fond sans boucle Python (masques NumPy), et `gpx3mf` produit directement un maillage indexé. [Fait] `dem-to-stl` utilise une double boucle Python pour les parois et une fusion de sommets par tri pour 3MF/OBJ. [Hypothèse — à vérifier par benchmark] Gain de temps de quelques secondes et de mémoire sur les gros maillages ; aucun effet sur la géométrie. Classe provisoire : C. Difficulté faible à moyenne, risque faible (le résultat doit rester identique triangle pour triangle, ce que vos tests savent déjà vérifier).

### Plafond de 4 itérations du comblement de pincements

[Fait] `fix_diagonal_pinches` s'arrête à 4 itérations ; Terrain2STL en autorise 12. [Fait, mesuré] Sur des masques aléatoires (stress) de 300×300, 4 itérations laissent jusqu'à 8 pincements, et jusqu'à 7 itérations sont nécessaires pour converger ; sur des disques bruités réalistes (0,5 à 5 % de pixels retirés) passés par `clean_mask`, **0 pincement restant** dans mes 30 tirages. [Opinion] Faiblesse théorique, quasi sans effet en pratique d'après ce test synthétique (les masques réels de côtes ou de trous LiDAR ne sont pas testés) ; correction triviale (boucler jusqu'à stabilité). Classe provisoire : C.

### Pincement diagonal : problème confirmé par d'autres

[Fait] TouchTerrain (option `clean_diags`, désactivée par défaut), `gpx3mf` (`remove_pinches`), Terrain2STL (`close_diagonals`) et `my3dmaps` traitent le même défaut. [Déduction] Votre correction est dans la norme du domaine, et votre test de comblement (au lieu de retrait de matière) est un choix défendable.

### Diagonale de la surface

[Fait] TouchTerrain, `gpx3mf` et Terrain2STL découpent chaque cellule avec une diagonale fixe ; la logique « anti-triangles fins » de TouchTerrain ne concerne que les parois et est commentée. [Déduction] Parmi les dépôts lus, aucun n'obtient une meilleure surface par un choix de diagonale : `dem-to-stl` est à égalité. [Hypothèse — à vérifier par benchmark] Choisir la diagonale selon le relief local réduirait les artefacts en escalier sur les arêtes vives, mais à 0,2 mm de pas l'effet visible est probablement faible. **Correction (rapport 2) :** cette dernière phrase était fausse pour le lot 1. `my3dmaps` choisit bien la diagonale selon le relief ; le benchmark est détaillé dans le rapport 2 (fiche 2).

### Raffinement adaptatif de `jkittner/dem-to-stl`

[Fait] Triangulation de Delaunay avec insertion de points là où le dénivelé par triangle dépasse un seuil ; elle **ajoute** des points à un maillage de base au lieu de simplifier, n'a pas de borne d'erreur, échantillonne le DEM par interpolation, et dépend d'Earth Engine, FastAPI et rioxarray. [Opinion] Semble plus sophistiqué, mais n'apporte rien à un maillage déjà calé sur la résolution d'impression. Classe provisoire : E. Sans licence détectée.

### « Frontière exacte » de TopoTwin

[Fait] Delaunay d'un nuage de points puis conservation des triangles dont le **centroïde** est dans le polygone. [Déduction] Le contour suit des arêtes de triangles, il n'est pas exactement sur le polygone, sauf si les sommets du contour sont insérés dans le nuage (non vérifié). Le dépôt prévoit une boucle de réparation, ce qui suggère que ce procédé produit des défauts de manifold. Classe provisoire : E. Sans licence détectée.

### rayshader

[Fait] `save_3dprint` appelle `write_stl`, fonction propre à rayshader qui relit les sommets de la scène `rgl` ; la simplification par `max_error` est déléguée au paquet `terrainmeshr`. **(Corrigé au rapport 3 : la version initiale attribuait l'écriture STL à `rgl`.)** [Déduction] L'archive ne contient pas l'algorithme de maillage : rien à reprendre. Classe provisoire : E.

## 4. Lots proposés pour l'analyse approfondie

- **Lot 1** (17 dépôts) : maillage, murs, fermeture du solide, nodata. C'est là que se jouent la qualité du STL et la robustesse.
- **Lot 2** (12 dépôts) : gros DEM, parallélisme, autres langages (C#, Rust, R, MATLAB), CNC. Lecture ciblée sur les algorithmes, pas sur l'application.
- **Lot 3** (13 dépôts) : sources nationales et traitement de données réelles ; utile pour votre fonction `--source`.
- Les dépôts « hors cœur », « redondants » et « sans code » ne seront relus que pour y chercher une idée précise ; je les signale dans le rapport final avec la raison de l'écart.

## 5. Table des 93 dépôts (triage automatique)

Colonnes : LOC = lignes de code source ; Licence = détectée dans `LICENSE*` puis dans les métadonnées de paquet / README (**détection partielle, à confirmer avant tout emprunt**) ; Tests = nombre de fichiers de test détectés ; Statut = ce que j'ai déjà lu.

| # | Dépôt | Langages (LOC) | Licence détectée | Tests | Catégorie | Statut |
|---|---|---|---|---|---|---|
| 1 | github.com/USACE-WRISES/terrain-to-stl | Python 9214, TS 7552, HTML 12 | ABSENTE | 7 | Lot 1 — maillage / robustesse | — |
| 2 | github.com/ChHarding/TouchTerrain_for_CAGEO | Python 6477, JS 3953, Notebook 2787 | LICENSE: «This project falls under GPL (v3) license rules. For details | 4 | Lot 1 — maillage / robustesse | lu partiellement (grid_tesselate.py : classe quad) |
| 3 | github.com/ansonl/DEM2STL | Python 5580, Shell 1848, Go 26 | AGPL (LICENSE) | 2 | Lot 1 — maillage / robustesse | — |
| 4 | github.com/jawensi/DEMto3D-QGIS-Plugin | Python 4560, TS 646 | ABSENTE | 0 | Lot 1 — maillage / robustesse | — |
| 5 | github.com/pingud98/my3dmaps | Python 4227, JS 382, HTML 103 | MIT (LICENSE) | 7 | Lot 1 — maillage / robustesse | — |
| 6 | github.com/kelsorj/Lidr-STL | Python 4589 | ABSENTE | 2 | Lot 1 — maillage / robustesse | — |
| 7 | github.com/fabianbachl/geostl | Python 3191 | MIT (LICENSE) | 13 | Lot 1 — maillage / robustesse | — |
| 8 | github.com/fgebhart/mapa | Python 2931, Notebook 217 | MIT (LICENSE) | 12 | Lot 1 — maillage / robustesse | — |
| 9 | github.com/dsatory/terrain-to-stl | Python 2488 | MIT (LICENSE) | 10 | Lot 1 — maillage / robustesse | — |
| 10 | github.com/g4goza/gpx3mf | Python 2325 | MIT (LICENSE) | 21 | Lot 1 — maillage / robustesse | lu partiellement (mesh.py) |
| 11 | codeberg.org/stroebjo/dem-to-stl | Python 722, Shell 158, JS 134 | AGPL (LICENSE) | 0 | Lot 1 — maillage / robustesse | — |
| 12 | github.com/silvioheinze/dtm2stl | Notebook 876, Python 66 | LICENSE: «dtm2stl Copyright (c) 2025 dtm2stl contributors Licensed und | 0 | Lot 1 — maillage / robustesse | — |
| 13 | github.com/timsoleskos/topoPrinter | Python 621 | MIT (LICENSE) | 4 | Lot 1 — maillage / robustesse | — |
| 14 | github.com/anoved/phstl | Python 367 | MIT (LICENSE) | 0 | Lot 1 — maillage / robustesse | — |
| 15 | github.com/mchristoffersen/dem2stl | Python 317 | GPL-3 (LICENSE) | 0 | Lot 1 — maillage / robustesse | — |
| 16 | codeberg.org/tok/geotiff2stl | Python 312 | pyproject.toml: license = {text = | 0 | Lot 1 — maillage / robustesse | — |
| 17 | github.com/fergujm2/topo-print-maker | Python 254 | ABSENTE | 0 | Lot 1 — maillage / robustesse | — |
| 18 | github.com/Akilth/MapLab3D | Matlab 76834 | AGPL (LICENSE) | 0 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 19 | github.com/dem-net/DEM.Net | C# 34074 | MIT modifiée, restriction commerciale (corrigé au rapport 3 ; détection automatique erronée) | 13 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 20 | github.com/domlysz/BlenderGIS | Python 33828 | GPL-3 (LICENSE) | 0 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 21 | github.com/tylermorganwall/rayshader | R 19315, C++ 2191 | DESCRIPTION: License: GPL-3 | 6 | Lot 2 — gros DEM / autres écosystèmes / CNC | survolé (save_3dprint) |
| 22 | github.com/DavidIngraham/contour_workbench | TS 10730, Rust 4210, JS 125 | GPL-3 (LICENSE) | 40 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 23 | github.com/twigley/terrology | Python 8948, JS 1088, HTML 282 | Apache-2.0 (LICENSE) | 11 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 24 | github.com/suchanek/satprint | Python 8668, JS 488, HTML 168 | MIT (LICENSE) | 11 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 25 | github.com/shanto462/HeightMapTo3dTerrain | Rust 2924, Python 260 | MIT (LICENSE) | 1 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 26 | github.com/tfoss/DEM_STL_carver | Python 1078 | ABSENTE | 0 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 27 | github.com/kevinstadler/dem2stl | R 254, HTML 221, JS 8 | MIT (LICENSE) | 0 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 28 | github.com/DomPrott/Topographic-3D-Printing | R 383 | MIT (LICENSE) | 0 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 29 | github.com/kilimnik/heightmap_to_stl | Rust 216 | GPL-3 (LICENSE) | 0 | Lot 2 — gros DEM / autres écosystèmes / CNC | — |
| 30 | github.com/fkhr79/nrw-terrain-print | Python 8286 | LICENSE: «# PolyForm Noncommercial License 1.0.0 <https://polyformproj | 8 | Lot 3 — sources nationales / données | — |
| 31 | github.com/stuart-mckay-dev/portland-topo-map | Python 4363 | Unlicense/CC0 (LICENSE-docs.md) | 1 | Lot 3 — sources nationales / données | — |
| 32 | github.com/aareblau/Geo3Dprint | Python 3829 | ABSENTE | 5 | Lot 3 — sources nationales / données | — |
| 33 | github.com/stretchyboy/geo-printing | Python 3084 | ABSENTE | 2 | Lot 3 — sources nationales / données | — |
| 34 | github.com/jeremiahgivens/TerrainToSTL | C# 2053, HTML 102 | ABSENTE | 2 | Lot 3 — sources nationales / données | — |
| 35 | github.com/tsukada-eri/japan-3d-terrain | Python 891 | MIT (LICENSE) | 0 | Lot 3 — sources nationales / données | — |
| 36 | github.com/Fluorax/dem2stl-relief | Python 802 | ABSENTE | 0 | Lot 3 — sources nationales / données | — |
| 37 | github.com/katie-jones/elevation-to-stl | Python 422 | MIT (LICENSE) | 0 | Lot 3 — sources nationales / données | — |
| 38 | github.com/wolfgangB33r/PeakPrinter | Python 302 | GPL-3 (LICENSE) | 0 | Lot 3 — sources nationales / données | — |
| 39 | github.com/heptaliane/jp_geo_generator | Python 263 | MIT (LICENSE) | 0 | Lot 3 — sources nationales / données | — |
| 40 | github.com/modestustr/turkeytopgraphyDEM2STL | Python 249 | MIT (LICENSE) | 0 | Lot 3 — sources nationales / données | — |
| 41 | github.com/Spanholz/DEMtoPrint | Notebook 232 | MIT (LICENSE) | 0 | Lot 3 — sources nationales / données | — |
| 42 | github.com/OpenFlightHub/convert-dem-to-stl | Python 98, Shell 66 | MIT (LICENSE) | 0 | Lot 3 — sources nationales / données | — |
| 43 | github.com/nitinkachhadiya3/Terrain2STL | Python 14647, TS 3591, JS 871 | ABSENTE | 29 | Hors cœur | lu (meshing.py, surface/processing/quality : survol) |
| 44 | github.com/Wenshuishi0528/topotile-studio | Python 12012, HTML 2927, Shell 9 | CC-BY (LICENSE) | 14 | Hors cœur | — |
| 45 | github.com/pavelzaitsau/gora-print-router | Python 13739 | AGPL (LICENSE) | 23 | Hors cœur | — |
| 46 | github.com/agrignon/MapMaker | TS 12412, HTML 13 | ABSENTE | 19 | Hors cœur | — |
| 47 | github.com/ginoclement/terrain | JS 5526, HTML 376 | ABSENTE | 6 | Hors cœur | — |
| 48 | github.com/irjudson/terraprint | Python 3882, HTML 1082, Shell 406 | ABSENTE | 3 | Hors cœur | — |
| 49 | github.com/ModelRift/terrain-to-3d | TS 4413, JS 23, HTML 15 | ABSENTE | 0 | Hors cœur | — |
| 50 | github.com/MattMullin/openmoon | JS 2247, Python 1786, HTML 214 | ABSENTE | 0 | Hors cœur | — |
| 51 | github.com/formscope/gps-terrain-stl | Python 3104, HTML 891 | ABSENTE | 0 | Hors cœur | — |
| 52 | github.com/andrewsheerin/topographic-maps | Python 2386, JS 1282, HTML 12 | ABSENTE | 7 | Hors cœur | — |
| 53 | github.com/cjorgensenmd/dem-to-3d-print-blender | Python 3584 | MIT (LICENSE) | 0 | Hors cœur | — |
| 54 | github.com/drewjobson/TopoTwin | Python 3569 | ABSENTE | 0 | Hors cœur | lu partiellement (mesh_builder.py) |
| 55 | github.com/jkittner/dem-to-stl | Python 2283, JS 930, HTML 330 | ABSENTE | 0 | Hors cœur | lu partiellement (converter.py) |
| 56 | github.com/hicatlab/DEM-Cropper | JS 2093, Python 826, HTML 258 | MIT (LICENSE.md) | 0 | Hors cœur | — |
| 57 | github.com/chris17453/terraprint3d | Python 2759 | README.md: MIT | 2 | Hors cœur | — |
| 58 | github.com/Lord0fBytes/gpx-terrain-studio | TS 2332, HTML 13 | ABSENTE | 11 | Hors cœur | — |
| 59 | github.com/avchaykin/3d-relief-map | JS 2042, HTML 153 | ABSENTE | 2 | Hors cœur | — |
| 60 | github.com/darrenoakey/map-generator | Python 2046 | README.md: licensed under [CC BY-NC 4.0](https://darren-static.waft.de | 8 | Hors cœur | — |
| 61 | github.com/niyoseris/3droute | Rust 834, JS 690, HTML 268 | MIT (LICENSE) | 0 | Hors cœur | — |
| 62 | github.com/jkoether/Terrain_Trails | Python 1668 | GPL-3 (LICENSE) | 0 | Hors cœur | — |
| 63 | github.com/ThatcherC/Terrain2STL | JS 536, C 407, HTML 274 | package.json: "license": "ISC | 1 | Hors cœur | — |
| 64 | github.com/JGGroot/ODTT | JS 1042, HTML 38 | ABSENTE | 1 | Hors cœur | — |
| 65 | github.com/alb0/3d-greenland | Notebook 440, Matlab 149 | Apache-2.0 (LICENSE.txt) | 0 | Hors cœur | — |
| 66 | github.com/stoccolma/terrainprint | HTML 578 | ABSENTE | 0 | Hors cœur | — |
| 67 | github.com/Hands3DLab/Generate-Terrain-with-GPX | Python 469 | README.md: MIT | 0 | Hors cœur | — |
| 68 | github.com/aparshin/srtm-to-stl | JS 202 | package.json: "license": "ISC | 0 | Hors cœur | — |
| 69 | github.com/aligundogdu/TerraStl | TS 145, JS 12 | README.md: licensed under the MIT License - see the https://www.google | 0 | Hors cœur | — |
| 70 | github.com/riggsd/ArcGIS-Pro-DEM-3D-Print-Toolbox | — | GPL-3 (LICENSE.txt) | 0 | Hors cœur | — |
| 71 | github.com/PaulZC/SRTM_to_STL | Python 1539 | MIT (LICENSE) | 0 | Redondant probable | — |
| 72 | github.com/spelfrog/3d-aachen | Python 1004 | MIT (LICENSE) | 0 | Redondant probable | — |
| 73 | github.com/coelacant1/ElevationDataTo3DModel | Python 580 | GPL-3 (LICENSE) | 0 | Redondant probable | — |
| 74 | github.com/r-barnes/DEMto3D | Python 547 | MIT (LICENSE) | 0 | Redondant probable | — |
| 75 | github.com/IdeaShakerLab/dem-to-stl | Python 542 | MIT (LICENSE) | 0 | Redondant probable | — |
| 76 | github.com/billyauhk/dem2stl | C 364, C++ 151 | ABSENTE | 0 | Redondant probable | — |
| 77 | github.com/lopoc/dem2stl | Python 504 | ABSENTE | 2 | Redondant probable | — |
| 78 | github.com/GeoffWilliams/dem2stl | Python 469 | LGPL (LICENSE) | 0 | Redondant probable | — |
| 79 | github.com/ScreamingHawk/geotiff-to-mesh | Python 454 | MIT (LICENSE) | 0 | Redondant probable | — |
| 80 | github.com/41pha1/Heightmap2STL | Python 440 | ABSENTE | 0 | Redondant probable | — |
| 81 | github.com/cvr/dem2stl | Python 364 | GPL-3 (LICENSE) | 0 | Redondant probable | — |
| 82 | github.com/oliver-rew/topo_tool | Python 287, Shell 6 | ABSENTE | 0 | Redondant probable | — |
| 83 | github.com/BrendanLuke15/DTM-2-STL | Python 203, Matlab 59 | MIT (LICENSE) | 0 | Redondant probable | — |
| 84 | github.com/Usuy-Leon/Topo3D-world | Python 242 | MIT (LICENSE) | 0 | Redondant probable | — |
| 85 | github.com/ewandennis/geotiff2stl | Python 162 | ABSENTE | 0 | Redondant probable | — |
| 86 | github.com/stefangluszek/GeoTiff2STL | Python 92 | ABSENTE | 0 | Redondant probable | — |
| 87 | github.com/F-Reid-McLain/macon-3d-print | — | ABSENTE | 0 | Sans code exploitable détecté | — |
| 88 | github.com/eatdeath86/terrain-relief-konfigurator | — | ABSENTE | 0 | Sans code exploitable détecté | — |
| 89 | github.com/leelakrishnane/terrain_map_3dprinting | — | ABSENTE | 0 | Sans code exploitable détecté | — |
| 90 | github.com/UnintelligibleMaker/PyTerrainModeler | — | ABSENTE | 0 | Sans code exploitable détecté | — |
| 91 | github.com/danny-kumpf/gpx-tracks-3d | — | ABSENTE | 0 | Sans code exploitable détecté | — |
| 92 | github.com/furrysalamander/mini-map-maker | — | ABSENTE | 0 | Sans code exploitable détecté | — |
| 93 | github.com/marcinbator/asc-to-3d-model-blender-script | — | ABSENTE | 0 | Sans code exploitable détecté | — |

## 6. Réserves sur ce triage

- Les catégories viennent des descriptions README et de signaux de code ; elles sont des **déductions** tant que le dépôt n'est pas lu.
- La détection de licence est partielle : un dépôt « ABSENTE » peut déclarer sa licence ailleurs ; sans licence explicite, le code est protégé par défaut et ne doit pas être réutilisé. `dem-to-stl` étant en GPL-3.0, les emprunts de code (hors scope pour l'instant) ne seraient envisageables qu'avec des licences compatibles (MIT, BSD, Apache-2.0, GPL-3) ; AGPL, PolyForm Noncommercial et CC-BY-NC sont problématiques. Ce n'est pas un avis juridique.
- L'extension `.ts` est comptée comme « TS » : elle peut désigner du TypeScript ou des fichiers de traduction Qt (cas probable du plugin QGIS).
- Plusieurs dépôts n'ont aucun code détecté (README seul, fichier non analysable, ou boîte à outils ArcGIS en `.pyt`, extension que mon scan ne lit pas).
- Les téléchargements correspondent à la branche par défaut au 9 octobre 2026, pas forcément à la version que vous avez examinée.