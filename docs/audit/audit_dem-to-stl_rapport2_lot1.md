# Audit comparatif dem-to-stl — Rapport 2 : lot 1 (maillage et robustesse)

Date : 9 octobre 2026. **Aucun fichier de votre projet n'a été modifié.** Les mesures ci-dessous ont été faites sur des copies de travail, avec des scripts écrits pour l'occasion (je peux vous les fournir).

Convention : **[Fait]** constaté dans le code ou mesuré ; **[Déduction]** ; **[Hypothèse — à vérifier par benchmark]** ; **[Opinion]**.

---

## 0. Résumé

1. **La meilleure idée du lot 1 est le fond plat fusionné de `my3dmaps`** : couvrir le fond par des rectangles, triangulés en éventail sans jonction en T. Pour un disque de 250 mm à 0,2 mm, elle ramène le fond de 2 454 368 à ≈ 6 500 triangles, soit **−50 % de triangles au total**, sans changer la géométrie. Classe B. Reste à valider sur de vrais slicers (triangles très allongés).
2. **Le choix de diagonale selon le relief (`my3dmaps`) apporte peu** : −3 % d'erreur de surface en moyenne sur un vrai DEM, rien sur des terrains synthétiques. Classe C. **Cela corrige le rapport 1.**
3. **Sur les données imparfaites, aucun dépôt du lot 1 ne fait mieux que `dem-to-stl`.** Sur la validation et les tests, vous êtes au-dessus de tout ce que j'ai lu dans ce lot.
4. **Trois idées séduisantes sont à écarter** : la scission de sommets aux coins en damier (non valable en STL), le maillage « 4 triangles par pixel », les triangulations adaptatives.
5. **Deux corrections triviales** : boucler les comblements de pincements jusqu'à stabilité (plafond de 4 itérations), et ajouter le volume signé dans `check_stl.py`.

## 1. Correction du rapport 1

Le rapport 1 écrivait qu'aucun dépôt lu ne choisissait la diagonale selon le relief. C'était vrai des cinq dépôts lus à ce moment-là, et **faux pour le lot 1** : `my3dmaps` (`mesh.py`) choisit, pour chaque cellule, la diagonale qui relie les deux coins de plus faible écart d'altitude. Le rapport 1 a été annoté en conséquence.

## 2. Dépôts du lot 1 : niveau de lecture et verdict

« Lu » = code de maillage lu en détail ; « partiel » = fonctions clés lues, le reste non ; « survol » = README, signatures et recherches ciblées. **Ce qui n'est pas lu n'est pas couvert par mes conclusions.**

| Dépôt | Licence (détectée) | Lecture | Apport pour dem-to-stl | Verdict |
|---|---|---|---|---|
| **my3dmaps** (pingud98) | MIT | `mesh.py` lu en détail | Fond/dessus fusionné par rectangles + éventails ; diagonale adaptative ; rapport de manifold | **Fiche 1, 2** |
| **geostl** (fabianbachl) | MIT | `mesh.py` lu ; `tiling`, `rectify` non lus | Architecture propre ; maillage indexé vectorisé ; tuiles raccordées | Fiche 3, 10 |
| **gpx3mf** (g4goza) | MIT | `mesh.py`, `raster.py` (signatures) | Maillage indexé vectorisé ; `remove_pinches` par retrait | Fiche 3, 5 |
| **mapa** (fgebhart) | MIT | `algorithm.py` lu | Fond en éventail (plaque rectangulaire) ; tests à STL de référence | Fiche 1, 6, 9 |
| **tok/geotiff2stl** | AGPL-3.0-or-later (en-tête SPDX) | `main.py` lu en entier | Fond en éventail, parois vectorisées | Fiche 1, 3 (idée seule) |
| **TouchTerrain** (ChHarding) | GPL-3 | `grid_tesselate.py` partiel (~6 500 lignes Python non lues) | NaN → cellules masquées ; diagonale fixe | Fiche 2, 5 |
| **Terrain2STL** (nitinkachhadiya3) | aucune | `meshing.py` lu (rapport 1) | Fond en éventail ; passes de pincements = 12 | Fiche 1, 3, 5 |
| **phstl** (anoved) | MIT | partiel | Écriture facette par facette | Fiche 8 |
| **topoPrinter** (timsoleskos) | MIT | `mesh_builder.py` lu | Rien (boucles Python, diagonale fixe) ; aplatissement des lacs | Fiche 11 |
| **jawensi/DEMto3D** | en-têtes GPL-2.0-or-later, pas de fichier LICENSE | partiel | Rien d'algorithmique ; modèle inversé (négatif) | Fiche 11 |
| **dsatory/terrain-to-stl** | MIT | README + docstrings | Hors cœur (bâtiments NYC, tuiles à pions) | Écarter |
| **USACE/terrain-to-stl** | **aucune** | liste de fonctions + README | Bandes de transition multirésolution (HEC-RAS) | Fiche 7 |
| **ansonl/DEM2STL** | AGPL | README | Workflow autour de TouchTerrain, rien de nouveau | Écarter |
| **stroebjo/dem-to-stl** | AGPL | fonction masque lue | Non étanche, zéro = absence de donnée | Écarter |
| **mchristoffersen/dem2stl** | GPL-3 | signatures + début | Rien | Écarter |
| **topo-print-maker** (fergujm2) | **aucune** | signaux de code | Lissage gaussien σ=1 avant maillage | Écarter |
| **dtm2stl** (silvioheinze) | EUPL | libs du notebook | Rien (numpy-stl, notebooks) | Écarter |
| **Lidr-STL** (kelsorj) | **aucune** | README | Hors cœur (ville LiDAR) | Écarter |

Une **réimplémentation indépendante** d'une idée reste possible quelle que soit la licence ; pour un emprunt de code, seuls MIT (my3dmaps, geostl, gpx3mf, mapa, phstl, topoPrinter, dsatory) et GPL-2+/GPL-3 sont compatibles avec votre GPL-3 ; AGPL, EUPL et « aucune licence » sont à traiter comme inspiration conceptuelle uniquement. Ce n'est pas un avis juridique.

## 3. Fiches comparatives

### Fiche 1 — Fond plat fusionné (rectangles + éventails) · Classe B

**Projet concurrent.** *my3dmaps* ([Fait] `sheet_mesh`, `_flat_rects`, `_fan_rects`) : les cellules dont les quatre coins ont exactement la même hauteur sont couvertes par des rectangles gloutons (ligne par ligne, extension vers le bas tant que la ligne suivante a la même étendue). Chaque rectangle est triangulé en éventail depuis son centre, en gardant tous les nœuds « requis » de son contour (nœuds du bord de l'ensemble fusionné et coins des rectangles voisins). Il n'y a donc pas de jonction en T. Cela fonctionne pour un masque quelconque. *Terrain2STL*, *mapa* et *tok* font la version simple (éventail unique depuis le centre) pour une plaque rectangulaire.

**Mon projet.** [Fait] Le fond est triangulé à la même densité que le dessus : 2 triangles par cellule.

**Différence.** [Fait] Même plan, même contour, même volume ; seule la triangulation change. [Mesure, estimation par simulation de l'algorithme, pas par une implémentation intégrée] :

| Cas | Fond actuel | Fond fusionné | Total hors parois |
|---|---|---|---|
| Ø250 mm, 0,2 mm (1 250², 1 227 184 cellules) | 2 454 368 | ≈ 6 464 (733 rectangles) | 4 908 736 → 2 460 832 (**−49,9 %**) |
| Ø120 mm, 0,2 mm | 565 584 | ≈ 3 100 | −49,7 % |
| Ø250 mm, 0,1 mm | 9 817 472 | ≈ 12 928 | 19,6 M → 9,8 M (−49,9 %) |

Le premier cas retrouve les 4,9 M de triangles de votre README. La couverture gloutonne en Python pur prend 1,8 s (0,2 mm) et 6,4 s (0,1 mm) ; non optimisé.

**Gain.** Important sur la taille du fichier, la RAM et le temps ; la garde RAM refuserait moins de cas. **Aucun effet sur la qualité de la surface.**

**Difficulté.** Moyenne : ≈ 60 lignes d'algorithme, plus l'intégration dans `build_mesh_from_arrays`, `index_mesh` (sommets centraux nouveaux) et les 3 formats d'export. [Déduction] `check_stl.py` reste valable : les sommets de contour du fond restent ceux des parois.

**Risque.** Moyen, pour trois hypothèses :
- [Hypothèse — à vérifier par benchmark] PrusaSlicer, Cura et OrcaSlicer acceptent sans réparation des triangles très allongés. Dans le cas Ø250 mm, un rectangle d'une cellule de haut sur ~1 000 cellules de large donne des triangles d'environ 100 mm × 0,1 mm (rapport ~1 000:1). `mapa` et `tok` produisent déjà de tels triangles sur leur fond, mais je n'ai trouvé dans leurs archives aucun test avec un slicer.
- [Hypothèse — à vérifier] Votre CAM CNC et vos visionneuses importent ce fond sans erreur.
- [Fait] Vos tests de propriétés (étanchéité, orientation, volume exact calculé à part) constituent le filet de sécurité naturel ; ils devraient passer sans modification.

**Verdict.** **Intégrer après prototype et test sur slicers.** Priorité n° 1 de ce lot.

### Fiche 2 — Diagonale adaptative · Classe C

**Projet concurrent.** *my3dmaps* choisit la diagonale a–c si |za−zc| ≤ |zb−zd|, sinon b–d. [Fait] *TouchTerrain*, *gpx3mf*, *Terrain2STL*, *geostl*, *tok* et *topoPrinter* utilisent une diagonale fixe ; la logique « anti-triangles fins » de TouchTerrain ne concerne que les parois et est commentée.

**Mon projet.** [Fait] Diagonale fixe.

**Mesure** (erreur quadratique moyenne entre la surface triangulée et le terrain échantillonné plus finement, toutes positions de la cellule ; ratios par rapport à une diagonale fixe de référence, 1,000) :
- **Terrains synthétiques isotropes** (3 types, pas 4 et 16) : **écart < 0,3 %** (ratios 0,997 à 1,002), la règle adaptative incluse.
- **Vrai DEM** (SRTM3 N44W070, Maine, 4 fenêtres de 577², pas 4 et 16, 8 cas) : adaptative −7,2 % à +1,0 % d'erreur, **moyenne −3,4 %** ; dans les 8 cas elle reste à ≤ 3,6 % de la meilleure diagonale fixe. Une diagonale fixe mal orientée par rapport aux crêtes dominantes peut être jusqu'à ~12 % moins bonne que l'autre.

**Différence.** [Déduction] Le gain existe quand le relief a une orientation dominante (crêtes parallèles), ce que les terrains synthétiques isotropes n'ont pas. La règle adaptative évite de dépendre de cette orientation sans la connaître.

**Gain.** Faible à moyen en erreur de surface (~3 %) ; [Hypothèse — à vérifier] invisible à 0,2 mm de pas, plus pertinent pour un pas d'impression grossier ou un DEM sous-échantillonné.

**Difficulté.** Faible (un `np.where` vectorisé). **Risque.** Faible à moyen : tous les triangles de surface changent, donc les références de test à régénérer ; les arêtes de bord et l'étanchéité sont inchangées.

**Limites de ma mesure.** Un seul DEM réel de 90 m (SRTM3), pas de LiDAR ; l'effet réel sur des LiDAR 1 m (IGN, GSI) est inconnu. **Verdict : étudier** (benchmark sur un de vos DEM réels avant toute décision).

### Fiche 3 — Parois et indexation vectorisées · Classe C

**Projet concurrent.** [Fait] *geostl*, *gpx3mf*, *tok* et *Terrain2STL* construisent les parois et le fond sans boucle Python ; *geostl* et *gpx3mf* produisent directement un maillage indexé (sommets partagés).

**Mon projet.** [Fait] Double boucle Python sur les arêtes de bord ; triangles en soupe float64 ; soudure des sommets par tri (`index_mesh`) pour 3MF/OBJ.

**Différence.** [Hypothèse — à vérifier par benchmark] gain de quelques secondes sur les gros maillages et, avec indexation native, une baisse notable de la RAM de la phase de construction ; aucun effet sur la géométrie.

**Gain** faible à moyen. **Difficulté** moyenne. **Risque** faible (résultat identique triangle pour triangle, vérifiable par vos tests). **Verdict : étudier, à faire en même temps que la fiche 1** (le fond est de toute façon réécrit).

### Fiche 4 — Fusion des zones plates du *dessus* · Classe D

**Projet concurrent.** *my3dmaps* (`merge_flat`) fusionne aussi les cellules du dessus dont les 4 coins sont **exactement** égaux.

**Différence.** [Déduction] Sur un DEM rééchantillonné par moyenne, des plateaux exactement plats sont rares, sauf mer écrêtée ou plans d'eau aplatis. [Hypothèse] gain nul sur vos reliefs terrestres standards. **Verdict : ne pas intégrer** pour l'instant (réexaminer si vous ajoutez un aplatissement des plans d'eau).

### Fiche 5 — Pincements en diagonale : trois stratégies · Classe C (boucle) / E (scission)

| Stratégie | Dépôt | Effet |
|---|---|---|
| **Combler** la cellule | **dem-to-stl** | Ajoute de la matière (≤ 1 cellule), continuité du relief |
| **Retirer** une cellule | gpx3mf, Terrain2STL (passes ≤ 12) | Enlève de la matière |
| **Scinder** les sommets | my3dmaps | Garde le contour exact, deux coins coïncidents |

- [Fait] *my3dmaps* duplique les sommets (blocs « alt ») aux coins en damier. Son contrôle de manifold compte les arêtes **par indice** de sommet.
- [Déduction] Un STL n'a pas d'indices : les deux sommets coïncidents sont ressoudés au chargement par les slicers, et l'arête verticale se retrouve partagée par 4 faces (non-manifold). La garantie ne tient que pour des formats indexés comme le 3MF. [Hypothèse — comportement des slicers à vérifier.] Pour du STL, **votre comblement est le bon choix**. Verdict scission : **ne pas intégrer** (classe E).
- [Fait, mesuré] Plafond de 4 itérations : sur des masques aléatoires de stress (300×300), il laisse jusqu'à 8 pincements, et jusqu'à 7 itérations sont nécessaires ; sur des disques bruités réalistes passés par `clean_mask`, 0 pincement restant (30 tirages). **Corriger en bouclant jusqu'à stabilité (avec un plafond de sécurité large)** : difficulté triviale, risque nul, classe C.

### Fiche 6 — Maillage « 4 triangles par pixel » (*mapa*) · Classe E

[Fait] Chaque pixel devient 4 triangles autour d'un centre ; les coins sont la moyenne de 4 pixels voisins (lissage implicite 2×2) ; le tout est écrit avec `numba`. [Déduction] Deux fois plus de triangles pour le même échantillonnage, aucune borne d'erreur, et une dépendance lourde. Votre grille est déjà calée sur le pas d'impression. **Ne pas intégrer.**

### Fiche 7 — Triangulations adaptatives et multirésolution · Classe E

*jkittner* (raffinement de Delaunay qui **ajoute** des triangles), *TopoTwin* (Delaunay avec masque par centroïde, rapport 1), *USACE/terrain-to-stl* (zones raffinées raccordées par des bandes de transition en éventail, [Fait] pensées pour des régions imposées par des données HEC-RAS), *rayshader* (simplification déléguée au paquet `terrainmeshr`). [Opinion] Aucune n'offre de borne d'erreur contrôlée ni de gain démontré sur un maillage déjà à la résolution d'impression. Les bandes de transition d'USACE sont l'idée la plus propre, mais répondent à un besoin que `dem-to-stl` n'a pas. **Ne pas intégrer.**

### Fiche 8 — Écriture en flux et RAM · Classe B conditionnelle

[Fait] *phstl* écrit facette par facette avec une mémoire constante, mais ne produit qu'une **surface** (pas de solide fermé), et saute le nodata facette par facette (surface trouée). [Fait] `dem-to-stl` construit tout le tableau de triangles (≈ 220 octets par triangle d'après votre propre mesure).

[Hypothèse — à vérifier par benchmark] Générer et écrire par bandes de lignes ramènerait la RAM à l'ordre des grilles (quelques octets par cellule). Le nombre de triangles de l'en-tête STL peut être calculé à l'avance ou réécrit à la fin (vous écrivez déjà un `.part`). Le 3MF et l'OBJ exigent un index global, donc ce serait un chemin à part.

**Difficulté** moyenne à importante. **Risque** moyen. **Verdict : étudier, mais seulement si la fiche 1 ne suffit pas pour vos tailles cibles** (elle divise déjà par deux la RAM de maillage).

### Fiche 9 — Validation et tests · Classe C / D

- [Fait] **Vos tests sont supérieurs à ceux du lot 1** : tests de propriétés sur blobs aléatoires, volume exact calculé indépendamment du code de maillage, relecture du STL par `check_stl.py`. Dans le lot 1 : *mapa* conserve des STL de référence (dossier `tests/stl` et script de régénération ; je n'ai pas lu les comparaisons elles-mêmes) ; *geostl* propose `is_watertight` de trimesh en option ; *my3dmaps* a `edge_manifold_report` et `signed_volume` ; *gpx3mf* a 21 fichiers de test détectés.
- **Idée à reprendre (D)** : le volume signé est déjà testé sous pytest, mais **n'est pas contrôlé par `check_stl.py` à l'exécution** (il vérifie les arêtes ouvertes, non-manifold et mal orientées). Un solide entièrement « retourné » mais cohérent passerait ce contrôle. Ajouter volume signé > 0 et comptage des triangles dégénérés : trivial.
- **Idée à reprendre temporairement (C)** : un STL de référence minimal, sur un petit DEM synthétique, protège contre un changement involontaire de triangulation pendant les fiches 1 à 3. C'est le cas où *mapa* a raison. À retirer ou à régénérer volontairement ensuite.
- **À ne pas faire (E)** : ajouter `trimesh` comme second avis. Votre contrôle indépendant fait déjà le travail, et la dépendance est lourde.

### Fiche 10 — Découpe en tuiles pour le plateau d'impression · Classe D

*geostl* (tuiles raccordées), *TouchTerrain*, *dsatory* (pions d'alignement), *my3dmaps*. [Fait] Ce tuilage est une **découpe pour l'imprimante**, pas une stratégie de mémoire. Hors du cœur de `dem-to-stl` (relief circulaire). [Opinion] À reconsidérer seulement si vous voulez imprimer des reliefs plus grands que le plateau.

### Fiche 11 — Fonctions de style · Classe D

Aplatissement des lacs (*topoPrinter*, *dsatory*), lissage gaussien σ=1 avant maillage (*topo-print-maker*), **modèle inversé / négatif** (*DEMto3D*, miroir et inversion de l'altitude). [Opinion] Aucune n'améliore la fidélité. Le **négatif** pourrait vous intéresser pour des moules ou des gabarits CNC : à vous de me dire si c'est un besoin réel.

## 4. Constats transversaux

**Qualité du maillage.** [Fait] Dans les dépôts dont j'ai lu le code de maillage, la surface repose sur une grille régulière à 2 triangles par cellule, sauf *mapa* (4 triangles par pixel) et *USACE* (zones raffinées). [Déduction] À même DEM et même pas, **aucun ne produit une surface supérieure à la vôtre**. Les différences portent sur le nombre de triangles (fiche 1) et, à la marge, sur l'orientation de la diagonale (fiche 2).

**Gros DEM.** [Fait] Dans le lot 1, le « tuilage » est une découpe d'impression (fiche 10) ; les mosaïques passent par `rasterio.merge` en mémoire (*tok*, *mapa*). [Déduction] Aucune stratégie de ce lot n'est supérieure à votre `gdal.Warp` vers une grille bornée par le pas d'impression. **Le vrai goulot est la mémoire du maillage de sortie** (fiches 1 et 8).

**Données imparfaites.** [Fait] *TouchTerrain* passe les NaN en cellules masquées et interpole les bords ; *geostl* refuse les valeurs non finies ; *tok* remplit avec l'altitude minimale ; *stroebjo* traite zéro comme absence de donnée ; *phstl* troue la surface sans fermer le solide ; *mapa* n'a aucune gestion. [Déduction] Votre combinaison (`FillNodata` GDAL sur 50 px, masque, ouverture morphologique, comblement des pincements, composante principale, bord lissé) est la plus complète du lot 1. Les lots 2 et 3 peuvent contredire cette conclusion.

**Architecture.** [Fait] *geostl* a la séparation la plus nette (`Region → Section → Mesher`, adaptateurs `ElevationSource.fetch(bbox)`, `TriangleMesh` indépendant de la bibliothèque d'export). [Opinion] Cela ressemble à votre `sources.py`. **Rien ne justifie une refonte.**

**Dépendances.** [Fait] `numpy-stl` (6 dépôts), `trimesh` (geostl en option, Terrain2STL), `numba` (mapa), `rasterio` (geostl, mapa, tok, dtm2stl). Vous : GDAL + NumPy + SciPy. [Déduction] Les fiches 1 à 3 se font en NumPy pur, **sans nouvelle dépendance** ; l'algorithme de *my3dmaps* est en Python/NumPy.

## 5. Mesures effectuées

1. **Pincements** : `fix_diagonal_pinches` de votre `mesh.py` appelé directement ; 60 masques aléatoires 300×300 (3 densités × 20) et 30 disques bruités (rayon 600, bruit 0,5 / 2 / 5 %) passés par `clean_mask`.
2. **Diagonale** : surface linéaire par morceaux comparée au terrain fin sur toute la cellule, quatre règles (fixe a–c, fixe b–d, adaptative, adaptative inverse). Terrains synthétiques par synthèse spectrale (2 049², trois types) et SRTM3 réel (577², 4 fenêtres).
3. **Fond fusionné** : couverture gloutonne par rectangles + comptage des nœuds de contour requis, sur un masque en disque ; **simulation de l'algorithme, pas une implémentation intégrée à votre code**.

**Non mesuré** : suite de tests de `dem-to-stl` (GDAL absent chez moi), temps réels de bout en bout, compatibilité slicers/CAM.

## 6. Réponses provisoires aux dix questions (lot 1 seulement)

1. **Meilleures idées algorithmiques** : *my3dmaps* (fond fusionné, diagonale adaptative).
2. **Meilleurs STL** : même géométrie que `dem-to-stl` ; la différence est la triangulation (−50 %), pas la fidélité.
3. **Gros DEM** : aucun dans ce lot.
4. **Données imparfaites** : `dem-to-stl`.
5. **Architecture** : *geostl* ≈ `dem-to-stl`.
6. **Meilleure idée découverte** : le fond fusionné.
7. **Meilleur rapport gain/coût** : boucle de pincements jusqu'à stabilité et volume signé dans `check_stl` (coût nul) ; ensuite le fond fusionné.
8. **Plus difficile mais intéressante** : écriture par bandes (fiche 8), seulement si nécessaire.
9. **Ne valent pas la complexité** : maillage 4 triangles, triangulations adaptatives, scission de sommets, `trimesh`, `numba`.
10. **Remise en cause de l'architecture** : **non**, rien dans ce lot.

## 7. Feuille de route provisoire (lot 1)

**Phase 1 — gains faciles** : boucle de pincements jusqu'à stabilité ; volume signé et triangles dégénérés dans `check_stl.py`.

**Phase 2 — améliorations importantes** : fond fusionné par rectangles (après prototype et test PrusaSlicer/Cura/Orca/CAM), avec vectorisation des parois dans le même chantier ; STL de référence minimal pendant la transition.

**Phase 3 — recherche** : diagonale adaptative (benchmark sur un de vos DEM LiDAR réels) ; écriture par bandes.

**À ne pas faire** : maillage 4 triangles par pixel, Delaunay adaptatif, scission de sommets, `trimesh`/`numba`, tuilage, fonctions de style.

## 8. Suite

Lot 2 (gros DEM, autres écosystèmes, CNC) : `DEM.Net`, `HeightMapTo3dTerrain`, `heightmap_to_stl`, `contour_workbench`, `BlenderGIS`, `satprint`, `terrology`, `MapLab3D`, `dem2stl` (R), `Topographic-3D-Printing` (R), `DEM_STL_carver`, et la lecture complète de `rayshader`.

Parties du lot 1 **non lues** et à lire si vous le souhaitez : `TouchTerrain` (hors `grid_tesselate.py`), `USACE/terrain-to-stl` (1 551 lignes), `geostl` (`tiling.py`, `rectify.py`).
