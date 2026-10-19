# Audit comparatif dem-to-stl — Rapport 5 : synthèse globale et feuille de route

Date : 10 octobre 2026. **Aucun fichier de votre projet n'a été modifié.** Ce document synthétise les rapports 1 à 4 ; les détails, mesures et réserves sont dans ces rapports (fiches 1 à 24).

Convention : **[Fait]** constaté dans le code ou mesuré ; **[Déduction]** ; **[Hypothèse — à vérifier par benchmark]** ; **[Opinion]**.

---

## 1. Couverture réelle de l'audit

| | Dépôts | Niveau |
|---|---|---|
| Lot 1 (maillage, robustesse) | 17 | code de maillage lu (en détail pour 5, partiel ou ciblé pour les autres) |
| Lot 2 (gros DEM, autres langages, CNC) | 12 | lecture ciblée (`satprint` : construction du maillage ; `hmterrain` : en-tête de l'algorithme et README ; `terrology` : fonction de découpe) |
| Lot 3 (sources nationales, données réelles) | 13 | ciblé |
| Hors lots, lus partiellement | 3 (Terrain2STL, jkittner, TopoTwin) | fonctions clés |
| **Triage seulement** (README, signaux de code, licence) | **48** | non lus : classés hors cœur (25), redondants probables (16), sans code exploitable (7) |

**Ce que cela veut dire.** Mes conclusions valent pour les 45 dépôts lus au moins en partie. Pour les 48 autres, « redondant » ou « hors cœur » est une **déduction** tirée de leurs descriptions : un de ces dépôts pourrait contenir une bonne idée. Les parties non lues de chaque lot sont listées en fin de rapports 2, 3 et 4.

**Ce que je n'ai pas fait.** Exécuter votre suite de tests (GDAL absent chez moi, sauf le noyau de rééchantillonnage via `rasterio`) ; tester les services réseau ; mesurer sur **vos** DEM. Toutes les mesures de maillage portent sur un DEM réel de 90 m (SRTM3, Maine) et sur des terrains synthétiques.

## 2. Classement des idées retenues

Gain / complexité / risque jugés **pour `dem-to-stl`**. « Source » indique d'où vient l'idée ; ce n'est pas toujours là que se trouve le code utilisable (voir section 5).

| Priorité | Idée | Projet source | Gain potentiel | Complexité | Risque | Verdict |
|---|---|---|---|---|---|---|
| 1 | **Rééchantillonner selon le sens** : `bilinear`/`cubic` quand le pixel de travail est plus fin que celui de la source (aujourd'hui `average` = `nearest` en agrandissement) | dem2stl-relief (idée) + mes mesures | **important** pour les sources grossières (pente ÷ ~2 sur SRTM ; 25–77 % de cellules plates évitées), nul sinon | faible | faible | **intégrer** (fiche 18) |
| 2 | **Fond plat fusionné** par rectangles + éventails sans jonction en T | my3dmaps | **important** (−50 % de triangles, RAM, taille, temps) ; aucun effet sur la surface | moyenne | moyen | **intégrer après prototype** testé sur slicers et CAM (rapport 2, fiche 1) |
| 3 | **Valider les téléchargements** (rattraper `IncompleteRead`, contrôler le contenu avant mise en cache) | nrw-terrain-print, Geo3Dprint (idées) | moyen (évite un cache empoisonné) | faible | faible | **intégrer** (fiche 19) |
| 4 | **Repli de résolution GSI** par couche (zoom ≥ 16) | japan-3d-terrain (idée) | moyen (Japon, pixel fin) | moyenne | faible | **intégrer** après le n° 1 (fiche 20) |
| 5 | **Boucler les pincements jusqu'à stabilité** | Terrain2STL (12 passes) | faible | triviale | nul | **intégrer** (rapport 2, fiche 5) |
| 6 | **Volume signé et triangles dégénérés dans `check_stl.py`** | my3dmaps | faible | triviale | nul | **intégrer** (rapport 2, fiche 9) |
| 7 | **Parois vectorisées, maillage indexé natif** | gpx3mf, geostl, satprint | faible à moyen (temps, RAM de construction) | moyenne | faible | **étudier**, avec le n° 2 (fiche 3) |
| 8 | **STL de référence minimal** pendant les chantiers 2 et 7 | mapa | faible (filet de sécurité temporaire) | faible | faible | **intégrer temporairement** (fiche 9) |
| 9 | **Diagonale adaptative** | my3dmaps | faible (−3 % d'erreur en moyenne sur SRTM, jusqu'à −7 %) | faible | faible à moyen | **étudier** sur un DEM LiDAR (fiche 2) |
| 10 | **Écriture par bandes** (tuiles à sommets partagés) | phstl, contour_workbench | conditionnel | moyenne à importante | moyen | **étudier**, seulement si le n° 2 ne suffit pas (fiche 8) |
| 11 | **Maillage adaptatif à erreur bornée** (RTIN à bord forcé) | hmterrain, contour_workbench | moyen en taille, **nul en fidélité** | très importante | important | **étudier**, phase 3 (rapport 3, fiche 12) |
| 12 | **Débit poli et `Retry-After`** | Geo3Dprint | faible | faible | faible | **étudier** si des 429 apparaissent (fiche 21) |
| 13 | **Sources NRW et Basse-Saxe** | nrw-terrain-print (documentation officielle) | selon vos projets | moyenne | faible | **à vous de décider** (fiche 24) |

## 3. Réponses aux dix questions

**1. Quel projet contient les meilleures idées algorithmiques ?**
`my3dmaps` (fond fusionné généralisé à un masque quelconque, diagonale selon le relief) pour votre cas ; `hmterrain` pour la rigueur (TIN à erreur bornée et testée).

**2. Quel projet produit potentiellement les meilleurs STL ?**
**Aucun de ceux que j'ai lus ne produit une surface plus juste que la vôtre à partir du même DEM** : tous reposent sur la même grille à 2 triangles par cellule (sauf `mapa`, 4 triangles). Les différences sont la taille (`my3dmaps`, `hmterrain`) et, pour les sources grossières, le rééchantillonnage (fiche 18, qui est un défaut de **votre** pipeline).

**3. Quel projet gère le mieux les gros DEM ?**
`dem-to-stl` : `gdal.Warp` fenêtré, `/vsicurl/`, grille bornée par le pas d'impression. Les autres assemblent les tuiles en mémoire ou plafonnent les entrées. Le vrai goulot chez vous est la mémoire du maillage de sortie (n° 2, 7 et 10).

**4. Quel projet gère le mieux les données imparfaites ?**
`dem-to-stl` (comblement GDAL, masque, ouverture, pincements, repli Copernicus recalé en altitude). Plusieurs dépôts remplissent les vides par 0 sans le signaler (rapport 4, section 3). **Deux faiblesses subsistent chez vous** : la validation du cache et le repli de résolution GSI (n° 3 et 4).

**5. Quel projet possède la meilleure architecture ?**
`geostl` (séparation `Region → Section → Mesher`, adaptateurs de sources) et `contour_workbench` (document d'architecture) ; votre structure plate est adaptée à un outil en ligne de commande. **Rien ne justifie une refonte.**

**6. Quelle est la meilleure idée découverte ?**
Le **fond fusionné par rectangles** (`my3dmaps`) pour l'efficacité. Pour la **qualité de surface**, la découverte la plus utile n'est pas dans un dépôt mais dans la comparaison : `average` n'interpole pas en agrandissement.

**7. Quelle amélioration apporterait le plus pour le moindre coût ?**
**Le choix du rééchantillonnage selon le sens (n° 1)** : [Opinion] quelques lignes, et seul changement qui améliore la surface pour un DEM donné. Viennent ensuite les trois corrections triviales (n° 3, 5, 6).

**8. Quelle amélioration serait la plus difficile mais potentiellement la plus intéressante ?**
Le **maillage adaptatif à erreur bornée sur emprise libre** (n° 11). Gain de taille attendu de l'ordre de 2 à 10 fois sur des données lissées (simulation RTIN : 10 à 53 % des triangles conservés pour un lissage σ=2 px), mais **sans amélioration de fidélité** et avec une complexité très importante.

**9. Quelles fonctionnalités ne valent pas la complexité supplémentaire ?**
Maillage « 4 triangles par pixel » (`mapa`), scission de sommets aux coins en damier (`my3dmaps`, non valable en STL), raffinement de Delaunay sans borne d'erreur (`jkittner`), « frontière exacte » par centroïde (`TopoTwin`), fusion planaire à tolérance (`Geo3Dprint`), booléen `manifold3d` pour découper (`terrology`), courbure terrestre (`dem2stl` en R), dégradation automatique de la résolution (`contour_workbench`), `trimesh` et `numba` comme dépendances, tuilage pour plateau d'impression, fonctions de style (lissage, lacs aplatis), manifeste JSON.

**10. Y a-t-il quelque chose de fondamentalement meilleur ailleurs qui remette en cause l'architecture ?**
**Non.** L'unique évolution structurelle envisageable est un maillage non régulier pour le dessus (n° 11), à n'engager qu'après un benchmark concluant sur vos données. Les n° 2, 7 et 10 touchent à la construction du maillage mais pas à l'architecture générale.

## 4. Feuille de route proposée

### Phase 1 — gains faciles (à faible risque)

1. **Rééchantillonnage selon le sens** avec réglage `--resample auto|average|bilinear|cubic`. À valider d'abord sur un DEM global (Copernicus ou SRTM) à petit rayon, en regardant `bilinear` contre `cubic` (dépassements près des falaises, hypothèse).
2. **Validation des téléchargements** : rattraper `http.client.HTTPException`, contrôler le type de fichier avant mise en cache, vérifier la longueur annoncée, et revalider un fichier déjà en cache (GeoTIFF/PNG lisible).
3. **Boucle des pincements** jusqu'à stabilité (avec plafond de sécurité large).
4. **`check_stl.py`** : volume signé > 0 et comptage des triangles dégénérés.

### Phase 2 — améliorations importantes

5. **Repli de résolution GSI** par couche (demander à `min(z, zmax)` et suréchantillonner avec la méthode du point 1), avec test sur une zone hors DEM1A.
6. **Fond fusionné** : prototype de la couverture gloutonne + éventails, avec tests de propriétés existants, **puis test sur PrusaSlicer, Cura, OrcaSlicer et votre CAM** (triangles très allongés, jusqu'à ~1 000:1). Dans le même chantier : parois vectorisées et indexation native. Un STL de référence minimal protège la transition.

### Phase 3 — recherche et expérimentation (benchmarks avant décision)

7. **Diagonale adaptative** : mesurer sur un DEM LiDAR (1 m) rééchantillonné à votre pas.
8. **Maillage adaptatif** : appliquer RTIN avec bord forcé à un de vos DEM déjà rééchantillonnés ; ne poursuivre que si le facteur de réduction dépasse 3 à une tolérance inférieure à l'épaisseur de couche.
9. **Écriture par bandes** : seulement si les tailles visées (par exemple Ø250 mm à 0,1 mm, ≈ 19,6 M de triangles sans fond fusionné) restent bloquées par la garde RAM après le point 6.
10. **Contrôle de débit** : seulement en cas de 429.

### À ne pas faire

Voir la réponse 9. En résumé : ne pas ajouter `trimesh`, `numba`, `manifold3d` ; ne pas remplacer la grille par un maillage « 4 triangles » ; ne pas scinder les sommets ; ne pas chercher un contour « exact » par booléen à 0,2 mm de pas.

## 5. Licences des idées retenues

| Idée | Source | Licence du dépôt | Ce que je recommande |
|---|---|---|---|
| Fond fusionné, diagonale, volume signé | my3dmaps | **MIT** | **Réimplémentation indépendante** à partir de la description ; un emprunt de code serait possible dans un projet GPL-3 avec reprise de l'avis MIT, mais inutile |
| Fond en éventail, passes de pincements | Terrain2STL | **aucune** | Idée seule ; ne rien copier |
| Rééchantillonnage selon le sens | dem2stl-relief | **aucune** | Idée seule (pratique courante, un commentaire de code) |
| Validation des téléchargements | nrw-terrain-print | **PolyForm Noncommercial** | Idée seule ; incompatible avec une reprise de code |
| Débit poli | Geo3Dprint | **aucune** | Idée seule |
| TIN à erreur bornée | hmterrain | **MIT** | Algorithme publié (Garland et Heckbert, 1995) ; implémenter depuis la publication |
| Repli GSI | japan-3d-terrain | **MIT** | Idée seule, rien à reprendre |

[Opinion] Aucune des phases 1 et 2 ne demande de reprendre du code : ce sont des idées simples, ou des corrections de votre propre code. Ce n'est pas un avis juridique.

## 6. Corrections apportées à mes rapports en cours d'audit

| Où | Ce qui était faux ou erroné | Correction |
|---|---|---|
| Rapport 1 | « Aucun dépôt lu ne choisit la diagonale selon le relief » | `my3dmaps` le fait (rapport 2) |
| Rapport 1 | Licence de DEM.Net « Apache-2.0 » | MIT modifiée avec restriction commerciale (rapport 3) |
| Rapports 1 et 2 | rayshader écrit le STL avec `rgl` ; simplification « déléguée » | `write_stl` propre à rayshader ; simplification par `terrainmeshr` (rapport 3) |
| Rapport 3 | Première simulation RTIN (97–100 % de triangles conservés) | Bogue détecté par un test de validation (plan incliné) ; chiffres corrigés |

## 7. Ce qu'il reste à trancher avec vous

1. **Quelles améliorations implémenter**, dans quel ordre (je propose la phase 1, puis le prototype de la phase 2).
2. **Vos tailles cibles** (diamètre, pas d'impression) : elles décident si le fond fusionné suffit ou si l'écriture par bandes devient utile.
3. **Un de vos DEM LiDAR et un de vos DEM grossiers** pour trancher les benchmarks de la phase 3 (diagonale, RTIN, `bilinear` contre `cubic`).
4. **Faut-il lire les 48 dépôts restants** plus à fond ? Je peux faire une lecture rapide ciblée (par exemple les 16 « redondants probables ») pour confirmer qu'ils n'apportent rien.
