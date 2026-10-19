# Audit comparatif dem-to-stl — Rapport 3 : lot 2 (gros DEM, autres écosystèmes, CNC)

Date : 9 octobre 2026. **Aucun fichier de votre projet n'a été modifié.** Les mesures ont été faites sur des copies de travail, avec des scripts écrits pour l'occasion.

Convention : **[Fait]** constaté dans le code ou mesuré ; **[Déduction]** ; **[Hypothèse — à vérifier par benchmark]** ; **[Opinion]**.

---

## 0. Résumé

1. **Le lot 2 contient la seule approche réellement différente rencontrée jusqu'ici : le TIN adaptatif à erreur bornée** (`hmterrain`, licence MIT). Contrairement au raffinement de `jkittner`, il garantit qu'aucun pixel n'est à plus d'une tolérance donnée du maillage. Ma simulation (RTIN, borne haute d'un TIN glouton) suggère de **conserver 10 à 53 % des triangles de dessus sur des données lissées (σ=2 px, tolérance 0,1 à 0,5 %)**, mais 90 % ou plus sur des données brutes à tolérance fine. Gain de taille uniquement, **pas de gain de fidélité**. Classe C ; à étudier après le fond fusionné du rapport 2.
2. **Gros DEM : aucun dépôt ne fait mieux que votre `gdal.Warp`.** Les seules idées utiles viennent de `contour_workbench` : une estimation mémoire avant construction avec dégradation automatique, et une construction par tuiles qui partagent leurs sommets de raccord.
3. **Données imparfaites : vous restez devant.** `DEM.Net` traite l'altitude 0 comme invalide dans sa fonction de filtrage, ce qui peut retirer le niveau de la mer.
4. **Le fond en éventail se confirme** : `satprint` est le quatrième dépôt à l'utiliser, après `Terrain2STL`, `mapa` et `tok`. Cela renforce la fiche 1 du rapport 2.
5. **Trois corrections de mes rapports précédents** (section 1), dont une licence mal détectée (`DEM.Net`) et une erreur sur `rayshader`.

## 1. Corrections des rapports précédents

| Point | Ce que j'avais écrit | Ce qui est vrai |
|---|---|---|
| **Licence de DEM.Net** (rapport 1, tableau) | Apache-2.0 | **MIT modifiée avec restriction commerciale** : usage personnel ou entreprises de moins de 100 000 $ de chiffre d'affaires annuel. Ce n'est pas une licence libre. Idées seulement, aucun emprunt de code. La détection automatique avait lu un fichier de licences tierces. |
| **rayshader, export STL** (rapport 1) | « appelle l'écriture STL de la bibliothèque `rgl` » | `save_3dprint` appelle **`write_stl` de rayshader**, qui relit les sommets de la scène `rgl` et écrit chaque triangle dans une boucle R. Le maillage lui-même vient de `rgl` ou du paquet `terrainmeshr`. |
| **rayshader, simplification** (rapports 1 et 2) | « déléguée à un autre paquet » | Exact : `plot_3d` appelle `terrainmeshr::triangulate_matrix(maxError, maxTriangles)`. |

Les rapports 1 et 2 ont été annotés en conséquence.

## 2. Dépôts du lot 2 : niveau de lecture et verdict

« Lu » = code de maillage lu en détail ; « partiel » = fonctions clés lues ; « survol » = README, signatures, recherches ciblées.

| Dépôt | Licence | Lecture | Apport pour dem-to-stl | Verdict |
|---|---|---|---|---|
| **hmterrain** (shanto462) | MIT | `tin.rs` (en-tête documenté), `mesh.rs` (signatures), README, tests | TIN glouton à erreur bornée, socle fermé, tests de propriétés | **Fiche 12** |
| **contour_workbench** (DavidIngraham) | GPL-3 | `HLD.md` lu, `lib.rs` (2 650 lignes) seulement cherché | Triangulation contrainte sur emprise libre ; estimation mémoire ; tuiles à sommets de raccord | Fiche 12, 13 |
| **terrology** (twigley) | Apache-2.0 | `_clip_mesh_to_polygon` lu ; reste survolé | Découpe exacte au polygone par booléen `manifold3d` | **Fiche 14** |
| **satprint** (suchanek) | MIT | `heightmap_to_mesh` lu | Maillage vectorisé, fond en éventail | Fiche 1 (renfort) |
| **DEM.Net** (C#) | **MIT modifiée, restriction commerciale** | `MeshService.cs`, triangulation, NoData, lecture GeoTiff : survol ciblé (34 000 lignes) | Rien à reprendre | Écarter |
| **rayshader** (R) | GPL-3 (`DESCRIPTION`) | `save_3dprint.R`, `write_stl.R`, appel `plot_3d` | Rien : export de la scène `rgl` | Écarter |
| **dem2stl** (kevinstadler, R) | MIT | `dem2stl.R` : survol | Forme arbitraire par masque NA ; courbure terrestre | Fiche 15 |
| **DEM_STL_carver** (tfoss) | **aucune** | survol (comblement de vides, lissage) | Rien d'original ; flux Fusion 360 | Écarter |
| **heightmap_to_stl** (kilimnik) | GPL-3 | survol (216 lignes) | Rien (grille complète, boucles) | Écarter |
| **BlenderGIS** | GPL-3 | options d'import DEM seulement (33 000 lignes) | Option « step » (sous-échantillonnage) ; comblement du nodata | Écarter |
| **MapLab3D** (MATLAB) | AGPL | recherche de mots-clés | Cartes OSM multicolores, hors cœur | Écarter |
| **Topographic-3D-Printing** (DomPrott, R) | MIT | README | Application Shiny, villes | Écarter |

## 3. Fiches comparatives

### Fiche 12 — TIN adaptatif à erreur bornée · Classe C (taille), difficulté très importante

**Projet concurrent.** *hmterrain* ([Fait] en-tête de `tin.rs`) implémente l'insertion gloutonne de Delaunay de Garland et Heckbert (1995) : on part de 2 triangles, on insère à chaque étape le pixel de plus grande erreur verticale, on rétablit Delaunay par échange d'arêtes, avec des prédicats entiers exacts. Le paramètre `--max-error` fixe la distance verticale maximale entre **tout pixel** et le maillage ; `--max-triangles` et `--max-vertices` donnent des budgets. Défaut : 0,1 % de l'amplitude du relief. `--base` ajoute parois et fond plan. [Fait, annoncé par le projet, non vérifié par moi] 3 à 45 % des pixels deviennent des sommets selon le relief.

*contour_workbench* ([Fait] `HLD.md`) fait de la triangulation **adaptative et contrainte** sur des emprises libres (spline, cercle à coins arrondis, polygone), avec détection mémoire préalable. C'est une preuve de faisabilité pour un contour non rectangulaire. Le projet est en Rust/WASM, GPL-3.

**Mon projet.** [Fait] Grille régulière, 2 triangles par cellule, murs en escalier ou lissés pour le cercle.

**Différence.** [Fait] `hmterrain` ne gère que des rectangles, et `rayon` n'y apparaît que dans `heightmap.rs` (chargement, validation) : la triangulation est séquentielle d'après mon grep. [Déduction] Contrairement à `jkittner` (rapport 1), l'erreur est **bornée** et vérifiée par test (`error_limit_is_respected_and_reduces_size`, `random_heightmaps_stay_valid_and_within_error`).

**Mesure** (simulation RTIN/Martini en Python, **borne haute** de ce qu'un TIN glouton garde ; part des triangles de la grille de pixels complète conservée ; erreur max exprimée en % de l'amplitude du relief ; 513×513) :

| Données | 0,05 % | 0,10 % | 0,25 % | 0,50 % |
|---|---|---|---|---|
| SRTM3 Maine brut, 4 fenêtres (entier en mètres, bruit de pixel) | 90–92 % | 90–92 % | 67–80 % | 44–69 % |
| SRTM3 lissé σ=1 px (2 fenêtres) | 82–84 % | 66–71 % | 35–42 % | 18–22 % |
| SRTM3 lissé σ=2 px (2 fenêtres) | 66–74 % | 44–53 % | 20–25 % | 10–13 % |
| fBm synthétique β=3,0 / 2,2 / 1,6 | 98 / 99,5 / 99,7 % | 96 / 99 / 99,5 % | 91 / 98 / 99 % | 83 / 95 / 97 % |

**Validation de ma simulation** : plan incliné → 2 triangles ; bruit blanc à tolérance 0 → 524 288 (grille complète). **Elle avait un bogue** que ces tests ont révélé : une première version donnait 97–100 %, chiffres que je n'ai pas retenus.

**Lecture.**
- [Déduction] Le gain dépend du lissage au niveau du pixel. Un DEM brut (quantifié, bruité) ne se simplifie presque pas. Un DEM moyenné se simplifie beaucoup.
- [Hypothèse — à vérifier par benchmark] Votre grille d'impression moyenne beaucoup de pixels source (par exemple des dizaines de pixels LiDAR par pixel imprimé) : vos données ressembleraient alors plutôt aux cas lissés, donc à une conservation de l'ordre de 20 à 53 % des triangles à 0,1–0,25 % de tolérance. À mesurer sur un de vos DEM rééchantillonnés.
- [Opinion] À 0,25 % du relief, une plaque de 40 mm de relief tolère 0,1 mm : de l'ordre de l'épaisseur d'une couche. Pour un usage CNC, [Hypothèse] votre CAM rééchantillonne probablement le STL en grille, donc la réduction ne change pas le résultat usiné.

**Gain.** Moyen (taille, RAM, temps de tranchage) ; **nul pour la fidélité**, qui devient bornée au lieu d'être exacte aux nœuds. **Difficulté** très importante : emprise non rectangulaire (murs en escalier ou polygone), cohérence parois/fond, 3 exports, tests. **Risque** important.

**Piste plus simple** ([Hypothèse], non prototypée) : un RTIN dont l'erreur est forcée à l'infini sur les nœuds du bord du masque, pour que le contour reste exactement celui de votre escalier. Il demande un raster de taille 2^k+1 (remplissage) et produit une triangulation sans fissure.

**Verdict : étudier, Phase 3, seulement après la fiche 1 du rapport 2** (qui donne déjà −50 % sans toucher au dessus) et si le benchmark sur vos données montre un facteur ≥ 3 à tolérance sous l'épaisseur de couche.

### Fiche 13 — Estimation mémoire avec dégradation automatique et tuiles de raccord · Classe D

**Projet concurrent.** *contour_workbench* ([Fait] `HLD.md`) estime le pic mémoire avant de construire ; au-delà de 72 % du budget, le « plein détail » automatique est remplacé par une cible d'erreur échantillonnée. Les faces à plein détail sont émises par tuiles de 128×128 cellules qui partagent leurs sommets de raccord. [Fait] Limites : 2 degrés par côté et 2 millions d'échantillons.

**Mon projet.** [Fait] Garde-fou RAM (≈ 220 octets par triangle) qui **refuse** ; `--dry-run`.

**Différence.** [Opinion] Votre refus explicite avec estimation est cohérent avec un outil en ligne de commande dont l'utilisateur maîtrise le pas. Une dégradation automatique ne se justifie que dans une application web sans contrôle de l'utilisateur. Les tuiles de raccord sont l'équivalent de la construction par bandes de la fiche 8 du rapport 2.

**Gain** faible. **Verdict : ne pas intégrer** la dégradation automatique ; l'idée des tuiles à sommets partagés est à retenir si vous implémentez l'écriture par bandes.

### Fiche 14 — Contour exact par booléen sur un prisme · Classe D

**Projet concurrent.** *terrology* (Apache-2.0, [Fait] `_clip_mesh_to_polygon`) intersecte le solide avec un prisme construit sur le polygone (multipolygone par union), via `manifold3d` appelé par `trimesh`. Le contour est **exactement** le polygone, sans escalier, et le résultat est un solide valide.

**Mon projet.** [Fait] Contour du masque en escalier (au pas d'impression) ; bord circulaire lissé en option ; polygone en escalier.

**Différence.** [Déduction] Seule approche rencontrée qui donne un contour exact quelle que soit la forme. Elle exige `trimesh`, `manifold3d`, `shapely`, et un booléen sur le maillage complet : [Hypothèse — à vérifier par benchmark] coûteux en temps et en RAM sur plusieurs millions de triangles. Le gain visible n'existe que si le pas est grossier par rapport à la précision recherchée : à 0,2 mm, l'escalier est de l'ordre du pas d'impression.

**Gain** faible à moyen. **Difficulté** importante, **dépendances** lourdes. **Verdict : ne pas intégrer** pour votre usage actuel. À reconsidérer seulement si vous découpez des polygones (frontières) à pas grossier.

### Fiche 15 — Courbure terrestre et formes arbitraires · Classe D

**Projet concurrent.** *dem2stl* (R, MIT, [Fait]) construit des modèles de forme arbitraire à partir d'un masque NA (seule la plus grande composante est retenue) et projette l'altitude sur un **sphéroïde terrestre**, avec une base courbe ou plane tangente.

**Mon projet.** [Fait] Reprojection UTM plane ; composante principale conservée en mode île.

**Différence.** [Déduction] La courbure n'a d'effet perceptible que pour de très grandes emprises : pour un diamètre de 250 km, elle fait environ 1,2 km de flèche (r²/2R), soit ≈ 1,2 mm à l'échelle 1:1 000 000. Vos reliefs de montagne ne sont pas concernés. **Verdict : ne pas intégrer** sauf projet de modèle à l'échelle d'un pays ou d'un continent.

### Fiche 16 — Traitement du nodata chez les autres · aucune idée à reprendre

| Dépôt | Traitement |
|---|---|
| **DEM.Net** | [Fait] `TriangulateHeightMapFiltered` marque invalides les altitudes absentes **ou égales à 0**, puis retire les régions invalides par remplissage. [Déduction] Le niveau de la mer disparaît. |
| **hmterrain** | [Fait] Rejette les valeurs non finies avec une erreur. |
| **BlenderGIS** | [Fait] Option « Fill nodata values » (interpolation) pour une texture de déplacement. |
| **DEM_STL_carver** | [Fait] Comblement de vides puis lissage gaussien σ=1 répété. |
| **dem2stl (R)** | [Fait] Comble les trous intérieurs en leur donnant l'altitude 0 pour éviter des parois supplémentaires. |
| **dem-to-stl** | `gdal.FillNodata` (50 px), masque des cellules restantes, ouverture morphologique, comblement des pincements, composante principale. |

[Opinion] Votre traitement est le plus complet de ce que j'ai lu.

### Fiche 17 — Validation et tests (lot 2)

[Fait] *hmterrain* a 40 tests marqués `#[test]`, dont des tests de propriétés (plan = 2 triangles, erreur zéro = chaque pixel reproduit, erreur limitée sur des cartes aléatoires) ; ses tests de socle définissent un contrôle d'arêtes appariées en sens opposés et un calcul de volume signé ([Déduction] utilisés ensemble dans `solid_base_is_closed_and_outward`, dont je n'ai pas lu le corps). *contour_workbench* (37 tests Rust) valide « fermé et orienté » et, avant livraison d'un 3MF, les indices, les parties OPC, la syntaxe XML, les unités, les identifiants d'objets. *satprint* (131) et *terrology* (238) testent largement mais sur l'application autour du maillage.

[Opinion] Aucune méthode supérieure à la vôtre. **Idée mineure (D)** : un contrôle de structure du 3MF (parties obligatoires, XML bien formé, unités, références de construction). [Hypothèse] Vos tests le font peut-être déjà : à vérifier dans `tests/`.

## 4. Constats transversaux (lots 1 et 2)

**Qualité du maillage.** [Fait] Le fond en éventail est utilisé par 4 dépôts (Terrain2STL, mapa, tok, satprint) pour des plaques rectangulaires. Seul *my3dmaps* l'étend à un masque quelconque. [Déduction] Pour votre cas (disque ou polygone), la fusion par rectangles reste la bonne généralisation.

**Gros DEM.** [Fait] Dans les lots 1 et 2, aucune stratégie de lecture par blocs d'un DEM n'est supérieure à `gdal.Warp` vers une grille bornée. *DEM.Net* lit les GeoTIFF par tuiles ou lignes de balayage via une bibliothèque TIFF ; *contour_workbench* limite à 2 millions d'échantillons ; *BlenderGIS* propose un « pas de pixels » ([Hypothèse] saut de pixels sans moyenne, non vérifié). [Opinion] Votre rééchantillonnage par moyenne sur la fenêtre évite le crénelage qu'un simple saut de pixels produirait.

**Performances.** [Fait] Aucun dépôt lu ne parallélise la construction du maillage. *hmterrain* ne parallélise que le chargement. [Hypothèse — à vérifier par benchmark] La construction par bandes de la fiche 8 (rapport 2) se prêterait à `multiprocessing`, mais le gain n'est pas évalué.

**Architecture.** [Fait] *contour_workbench* sépare un cœur géométrique déterministe, une frontière WASM et une application, avec un document d'architecture (`HLD.md`). [Opinion] Bon modèle de documentation, mais votre structure plate à 9 fichiers suffit pour un outil en ligne de commande.

**Dépendances et portabilité.**

| Dépôt | Dépendances notables | Remarque |
|---|---|---|
| hmterrain | Rust ≥ 1.88, `rayon` | Binaires multiplateformes ; PNG/TIFF/JPEG seulement, **pas de géoréférencement** |
| contour_workbench | Rust/WASM, `manifold-3d` 3.5.4 | Navigateur ; GeoTIFF géographiques seulement |
| terrology | `trimesh`, `manifold3d`, `shapely`, `geopandas` | Lourd |
| satprint | Python ≥ 3.12 | |
| DEM.Net | .NET | Licence restrictive |
| vous | GDAL, NumPy, SciPy | Une seule dépendance lourde (GDAL) |

**Licences.** MIT : hmterrain, satprint, dem2stl (R), my3dmaps, geostl, gpx3mf, mapa, phstl. Apache-2.0 : terrology. GPL-3 : contour_workbench, BlenderGIS, heightmap_to_stl, rayshader (compatibles avec votre GPL-3 pour un éventuel emprunt, mais en Rust, R ou C++, donc surtout conceptuels). **Restrictive : DEM.Net.** AGPL : MapLab3D. **Absente : DEM_STL_carver.**

## 5. Réponses provisoires aux dix questions (lots 1 et 2)

1. **Meilleures idées algorithmiques** : *my3dmaps* (fond fusionné) et *hmterrain* (TIN à erreur bornée).
2. **Meilleurs STL potentiels** : à même DEM, aucun ne produit une surface plus fidèle que la grille ; *hmterrain* produit des maillages plus petits à erreur bornée, pas plus justes.
3. **Gros DEM** : `dem-to-stl` (GDAL), à égalité avec rien d'autre dans ces deux lots.
4. **Données imparfaites** : `dem-to-stl`.
5. **Architecture** : *geostl* et *contour_workbench* (par la documentation) ; pas de refonte à faire.
6. **Meilleure idée découverte** : le fond fusionné par rectangles (*my3dmaps*).
7. **Meilleur gain pour le moindre coût** : boucle des pincements jusqu'à stabilité, volume signé dans `check_stl`, puis fond fusionné.
8. **Plus difficile mais potentiellement intéressante** : maillage adaptatif à erreur bornée sur emprise libre (fiche 12).
9. **Ne valent pas la complexité** : booléen de découpe (`manifold3d`), courbure terrestre, dégradation automatique, TIN glouton écrit à la main avant benchmark.
10. **Remise en cause de l'architecture** : **non**. Si le benchmark de la fiche 12 est concluant, la seule évolution structurelle serait un maillage non régulier pour le dessus.

## 6. Feuille de route provisoire (lots 1 et 2)

**Phase 1 — gains faciles** : boucle des pincements jusqu'à stabilité ; volume signé et triangles dégénérés dans `check_stl.py`.

**Phase 2 — améliorations importantes** : fond fusionné par rectangles (prototype + test PrusaSlicer/Cura/Orca/CAM), parois vectorisées dans le même chantier, STL de référence minimal pendant la transition.

**Phase 3 — recherche** : benchmark de la fiche 12 sur un de vos DEM (RTIN avec bord forcé) ; diagonale adaptative (rapport 2) ; écriture par bandes avec tuiles à sommets partagés.

**À ne pas faire** : triangulation adaptative sans borne d'erreur, maillage 4 triangles par pixel, scission de sommets, booléens `manifold3d`/`trimesh`, courbure terrestre, dégradation automatique de la résolution.

## 7. Mesures effectuées dans ce rapport

1. **RTIN** : implémentation de l'algorithme Right-Triangulated Irregular Network (erreur propagée du plus petit au plus grand triangle, comptage des feuilles pour une tolérance). Validée sur plan incliné et bruit blanc ; bogue corrigé. Appliquée à SRTM3 Maine (4 fenêtres), SRTM lissé, fBm synthétique, et SRTM bruité. **C'est une borne haute pour un TIN glouton ; je n'ai pas exécuté `hmterrain` (Rust absent de mon environnement).**
2. **Aucune mesure sur vos propres DEM** : elles seraient plus décisives que tout ce qui précède.

## 8. Suite

Lot 3 (sources nationales et données réelles), 13 dépôts : `nrw-terrain-print`, `Geo3Dprint`, `jp_geo_generator`, `japan-3d-terrain`, `dem2stl-relief`, `portland-topo-map`, `geo-printing`, `DEMtoPrint`, `TerrainToSTL`, `PeakPrinter`, `elevation-to-stl`, `convert-dem-to-stl`, `turkeytopgraphyDEM2STL`.

Parties du lot 2 **non lues** : `terrain-core/src/lib.rs` de contour_workbench (2 650 lignes), `mesh.rs` et `tin.rs` de hmterrain au-delà de l'en-tête (ma lecture de l'algorithme vient du commentaire de tête et des tests), `ElevationService.cs` de DEM.Net.
