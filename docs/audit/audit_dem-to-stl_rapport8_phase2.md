# Audit comparatif dem-to-stl — Rapport 8 : phase 2

Date : 11 octobre 2026. Comme pour le rapport 7, les modifications sont faites **sur une copie** de votre projet : vous recevez des patchs et un zip, votre dépôt n'est pas touché. Les fichiers de la phase 1 déjà livrés restent valables ; les nouveaux fichiers cumulent les deux phases.

Convention : **[Fait]** constaté dans le code ou mesuré ; **[Déduction]** ; **[Hypothèse — à vérifier par benchmark]** ; **[Opinion]**.

---

## 0. Résumé

- **Suite de tests : 365 réussis + 1 ignoré** (165 à l'origine, 316 après la phase 1, **+49 en phase 2**). Le patch cumulatif s'applique à votre archive d'origine et la suite passe dans le zip décompressé.
- **Repli de résolution GSI : fait.** Une tuile absente de la couche fine est remplacée par la tuile des couches suivantes (zoom 15, puis 14), interpolée. Avant, la zone entière pouvait ressortir « aucune tuile ».
- **Parois vectorisées : fait, et c'est le gain de temps le plus net** : à Ø250 mm / 0,2 mm, **STL 5,2 s → 2,6 s**, **3MF 23,1 s → 15,4 s**, avec des **fichiers identiques octet pour octet** à ceux de votre code d'origine (STL, 3MF hors date, OBJ).
- **Non-régression : faite.** Empreintes SHA-256 de sorties calculées avec votre code d'origine : le maillage par défaut n'a pas changé d'un bit.
- **Indexation native : non faite, par décision mesurée** (section 4).
- **Le fond fusionné reste un prototype** (`--bottom merged`), en attente de vos essais dans votre trancheur et votre FAO. **Non testé** : réseau GSI réel, vos DEM.

## 1. Modifications de la phase 2

9 fichiers, +509 lignes, −69 (dont 353 lignes de tests dans 3 nouveaux fichiers).

| Fichier | Changement |
|---|---|
| `sources.py` | `fetch_gsi` renvoie une liste de mosaïques (la plus grossière d'abord) ; `_gsi_collect`, `_gsi_write_mosaic` |
| `raster.py` | `source_pixel_size_m(reduce="min"|"max")` ; `warp_window` utilise la source la plus grossière |
| `mesh.py` | `build_plain_walls` : parois sans boucle Python |
| `mesh_export.py` | `_positional_lines` : chaque valeur distincte formatée une seule fois |
| `README.md`, `README.fr.md` | repli GSI et règle « la plus grossière décide » |
| `tests/test_gsi_fallback.py` | 10 tests |
| `tests/test_walls_vectorized.py` | 26 tests |
| `tests/test_export_formatting.py` | 13 tests |

## 2. Repli de résolution GSI

### 2.1 Le problème (rapport 4, fiche 20)

[Fait] À un zoom de travail de 16 ou 17 (pixel de travail ≲ 3 m), `_gsi_tile` n'interrogeait que `dem1a_png`. Si cette couche n'avait pas la tuile, aucune couche plus grossière n'était essayée : trou dans la mosaïque, ou `NoCoverage` pour toute la zone puis repli Copernicus à 30 m.

### 2.2 Le changement

- La mosaïque fine est construite comme avant. Pour les tuiles absentes, `fetch_gsi` demande la **tuile ancêtre** aux zooms plus grossiers (15 puis 14), en reprenant `_gsi_tile` (mêmes couches, même cache, mêmes contrôles de contenu que pour la phase 1).
- Chaque zoom de repli produit **sa propre mosaïque GeoTIFF** ; la fonction renvoie une liste, **de la plus grossière à la plus fine**. `gdal.Warp` les superpose dans l'ordre : la plus fine l'emporte là où elle a des données. Aucun code d'interpolation maison, aucun raccord de tuiles à gérer.
- Aucune requête de repli quand toutes les tuiles fines existent (testé). Au zoom 14 il n'y a pas de repli : toutes les couches y sont déjà essayées.
- **Conséquence pour `raster.py`** : avec des mosaïques de résolutions différentes, choisir le rééchantillonnage d'après la source la **plus fine** aurait répliqué par blocs la mosaïque grossière. `warp_window` décide donc d'après la source la **plus grossière** (`reduce="max"`). Si une source est agrandie, `bilinear` pour toutes.

### 2.3 Mesures

- [Fait] **Réseau simulé, plan incliné** (défini en mètres Web Mercator, identique à tout zoom), DEM1A absent, DEM5A au zoom 15 : **l'altitude restituée par la mosaïque de repli puis `warp_window` écarte de la vérité de moins de 0,05 m** ; avec `average` (l'ancien comportement du rééchantillonnage) l'écart dépasse 0,5 m (terrasses de 4,8 m sur une pente de 0,3 à 0,5).
- [Fait] **Bout en bout par la ligne de commande** (`--source gsi`, pixel de 1 m, aucune tuile DEM1A) : le journal indique « repli de résolution au zoom 15 » et « ré-échantillonnage bilinear », `check_stl` valide tout, le dénivelé du STL correspond à celui attendu à 2 % près. [Fait] Avec le code de la phase 1, la même commande échoue (« aucune tuile »).
- [Fait] Les **7 de ces 10 tests échouent avec le code de la phase 1** ; les 3 autres vérifient des comportements inchangés (aucune requête de repli si tout existe, `NoCoverage` si rien n'existe, pas de repli au zoom 14).

### 2.4 Effet de bord mesuré de la règle « la plus grossière décide »

[Fait] Quand `--resample auto` choisit `bilinear` à cause d'une source agrandie, **les sources fines de la même liste sont réduites avec `bilinear` au lieu de `average`**. Mesure (source fine de 1 m réduite à 5 m, comparée à la moyenne exacte par blocs) : écart de **7,8 %** de l'écart-type du relief sur un relief rugueux synthétique, **2,8 %** sur un relief doux (`cubic` : 5,7 % et 1,9 %). C'est faible, mais ce n'est pas nul.

[Déduction] Ce cas ne concerne que les **listes mêlant des résolutions** et seulement quand le pixel de travail est entre les deux : mosaïques GSI avec repli, et listes USGS qui mêlent 1 m et 1/3″ (≈10 m). Pour des sources homogènes, le comportement est inchangé. [Opinion] Le compromis est acceptable : la source grossière y gagne une interpolation à la place de blocs. Pour revenir à l'ancien choix : `--resample average`.

### 2.5 Limites

- [Fait] **Je n'ai pas testé le réseau GSI réel** (aucun accès). La logique est testée avec des tuiles simulées.
- **Couverture de DEM1A : non vérifiée.** Ma recherche n'a trouvé que des mentions indirectes de différences selon les régions. Le code n'en dépend pas (il traite toute tuile absente) ; j'ai **retiré du README** une phrase affirmant une couverture partielle, que j'avais d'abord écrite.
- [Hypothèse — à vérifier] Les valeurs d'une tuile de repli (couche DEM5A) et celles d'une tuile fine voisine (DEM1A) peuvent présenter un petit décalage à la jonction : non mesuré.

## 3. Parois vectorisées

### 3.1 Ce que j'ai mesuré avant de coder

[Fait] Profil ligne à ligne de `build_mesh_from_arrays` à Ø250 mm : les deux boucles Python des parois parcourent **tous** les nœuds de la grille (1,56 M × 2 passes) pour n'en retenir qu'environ 10 000 sur le bord ; elles représentent **environ 90 % du temps de cette fonction** sous le profileur (qui gonfle les chiffres : 8,8 s mesurés pour ~4,5 s réels). Mon hypothèse du rapport 2 (« quelques secondes ») était donc **sous-estimée**.

### 3.2 Le changement

`mesh.build_plain_walls` : mêmes arêtes (une cellule incluse contre une cellule exclue ou hors grille), même orientation, **même ordre** des triangles que la double boucle (arêtes verticales ligne par ligne, puis horizontales colonne par colonne).

### 3.3 Preuves d'équivalence

- [Fait] **Tableau contre tableau** (`np.array_equal`) avec l'ancienne boucle, conservée dans les tests comme implémentation de référence : 20 masques aléatoires, plus cinq cas limites (une cellule, grille pleine, bande contre le bord, damier, trou central) et le masque vide.
- [Fait] **Bit à bit avec votre code d'origine** : j'ai chargé `mesh.py` du commit de base et comparé les sorties de `build_mesh_from_arrays` sur 12 masques aléatoires et un disque avec `--smooth-wall` : identiques, même ordre.
- [Fait] **Non-régression permanente** : le fichier `test_walls_vectorized.py` contient des empreintes SHA-256 de sorties du code d'origine sur 4 entrées entièrement déterministes (disque, anneau, forme en L, disque avec bord lissé ; entiers et multiples exacts, sans générateur aléatoire, donc stables d'une version de NumPy à l'autre). C'est le « STL de référence minimal » de la feuille de route.
- [Fait] **Test de mutation** : inverser l'orientation des parois fait échouer **24 de ces 26 tests** (les 2 qui passent encore sont les tests de bords de masque vide et de masque sans paroi extérieure, indépendants de l'orientation).
- [Fait] **Fichiers à Ø250 mm** : STL identique octet pour octet à celui du code d'origine (`cmp`).

## 4. Formatage de l'export 3MF, et indexation native

### 4.1 Ce que j'ai mesuré

[Fait] 3MF à Ø250 mm, grille, code d'origine : **23,1 s**, contre 5,2 s pour le STL. Profil de l'écriture du 3MF (4,9 M de triangles) : **16 s** dans 7,4 M d'appels à `_positional` (formatage valeur par valeur), 3,6 s de compression, 2,4 s de lignes de faces, **3,5 s** de soudure des sommets (`index_mesh`).

### 4.2 Décision sur l'indexation native

**Non faite.** [Fait] `index_mesh` ne pèse que 3,5 à 3,8 s sur ~21 s de surcoût du 3MF : une indexation native gagnerait au mieux 16 % de l'export 3MF et exigerait de refondre la construction du maillage. [Opinion] Mauvais rapport bénéfice/complexité, tel que le demandait votre cahier des charges.

### 4.3 Ce que j'ai trouvé à la place

[Fait] Dans `_vertex_lines`, le test `"e" in text` censé détecter une notation à exposant est **toujours vrai pour le 3MF**, parce que le gabarit `<vertex …>` contient la lettre « e ». La voie lente (sans exposant) est donc toujours utilisée pour le 3MF, alors que le commentaire laisse croire qu'elle est rare. Cette voie produit les valeurs les plus courtes (« 12.4 » plutôt que « 12.3999996 »).

**Ce que j'ai fait** : **conserver exactement ce texte** et accélérer la voie lente en ne formatant qu'une fois chaque valeur distincte (sur une grille régulière, x et y n'ont que quelques centaines de valeurs distinctes pour des millions de sommets). Le comportement est documenté dans la fonction. `-0.0` et `0.0` restent distincts (unicité sur les bits).

[Fait] **Le texte produit est identique caractère pour caractère** : tests (13) contre l'ancien formatage sur des tableaux aléatoires, une grille et des valeurs particulières (zéros signés, 1e-30, 3,4e38), pour les deux gabarits ; à Ø250 mm, **XML du modèle 3MF identique octet pour octet (333 Mo, hors date) et OBJ identique octet pour octet** à ceux du code d'origine.

## 5. Mesures d'ensemble

Ø250 mm, pas 0,2 mm, DEM synthétique de 3 001² pixels, **un seul cœur**, fichiers déjà en cache disque (la première exécution à froid est plus lente : 8,3 s pour le STL d'origine).

| Export | Code d'origine | Phase 2, fond `grid` | Phase 2, fond `merged` |
|---|---|---|---|
| STL | 5,2 s | **2,5 à 3,1 s** | 2,2 à 2,4 s |
| 3MF | 23,1 s | **15,4 s** | 10,5 à 11,4 s |
| OBJ | 12,6 s | **10,0 à 10,3 s** | 6,6 à 6,8 s |
| Pic de RAM (STL) | 1 102 Mo | 1 100 Mo | 645 à 687 Mo |

[Fait] La RAM de pointe n'est pas améliorée par les parois (le dessus et le fond dominent) : c'est le fond fusionné qui la réduit. [Opinion] Les gains de temps viennent de la suppression de boucles Python ; les mesures sont sur un relief synthétique : vos DEM changeraient les chiffres, pas la nature du gain.

## 6. Limites

1. **Réseau GSI réel non testé** ; couverture DEM1A non vérifiée (section 2.5).
2. **Fond fusionné toujours en prototype** ; Cura, OrcaSlicer, Bambu Studio et votre FAO non testés (rapport 7).
3. **Un seul cœur, Linux, GDAL 3.8.4** ; aucun test sous Windows.
4. La règle « la plus grossière décide » change le rééchantillonnage des listes à résolutions mixtes (§2.4).
5. Les empreintes de non-régression valent pour le **fond `grid`** ; si vous changez le défaut vers `merged`, elles devront être régénérées volontairement.

## 7. Reprendre ces changements

```fish
cd dem-to-stl
patch -p1 < ../dem-to-stl_phases1-2_cumul.patch   # sur l'archive d'origine ; ou décompresser le zip
python -m pytest -q
python circle_dem_to_stl.py --source gsi --lat 35.36 --lon 138.73 --radius 150 \
  --diameter-mm 60 --output fuji.stl
```

Si vous avez déjà appliqué le patch de la phase 1 : `dem-to-stl_phase2_incremental.patch` ne contient que la phase 2.

## 8. Ce qui reste de la feuille de route

| Phase | Élément | État |
|---|---|---|
| 1 | Rééchantillonnage, téléchargements, pincements, `check_stl` | fait (rapport 7) |
| 2 | Fond fusionné | **prototype fait** ; à tester dans vos trancheurs |
| 2 | Repli de résolution GSI | **fait** |
| 2 | Parois vectorisées | **fait** |
| 2 | Indexation native | **écartée** (mesure) |
| 2 | STL de référence minimal | **fait** (empreintes) |
| 3 | Diagonale adaptative (benchmark LiDAR) | non fait |
| 3 | TIN adaptatif (RTIN à bord forcé), écriture par bandes | non fait |
| 3 | Débit poli / `Retry-After` | non fait |
