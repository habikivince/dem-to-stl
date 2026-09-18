# dem-to-stl

Découpe un modèle numérique de terrain (DEM/MNT, GeoTIFF projeté en
mètres) selon un cercle, et exporte un maillage **STL solide et
étanche** (relief + paroi verticale + fond plat), prêt pour
l'impression 3D ou l'usinage CNC — sans trou, sans paroi non-manifold,
échelle physique précise.

## Origine du projet

Ce dépôt est le fruit d'une longue conversation entre un amateur de
travail du bois et de CNC (aucune formation en géomatique ni en
développement logiciel) et Claude (Anthropic), pour produire des
reliefs imprimables à partir de vraies données topographiques du monde
entier — le Mont Fuji, l'Everest, Monument Valley, La Réunion, le
Cervin, entre autres.

**La quasi-totalité du code a été écrite par l'IA sous la direction et
la validation d'un humain** : chaque fonctionnalité est née d'un besoin
concret rencontré en cours de route (un maillage percé, des parois trop
fines pour l'impression, des artefacts de données en pleine mer...),
implémentée par Claude, puis testée sur des données réelles avant
d'être acceptée. Les corrections les moins évidentes (pincement
diagonal, séparation terre/mer par composante connexe plutôt que par
seuil d'altitude) sont documentées ci-dessous et dans les commentaires
du code, précisément parce qu'elles ne sont pas intuitives et qu'on ne
les découvre qu'en pratique.

Aucune garantie de qualité "production" au sens habituel du terme —
c'est un outil de passionné, testé sur des dizaines de cas réels par
un seul utilisateur, pas audité par une communauté. Les retours,
corrections et *pull requests* sont bienvenus.

## Fonctionnalités

- Découpe circulaire d'un DEM autour d'un centre (X, Y) et d'un rayon,
  avec mise à l'échelle physique précise vers un diamètre imprimé cible.
- Maillage garanti étanche (vérifié : chaque arête partagée par
  exactement 2 triangles) et normales cohérentes (volume signé positif).
- Élimination automatique des parois trop fines pour l'impression FDM
  (ouverture morphologique) et des pincements en diagonale
  (sommets non-manifold invisibles à un simple test d'étanchéité).
- **Mode île** (`--sea-level`) : pour un DEM où la donnée en mer est
  absente ou aberrante (cas fréquent des LiDAR topographiques), les
  pixels invalides sont remplacés par un plat au niveau de la mer
  plutôt qu'exclus — le relief réel émerge d'un disque toujours plein.
  Sépare la vraie terre du bruit de surface d'eau (vagues classées par
  erreur comme "sol") par composante connexe **et** seuil d'altitude —
  la connectivité seule ne suffit pas quand ce bruit touche le vrai
  littoral.
- Barres de progression (natives GDAL pour le ré-échantillonnage,
  maison pour les boucles Python pures) — utile sur les gros volumes
  (La Réunion : ~40 Go de tuiles sources).
- Diamètre et épaisseur de socle demandés interactivement si omis en
  ligne de commande.

## Installation

```bash
pip install -r requirements.txt
```

GDAL nécessite les bibliothèques système correspondantes (`libgdal-dev`
sur Debian/Ubuntu, `gdal` sur Arch) — si `pip install gdal` échoue,
installe GDAL via le gestionnaire de paquets de ton système plutôt que
pip seul.

## Usage

```bash
python circle_dem_to_stl.py \
    --input mont_fuji_dem1a.tif \
    --output fuji_disc.stl \
    --cx 293526 --cy 3915399 --radius 9300 \
    --diameter-mm 250 --vexag 1.0
```

Sans `--diameter-mm` ni `--base-mm`, le script les demande de façon
interactive au lancement.

### Mode île (relief entouré d'eau sans donnée fiable)

```bash
python circle_dem_to_stl.py \
    --input reunion_mosaic.vrt \
    --output reunion_disc.stl \
    --cx 347000 --cy 7664000 --radius 38000 \
    --diameter-mm 250 --vexag 3.0 \
    --sea-level 0
```

### Toutes les options

| Option | Description | Défaut |
|---|---|---|
| `--input` | GeoTIFF ou VRT source, CRS projeté en mètres | requis |
| `--output` | Fichier STL de sortie | requis |
| `--cx`, `--cy` | Centre du cercle (unités du CRS source) | requis |
| `--radius` | Rayon du cercle, en mètres | requis |
| `--diameter-mm` | Diamètre final imprimé, en mm | demandé si omis (250) |
| `--base-mm` | Épaisseur du socle plat, en mm | demandé si omis (3) |
| `--vexag` | Exagération verticale | 1.0 |
| `--print-spacing-mm` | Résolution cible du maillage, en mm | 0.2 |
| `--sea-level` | Active le mode île à cette altitude (m) | désactivé |
| `--min-elevation` | Seuil d'exclusion des gouffres d'interpolation, en m | -50 |
| `--land-threshold` | Seuil terre/bruit d'eau au-dessus de `--sea-level`, en m | 1.0 |

## `island_dem_to_stl.py`

Variante conservée pour le cas où l'on veut suivre le **contour réel du
littoral** (silhouette exacte, pas un disque) plutôt qu'un cercle avec
mer plate autour — moins utilisée en pratique au fil du projet
(le rendu en dents de scie du littoral brut s'est révélé moins
satisfaisant visuellement que le disque avec mer plate), mais gardée
disponible.

## Sources de données compatibles

Testé avec succès sur, entre autres : GSI DEM1A (Japon, XML/JPGIS —
nécessite un pré-traitement, non inclus ici), IGN LiDAR HD (France),
USGS 3DEP (États-Unis), swissALTI3D (Suisse), Kartverket LiDAR
(Norvège), Copernicus GLO-30 (mondial). Toujours vérifier le CRS du
fichier source (`gdalinfo`) avant utilisation — reprojeter avec
`gdalwarp` si nécessaire, le script suppose un CRS déjà projeté en
mètres.

## Tests

```bash
pytest tests/
```

Les tests reproduisent les scénarios synthétiques utilisés pendant le
développement (île conique, artefact isolé, bruit de surface d'eau
raccordé au littoral, motif en damier) et vérifient l'étanchéité, la
cohérence des normales, et le comportement attendu de chaque option.

## Licence

GNU GPLv3 — voir [LICENSE](LICENSE). Toute version modifiée ou dérivée
doit rester open source sous la même licence.
