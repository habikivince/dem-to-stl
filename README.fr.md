# dem-to-stl

*[English version](README.md)*

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

## Galerie

<table>
<tr>
<td><img src="screenshots/mont_fuji_5.png" width="400" alt="Disque du Mont Fuji, détail du cratère"></td>
<td><img src="screenshots/matterhorn_2.png" width="400" alt="Disque du Cervin"></td>
</tr>
<tr>
<td align="center">Mont Fuji — cratère sommital</td>
<td align="center">Cervin</td>
</tr>
</table>

D'autres angles pour chaque site sont disponibles dans [`screenshots/`](screenshots/).

## Installation

```bash
pip install -r requirements.txt
```

GDAL nécessite les bibliothèques système correspondantes (`libgdal-dev`
sur Debian/Ubuntu, `gdal` sur Arch) — si `pip install gdal` échoue,
installe GDAL via le gestionnaire de paquets de ton système plutôt que
pip seul.

Un tutoriel complet, pas à pas — de la prise en main rapide sans
téléchargement jusqu'au traitement d'un vrai relief — est disponible
dans [TUTORIAL.md](TUTORIAL.md).

## Usage

```bash
python circle_dem_to_stl.py \
    --input mont_fuji_dem1a.tif \
    --output fuji_disc.stl \
    --cx 293526 --cy 3915399 --radius 9300 \
    --diameter-mm 250 --vexag 1.0
```

Ou directement en coordonnées géographiques (WGS84), que la source soit
déjà projetée ou non — si elle est en CRS géographique, une zone UTM
est déterminée et appliquée automatiquement, sans reprojection manuelle
préalable :

```bash
python circle_dem_to_stl.py \
    --input copernicus_glo30.tif \
    --output relief_disc.stl \
    --lat 35.3606 --lon 138.7274 --radius 9300 \
    --diameter-mm 250 --vexag 1.0
```

### Téléchargement Copernicus DEM intégré

Pas de fichier source sous la main ? `--download-copernicus` récupère
automatiquement les tuiles GLO-30 (30 m, par défaut) ou GLO-90 (90 m)
nécessaires depuis le compartiment S3 public de Copernicus DEM (accès
direct, sans compte ni clé API) :

```bash
python circle_dem_to_stl.py \
    --download-copernicus --cache-dir ./cache \
    --lat 35.3606 --lon 138.7274 --radius 9300 \
    --output relief_disc.stl --diameter-mm 250 --vexag 1.0
```

Les tuiles sont mises en cache dans `--cache-dir` (une par degré carré,
réutilisées d'un lancement à l'autre) et assemblées automatiquement en
mosaïque. Requiert `--lat`/`--lon` (pas `--cx`/`--cy`, puisqu'il faut
des coordonnées géographiques pour déterminer les tuiles). Les tuiles
purement océaniques (absentes du compartiment par convention) sont
ignorées sans faire échouer le téléchargement.

Sans `--diameter-mm`, `--base-mm` ni `--vexag`, le script les demande de
façon interactive au lancement.

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
| `--input` | GeoTIFF ou VRT source, CRS géographique ou projeté | requis |
| `--output` | Fichier STL de sortie | requis |
| `--cx`, `--cy` | Centre du cercle, dans le CRS du fichier source | voir `--lat`/`--lon` |
| `--lat`, `--lon` | Centre du cercle en WGS84 — alternative à `--cx`/`--cy` | voir `--cx`/`--cy` |
| `--radius` | Rayon du cercle, en mètres | requis |
| `--diameter-mm` | Diamètre final imprimé, en mm | demandé si omis (250) |
| `--base-mm` | Épaisseur du socle plat, en mm | demandé si omis (3) |
| `--vexag` | Exagération verticale | demandée si omise (1.0) |
| `--print-spacing-mm` | Résolution cible du maillage, en mm | 0.2 |
| `--sea-level` | Active le mode île à cette altitude (m) | désactivé |
| `--min-elevation` | Seuil d'exclusion des gouffres d'interpolation, en m | -50 |
| `--land-threshold` | Seuil terre/bruit d'eau au-dessus de `--sea-level`, en m | 1.0 |
| `--download-copernicus` | Télécharge les tuiles Copernicus DEM nécessaires au lieu de fournir `--input` | désactivé |
| `--copernicus-product` | Résolution Copernicus à télécharger : `30` (GLO-30) ou `90` (GLO-90) | 30 |
| `--cache-dir` | Dossier de cache des tuiles Copernicus téléchargées | `./copernicus_cache` |

Fournis soit `--cx`/`--cy`, soit `--lat`/`--lon` — pas les deux.

## Architecture

```
circle_dem_to_stl.py   # argparse, prompts interactifs, orchestration
raster.py              # GDAL : ouverture, détection/reprojection CRS, ré-échantillonnage
mesh.py                # numpy/scipy pur : nettoyage du masque, construction du maillage
stl_io.py              # écriture STL binaire
download.py            # téléchargement Copernicus DEM (compartiment S3 public)
tests/
```

`mesh.py` ne dépend jamais de GDAL ni du disque — testable avec de
simples tableaux numpy, ce qui permet de vérifier étanchéité et
cohérence des normales sans jamais ouvrir un fichier.

## Sources de données compatibles

Testé avec succès sur, entre autres : GSI DEM1A (Japon, XML/JPGIS —
nécessite un pré-traitement, non inclus ici), IGN LiDAR HD (France),
USGS 3DEP (États-Unis), swissALTI3D (Suisse), Kartverket LiDAR
(Norvège), Copernicus GLO-30 (mondial). Le CRS géographique ou projeté
est détecté automatiquement (`gdalinfo` reste utile pour vérifier ce
que contient un fichier avant de s'en servir, mais la reprojection
manuelle préalable n'est plus nécessaire).

## Tests

```bash
pytest tests/
```

Les tests reproduisent les scénarios synthétiques utilisés pendant le
développement (île conique, artefact isolé, bruit de surface d'eau
raccordé au littoral, motif en damier) et vérifient l'étanchéité, la
cohérence des normales, et le comportement attendu de chaque option.

## Historique des changements notables

- **Restructuration modulaire** (`raster.py`/`mesh.py`/`stl_io.py`) : élimine la duplication qui existait entre l'ancien `circle_dem_to_stl.py` monolithique et `island_dem_to_stl.py`.
- **Reprojection CRS automatique** : plus besoin de `gdalwarp` manuel avant chaque nouveau site quand la source est en coordonnées géographiques.
- **`--lat`/`--lon`** en alternative à `--cx`/`--cy`.
- **`island_dem_to_stl.py` retiré** : son cas d'usage (relief entouré d'eau) est couvert par `--sea-level`, dont l'approche (cercle à mer plate) s'est révélée préférable en pratique à la silhouette exacte du littoral.

## Licence

GNU GPLv3 — voir [LICENSE](LICENSE). Toute version modifiée ou dérivée
doit rester open source sous la même licence.
