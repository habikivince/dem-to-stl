# Tutoriel — De la donnée topographique au STL imprimable

Ce tutoriel suppose que tu as déjà installé les dépendances
(`pip install -r requirements.txt`, voir le README) et que tu es dans
le dossier du dépôt. Les commandes sont données en syntaxe **bash**
(la plus répandue) — adapte si besoin pour un autre shell (fish, zsh...).

## Partie 1 — Prise en main rapide (sans télécharger de données)

Avant de partir sur un vrai relief, vérifions que tout fonctionne avec
un DEM synthétique — une simple colline conique, générée en quelques
lignes, pas besoin d'attendre un téléchargement.

```bash
python3 -c "
from osgeo import gdal, osr
import numpy as np
gdal.UseExceptions()

size = 200
cx0, cy0 = 500000.0, 5000000.0
xs = np.arange(size) + cx0 - size/2
ys = np.arange(size) + cy0 - size/2
X, Y = np.meshgrid(xs, ys[::-1])
dist = np.sqrt((X-cx0)**2 + (Y-cy0)**2)
elev = np.where(dist <= 80, 120 - dist*1.2, -9999.0).astype(np.float32)

driver = gdal.GetDriverByName('GTiff')
ds = driver.Create('colline_test.tif', size, size, 1, gdal.GDT_Float32)
ds.SetGeoTransform([cx0-size/2, 1.0, 0, cy0+size/2, 0, -1.0])
srs = osr.SpatialReference(); srs.ImportFromEPSG(32631)
ds.SetProjection(srs.ExportToWkt())
band = ds.GetRasterBand(1)
band.SetNoDataValue(-9999.0)
band.WriteArray(elev)
"
```

Ça crée `colline_test.tif` : une colline de 120 m de haut, 80 m de
rayon, centrée sur (500000, 5000000) en UTM 31N (EPSG:32631) — un CRS
projeté en mètres choisi arbitrairement pour l'exemple.

Lance maintenant l'outil dessus :

```bash
python3 circle_dem_to_stl.py \
    --input colline_test.tif --output colline_disc.stl \
    --cx 500000 --cy 5000000 --radius 80 \
    --diameter-mm 60 --base-mm 3 --vexag 1.0
```

(`--lat`/`--lon` en WGS84 fonctionnent aussi comme alternative à
`--cx`/`--cy` — pratique si tu ne connais que les coordonnées
géographiques du site, y compris quand la source est déjà projetée.)

Tu devrais voir défiler la résolution de ré-échantillonnage, le
nombre de cellules incluses, puis la construction des parois et
l'écriture du STL (`colline_disc.stl`, ~280 000 triangles pour cet
exemple). Ouvre-le dans FreeCAD, Blender, ou ton visualiseur STL
habituel : tu dois voir un disque solide avec une colline conique au
centre — pas de trou, pas de fissure.

**Ce que fait chaque paramètre ici** :
- `--cx`/`--cy` : le centre du cercle, dans le CRS du fichier source (ici, UTM 31N).
- `--radius 80` : rayon de 80 m en réalité — pile la taille de la colline.
- `--diameter-mm 60` : la pièce imprimée fera 60 mm de diamètre.
- `--base-mm 3` : socle plat de 3 mm sous le relief.
- `--vexag 1.0` : pas d'exagération verticale — le relief est déjà bien marqué à cette échelle.

## Partie 2 — Passer à un vrai relief

C'est là que la majorité du travail se situe en pratique — pas dans
l'outil lui-même, mais dans la préparation des données en amont.
Voici la méthode complète, telle qu'elle a émergé au fil de
nombreux sites traités (Mont Fuji, Everest, Cervin, Monument Valley,
La Réunion, entre autres).

### 2.1 — Trouver une source de données

Selon la zone du monde, la meilleure source varie beaucoup :

| Zone | Source typique | Résolution |
|---|---|---|
| Japon | GSI DEM1A/5A (XML/JPGIS) | 1-5 m |
| France | IGN LiDAR HD | 0,5-1 m |
| États-Unis | USGS 3DEP | 1-10 m |
| Suisse | swissALTI3D | 0,5-2 m |
| Norvège | Kartverket Høydedata | 1 m |
| Reste du monde | Copernicus GLO-30 | 30 m |
| Planétaire (Lune, Mars...) | USGS Astrogeology | variable |

Règle générale : plus la zone d'intérêt est petite par rapport à la
résolution, plus le manque de détail se voit. Un pic isolé avec un
rayon de 1-2 km mérite du 1-10 m ; un massif entier avec un rayon de
plusieurs dizaines de km s'accommode très bien de 30 m.

### 2.2 — Vérifier le CRS (utile à connaître, plus bloquant)

```bash
gdalinfo mon_fichier.tif | grep -A 15 "Coordinate System is"
```

Depuis la reprojection automatique, ce n'est plus une étape obligatoire
avant de lancer le script — si le fichier est en CRS géographique
(degrés, souvent WGS84 ou NAD83), `circle_dem_to_stl.py` détecte le cas
et reprojette lui-même à la volée vers la zone UTM appropriée pendant
le ré-échantillonnage. `gdalinfo` reste utile pour vérifier ce que
contient un fichier, repérer le nodata déclaré, ou diagnostiquer un
problème — mais plus pour reprojeter manuellement au préalable.

### 2.3 — Donner le point d'intérêt

Deux façons, au choix :

```bash
# directement en coordonnées géographiques (le plus simple)
python3 circle_dem_to_stl.py --input mon_fichier.tif --output sortie.stl \
    --lat <latitude> --lon <longitude> --radius <rayon en m>

# ou dans le CRS déjà projeté du fichier, si tu l'as déjà sous cette forme
python3 circle_dem_to_stl.py --input mon_fichier.tif --output sortie.stl \
    --cx <X> --cy <Y> --radius <rayon en m>
```

Si tu dois malgré tout convertir un point à la main pour une autre
raison (vérifier une distance, par exemple) :

```bash
python3 -c "
from osgeo import osr
osr.UseExceptions()
src = osr.SpatialReference(); src.ImportFromEPSG(4326)  # WGS84, à adapter
dst = osr.SpatialReference(); dst.ImportFromEPSG(<code du fichier>)
src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
ct = osr.CoordinateTransformation(src, dst)
lon, lat = <longitude>, <latitude>
x, y, z = ct.TransformPoint(lon, lat)
print(x, y)
"
```

### 2.4 — Vérifier que le point tombe bien dans les données, avec assez de marge

```bash
gdalinfo mon_fichier.tif | grep -E "Corner Coordinates" -A 5
```

Compare aux coordonnées du point d'intérêt : la marge la plus faible
dans les 4 directions plafonne le rayon utilisable. Si le point est
trop près d'un bord, il faut une tuile voisine supplémentaire — le
plus fiable est de fusionner plusieurs tuiles avant de continuer :

```bash
gdalbuildvrt mosaique.vrt tuile1.tif tuile2.tif tuile3.tif
```

(`.vrt` fonctionne directement en `--input`, pas besoin de matérialiser
un GeoTIFF fusionné pour un usage simple.)

### 2.5 — Si le point exact n'est pas connu avec certitude

Pour un sommet dont les coordonnées précises varient selon les
sources, le plus fiable est de chercher le maximum réel dans les
données elles-mêmes plutôt que de trancher entre des sources
contradictoires — en lisant une fenêtre ciblée, jamais la tuile
entière (risque réel de saturation mémoire sur de grandes tuiles) :

```bash
python3 -c "
from osgeo import gdal
import numpy as np
gdal.UseExceptions()

ds = gdal.Open('mon_fichier.tif')
gt = ds.GetGeoTransform()
band = ds.GetRasterBand(1)
nodata = band.GetNoDataValue()

tx, ty = <estimation approximative X, Y>
radius = 2000  # mètres, à ajuster large

px = max(0, int((tx - radius - gt[0]) / gt[1]))
py = max(0, int((ty + radius - gt[3]) / gt[5]))
size_px = int(2*radius / abs(gt[1])) + 2
size_px_x = min(size_px, ds.RasterXSize - px)
size_px_y = min(size_px, ds.RasterYSize - py)

arr = band.ReadAsArray(px, py, size_px_x, size_px_y)
xs = gt[0] + (px + np.arange(size_px_x) + 0.5) * gt[1]
ys = gt[3] + (py + np.arange(size_px_y) + 0.5) * gt[5]
X, Y = np.meshgrid(xs, ys)
dist2 = (X-tx)**2 + (Y-ty)**2
valid = (arr != nodata) if nodata is not None else np.ones_like(arr, dtype=bool)
mask = valid & (dist2 <= radius**2)
masked = np.where(mask, arr, -np.inf)
idx = np.unravel_index(np.argmax(masked), masked.shape)
print('Altitude max :', arr[idx])
print('Coordonnées :', xs[idx[1]], ys[idx[0]])
"
```

### 2.6 — Lancer la génération

```bash
python3 circle_dem_to_stl.py \
    --input mosaique.vrt --output mon_relief.stl \
    --cx <X> --cy <Y> --radius <rayon en m>
```

Sans `--diameter-mm` ni `--base-mm` ni `--vexag`, l'outil les demande
interactivement. Pour un enchaînement scripté sans intervention
manuelle, passe-les tous en ligne de commande.

**Choisir le rayon** : dépend de ce que tu veux capturer — juste un
sommet, ou le sommet et ses abords immédiats jusqu'à la selle qui le
relie à un voisin. Il n'y a pas de règle universelle ; regarder une
carte topographique de la zone pour repérer les cols/selles voisins
reste la meilleure approche.

**Choisir `--vexag`** : 1.0 (aucune exagération) convient aux reliefs
déjà spectaculaires (haute montagne). Pour un relief doux (colline,
plateau) vu sur un grand rayon, une exagération de 2-5× fait
généralement ressortir la forme sans la dénaturer.

## Partie 3 — Mode île (`--sea-level`)

Si le relief est entouré d'eau (île, presqu'île) et que les données
en mer sont absentes ou aberrantes — cas fréquent des LiDAR
topographiques, qui ne mesurent pas correctement l'eau — utilise
`--sea-level` :

```bash
python3 circle_dem_to_stl.py \
    --input mosaique.vrt --output ile.stl \
    --cx <X> --cy <Y> --radius <rayon large, incluant la mer autour> \
    --sea-level 0
```

Le cercle est alors toujours plein : les pixels invalides ou
aberrants sont remplacés par un plat au niveau de la mer plutôt
qu'exclus, et le relief réel émerge de ce plateau.

Deux réglages fins, à ajuster seulement si le résultat montre des
artefacts :
- `--min-elevation` (défaut -50 m) : seuil sous lequel une valeur
  hors nodata déclaré est traitée comme un gouffre d'interpolation.
- `--land-threshold` (défaut 1 m au-dessus de `--sea-level`) : seuil
  séparant la vraie terre du bruit de surface d'eau (vagues classées
  par erreur comme "sol" par le LiDAR) — ce bruit touche souvent
  directement le littoral, la connectivité seule ne suffit pas à
  l'en distinguer.

## Partie 4 — Vérifier avant d'imprimer

L'outil garantit un maillage étanche et des normales cohérentes par
construction (voir `tests/`), mais deux points restent à ta charge :

1. **Parois trop fines pour ta buse** : peuvent apparaître en bord de
   disque (quantification de la découpe circulaire à la résolution de
   ré-échantillonnage) ou sur un relief à pans très raides. Un
   vérificateur de maillage dédié (par ex. basé sur `trimesh`) peut
   détecter ça automatiquement avant impression.
2. **Surplombs réels** : un relief avec de vraies parois verticales ou
   en dévers (falaise, crête acérée) reste un vrai surplomb dans le
   STL — ce n'est pas un défaut du maillage, juste une question
   d'orientation d'impression ou de supports à prévoir.

## Dépannage rapide

| Symptôme | Cause probable | Solution |
|---|---|---|
| `ModuleNotFoundError: No module named 'scipy'` | dépendance manquante | `pip install scipy --break-system-packages` |
| `_ARRAY_API not found` / erreur numpy au moment d'écrire le raster | conflit numpy 2.x / bindings GDAL (peut survenir après une simple installation de `requirements.txt`, qui ne plafonne pas la version de numpy) | `pip install "numpy<2" --break-system-packages --force-reinstall` |
| `OSError: [Errno 122] Disk quota exceeded` ou processus tué (OOM) | `/tmp` trop petit, ou lecture d'une tuile entière en mémoire | rediriger `TMPDIR` vers un disque avec de la place ; lire par fenêtre plutôt que la tuile entière (voir §2.5) |
| Le cercle produit semble décalé par rapport au relief attendu | mauvais point de centre, ou confusion NW/SW dans le nommage des tuiles | revérifier la conversion de coordonnées (§2.3) et le sens du nommage des tuiles de la source utilisée |
| `RuntimeError: <fichier>: No such file or directory` | étape de fusion/téléchargement oubliée avant de lancer le script | vérifier `ls` avant de relancer |
