# Tutorial — From topographic data to a printable STL

*[Version française](TUTORIAL.fr.md)*

This tutorial assumes you've already installed the dependencies
(`pip install -r requirements.txt`, see the README) and that you're in
the repository folder. Commands are given in **bash** syntax (the most
common) — adapt as needed for another shell (fish, zsh...).

## Part 1 — Quick start (no data download required)

Before tackling a real relief, let's check that everything works with
a synthetic DEM — a simple conical hill, generated in a few lines, no
need to wait on a download.

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
ds = driver.Create('hill_test.tif', size, size, 1, gdal.GDT_Float32)
ds.SetGeoTransform([cx0-size/2, 1.0, 0, cy0+size/2, 0, -1.0])
srs = osr.SpatialReference(); srs.ImportFromEPSG(32631)
ds.SetProjection(srs.ExportToWkt())
band = ds.GetRasterBand(1)
band.SetNoDataValue(-9999.0)
band.WriteArray(elev)
"
```

This creates `hill_test.tif`: a hill 120 m tall, 80 m in radius,
centered on (500000, 5000000) in UTM 31N (EPSG:32631) — a projected
CRS in meters chosen arbitrarily for this example.

Now run the tool on it:

```bash
python3 circle_dem_to_stl.py \
    --input hill_test.tif --output hill_disc.stl \
    --cx 500000 --cy 5000000 --radius 80 \
    --diameter-mm 60 --base-mm 3 --vexag 1.0
```

(`--lat`/`--lon` in WGS84 also work as an alternative to `--cx`/`--cy`
— handy if you only know the site's geographic coordinates, even when
the source is already projected.)

You should see the resampling resolution, the number of included
cells, then wall construction and STL writing scroll by
(`hill_disc.stl`, ~280,000 triangles for this example). Open it in
FreeCAD, Blender, or your usual STL viewer: you should see a solid
disc with a conical hill in the center — no holes, no cracks.

**What each parameter does here**:
- `--cx`/`--cy`: the circle's center, in the source file's CRS (here, UTM 31N).
- `--radius 80`: an 80 m radius in reality — exactly the hill's size.
- `--diameter-mm 60`: the printed piece will be 60 mm in diameter.
- `--base-mm 3`: a 3 mm flat base under the relief.
- `--vexag 1.0`: no vertical exaggeration — the relief is already well marked at this scale.

## Part 2 — Moving to a real site

This is where most of the actual work lies in practice — not in the
tool itself, but in preparing the data beforehand. Here is the full
method, as it emerged over the course of many sites processed (Mount
Fuji, Everest, the Matterhorn, Monument Valley, La Réunion, among
others).

### 2.1 — Finding a data source

Depending on the region of the world, the best source varies a lot:

| Region | Typical source | Resolution |
|---|---|---|
| Japan | GSI DEM1A/5A (XML/JPGIS) | 1-5 m |
| France | IGN LiDAR HD | 0.5-1 m |
| United States | USGS 3DEP | 1-10 m |
| Switzerland | swissALTI3D | 0.5-2 m |
| Norway | Kartverket Høydedata | 1 m |
| Rest of the world | Copernicus GLO-30 | 30 m |
| Planetary (Moon, Mars...) | USGS Astrogeology | variable |

General rule: the smaller the area of interest relative to the
resolution, the more the lack of detail shows. An isolated peak with
a 1-2 km radius deserves 1-10 m data; an entire massif with a radius
of several dozen km works perfectly well with 30 m.

### 2.2 — Checking the CRS (useful to know, no longer a blocker)

```bash
gdalinfo my_file.tif | grep -A 15 "Coordinate System is"
```

Since automatic reprojection was added, this is no longer a mandatory
step before running the script — if the file is in a geographic CRS
(degrees, often WGS84 or NAD83), `circle_dem_to_stl.py` detects this
and reprojects on the fly to the appropriate UTM zone during
resampling. `gdalinfo` remains useful for checking what a file
actually contains, spotting the declared nodata value, or
troubleshooting an issue — just no longer for reprojecting manually
beforehand.

### 2.3 — Providing the point of interest

Two options, your choice:

```bash
# directly in geographic coordinates (the simplest)
python3 circle_dem_to_stl.py --input my_file.tif --output output.stl \
    --lat <latitude> --lon <longitude> --radius <radius in m>

# or in the file's already-projected CRS, if you already have it in that form
python3 circle_dem_to_stl.py --input my_file.tif --output output.stl \
    --cx <X> --cy <Y> --radius <radius in m>
```

If you still need to convert a point by hand for some other reason
(checking a distance, for instance):

```bash
python3 -c "
from osgeo import osr
osr.UseExceptions()
src = osr.SpatialReference(); src.ImportFromEPSG(4326)  # WGS84, adapt as needed
dst = osr.SpatialReference(); dst.ImportFromEPSG(<file's code>)
src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
ct = osr.CoordinateTransformation(src, dst)
lon, lat = <longitude>, <latitude>
x, y, z = ct.TransformPoint(lon, lat)
print(x, y)
"
```

### 2.4 — Checking that the point actually falls within the data, with enough margin

```bash
gdalinfo my_file.tif | grep -E "Corner Coordinates" -A 5
```

Compare against the point of interest's coordinates: the smallest
margin in the 4 directions caps the usable radius. If the point is
too close to an edge, you'll need an additional neighboring tile — the
most reliable approach is to merge several tiles before continuing:

```bash
gdalbuildvrt mosaic.vrt tile1.tif tile2.tif tile3.tif
```

(`.vrt` works directly as `--input` — no need to materialize a merged
GeoTIFF for simple use.)

### 2.5 — If the exact point isn't known with certainty

For a summit whose precise coordinates vary between sources, the most
reliable approach is to search for the actual maximum within the data
itself rather than arbitrating between conflicting sources — by
reading a targeted window, never the whole tile (a real risk of
running out of memory on large tiles):

```bash
python3 -c "
from osgeo import gdal
import numpy as np
gdal.UseExceptions()

ds = gdal.Open('my_file.tif')
gt = ds.GetGeoTransform()
band = ds.GetRasterBand(1)
nodata = band.GetNoDataValue()

tx, ty = <rough estimate X, Y>
radius = 2000  # meters, err on the large side

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
print('Max elevation:', arr[idx])
print('Coordinates:', xs[idx[1]], ys[idx[0]])
"
```

### 2.6 — Running the generation

```bash
python3 circle_dem_to_stl.py \
    --input mosaic.vrt --output my_relief.stl \
    --cx <X> --cy <Y> --radius <radius in m>
```

Without `--diameter-mm`, `--base-mm`, or `--vexag`, the tool prompts
for them interactively. For a scripted sequence with no manual
intervention, pass all of them on the command line.

**Choosing the radius**: depends on what you want to capture — just a
summit, or the summit and its immediate surroundings down to the col
connecting it to a neighbor. There's no universal rule; looking at a
topographic map of the area to spot nearby cols/saddles remains the
best approach.

**Choosing `--vexag`**: 1.0 (no exaggeration) suits reliefs that are
already spectacular (high mountains). For a gentle relief (hill,
plateau) viewed over a large radius, an exaggeration of 2-5× generally
brings out the shape without distorting it.

## Part 3 — Island mode (`--sea-level`)

If the relief is surrounded by water (island, peninsula) and the sea
data is missing or unreliable — a common case with topographic
LiDAR, which doesn't measure water correctly — use `--sea-level`:

```bash
python3 circle_dem_to_stl.py \
    --input mosaic.vrt --output island.stl \
    --cx <X> --cy <Y> --radius <large radius, including the surrounding sea> \
    --sea-level 0
```

The circle is then always full: invalid or aberrant pixels are
replaced with a flat plane at sea level rather than excluded, and the
real terrain emerges from that plateau.

Two fine-tuning settings, to adjust only if the result shows
artifacts:
- `--min-elevation` (default -50 m): threshold below which a value
  outside the declared nodata is treated as an interpolation chasm.
- `--land-threshold` (default 1 m above `--sea-level`): threshold
  separating true land from water-surface noise (waves misclassified
  as "ground" by the LiDAR) — this noise often touches the real
  shoreline directly, so connectivity alone isn't enough to
  distinguish it.

## Part 4 — Checking before printing

The tool guarantees a watertight mesh and consistent normals by
construction (see `tests/`), but two things remain your
responsibility:

1. **Walls too thin for your nozzle**: can appear at the edge of the
   disc (quantization of the circular clip at the resampling
   resolution) or on a relief with very steep faces. A dedicated mesh
   checker (e.g. based on `trimesh`) can detect this automatically
   before printing.
2. **Real overhangs**: a relief with genuine vertical or overhanging
   walls (cliff, sharp ridge) remains a real overhang in the STL —
   that's not a mesh defect, just a matter of print orientation or
   supports to plan for.

## Quick troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `ModuleNotFoundError: No module named 'scipy'` | missing dependency | `pip install scipy --break-system-packages` |
| `_ARRAY_API not found` / numpy error when writing the raster | numpy 2.x / GDAL bindings conflict (can occur after a plain `requirements.txt` install, which doesn't cap the numpy version) | `pip install "numpy<2" --break-system-packages --force-reinstall` |
| `OSError: [Errno 122] Disk quota exceeded` or process killed (OOM) | `/tmp` too small, or a whole tile read into memory | redirect `TMPDIR` to a disk with room; read by window rather than the whole tile (see §2.5) |
| The produced circle seems offset from the expected relief | wrong center point, or NW/SW mix-up in the tile naming | double-check the coordinate conversion (§2.3) and the tile-naming convention of the source used |
| `RuntimeError: <file>: No such file or directory` | a merge/download step was skipped before running the script | check `ls` before running again |
