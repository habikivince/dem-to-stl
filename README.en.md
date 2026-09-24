# dem-to-stl
 
*[Version française](README.fr.md)*
 
Clips a Digital Elevation Model (DEM/GeoTIFF, projected or geographic
CRS) to a circle and exports a **solid, watertight STL mesh** (terrain
surface + vertical wall + flat base), ready for 3D printing or CNC
machining — no holes, no non-manifold geometry, precise physical
scale.
 
## Project origin
 
This repository grew out of a long conversation between a woodworking
and CNC hobbyist (no background in geomatics or software development)
and Claude (Anthropic), aimed at producing printable relief models
from real-world topographic data from around the world — Mount Fuji,
Everest, Monument Valley, La Réunion, the Matterhorn, among others.
 
**Nearly all of the code was written by the AI under human direction
and validation**: every feature grew out of a concrete problem
encountered along the way (a mesh with a hole in it, walls too thin to
print, data artifacts in open water...), implemented by Claude, then
tested against real data before being accepted. The less obvious
fixes (diagonal pinch removal, land/sea separation by connected
component rather than a simple elevation threshold) are documented
below and in the code comments, precisely because they aren't
intuitive and only surface in practice.
 
No "production quality" guarantee in the usual sense — this is a
hobbyist's tool, tested on dozens of real cases by a single user, not
community-audited. Feedback, corrections, and pull requests are
welcome.
 
## Features
 
- Circular clipping of a DEM around a center point and radius, with
  precise physical scaling to a target printed diameter.
- Guaranteed watertight mesh (verified: every edge shared by exactly
  2 triangles) and consistent normals (positive signed volume).
- Automatic removal of walls too thin to print in FDM (morphological
  opening) and of diagonal pinches (non-manifold vertices invisible to
  a simple watertightness test).
- **Island mode** (`--sea-level`): for a DEM where sea data is missing
  or unreliable (a common case with topographic LiDAR), invalid pixels
  are replaced with a flat plane at sea level instead of being
  excluded — the real terrain emerges from a disc that's always full.
  Separates true land from water-surface noise (waves misclassified as
  "ground") using both connected components **and** an elevation
  threshold — connectivity alone isn't enough when this noise touches
  the real shoreline.
- Automatic CRS detection and reprojection: no manual `gdalwarp` step
  needed beforehand, even for a geographic-CRS source.
- Built-in Copernicus DEM download (`--download-copernicus`): fetches
  the needed GLO-30/GLO-90 tiles directly from the public S3 bucket,
  no account or API key required.
- Progress bars (native GDAL for resampling, custom for the pure
  Python loops) — useful on large volumes (La Réunion: ~40 GB of
  source tiles).
- Diameter, base thickness, and vertical exaggeration are prompted for
  interactively if omitted from the command line.
## Installation
 
```bash
pip install -r requirements.txt
```
 
GDAL requires the matching system libraries (`libgdal-dev` on
Debian/Ubuntu, `gdal` on Arch) — if `pip install gdal` fails, install
GDAL through your system's package manager rather than pip alone.
 
A complete, step-by-step tutorial — from an immediate hands-on demo
with no download required, to processing a real site — is available
in [TUTORIAL.md](TUTORIAL.md).
 
## Usage
 
```bash
python circle_dem_to_stl.py \
    --input mont_fuji_dem1a.tif \
    --output fuji_disc.stl \
    --cx 293526 --cy 3915399 --radius 9300 \
    --diameter-mm 250 --vexag 1.0
```
 
Or directly with geographic coordinates (WGS84), whether the source is
already projected or not — if it's in a geographic CRS, a UTM zone is
determined and applied automatically, with no manual reprojection
beforehand:
 
```bash
python circle_dem_to_stl.py \
    --input copernicus_glo30.tif \
    --output relief_disc.stl \
    --lat 35.3606 --lon 138.7274 --radius 9300 \
    --diameter-mm 250 --vexag 1.0
```
 
### Built-in Copernicus DEM download
 
No source file on hand? `--download-copernicus` automatically fetches
the needed GLO-30 (30 m, default) or GLO-90 (90 m) tiles from the
Copernicus DEM public S3 bucket (direct access, no account or API key):
 
```bash
python circle_dem_to_stl.py \
    --download-copernicus --cache-dir ./cache \
    --lat 35.3606 --lon 138.7274 --radius 9300 \
    --output relief_disc.stl --diameter-mm 250 --vexag 1.0
```
 
Tiles are cached in `--cache-dir` (one per square degree, reused
across runs) and mosaicked automatically. Requires `--lat`/`--lon`
(not `--cx`/`--cy`, since geographic coordinates are needed to
determine which tiles to fetch). Purely oceanic tiles (absent from the
bucket by convention) are skipped without failing the download.
 
### Island mode (relief surrounded by unreliable water data)
 
```bash
python circle_dem_to_stl.py \
    --input reunion_mosaic.vrt \
    --output reunion_disc.stl \
    --cx 347000 --cy 7664000 --radius 38000 \
    --diameter-mm 250 --vexag 3.0 \
    --sea-level 0
```
 
### All options
 
| Option | Description | Default |
|---|---|---|
| `--input` | Source GeoTIFF or VRT, geographic or projected CRS | omitted if `--download-copernicus` is used |
| `--output` | Output STL file | required |
| `--cx`, `--cy` | Circle center, in the source file's CRS | see `--lat`/`--lon` |
| `--lat`, `--lon` | Circle center in WGS84 — alternative to `--cx`/`--cy` | see `--cx`/`--cy` |
| `--radius` | Circle radius, in meters | required |
| `--diameter-mm` | Final printed diameter, in mm | prompted if omitted (250) |
| `--base-mm` | Flat base thickness, in mm | prompted if omitted (3) |
| `--vexag` | Vertical exaggeration | prompted if omitted (1.0) |
| `--print-spacing-mm` | Target mesh resolution, in mm | 0.2 |
| `--sea-level` | Enables island mode at this elevation (m) | disabled |
| `--min-elevation` | Threshold for excluding interpolation artifacts, in m | -50 |
| `--land-threshold` | Land/water-noise threshold above `--sea-level`, in m | 1.0 |
| `--download-copernicus` | Downloads the needed Copernicus DEM tiles instead of providing `--input` | disabled |
| `--copernicus-product` | Copernicus resolution to download: `30` (GLO-30) or `90` (GLO-90) | 30 |
| `--cache-dir` | Cache folder for downloaded Copernicus tiles | `./copernicus_cache` |
 
Provide either `--cx`/`--cy` or `--lat`/`--lon` — not both.
 
## Architecture
 
```
circle_dem_to_stl.py   # argparse, interactive prompts, orchestration
raster.py              # GDAL: opening, CRS detection/reprojection, resampling
mesh.py                # pure numpy/scipy: mask cleanup, mesh construction
stl_io.py              # binary STL writing
download.py            # Copernicus DEM download (public S3 bucket)
tests/
```
 
`mesh.py` never depends on GDAL or the filesystem — testable with
plain numpy arrays, which lets watertightness and normal consistency
be verified without ever opening a file.
 
## Compatible data sources
 
Successfully tested on, among others: GSI DEM1A (Japan, XML/JPGIS —
requires pre-processing, not included here), IGN LiDAR HD (France),
USGS 3DEP (United States), swissALTI3D (Switzerland), Kartverket LiDAR
(Norway), Copernicus GLO-30 (worldwide). Geographic or projected CRS
is detected automatically (`gdalinfo` remains useful for checking what
a file actually contains, but manual reprojection beforehand is no
longer necessary).
 
## Tests
 
```bash
pytest tests/
```
 
The tests reproduce the synthetic scenarios used during development
(conical island, isolated artifact, water-surface noise touching the
shoreline, checkerboard pattern) and verify watertightness, normal
consistency, and the expected behavior of each option.
 
## Notable changes
 
- **Modular restructuring** (`raster.py`/`mesh.py`/`stl_io.py`):
  eliminates the duplication that used to exist between the old
  monolithic `circle_dem_to_stl.py` and `island_dem_to_stl.py`.
- **Automatic CRS reprojection**: no more manual `gdalwarp` before
  each new site when the source is in geographic coordinates.
- **`--lat`/`--lon`** as an alternative to `--cx`/`--cy`.
- **Built-in Copernicus DEM download** (`--download-copernicus`).
- **`island_dem_to_stl.py` removed**: its use case (relief surrounded
  by water) is covered by `--sea-level`, whose approach (circle with a
  flat sea) proved preferable in practice to following the exact
  shoreline.
## License
 
GNU GPLv3 — see [LICENSE](LICENSE). Any modified or derivative version
must remain open source under the same license.
