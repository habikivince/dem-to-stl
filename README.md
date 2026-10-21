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
- Regular polygons (3 to 360 sides) inscribed in or circumscribed about a reference circle, with sharp
  corners or exact tangent-arc rounding (`--sides`, see below).
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
## Gallery

<table>
<tr>
<td><img src="screenshots/mont_fuji_5.png" width="400" alt="Mount Fuji disc, crater detail"></td>
<td><img src="screenshots/matterhorn_2.png" width="400" alt="Matterhorn disc"></td>
</tr>
<tr>
<td align="center">Mount Fuji — summit crater</td>
<td align="center">Matterhorn</td>
</tr>
</table>

More angles for each site are available in [`screenshots/`](screenshots/).

## Installation
 
```bash
pip install -r requirements.txt
```
 
GDAL requires the matching system libraries (`libgdal-dev` on
Debian/Ubuntu, `gdal` on Arch) — if `pip install gdal` fails, install
GDAL through your system's package manager rather than pip alone.
 
**Optional dependencies** (`requirements-optional.txt`): `laspy` for LAS point
clouds and `lazrs` for compressed `.laz` files; `netCDF4` and `h5py` are used
by the tests only. Without them everything else works, and the program says
what is missing when you give it a `.las`/`.laz` file. On Arch (where `pip`
refuses to touch the system environment), use a virtual environment that
still sees the system GDAL:

```fish
python -m venv --system-site-packages .venv
source .venv/bin/activate.fish
pip install -r requirements.txt -r requirements-optional.txt
```

(Arch package names for laspy/lazrs were not checked here; `pip` in a venv is
the route that was tested.)

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
 
### Non-linear exaggeration (`--vexag-curve`)

```bash
python circle_dem_to_stl.py --source auto --lat 27.9881 --lon 86.925 --radius 20000 \
    --diameter-mm 200 --vexag-curve "0:0,4000:4000,8848:6000"
```

Instead of a single factor, a table *real altitude → equivalent altitude* (metres, decimal point),
applied before the scale and `--vexag` (which still multiplies). In the example, relief is true to scale
up to 4000 m, then the next 4848 m count for only 2000 m. Interpolation is monotone PCHIP: no dip or spike
is invented between two points, the highest point stays the highest, and beyond the given points the slope of
the first or last segment continues. Two equal equivalent altitudes clip the relief (a flat). The log warns
that this is **not a faithful replica** (relative heights are altered) and prints the average exaggeration per
segment; the file name gets a `_courbe` suffix. With a curve, `check_stl.py --relief-m` must be given the
*equivalent* relief (`relief.relief_equivalent_m` in the `--report` JSON).

### Base altitude (`--base-elevation`)

```bash
python circle_dem_to_stl.py --input massif.tif --source-crs EPSG:2154 --cx 970000 --cy 6500000 --radius 8000 \
    --diameter-mm 200 --vexag 1 --base-mm 8 --base-elevation 1000
```

By default the top of the base sits at the lowest point of the piece, so two pieces cut from the same range do not
share the same zero. `--base-elevation M` fixes the altitude (metres) that sits at the top of the base: everything
lower is flattened to it (the log gives the share of nodes concerned), everything higher is kept, **including
summits that are not connected to the main one**. A value below the lowest point thickens the base by the
difference. A value at or above the highest point is refused. The file name gets a `_ref<M>m` suffix and the
`--report` JSON a `model.base_elevation_m` field. It cannot be combined with `--sea-level`: island mode already
sets the base at sea level, but keeps only the largest land component and flattens any other summit (measured: a
1300 m summit separate from the main massif is erased by `--sea-level 1000`, kept by `--base-elevation 1000`).

### All options
 
| Option | Description | Default |
|---|---|---|
| `--input` | One or more source files: GeoTIFF/VRT, ESRI ASCII Grid (`.asc`), SRTM (`.hgt`), GeoPackage raster (`.gpkg`), NetCDF (`.nc`), HDF5 (`.h5`/`.hdf5`, NetCDF-4/CF only), LAS/LAZ point cloud — see [Input formats](#input-formats) | omitted if `--download-copernicus` is used |
| `--layer` | GeoPackage raster table, NetCDF variable or HDF5 dataset to read (required when the file holds several candidates) | none |
| `--band` | Band to read: multiband file, or NetCDF/HDF5 slice (time, level...) — never chosen for you | none |
| `--source-crs` | CRS of a file that declares none (e.g. `EPSG:2056`); refused if the file already declares one | none |
| `--vertical-unit` | `m`, `ft` or `us-ft`: unit of the altitudes, replaces the declared (or missing) one | declared unit, else metres + warning |
| `--lidar-select` | LAS/LAZ: `ground` (class 2 mean per cell), `lowest` or `highest` per cell | `ground` |
| `--lidar-cell` | LAS/LAZ: cell size of the gridded point cloud, in m | working pixel |
| `--lidar-interp` | LAS/LAZ: hole filling `linear`, `nearest` or `none` | `linear` |
| `--lidar-max-gap` | LAS/LAZ: maximum fill distance, in m | 50 cells |
| `--lidar-spike-m` | LAS/LAZ: tolerance for isolated spike/pit cells, in m (`0` = off) | 10 |
| `--output` | Output STL file | required |
| `--cx`, `--cy` | Circle center, in the source file's CRS | see `--lat`/`--lon` |
| `--lat`, `--lon` | Circle center in WGS84 — alternative to `--cx`/`--cy` | see `--cx`/`--cy` |
| `--radius` | Circle radius (reference circle with `--sides`), in meters | required |
| `--diameter-mm` | Final printed diameter, in mm | prompted if omitted (250) |
| `--base-mm` | Flat base thickness, in mm | prompted if omitted (3) |
| `--vexag` | Vertical exaggeration | prompted if omitted (1.0) |
| `--vexag-curve` | Non-linear exaggeration: `alt:equivalent,...` points (m), PCHIP, applied before the scale and `--vexag`; not a faithful replica | disabled |
| `--base-elevation` | Altitude (m) at the top of the base: lower is flattened, higher is kept (default: lowest point); excludes `--sea-level` | lowest point |
| `--print-spacing-mm` | Target mesh resolution, in mm | 0.2 |
| `--resample` | DEM resampling: `auto`, `average`, `bilinear`, `cubic`, `cubicspline` | `auto` |
| `--top` | Top-surface triangulation: `grid`, `adaptive` (prototype) | `grid` |
| `--top-tol-mm` | With `--top adaptive`: vertical tolerance, in print mm | 0.05 |
| `--diagonal` | Cell diagonal of the top surface: `fixed`, `adaptive` (prototype) | `fixed` |
| `--bottom` | Flat-bottom triangulation: `grid`, `merged` (prototype), `merged-rows` (comparison) | `grid` |
| `--bottom-max-cells` | With `merged-rows`: maximum rectangle side, in cells (0 = unlimited) | 0 |
| `--sides` | Regular N-gon (3-360) instead of the disc (see below) | disabled |
| `--circle-mode` | With `--sides`: `inscribed` (vertices on the reference circle) or `circumscribed` (sides tangent) | required with `--sides` |
| `--corner-radius` | With `--sides`: corner rounding radius, in printed mm | 0 (sharp) |
| `--ref-outline` | With `--sides`, `--circle-mode inscribed` and `--corner-radius`: what the reference circle describes, `sharp` (the polygon before rounding) or `rounded` (tangent to the rounded corners) | `sharp` |
| `--rotation` | With `--sides`: rotation in degrees, counter-clockwise | 0 |
| `--corner-tol-mm` | With `--sides` and `--corner-radius`: maximum arc sagitta, in printed mm | 0.005 |
| `--corner-segments` | With `--sides` and `--corner-radius`: segments per arc (replaces `--corner-tol-mm`) | computed |
| `--sea-level` | Enables island mode at this elevation (m) | disabled |
| `--min-elevation` | Threshold for excluding interpolation artifacts, in m | -50 |
| `--land-threshold` | Land/water-noise threshold above `--sea-level`, in m | 1.0 |
| `--extrapolate-edges` | Extends the relief by up to 50 working pixels beyond the data (former gap filling); invented cells are logged | disabled |
| `--download-copernicus` | Downloads the needed Copernicus DEM tiles instead of providing `--input` | disabled |
| `--copernicus-product` | Copernicus resolution to download: `30` (GLO-30) or `90` (GLO-90) | 30 |
| `--cache-dir` | Cache folder for downloaded Copernicus tiles | `./copernicus_cache` |
| `--preview` | Shaded PNG preview of the cut area, written before meshing; optional file name | disabled |
| `--report` | JSON report next to the model (source, attribution, scale, relief, triangles, command line); optional file name | disabled |
 
Provide either `--cx`/`--cy` or `--lat`/`--lon` — not both.
 
## Architecture
 
```
circle_dem_to_stl.py   # argparse, interactive prompts, orchestration
formats.py             # input formats: detection, per-format reading rules, normalisation (CRS, units, NoData)
lidar.py               # LAS/LAZ point clouds: chunked reading, class filtering, rasterisation, hole filling
raster.py              # GDAL: opening, CRS detection/reprojection, resampling
mesh.py                # pure numpy/scipy: mask cleanup, mesh construction
regular_polygon.py     # pure numpy: regular polygons (--sides), exact membership, projection, arcs
vexag_curve.py         # non-linear exaggeration (--vexag-curve): altitude -> equivalent altitude, monotone PCHIP
stl_io.py              # binary STL writing
preview.py             # shaded PNG preview (--preview): hillshade and PNG writer, no GDAL
report.py              # printed relief height, indicative exaggeration, JSON report (--report); no GDAL
download.py            # Copernicus DEM download (public S3 bucket)
tests/
```
 
Every source goes through the same chain: **detection → reader → normalisation
(altitudes in metres, CRS known, NoData marked) → geospatial processing
(reprojection, clipping, resampling) → common meshing engine → export**. The
meshing engine (`mesh.py`) knows nothing about file formats. A point cloud is
not forced through a raster reader: `lidar.py` rasterises it directly onto the
working grid.

`mesh.py` never depends on GDAL or the filesystem — testable with
plain numpy arrays, which lets watertightness and normal consistency
be verified without ever opening a file.
 
## Input formats

`--input` takes one or more files; the format is recognised from the extension
and, for point clouds, from the file signature. Whatever the format, the program
never guesses silently: when a choice is ambiguous (several bands, tables or
variables, a missing CRS, an unclassified point cloud) it **stops and says what to
add**. The three levels below are what was actually tested.

| Format | Extensions | Level | What you must know |
|---|---|---|---|
| GeoTIFF, VRT | `.tif` `.tiff` `.vrt` | **fully supported** | Any compression GDAL has (DEFLATE, LZW, PackBits, tiles, BigTIFF; ZSTD and LERC were tested because this GDAL has them). Int16/UInt16/Int32/Float32/Float64. One band, or `--band N`. |
| ESRI ASCII Grid | `.asc` | **fully supported** | Header validated (malformed or truncated files are refused, never read blindly). CRS from the `.prj` next to it, otherwise `--source-crs`. |
| SRTM | `.hgt` | **fully supported** | Name must look like `N46E007.hgt` (it carries the position). SRTM1 (3601²) and SRTM3 (1201²). Voids are NoData. Give neighbouring tiles together to assemble them. |
| GeoPackage (raster) | `.gpkg` | **needs a selection** | Only **raster elevation tables** (one band Int16/UInt16/Float32). With 2+ tables, `--layer TABLE` is required. Vector-only packages and image tiles (RGB, palette) are refused: *not every GeoPackage is a DEM*. |
| NetCDF | `.nc` | **needs a selection** | Regular spatial axes (lat/lon or x/y) required. Several variables → `--layer NAME` unless exactly one declares itself an elevation. Extra dimensions (time, level) → `--band N`. |
| HDF5 | `.h5` `.hdf5` | **partial** | Only HDF5 that follows the NetCDF-4/CF model with spatial axes (read like NetCDF). A generic HDF5 (bare datasets) is **reported as unsupported**, with the datasets it contains. |
| LAS, LAZ | `.las` `.laz` | **supported, own pipeline** | Needs `laspy` (+`lazrs` for `.laz`). Ground class 2 by default; an unclassified cloud needs an explicit `--lidar-select`. See [Point clouds](#point-clouds-lasllaz). |
| XYZ text, **regular grid** | `.xyz` | **supported** (GDAL "ASCII Gridded XYZ") | One `x y z` line per cell on a regular grid, rows in either order; a missing point or row is NoData (never zero; a real 0.000 m value stays a value). The file carries no CRS: `--source-crs` is required. A list of scattered points is **not** a grid and is refused ("Ungridded dataset", or "not recognized as a supported file format" when the file has a header line or extra columns); grid it first (e.g. `gdal_grid`), knowing that interpolation invents values between the points. |
| Any other raster GDAL opens | | read through GDAL (not specifically tested) | Same rules for bands, CRS, units and NoData. |
| **Not supported** | `.zip`/`.hgt.zip` archives, E57/PLY clouds, scattered XYZ point lists, curvilinear or irregular NetCDF grids, generic HDF5 | refused with a message | Unzip / convert first (`gdal_translate`, PDAL). |

### Rules shared by all formats

* **Altitudes are converted to metres** from the declared unit (`feet`, `US survey
  foot`, `cm`, `km`...). The unit of a band, a NetCDF variable, or the vertical
  part of a compound CRS is used; `--vertical-unit m|ft|us-ft` overrides it. If
  none is declared, **metres are assumed and a warning says so**; a non-length unit
  (e.g. `K`) is refused. A declared `scale`/`offset` (packed integers) is applied;
  NetCDF `positive = down` (depths) flips the sign.
* **Coordinates are never assumed.** A file without CRS is refused until you give
  `--source-crs EPSG:xxxx` (which is itself refused if the file already declares
  one). Degrees vs metres is checked against the extent (a "geographic" CRS with
  coordinates like 2 600 000 is refused; a metric CRS with degree-sized
  coordinates triggers a warning). NetCDF with `degrees_north/east` axes but no datum
  is read as WGS84 **with an explicit warning**. A source in feet is reprojected to
  UTM metres; `--radius` is always in metres.
* **No silent vertical datum change.** A compound CRS (e.g. EPSG:2056+5728) has its
  vertical part ignored for positioning and the log says no geoid/ellipsoid
  conversion was applied.
* **NoData is never zero, and the relief is not extended beyond the data.** A declared
  NoData, NaN/inf, SRTM voids and the area outside the file's footprint all stay NoData.
  Only holes are filled (the number of filled cells is logged): gaps surrounded by data
  or crossing the window from side to side (a river, a seam between tiles), at most 50
  pixels from valid data. Everything else is left out of the model: the area beyond the
  edge of the file, the sea of a source that leaves it empty, any larger gap. See
  [Gaps and the edge of the data](#gaps-and-the-edge-of-the-data). Integer sources are
  resampled in floating point.
* **Nothing is silently downsampled**, and the log prints one `entrée ...` line per
  file: format, size, type, CRS, unit, and whether a band/variable/layer was chosen.

### Examples

```bash
# Multiband GeoTIFF: choose the band
python circle_dem_to_stl.py --input dem.tif --band 2 --cx 2601000 --cy 1201000 --radius 300 --diameter-mm 60 --vexag 1 --output out.stl
# Altitudes in feet that the file does not declare
python circle_dem_to_stl.py --input dem.tif --vertical-unit ft ...
# ASC without .prj
python circle_dem_to_stl.py --input terrain.asc --source-crs EPSG:2056 ...
# SRTM: two neighbouring tiles assembled, centre on the seam
python circle_dem_to_stl.py --input N46E007.hgt N46E008.hgt --lat 46.5 --lon 8.0 --radius 3000 ...
# GeoPackage with several raster tables
python circle_dem_to_stl.py --input data.gpkg --layer mnt_2m ...
# NetCDF: variable and time step
python circle_dem_to_stl.py --input grid.nc --layer elevation --band 3 --lat 46.1 --lon 7.15 --radius 3000 ...
```

### Point clouds (LAS/LAZ)

A point cloud has no grid, no NoData and no native resolution, so it is **not**
pushed through a raster reader: `lidar.py` reads it in chunks of 2 million points
(constant memory) and rasterises it **directly on the working grid** (same origin,
step and CRS), so there is no extra resampling unless you ask for `--lidar-cell`.

1. Points outside the requested window are dropped first; the file CRS (WKT record,
   extended record, or GeoTIFF keys) is reprojected chunk by chunk when it differs
   from the working CRS.
2. *Withheld* points and ASPRS noise classes 7 and 18 are always excluded. With
   `--lidar-select ground` (default) only class 2 is kept and each cell takes the
   **mean** of its ground points; `lowest`/`highest` take the min/max of all kept
   points. An **unclassified** cloud with `ground` is refused, with the classes
   present and the alternatives; `lowest` is only an approximation of the ground
   (too high under dense canopy or buildings) and the log says so.
3. Exact duplicate points (same x, y, z) are removed within each chunk.
4. Isolated aberrant cells (a spike or pit exceeding all ≥5 valid neighbours by more
   than `--lidar-spike-m` plus their own spread) are dropped; slopes and cliffs are not.
5. Empty cells are filled by `--lidar-interp`: `linear` (Delaunay triangulation of the
   valid cells around the holes; never extrapolates outside the cloud), `nearest`, or
   `none`; only up to `--lidar-max-gap` (default 50 cells). The log reports how many
   cells were filled and how many remain NoData, and warns when fewer than half of
   the cells contain a point (the filled detail would be invented).
6. `--lidar-cell M` (default: the working pixel) sets the cell of the gridded cloud;
   if it differs from the working pixel the grid is resampled **once** (bilinear if
   coarser, average if finer) and the log says so.

Altitude units come from the CRS vertical unit (or `--vertical-unit`); without
one the horizontal unit is assumed, **with a warning**. A CRS-less cloud needs
`--source-crs`. Several files must share one CRS. Tested on LAS 1.2/1.3/1.4
(point formats 0, 1, 3, 4, 6, 8) and LAZ.

### Common errors and fixes

| Message (start) | Cause | Fix |
|---|---|---|
| `... contient N bandes ... --band N` | multiband file or NetCDF/HDF5 time dimension | pick the band (the list shows each band's time/unit) |
| `... aucun système de coordonnées ... --source-crs` | no CRS in the file | `--source-crs EPSG:xxxx` |
| `... déclare déjà un CRS` | `--source-crs` given for a georeferenced file | drop `--source-crs` |
| `... impossible en degrés` | geographic CRS declared for metre coordinates | fix the file's CRS or give the right `--source-crs` |
| `... n'est pas une unité de longueur` | variable in K, %, ... | not an elevation; or force `--vertical-unit` if it really is |
| `fichier ASC tronqué` / `en-tête ASC mal formé` | cut-off or malformed `.asc` | re-download/regenerate the file |
| `tuile SRTM ... le nom doit être ...` | `.hgt` renamed or in an archive | restore `N46E007.hgt`-style names, unzip |
| `... N tables raster, choisis-en une avec --layer` | several raster tables | `--layer TABLE` |
| `... ne se déclare pas comme une altitude ... --layer` | NetCDF variable(s) not identified as elevation | `--layer NAME` if you know it is |
| `HDF5 non géoréférencé` | generic HDF5 without axes | convert to GeoTIFF/NetCDF with `gdal_translate` first |
| `aucun point classé « sol »` | unclassified cloud | classify it (PDAL `filters.smrf`, LAStools `lasground`) or `--lidar-select lowest`/`highest` |
| `la lecture des fichiers LAS/LAZ requiert laspy` / `compressé (LAZ) ... aucun décompresseur` | missing optional dependency | `pip install laspy lazrs` (in a venv) |
| `aucun point du nuage dans la zone demandée` | centre/radius outside the cloud | check `--cx/--cy` and the CRS of the cloud |

### Known limits

* HDF5: only the NetCDF-4/CF layout; no generic HDF5, no HDF-EOS.
* NetCDF: one variable and one slice per run; no curvilinear or irregular grids.
* GeoPackage: raster elevation tables only. GDAL applies the table's own scale/offset
  and declared precision: in our test, a table declaring a precision of 1 returned values
  rounded to the metre (that is GDAL reading the file as declared).
* SRTM: `.hgt` files only (no `.zip`); SRTM1 and SRTM3 tiles are not mixed in one assembly (warning).
* The vertical datum (geoid vs ellipsoid) is never converted; an undeclared unit is assumed to be metres.
* LAS/LAZ: the whole file is streamed once (time grows with the number of points,
  not with the size of the window); exact duplicates are removed only within a
  2-million-point chunk; `lowest` approximates the ground; `linear` filling can
  bridge a concavity inside the convex hull of the surrounding cells.
* Pre-existing behaviour that is now logged: `--min-elevation` (default −50 m)
  silently dropped lower terrain before; the number of excluded cells is now printed.


## Compatible data sources
 
Successfully tested on, among others: GSI DEM1A (Japan, XML/JPGIS —
requires pre-processing, not included here), IGN LiDAR HD (France),
USGS 3DEP (United States), swissALTI3D (Switzerland), Kartverket LiDAR
(Norway), Copernicus GLO-30 (worldwide). Geographic or projected CRS
is detected automatically (`gdalinfo` remains useful for checking what
a file actually contains, but manual reprojection beforehand is no
longer necessary).
 
## Automatic national sources (`--source`)

With `--lat/--lon`, the DEM can be fetched automatically instead of
passing `--input`:

```
python circle_dem_to_stl.py --source auto --lat 46.0 --lon 7.6 --radius 5000 --output out.stl
```

| `--source` | Provider | Access method |
|---|---|---|
| `swisstopo` | swissALTI3D (Switzerland) | STAC API, 0.5 m / 2 m GeoTIFF tiles |
| `kartverket` | Høydedata DTM1 (Norway) | WCS, resampled to the working resolution |
| `usgs` | 3DEP (United States) | TNM Access API, partial reads via `/vsicurl/` |
| `gsi` | GSI elevation tiles (Japan) | PNG tiles (DEM1A → DEM5A/5B → DEM10B), Web Mercator |
| `ign` | LiDAR HD MNT 0.5 m, fallback RGE ALTI 1 m (France, La Réunion) | WMS-R, raw float GeoTIFF; rendered/hillshade answers are rejected |
| `austria` | BEV airborne-laser DTM 1 m (Austria) | INSPIRE ATOM feed, 50 km Cloud-Optimized GeoTIFF tiles (EPSG:3035), partial reads via `/vsicurl/`; **explicit only**: `auto` never tries it |
| `copernicus` | GLO-30/90 (worldwide) | same as `--download-copernicus` |
| `auto` | national source covering the point, otherwise Copernicus | — |

The resolution fetched follows the working pixel size (derived from
`--diameter-mm`, `--radius` and `--print-spacing-mm`), so no 1 m data is
downloaded when a 15 m grid is enough. Data is cached under `--cache-dir`.

**Status.** These providers were written from official documentation but
could not be tested against the live servers in the development
environment (tests use simulated HTTP). USGS dataset names and the
Kartverket coverage name in particular should be confirmed on first real
use. Kartverket and IGN (WCS/WMS servers whose exact request syntax could not be confirmed) try several request variants and log the one accepted; IGN LiDAR HD coverage of France is still incomplete (a warning reports the share of the disc without data). GSI's FGD XML
download requires a login and is not automated; the tile service is used
instead (interpolated values, not the raw JPGIS mesh).

`austria` (BEV): licence (CC BY 4.0), 1 m resolution and the ATOM download
service come from the provider's official record; the feed address, the title
format of its entries and the idea that the COGs declare an unusable CRS come
from the open-source `geostl` project. However, the two tiles examined (E = 4,500,000
and 4,550,000 m) declare a correct CRS ("ETRS89-extended / LAEA Europe", EPSG:3035,
read with GDAL 3.13.3), so the small local VRT that wraps each tile and forces
EPSG:3035 is probably not necessary; it is kept as a precaution. **It has been run against the real BEV server at two places
only** (October 2026, series of 15.09.2025, both in high mountains): the
Grossglockner (one tile) and an area straddling the boundary between two tiles
(E = 4,550,000 m, vertical seam); in both cases the STL is watertight
(`check_stl.py`) and the seam is invisible on the preview. Not yet tried: a flat
area in the east, Vorarlberg, a four-tile corner, a radius of several tens of km.
This is why `auto` does not use this source: elsewhere, try `--source austria`
first and report any failure. Each tile is 50 km wide and is never downloaded
whole. The working CRS is the source CRS (EPSG:3035); its scale
distortion over Austria is at most about 0.08 % (measured at six points with GDAL
3.8.4), i.e. 0.15 mm on a 200 mm model. The BEV height accuracy is generally
±0.5 m, more in high mountains (provider's record).

### Place names (`--place`)

Instead of coordinates, give a place name; it is resolved with Nominatim
(OpenStreetMap) and, without `--input`/`--source`, the DEM source is chosen
automatically:

```
python circle_dem_to_stl.py --place "Mont Fuji" --radius 9300 --diameter-mm 250 --output fuji.stl
```

The first result (Nominatim's relevance order) is used; the other candidates
are listed in the log and `--place-pick N` selects another one. Administrative
areas and islands resolve to their centre, not to a summit: the log warns, and
`--lat/--lon` stays available for exact points. The public Nominatim server
allows 1 request per second and expects an identifying User-Agent (override with
`DEM2STL_USER_AGENT`, ideally with a contact); results are cached under
`--cache-dir`. Geocoding data © OpenStreetMap contributors (ODbL).

### Resolution and quality of automatic sources

The mesh is limited by `--print-spacing-mm`: the DEM is resampled (`average` when the
working grid is coarser than the source, see `--resample` below) to
`print-spacing / scale` metres per pixel (shown as "résolution de ré-échantillonnage"
in the log, e.g. 8 m for a 150 mm print of a 3 km radius). Sources therefore fetch
only the resolution needed for that grid (swisstopo 2 m instead of 0.5 m, GSI tiles
at the matching zoom, WCS/WMS requests at 2x the working resolution), which is why
cached files are much lighter than full-resolution downloads. Averaging 0.5 m or
2 m data down to 8 m gives nearly the same grid.

* `--full-res` fetches each source's finest resolution (heavy downloads, same mesh
  size). Kartverket and IGN otherwise resample server-side (nearest neighbour) at
  twice the working resolution before averaging; `--full-res` removes that shortcut.
* GSI resolution fallback: at fine working resolutions (zoom 16-17) only DEM1A is queried. A
  tile it does not have used to leave a hole (or "no tiles" for the whole area). It is now replaced by the next layers at their own zoom (DEM5A/5B at zoom 15,
  then DEM10B at zoom 14); those coarser mosaics are passed to the resampling step, which
  interpolates them (bilinear) under the finer one. The log says "repli de résolution". With
  several sources of different resolutions, the coarsest one decides the `--resample auto` choice.
* GSI elevation tiles are derived from the DEM (interpolated, 1 cm steps), not the
  raw JPGIS mesh; the number of tiles is capped, so very large areas use a coarser zoom.
* `--resample auto` (default) uses `average` when the working grid is coarser than the
  source (averaging is anti-aliasing) and `bilinear` when it is finer, i.e. when the DEM
  is enlarged (for example a 30 m global DEM on a small radius, or the Copernicus
  fallback). GDAL's `average` then behaves like nearest neighbour: the surface becomes
  flat plateaus separated by steps (in a synthetic test, 25 to 77 % of the cells had
  exactly zero slope at 4x to 16x enlargement). The log says "ré-échantillonnage bilinear"
  when this applies. `cubic` can overshoot near cliffs (about 7 % of a step, i.e. 71 m on
  a 1000 m cliff); `cubicspline` smooths without overshooting; `average` forces the
  previous behaviour.
* `--top adaptive` (prototype, default stays `grid`) replaces the top surface with a bounded-error
  RTIN mesh (right-triangulated irregular network): very few triangles where the relief is smooth,
  full resolution along the outline and along the 1024-cell block seams, so the walls and the
  bottom are unchanged and the solid stays watertight. The tolerance is nominal: measured against the
  DEM nodes on four real DEMs, the largest deviation is 1.5 to 1.8 times `--top-tol-mm` (to guarantee
  a given maximum error, use half of it). The reduction depends on how much real detail the working
  grid carries: a 5 m LiDAR-derived DEM printed at 250 mm / 0.2 mm goes from 2.5 M to 144 000 triangles
  at 0.05 mm (-94 %); noise-like terrain loses almost nothing. Memory is bounded by the block size.
  `--diagonal adaptive` is meaningless with it (refused).
* `--diagonal adaptive` (prototype, default stays `fixed`) cuts each top cell along the diagonal that
  joins the two corners with the smaller height difference. Same number of triangles, same
  vertices, walls and bottom: only the surface changes. Measured at the cell centres against the
  DEM itself on six real DEMs (5 m to 240 m pixels), the mean error drops by 4 % (grid step 4) to
  7 % (step 8) compared with the fixed diagonal, by 13 to 25 % on the 5 m LiDAR-derived DEM, and
  rises by 2.5 % at worst (the 240 m DEM). At the default 0.2 mm spacing the difference is very
  small; it matters for coarse spacings relative to the terrain detail.
* `--bottom merged` (prototype, default stays `grid`) triangulates the flat bottom with the squares
  of a balanced quadtree, each one as a fan that keeps every required boundary node, so there
  is no T-junction and the solid stays watertight. The geometry is identical (same plane, same
  outline, same volume: the tests compare them to 1e-12) but the file has about half the
  triangles (the bottom drops by 96 to 99 %), a lower RAM peak and a smaller file. The
  triangles stay compact (aspect ratio at most 4 in the tests). `merged-rows` uses greedy
  rectangles instead (the method of the my3dmaps project): even fewer triangles, but extremely
  elongated ones (aspect ratio above 10^5 measured on a disc), kept for comparison only.
  Check the result in your slicer and CAM before making it the default.
* HTTP 429 and 503 answers with a `Retry-After` header (seconds or HTTP date) are waited for, at
  least that long, capped at 60 s.
* Downloads are checked: a truncated answer is retried, and an error page received with
  HTTP 200 (XML/HTML instead of a GeoTIFF or PNG) is neither cached nor reused if an
  earlier run stored it (the invalid cached file is downloaded again).
* `compare_dems.py A B --lat .. --lon .. --radius .. --pixel-size ..` measures the
  difference between two DEMs on the working grid (mean offset, RMS, 95th percentile)
  to check an automatic source against manually downloaded data.

### Gaps and the edge of the data

Cells of the area without data are treated in two ways:

* **Holes** (voids, a missing strip between two tiles, a river left empty by the source):
  filled by inverse-distance interpolation (`gdal.FillNodata`), as long as they are
  surrounded by data or cross the window from side to side, and no cell is more than 50
  working pixels from valid data. The number of filled cells is logged.
* **Everything else** stays empty: the area beyond the edge of the file, the sea of a source
  that leaves it as NoData, a large missing zone. A disc cut by
  the edge of the file is therefore truncated instead of completed, which `--preview` shows;
  with `--sea-level` these cells become sea. The log counts the cells of the area that the
  former filling would have extended from the edge of the data.

The former filling extended the relief up to 50 working pixels beyond the data. At the
default `--print-spacing-mm 0.2` that is 10 mm on the printed part, whatever the radius, and
it was invented relief: on a cliff coast, for instance, a plateau at the height of the cliff
advanced into the sea. `--extrapolate-edges` brings it back (the log then says how many cells
were extrapolated). With `--sides`, a shape that overshoots the data is refused, as before;
`--extrapolate-edges` accepts it when the overshoot is within 50 pixels. The option has no
effect on LAS/LAZ input, which never extends beyond the points (see `--lidar-interp`).
The `--report` JSON counts these cells in `source.no_data`.

### Filling gaps with Copernicus (`--fill-with-copernicus`)

National sources stop at their border (e.g. swissALTI3D at the Italian side of the
Matterhorn) and some are incomplete (IGN LiDAR HD). With `--fill-with-copernicus`
(needs `--lat/--lon` or `--place`; works with `--source` or `--input`), parts of the
disc without data are filled from Copernicus GLO-30, which is downloaded only if
there are gaps. Copernicus is first shifted in altitude by the median difference
measured on the area both rasters cover (logged, with its interquartile range); the
fill is abandoned if that shift exceeds 50 m. Limits: Copernicus is a surface model
(trees, buildings) at 30 m, the national data are usually bare-earth, so the
junction shows a resolution and character change; the log reports the share of the
disc filled.

### Regular polygons (`--sides`)

Instead of a disc, a regular N-gon (3 to 360 sides) tied to a **reference circle**, with optionally rounded
corners:

```
python circle_dem_to_stl.py --input dem.tif --cx 2601000 --cy 1201000 --radius 300 --diameter-mm 60 \
    --base-mm 3 --sides 6 --circle-mode circumscribed --corner-radius 2 --output hexagon.stl
```

* The reference circle is the one given by `--radius` (ground metres) and `--diameter-mm` (its printed diameter):
  same scale as a disc. `--circle-mode` is mandatory: `inscribed` (vertices lie on the circle) or `circumscribed`
  (sides are tangent to the circle). In circumscribed mode the part **exceeds** `--diameter-mm` (by up to a factor
  1/cos(π/N) at the vertices): the log and `--dry-run` print the real printed size (`pièce : W x H mm`).
* `--corner-radius R` rounds the corners with circular arcs tangent to the sides, `R` in **printed mm** (0 or
  omitted: sharp corners). Maximum: the apothem of the part; at that value the outline becomes a circle.
* With rounded corners in `inscribed` mode, `--ref-outline` chooses **what the reference circle describes**:
  * `sharp` (default): the polygon *before* rounding. Its vertices lie on the circle, then the rounding cuts
    them back, so the rounded part sits slightly inside the circle and no longer touches it (extreme radius
    of the outline: (a − R)/cos(π/N) + R = ref − R·(1/cos(π/N) − 1), printed in the log).
  * `rounded`: the **rounded outline**. The reference circle is tangent to the circle of every corner
    rounding: it touches the top of each arc and never exceeds it, so the part fits in `--diameter-mm`
    exactly. The underlying sharp polygon is then larger (its vertices lie beyond the circle). The centres of
    the roundings are at (reference radius − R) from the centre; the maximum `--corner-radius` becomes the
    reference radius (the outline is then the reference circle itself). File names get the suffix
    `_cercle-tangent`, so the two readings never overwrite each other.
  In `circumscribed` mode nothing changes between the two: rounding does not move the sides, which stay
  tangent to the circle (the option is accepted and the log says it has no effect).
* Default orientation: a vertex pointing north for odd N, a side facing north for even N. `--rotation` turns it
  by that many degrees, counter-clockwise.
* The edge is **exact**: straight sides (no staircase), vertices and tangent points on the theoretical outline,
  vertical wall. Arcs are split so that the chord-to-arc gap never exceeds `--corner-tol-mm` (default 0.005 mm
  printed; vertices lie on the arc, chords inside it), or into `--corner-segments` segments per arc if you
  prefer to force that number. The log gives the segment count and the resulting gap.
* Incompatible with `--polygon`, `--clip-to-place` and `--smooth-wall`. Like `--smooth-wall`, the mask must be
  a single piece without holes: use `--sea-level` or `--fill-with-copernicus` for areas without data.
  Downloaded data (`--source`, `--download-copernicus`) covers the disc that encloses the whole part.
* `check_stl.py` accepts the same `--sides`, `--circle-mode`, `--corner-radius` (mm), `--ref-outline` and `--rotation` to check
  the exact footprint of the part.

### Cutting by polygon (`--polygon`)

Instead of a disc, model only the inside of a contour (island, border, commune...).
The contour can be any vector file OGR reads (GeoJSON, GPKG, shapefile), in any CRS:

```
python circle_dem_to_stl.py --polygon corse.geojson --polygon-where "NAME='Corse'" \
    --diameter-mm 200 --vexag 2 --base-mm 3 --output corse.stl
```

* The centre and size come from the contour: `--radius`, `--lat/--lon`, `--cx/--cy`,
  `--place` and `--smooth-wall` are not allowed with `--polygon`. Without `--input` or
  `--source`, the DEM source is chosen automatically (`--source auto`).
* `--diameter-mm` is the **largest printed dimension** of the contour.
* All polygons of the file (or those selected by `--polygon-where`, an OGR attribute
  filter) are merged; holes and disjoint parts (archipelago) are kept.
* The wall follows the contour as a staircase at the working resolution. A tip thinner
  than one grid step is absent from the mesh, so the printed dimension can be a few
  steps shorter than requested for very pointed shapes.
* `--fill-with-copernicus` and `--sea-level` work with a contour; polygons crossing the
  antimeridian are not handled.

### Checking an STL (`check_stl.py`)

```
python check_stl.py cervin.stl --diameter-mm 150 --base-mm 3 --radius 3000 --relief-m 2247
```

Re-reads the binary STL (independently of the meshing code) and checks the largest
horizontal dimension, the base thickness, the elevation range implied by the scale
(with `--radius`, optionally against `--relief-m`, tolerance 1 %) and watertightness
(open, non-manifold or badly oriented edges), the sign of the signed volume (an entirely
inside-out solid is consistent edge by edge and only this sign reveals it) and degenerate
triangles. Exit code 1 on failure; `--tol-mm`
adjusts the dimension tolerance, `--no-watertight` skips the edge test on huge files.

Two further checks run with the watertightness test (`--no-topology` skips them):

* **Shells and Euler characteristic.** Shells are the parts connected through a shared edge;
  for each one the tool reports V − E + F (2 for a closed surface without a handle, 0 for a
  ring-shaped solid, and so on). Several islands or a hole in the data legitimately give several
  shells or a handle, so without an expectation this is only described; two shells touching at
  a single vertex are always refused. `--expect-shells N` turns it into a check (N shells, each
  with V − E + F = 2).
* **Inverted triangles.** A triangle moved across its neighbours (for instance by a contour
  projection) keeps its edges consistent, so neither the edge test nor the signed volume sees
  it; it is reported when its normal opposes those of all three neighbours (dot product < −0.5).
  A sharp ridge or a cliff edge is not an inversion. Neither check detects self-intersections.

### Place outlines (`--place ... --clip-to-place`)

`--clip-to-place` cuts along the OpenStreetMap outline of the place (island, commune,
country) instead of a circle; the outline is requested from Nominatim
(`polygon_geojson`, simplified server-side to about 30 m) and then behaves exactly like
`--polygon`. A result that is only a point (a summit, say) has no outline: the tool
says so, and `--place-pick` lets you choose another candidate.

```
python circle_dem_to_stl.py --place "La Réunion" --clip-to-place --diameter-mm 200 --output reunion.stl
```

## Data sources, licences and attribution

The tool only downloads data; the conditions below matter when you **share or sell** a
model built from them. Wording checked against each provider's own pages in October
2026; this is not legal advice, so check the current terms before publishing.

| Source | Terms | Attribution |
|---|---|---|
| swisstopo (swissALTI3D) | Open Government Data: free use, including commercial | Mandatory: "© swisstopo" or "Federal Office of Topography swisstopo" |
| IGN (LiDAR HD, RGE ALTI) | Licence Ouverte Etalab 2.0 (open data) | Source "IGN"; RGE ALTI is less accurate on steep terrain (about 7 m vertical, per the Earth Engine catalogue) |
| Kartverket (Høydedata DTM1) | CC BY 4.0 | Credit Kartverket |
| BEV (Austria, airborne-laser DTM 1 m) | CC BY 4.0 (stated on the provider's record) | "© BEV – Bundesamt für Eich- und Vermessungswesen"; the record prescribes no credit line, and CC BY also asks you to indicate changes |
| USGS 3DEP | Public domain | Credit requested ("U.S. Geological Survey, 3D Elevation Program") |
| GSI Japan (elevation tiles) | GSI content terms (PDL 1.0) | State the source (出典:国土地理院ウェブサイト) **and** that the data was processed, e.g. "地理院タイル (標高タイル(基盤地図情報数値標高モデル))を加工して作成"; never present it as made by GSI |
| Copernicus DEM GLO-30 | Free licence | "produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved" when adapted; the licence also asks for a no-liability notice when distributing |
| OpenStreetMap / Nominatim (`--place`) | ODbL; public server limited to 1 request/second with an identifying User-Agent | "© OpenStreetMap contributors" |

### Safety checks and speed

* Parameters are validated before any download or computation: `--radius`, `--diameter-mm`,
  `--base-mm` and `--print-spacing-mm` must be strictly positive, `--vexag` must not be
  negative (0 only warns: flat disc), latitude/longitude must be in range. Invalid or missing
  input rasters, unreachable data sources and unwritable outputs give a one-line
  `Erreur : ...` instead of a traceback (`--debug` shows the full trace).
* The output folder is created if needed and checked first; the STL is written under a
  temporary name and renamed at the end, so a failure never leaves a truncated file.
* Before computing, the tool prints an estimate (grid, triangles, STL size, RAM). If the
  estimated RAM exceeds 80 % of the available RAM (Linux), it stops; `--force` overrides.
  The constants come from a measured case (4.9 M triangles: about 1.05 GB peak); beyond it
  the estimate is an extrapolation.
* The STL writer is vectorised: a 4.9 M-triangle model now takes about 19 s end to end on
  the development machine, versus several minutes before. Vertices are identical to the old
  writer; normals match to float32 rounding.

### Output folder and formats (STL, 3MF, OBJ)

Generated files go to `output/` (created if needed; change it with `--output-dir`).
`--output name.stl` (no folder) is placed there; a path with a folder is used as given.
`--output` is optional: the name is then derived from the place or contour and the
diameter (`mont-fuji_250mm.stl`) and never replaces an existing file (`_2`, `_3`...).

`--format stl|3mf|obj`, a comma-separated list, or `all` writes several formats from the
same mesh; the extension of `--output` selects the format when `--format` is absent.
A name with an unknown extension is still written as an STL under exactly that name.

| Format | Role | Content |
|---|---|---|
| STL | compatibility and manufacturing | triangles and normals; no units (millimetres by convention) |
| 3MF | 3D printing | indexed mesh, **unit declared (millimetre)**, metadata (title, description, data-source attribution, date, application), about 4x smaller |
| OBJ | exchange and visualisation | vertices and faces only, millimetres stated in a comment; no normals, no `.mtl` (a DEM relief has no material to keep) |

All three contain the same triangles, in the same order, with the same orientation
(outward normals) and the same float32 coordinates: `check_stl.py` reads any of them, and
the test suite compares them. The 3MF is a real 3MF package written with the standard
library (ZIP/OPC, core specification); the official `lib3mf` is used by the tests, when
installed, as an independent validator. Not checked here: opening in PrusaSlicer itself.
Measured on 4.9 M triangles: STL 245 MB, 3MF 61 MB, OBJ 192 MB, about 40 s for the three.
Blender and some viewers import OBJ with Y up: rotate if needed (the file is Z up).

### Relief height and JSON report (`--report`)

Before meshing, the log states the printed relief: altitudes, vertical range, height above the
base in mm with and without the vertical exaggeration, and share of the diameter. When the
exaggeration is far from a rough reference (relief of about 22.5 % of the diameter, exaggeration
between 1 and 8), an indicative `--vexag` is printed; it is never applied, the use (cutter, board
thickness, intended look) decides. The relief needs the data, so `--dry-run` cannot print it.

`--report` (optionally `--report file.json`; by default next to the model, same name in `.json`)
writes a JSON file with the source and its attribution text, scale, relief, exaggeration, triangle
counts, mesh options, the cells without data (`source.no_data`) and the command line. It is written after the models; if it cannot be written,
a warning is printed and the models are kept.

### Shaded preview (`--preview`)

`--preview` (optionally `--preview file.png`; by default next to the model, same name in `.png`)
writes a hillshade PNG of the area that is actually cut, **before** meshing, so a wrong centre,
holes in the data or a flat-looking relief show up without waiting for the mesh. North is up,
light comes from the north-west, the shading is computed on the printed relief (altitudes times
`--vexag`), and everything outside the piece is transparent. Images wider than 1600 px are averaged
down. Like `--report`, it needs the data, so `--dry-run` only lists the file it would write.

### Previewing and smooth walls

* `--dry-run` validates the parameters, prints the estimate (grid, triangles, size of each
  requested format, RAM) and the files that would be written, then stops: nothing is
  downloaded, computed or written.
* `--smooth-wall` needs a mask with a single outline. When the disc contains a hole (a large
  area without data) or several parts, the tool now stops with an explanation instead of
  producing an open mesh; use `--sea-level` if the hole is sea, or drop `--smooth-wall`.
  Triangles without area (a contour vertex exactly on the circle) are removed.

## Tests
 
```bash
pytest tests/
```

The input-format tests build tiny synthetic files (a tilted plane, so the
expected altitude of every cell is known exactly) and check values, size,
orientation, NoData, corruption and units for each format, then run the whole
command line and measure the STL that comes out (watertight, outward faces,
relief matching the known slope, summit on the expected side):
`tests/test_input_formats.py` (GeoTIFF, ASC, SRTM), `tests/test_input_containers.py`
(GeoPackage, NetCDF, HDF5), `tests/test_lidar.py` (LAS/LAZ), and
`tests/test_input_formats_cli.py`. LAS/LAZ tests are skipped if `laspy` is
missing, and the NetCDF/HDF5 ones if `netCDF4`/`h5py` are (the reason is
printed). Run one module with e.g. `pytest tests/test_lidar.py -q`.
 
The tests reproduce the synthetic scenarios used during development
(conical island, isolated artifact, water-surface noise touching the
shoreline, checkerboard pattern) and verify watertightness, normal
consistency, and the expected behavior of each option.
 
## Notable changes
 
- **No relief beyond the data. Behaviour change**: `gdal.FillNodata` used to extend the relief up to 50
  working pixels beyond the data (10 mm on the printed part at the default `--print-spacing-mm`), outside
  the file's footprint and around a large empty zone (a sea left as NoData). Only holes are filled now;
  `--extrapolate-edges` brings the former behaviour back, and the log and the `--report` JSON count the
  cells concerned (`tests/test_fill_edges.py`). A disc cut by the edge of the file is truncated, and a
  `--sides` shape that overshoots the data is refused unless `--extrapolate-edges` is given.
- **`--preview`**: shaded PNG preview of the cut area, written before meshing
  (`tests/test_preview.py`, `tests/test_preview_cli.py`).
- **`check_stl.py`: shells, Euler characteristic and inverted triangles** (`--expect-shells`,
  `--no-topology`), and **printed relief height in the log plus `--report`** (JSON next to the
  model). Tests: `tests/test_topology_checks.py`, `tests/test_report.py`, `tests/test_report_cli.py`.
- **`--vexag-curve`**: optional non-linear exaggeration (altitude → equivalent altitude, monotone PCHIP),
  announced as not faithful (`tests/test_vexag_curve.py`, `tests/test_vexag_curve_cli.py`).
- **Modular restructuring** (`raster.py`/`mesh.py`/`stl_io.py`):
  eliminates the duplication that used to exist between the old
  monolithic `circle_dem_to_stl.py` and `island_dem_to_stl.py`.
- **Automatic CRS reprojection**: no more manual `gdalwarp` before
  each new site when the source is in geographic coordinates.
- **`--lat`/`--lon`** as an alternative to `--cx`/`--cy`.
- **Built-in Copernicus DEM download** (`--download-copernicus`).
- **More input formats**: ESRI ASCII Grid, SRTM `.hgt`, GeoPackage raster, NetCDF,
  NetCDF-4-style HDF5 and LAS/LAZ point clouds, with explicit rules instead of
  guesses (bands, layers, CRS, units, NoData). GeoTIFFs in Float32/Float64 give
  byte-identical STL files to before. **Behaviour change**: integer DEMs
  (e.g. Int16) used to have their *resampled* elevations rounded to whole
  metres; they are now resampled in floating point.
- **`island_dem_to_stl.py` removed**: its use case (relief surrounded
  by water) is covered by `--sea-level`, whose approach (circle with a
  flat sea) proved preferable in practice to following the exact
  shoreline.
## License
 
GNU GPLv3 — see [LICENSE](LICENSE). Any modified or derivative version
must remain open source under the same license.
