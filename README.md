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
| `--radius` | Circle radius (reference circle with `--sides`), in meters | required |
| `--diameter-mm` | Final printed diameter, in mm | prompted if omitted (250) |
| `--base-mm` | Flat base thickness, in mm | prompted if omitted (3) |
| `--vexag` | Vertical exaggeration | prompted if omitted (1.0) |
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
| `--rotation` | With `--sides`: rotation in degrees, counter-clockwise | 0 |
| `--corner-tol-mm` | With `--sides` and `--corner-radius`: maximum arc sagitta, in printed mm | 0.005 |
| `--corner-segments` | With `--sides` and `--corner-radius`: segments per arc (replaces `--corner-tol-mm`) | computed |
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
regular_polygon.py     # pure numpy: regular polygons (--sides), exact membership, projection, arcs
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
  The reference circle always describes the sharp polygon, even when rounded, so a rounded part is slightly
  smaller at the corners (extreme radius of the outline: (a − R)/cos(π/N) + R, printed in the log).
* Default orientation: a vertex pointing north for odd N, a side facing north for even N. `--rotation` turns it
  by that many degrees, counter-clockwise.
* The edge is **exact**: straight sides (no staircase), vertices and tangent points on the theoretical outline,
  vertical wall. Arcs are split so that the chord-to-arc gap never exceeds `--corner-tol-mm` (default 0.005 mm
  printed; vertices lie on the arc, chords inside it), or into `--corner-segments` segments per arc if you
  prefer to force that number. The log gives the segment count and the resulting gap.
* Incompatible with `--polygon`, `--clip-to-place` and `--smooth-wall`. Like `--smooth-wall`, the mask must be
  a single piece without holes: use `--sea-level` or `--fill-with-copernicus` for areas without data.
  Downloaded data (`--source`, `--download-copernicus`) covers the disc that encloses the whole part.
* `check_stl.py` accepts the same `--sides`, `--circle-mode`, `--corner-radius` (mm) and `--rotation` to check
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
