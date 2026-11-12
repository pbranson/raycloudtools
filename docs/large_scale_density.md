# Volumetric stem density for large forests: domain decomposition and output

Design note for running the vertical-cylinder density workflow over forests of order 300 x 300 m scanned with many
billions of rays, producing a single 3D gridded output for hydrodynamic models.

**Status:** implemented as `raysplit ... tiles`, `raydensity` and `scripts/raydensity_tiles.sh`. Tested on synthetic
scenes (a 6 m plot in 3 x 3 tiles gave the same cells as a single run; with `--ground`, layers on a 10% slope were
within 0.5% of flat ground over the plot and 3.5% per cell), but not yet at the scale of billions of rays. The
estimator, its accuracy and its known biases are described in
[`pneumatophore_density_report.md`](pneumatophore_density_report.md).

```console
raydensity_tiles.sh forest.ply forest.zarr 400000,6100000 25 12,12 1,4.3 32 -- --ground forest_mesh.ply \
  --heights 0,0.6 --stem_diameter 0.01 --footprint 0.01,0.003,0.05 \
  --self_shadow_correction --hybrid --min_observed_fraction 0.9
```

The last three options are the recommended estimate (report, section 10): the hybrid, with cells observed less than
0.9 masked. Use the measured stem diameter and a footprint calibrated with a rod test (report, section 9.2), and
drop `--hybrid` where range or pose errors exceed about 1-2 cm.

Here `1,4.3` is the absolute height slab of rays to keep (required): for ground between 1.2 and 3.5 m, from the lowest
ground minus 0.2 m to the highest ground plus the top of `--heights` (0.6 m) plus 0.2 m. Keep it small: a tall slab
keeps the canopy and sky rays, which multiplies the size of the tile files. The script exits non-zero if the split
or any tile fails.

```python
import xarray as xr
cells = xr.open_zarr("forest.zarr", consolidated=False)   # density, solid_fraction, ... on (z, y, x)
```

## Why the current interface doesn't scale

| Issue | Detail |
|---|---|
| One dense in-memory grid | `DensityGrid` holds 12 bytes per voxel (+4 for the occlusion correction, +4 temporary). 300 x 300 m at 5 cm, 0.6 m high plus margins, is 6002 x 6002 x 14 = 5.0e8 voxels, about 10 GB (plus 0.3 GB of column peaks). At 2.5 cm it is 12000^2 x 26 = 3.7e9 voxels, which exceeds the int voxel index (2.1e9), so the grid can't be built at all. `rayrender` is for plots, not whole domains. |
| Single-threaded walk | `calculateDensities` walks every ray through the grid in one thread: roughly 1e6-1e7 rays/s, so hours for 1e10 rays. |
| Outputs are images and CSV | The rendered image is a column integral, `--profile` is plot-mean, and `--cells` is a CSV of 540k rows for 1 m cells. There is no volumetric, georeferenced, chunked output. |
| Absolute heights | Layers are horizontal in the cloud's frame. Over 300 m, ground relief is larger than the 0.1 m layers, so heights must be relative to the ground. |

## Proposed workflow

```
ground model  ->  ray binning into tiles (+ halo)  ->  per-tile density (parallel)  ->  chunked store (no merge step)
```

### 1. Global grid and tiles

* A global origin and voxel width w (e.g. 5 cm) define the voxel lattice. Output cells (e.g. 1 x 1 x 0.1 m) and
  tiles (e.g. 25 x 25 m) are integer multiples of w, aligned to the same origin, so every voxel belongs to exactly
  one cell and one tile, with no resampling at tile edges.
* Each tile is processed with its bounds expanded by a halo h >= max(pool_radius, occlusion_radius) + 2w (the pooled
  neighbourhood, or the occlusion neighbourhood when `--occlusion_radius` is larger, plus two voxels for the
  neighbour prior), 0.15 m by default (`raydensity` rejects a smaller `--halo`). Only the tile's
  core is output. With the halo, a tile's core voxels receive exactly the same rays, in the same order, as in a
  single global run, so the result is identical to float rounding: tiles have no seams and need no blending.
* The pooled neighbourhood is (2k+1) x (2k+1) voxels with k = round(pool_radius / w), i.e. +-(k + 1/2) w: +-0.075 m
  for the default pool_radius 0.05 at 5 cm voxels. A wider window (0.1: +-0.125 m, halo 0.2 m) reduces the
  self-shadowing bias of `density` (about +$`0.47\,d_\mathrm{eff}`$ / window width: +6-8% at the default, +2-4% with 0.1, for 1 cm
  stems), but reads lower where stems are patchy or clustered at a smaller scale (cable roots: 0.87 with the default,
  0.77 with 0.1). `--self_shadow_correction` removes the self-shadowing bias for uniform stands, but is off by
  default, as it makes patches and cable roots read low (report, sections 8.1 and 10). It needs `--stem_diameter`
  and the pooled estimator (the default), and `raydensity` refuses it otherwise. It doesn't change the halo. For
  clustered stands whose lower layers are well observed, the experimental single-voxel configuration (`--pool_radius
  0 --self_shadow_correction --self_shadow_second_order --occlusion_radius 0.05 --ratio_bias_correction
  --no_neighbour_priors`; report, section 8.1.2) resolves root lines and patches; there the halo is set by the
  occlusion radius, 0.15 m at 5 cm voxels, as by default.
* `--hybrid` (with `--self_shadow_correction`; off by default, and not with the four single-voxel options, which it
  sets itself) computes both estimates in the same pass and gives each cell the single-voxel one where that
  observes at least 0.9 of the cell's voxels, otherwise the default-window one (report, section 8.1.3). It is the
  recommended configuration where range and pose errors are about 1 cm or less; with larger errors use
  `--self_shadow_correction` alone (report, section 10). Its single-voxel part takes the occlusion response from the
  pool_radius neighbourhood, so the default halo covers it, and a tiled run equals a single run. It holds a second
  copy of the voxel grid while estimating, so a tile needs about twice the memory below, and takes about 5% longer
  than `--self_shadow_correction` alone.
* `raydensity`'s voxel grid covers the halo (rounded up to whole voxels) plus one more voxel. The split (step 2)
  must cover that grid with some margin, so that the ray ends it clips lie outside the grid: split with a halo 1.5
  voxels larger, 0.225 m by default (`scripts/raydensity_tiles.sh` does this). A split halo equal to the
  `raydensity` halo puts clipped ray ends on a voxel face of the grid; this used to change a few edge voxels (by up
  to 0.4%), and now also matches a single run (`RayDensityTiles` test), but the larger split halo is still advised.
* A 25 m tile at 5 cm with the default halo is 508 x 508 x nz voxels, with nz = (relief + 0.6 m)/w + 4: 26 voxels
  for 0.5 m of relief within the tile (about 110-135 MB in total), 56 for 2 m (about 290 MB). So 32 tiles run
  concurrently in 3-10 GB. Use 25 m tiles (12.5 m in steep terrain).

### 2. Ray binning (one streaming pass)

Each ray must go to every tile it crosses, not only the tile of its end point, otherwise tiles lose the rays that
pass through them and are biased high.

* `raysplit grid wx,wy,0 overlap` already does this (it tests each tile in the ray's bounding box and clips the
  ray to it, marking rays whose end is clipped as unbounded). But it also clips the **start**, which loses the
  scanner range that the beam footprint correction needs (the footprint, and so the effective stem width, depends
  on range from the scanner). The density tiles should keep the original start, and only clip the end (with zero
  alpha when the end is outside the tile); `DensityGrid` clips to its own bounds and computes the range offset
  from the original start.
* Before tiling, discard rays (or ray portions) outside the height slab of interest (ground to the top of the
  stems, plus margin). Hand-held scans are dominated by long rays into the canopy and sky, so this removes most of
  the data volume. Rays crossing the slab are short (a few metres), so the duplication across tiles is small.
* The pass is I/O bound: about 36 bytes per ray in .ply (48 with `RAYLIB_DOUBLE_RAYS`), so 1e10 rays is about
  360-480 GB read once. With the domain given (`raysplit ... tiles x0,y0 w halo zmin,zmax nx,ny`, as the script
  does), the cloud is read once per 256 tiles; without it, an extra pass first bounds the tiles by the rays. When
  there are more tiles than open files, `raysplit` makes several passes (12 x 12 = 144 tiles is one pass); a
  two-level split (100 m blocks, then tiles) keeps it to two.
* SLAM data is time ordered, so consecutive rays are spatially local. Writing per-tile buffers of a few MB keeps
  the writes sequential.

### 3. Height above ground

* Don't flatten the cloud with `raytranslate subtract ground_mesh.ply`: it shifts each **whole ray** by the ground
  height at its end point, so the rest of the ray is misplaced vertically by the slope times its horizontal distance
  from its end. Hand-held rays cross several metres of the stem layer, so on a 10% slope this put the rays' paths
  tens of centimetres out, and gave 0.43-0.71 of the flat-ground $`N d`$ in simulation (only this slope has been
  tested).
* Instead (as implemented), keep absolute heights in the voxel grid (each tile covers its local relief), and assign
  each voxel to a height-above-ground layer by the height of its centre above the ground beneath the column's centre
  (`--ground`, from the mesh at voxel resolution). Layers are therefore misregistered by up to half a voxel (2.5 cm
  at 5 cm voxels) where the ground height isn't a multiple of the voxel, so keep voxels at most half the layer
  height. The split's height slab must cover the ground relief plus `--heights`, with a margin of a few voxels. The
  ground mesh comes from `rayextract terrain` (itself run per tile, as in `scripts/rayextract_trees_large.sh`).
* Ground returns must be removed (zero alpha) before density estimation, using the ground model, as now.

### 4. Per-tile density (embarrassingly parallel)

A dedicated tool, `raydensity`, rather than more options on `rayrender`, which is an image renderer:

```
raydensity tile_11_7.ply forest.zarr --origin x0,y0 --tile_width 25 --tiles 12,12 --tile 11,7
    --voxel 0.05 --cell_size 1,1,0.1 --heights 0,0.6 --ground ground_mesh.ply
    --stem_diameter 0.01 --footprint 0.01,0.003,0.05 [--voxels]
```

`--origin`, `--tile_width`, `--tiles` and `--tile` are required (tile indices run from 0 to tiles - 1). Its defaults
differ from `rayrender`'s: `--density_model vertical_cylinder`, `--estimator pooled`, and the flat top model off
(`--flat_top` enables it).

It runs the existing `DensityGrid` pipeline (calculateDensities, applyOverlapCorrection, addNeighbourPriors) on
the tile plus halo, aggregates the core voxels into cells by height above ground, and writes its chunk of the
output store. Tiles are independent: run with a job array (SLURM), `xargs -P`, or threads. With `--ground`, only the
mesh triangles overlapping the tile and its halo are used, so one mesh for the whole domain costs each tile little.

### 5. Output: a chunked store with one chunk per tile

* Use Zarr (v2): a directory holding JSON metadata (`.zarray`, `.zattrs`) and one binary file per chunk. Without
  compression, a chunk is just raw little-endian float32 data, so it is writable from C++ with no dependencies
  (zlib compression is a small, optional step). It is read directly by xarray, zarr-python, GDAL and many
  visualisation tools, and converts to NetCDF or Cloud Optimised GeoTIFF.
* Chunk = one tile's core, so each tile writes only its own chunk files, and there is no merge step and no
  locking. Every tile writes the store's metadata (from the global extent and the options, so identical for every
  tile) before its estimate, by writing a temporary file and renaming it, so concurrent writers don't conflict. If
  the store already has different metadata (a tile run with other options or another domain), the tile fails
  without writing, so a store can't mix parameters. Each tile also writes `tiles/i_j.json` (its command line, input
  and halo), which Zarr readers ignore.
* Variables (as implemented), on dimensions (z, y, x) for 1 x 1 x 0.1 m cells (300 x 300 x 6 per variable, small):
  `density` ($`N d`$, the mean of the observed voxel densities: the main output, which reads high by the self-shadowing
  bias of the report's section 8.1 unless `--self_shadow_correction` is given), `pooled_density` (total corrected
  hits over total weighted path of the observed voxels, i.e. the path-weighted mean of the voxel densities: within 3% for stems uniform within the cell, but about
  half the truth where stems are patchy within the cell, and low where the density falls with height),
  `patchiness_ratio` (density / pooled_density: 1.03-1.09 for random stems in the lowest 0.2 m, 1.6-2.4 for patches
  and cable roots), `solid_fraction` ($`N \pi d^2 / 4`$ from `density`), `hits`, `weighted_path` and `observed_fraction`
  (the fraction of the cell's voxels with more than 2 rays after the neighbour top-up, or without it with
  `--no_neighbour_priors`, and with `--hybrid` that of the default-window estimate; cells below 0.9 read low) and,
  with `--hybrid`, `hybrid_weight` (the fraction of the cell's observed voxels taking the single-voxel estimate: 0
  where the cell takes the default-window estimate, otherwise at least 0.9) (float32; unobserved cells are NaN, with
  `observed_fraction` 0), plus coordinates `x`,
  `y`, `z` (cell centres, float64, so that georeferenced coordinates are exact). With `--voxels`, also
  `voxel_density` on (zv, yv, xv) as heights above ground (6000 x 6000 x 12 float32, 1.7 GB, for 3D inspection). Not
  yet output: the mean $`d_\mathrm{eff}/d`$, and a flag for cells at the solid fraction limit.
* Attributes (`.zattrs` of the group): title, quantity, density_model, origin, tile_width, num_tiles,
  voxel_width, cell_size, height_range, heights_above_ground, voxels, stem_diameter, footprint, return_range,
  estimator, pool_radius, min_gap_ratio, min_sin_zenith, max_sin_zenith, occlusion_correction and
  flat_top_compensation, plus entries written only when the corresponding option is set, so that stores written
  without them are unchanged: `"self_shadow_correction": true`, `"self_shadow_second_order": true`,
  `"occlusion_radius": <m>`, `"ratio_bias_correction": true`, `"neighbour_priors": false` (for
  `--no_neighbour_priors`), `"hybrid": true` and `"min_observed_fraction": <f>`. Tiles run with different options can't share a store. Not yet written: the CRS
  (the coordinates are in the cloud's frame) and the software version. The per-tile command lines are in `tiles/i_j.json`.
* Mask: `--min_observed_fraction f` (default 0, off; 0.9 recommended) sets `density`, `pooled_density`,
  `patchiness_ratio` and `solid_fraction` to NaN in cells observed less than f, and keeps `observed_fraction`,
  `hits`, `weighted_path` and `hybrid_weight`, so masked cells can still be inspected and aggregated. With `--ground`
  the layers are heights above the ground, so a low layer isn't blanked for lying partly below it, as it would be
  with absolute heights over uneven ground.
* Because `density` is a mean of voxel values, and `observed_fraction`, `hits` and `weighted_path` are also stored,
  coarser cells (e.g. 5 m for a larger model domain) can be aggregated from the store: `pooled_density` from the
  summed hits and paths, and `density` weighted by the number of observed voxels (`observed_fraction` times the
  voxels per cell). When the cell size is a multiple of the voxel width every cell has the same number of voxels
  (also with `--ground`), so weighting by `observed_fraction` is exact; otherwise the count varies with the ground
  height and isn't stored.

### 6. Validation

* Tiled versus single run: a 300 stems/m^2 synthetic scene (`raycreate pneumatophores`), 3 x 3 tiles of 2 m via
  `scripts/raydensity_tiles.sh`, gave cells identical to a single 6 m tile (`density`, `pooled_density`, `hits`,
  `weighted_path`, `observed_fraction`, `patchiness_ratio`). The test suite (`RayDensityTiles`) also splits with a
  halo equal to `raydensity`'s and checks agreement to 1e-4.
* Heights above ground: the scene sheared by a 10% x 5% plane and estimated with `--ground` matched the flat scene
  to within 0.5% per layer over the plot (4 scenes; `pooled_density` within 1.6%), and 3.5% per cell; absolute
  heights gave 0.72-0.95 of it.
* Still to do: the cell study scenes run through the tiled path, and a run at scale.

## Remaining work

* Run at scale: the split is a single-threaded streaming pass (I/O bound), and each tile is single threaded. A
  SLURM job array version of `scripts/raydensity_tiles.sh` is straightforward (one tile per array index).
* A CRS attribute, e.g. from a `--georeference` projection file.
* Optional zlib compression of chunks (the writer is uncompressed, which xarray reads directly).
* A consolidated metadata file (`.zmetadata`) to speed up opening stores with many arrays over a network.
