# Volumetric stem density for large forests: domain decomposition and output

Design note for running the vertical-cylinder density workflow over forests of order 300 x 300 m scanned with many
billions of rays, producing a single 3D gridded output for hydrodynamic models.

**Status:** implemented as `raysplit ... tiles`, `raydensity` and `scripts/raydensity_tiles.sh`. Tested on synthetic
scenes (tiled and single runs agree to float rounding; ground-relative layers on a 10% slope match flat ground to
within 1%), but not yet at the scale of billions of rays.

```console
raydensity_tiles.sh forest.ply forest.zarr 400000,6100000 25 12,12 32 -- --ground forest_mesh.ply \
  --heights 0,0.6 --stem_diameter 0.01 --footprint 0.01,0.003,0.05
```

```python
import xarray as xr
cells = xr.open_zarr("forest.zarr", consolidated=False)   # density, solid_fraction, ... on (z, y, x)
```

## Why the current interface doesn't scale

| Issue | Detail |
|---|---|
| One dense in-memory grid | `DensityGrid` holds 12 bytes per voxel (+4 for the occlusion correction, +4 temporary). 300 x 300 x 0.6 m at 5 cm is 6000 x 6000 x 12 = 4.3e8 voxels, about 9 GB. At 2.5 cm it is 8x more. |
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
* Each tile is processed with its bounds expanded by a halo h >= pool_radius + w (the pooled neighbourhood plus the
  one-voxel neighbour prior), e.g. 0.1 m. Only the tile's core is output. With the halo, a tile's core voxels
  receive exactly the same rays, in the same order, as in a single global run, so the result is identical: tiles
  have no seams and need no blending.
* A 25 m tile at 5 cm is 500 x 500 x (relief + 0.6 m)/w voxels: tens of MB, so many tiles run concurrently.

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
* The pass is I/O bound: about 36 bytes per ray in .ply, so 1e10 rays is about 360 GB read once. When there are
  more tiles than open files, `raysplit` already makes several passes; a two-level split (100 m blocks, then
  tiles) keeps it to two.
* SLAM data is time ordered, so consecutive rays are spatially local. Writing per-tile buffers of a few MB keeps
  the writes sequential.

### 3. Height above ground

* Don't flatten the cloud with `raytranslate subtract ground_mesh.ply`: it shifts each **whole ray** by the ground
  height at its end point, so the rest of the ray is misplaced vertically by the slope times its horizontal
  distance from its end. Hand-held rays cross several metres of the stem layer, so on a 10% slope this put the
  rays' paths tens of centimetres out, and underestimated N d by 30-50% in simulation (only this slope has been tested).
* Instead (as implemented), keep absolute heights in the voxel grid (each tile covers its local relief), and assign
  each voxel to a height-above-ground layer using the ground height beneath its centre (`--ground`, from the mesh at
  voxel resolution). The ground mesh comes from `rayextract terrain` (itself run per tile, as in
  `scripts/rayextract_trees_large.sh`).
* Ground returns must be removed (zero alpha) before density estimation, using the ground model, as now.

### 4. Per-tile density (embarrassingly parallel)

A dedicated tool, `raydensity`, rather than more options on `rayrender`, which is an image renderer:

```
raydensity tile_12_7.ply forest.zarr --origin x0,y0 --tile_width 25 --tiles 12,12 --tile 12,7 --halo 0.15
    --voxel 0.05 --cell 1,1,0.1 --heights 0,0.6 --ground ground_mesh.ply
    --stem_diameter 0.01 --footprint 0.01,0.003,0.05 [--voxels]
```

It runs the existing `DensityGrid` pipeline (calculateDensities, applyOverlapCorrection, addNeighbourPriors) on
the tile plus halo, aggregates the core voxels into cells by height above ground, and writes its chunk of the
output store. Tiles are independent: run with a job array (SLURM), `xargs -P`, or threads.

### 5. Output: a chunked store with one chunk per tile

* Use Zarr (v2): a directory holding JSON metadata (`.zarray`, `.zattrs`) and one binary file per chunk. Without
  compression, a chunk is just raw little-endian float32 data, so it is writable from C++ with no dependencies
  (zlib compression is a small, optional step). It is read directly by xarray, zarr-python, GDAL and many
  visualisation tools, and converts to NetCDF or Cloud Optimised GeoTIFF.
* Chunk = one tile's core, so each tile writes only its own chunk files, and there is no merge step and no
  locking. The metadata is written once, up front, from the global extent.
* Variables (as implemented), on dimensions (z, y, x) for 1 x 1 x 0.1 m cells (300 x 300 x 6 per variable, small):
  `density` (N d, mean of voxel estimates), `single_voxel_density` (total hits over total weighted path, as if
  the cell were one voxel; its ratio to `density` indicates patchiness), `solid_fraction` (N pi d^2 / 4), `hits`,
  `weighted_path` and `observed_fraction`, plus coordinates `x`, `y`, `z` (cell centres). With `--voxels`, also
  `voxel_density` on (zv, yv, xv) as heights above ground (6000 x 6000 x 12 float32, 1.7 GB, for 3D inspection).
  Not yet output: the mean d_eff / d, and a flag for cells at the solid fraction limit.
* Attributes: CRS and origin (as for `--georeference`), voxel and cell sizes, stem diameter, footprint, estimator,
  software version and the command line, so results are reproducible.
* Because the cell values are means of voxel values, and the sufficient statistics (sums and counts) are also
  stored, any coarser aggregation (e.g. 5 m cells for a larger model domain) can be done exactly from the store.

### 6. Validation

* Tiled versus single-run equality on a synthetic scene (`raycreate pneumatophores`), which the halo rule should
  make exact.
* The existing synthetic scenes and `scripts/pneumatophore_cell_study.py`, run through the tiled path.

## Remaining work

* Run at scale: the split is a single-threaded streaming pass (I/O bound), and each tile is single threaded. A
  SLURM job array version of `scripts/raydensity_tiles.sh` is straightforward (one tile per array index).
* Optional zlib compression of chunks (the writer is uncompressed, which xarray reads directly).
* A consolidated metadata file (`.zmetadata`) to speed up opening stores with many arrays over a network.
