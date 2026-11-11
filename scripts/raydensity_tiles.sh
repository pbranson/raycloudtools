#!/bin/bash
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Volumetric stem density over a large domain: splits the ray cloud into tiles (each with every ray crossing it,
# plus a halo), then runs raydensity on each tile in parallel, writing one chunk per tile of a single Zarr store.
#
# usage: raydensity_tiles.sh cloud.ply store.zarr x0,y0 tile_width nx,ny [jobs] [-- raydensity options]
#   x0,y0       domain origin (m). Tile (i, j) covers [x0 + i w, x0 + (i+1) w) and similarly in y
#   tile_width  tile width (m), a multiple of the voxel and cell widths
#   nx,ny       number of tiles in the domain
#   jobs        number of tiles to process at once (default: number of processors)
# Options after -- are passed to raydensity, e.g.
#   -- --ground ground_mesh.ply --heights 0,0.6 --stem_diameter 0.01 --footprint 0.01,0.003,0.05
# Environment variables: HALO (default 0.15, which suits the default 0.05 m voxels and pool radius), and HEIGHTS
# (the absolute height range of rays to keep, default -10,50: set it to just cover the ground and stems).
#
# For example, a 300 x 300 m forest in 25 m tiles:
#   raydensity_tiles.sh forest.ply forest.zarr 400000,6100000 25 12,12 32 -- --ground forest_mesh.ply --stem_diameter 0.01
set -e
if [ $# -lt 5 ]; then
  sed -n '6,22p' "$0"
  exit 1
fi
cloud=$1
store=$2
origin=$3
width=$4
tiles=$5
shift 5
jobs=$(nproc 2>/dev/null || echo 4)
if [ $# -gt 0 ] && [ "$1" != "--" ]; then
  jobs=$1
  shift
fi
if [ "$1" == "--" ]; then
  shift
fi
halo=${HALO:-0.15}
heights=${HEIGHTS:--10,50}
nx=${tiles%,*}
ny=${tiles#*,}
stub=${cloud%.ply}
workdir=$(mktemp -d "${stub}_tiles.XXXXXX")

# 1. split into tiles, keeping ray start points, within the domain only
# (a link, not a copy, names the tile files tile_i_j.ply without duplicating a large cloud)
ln -s "$(cd "$(dirname "$cloud")" && pwd)/$(basename "$cloud")" "$workdir/tile.ply"
(cd "$workdir" && raysplit tile.ply tiles "$origin" "$width" "$halo" "$heights" "$tiles")
rm "$workdir/tile.ply"

# 2. estimate each tile in parallel, writing its chunk of the store
export RD_STORE=$(cd "$(dirname "$store")" && pwd)/$(basename "$store")
export RD_DIR=$workdir RD_ORIGIN=$origin RD_WIDTH=$width RD_TILES=$tiles RD_HALO=$halo
export RD_OPTIONS="$*"
run_tile() {
  local i=$1 j=$2
  # shellcheck disable=SC2086
  if ! raydensity "$RD_DIR/tile_${i}_${j}.ply" "$RD_STORE" --origin "$RD_ORIGIN" --tile_width "$RD_WIDTH" \
    --tiles "$RD_TILES" --tile "$i,$j" --halo "$RD_HALO" $RD_OPTIONS > "$RD_DIR/tile_${i}_${j}.log" 2>&1; then
    echo "tile $i,$j failed: see $RD_DIR/tile_${i}_${j}.log"
  fi
}
export -f run_tile
for ((i = 0; i < nx; i++)); do
  for ((j = 0; j < ny; j++)); do
    if [ -f "$workdir/tile_${i}_${j}.ply" ]; then
      echo "$i $j"
    fi
  done
done | xargs -P "$jobs" -L 1 bash -c 'run_tile $0 $1'
echo "tiles written to $RD_STORE (tile clouds and logs in $workdir)"
