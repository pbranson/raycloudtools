#!/bin/bash
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
#
# Volumetric stem density over a large domain: splits the ray cloud into tiles (each with every ray crossing it,
# plus a halo), then runs raydensity on each tile in parallel, writing one chunk per tile of a single Zarr store.
#
# usage: raydensity_tiles.sh cloud.ply store.zarr x0,y0 tile_width nx,ny zmin,zmax [jobs] [-- raydensity options]
#   x0,y0       domain origin (m). Tile (i, j) covers [x0 + i w, x0 + (i+1) w) and similarly in y
#   tile_width  tile width (m), a multiple of the voxel and cell widths
#   nx,ny       number of tiles in the domain
#   zmin,zmax   absolute height range of the rays to keep (m): just cover the ground and the stems, e.g. from the
#               lowest ground minus 0.2 m to the highest ground plus the top of --heights plus 0.2 m. A tall
#               range keeps the canopy and sky rays, which multiplies the size of the tile files
#   jobs        number of tiles to process at once (default: number of processors)
# Options after -- are passed to raydensity, e.g.
#   -- --ground ground_mesh.ply --heights 0,0.6 --stem_diameter 0.01 --footprint 0.01,0.003,0.05
# The halo is the one raydensity uses (--halo if given, otherwise --pool_radius + 2 --voxel, by default 0.15 m).
# Exits non-zero if the split or any tile fails. Tile clouds and logs are kept in cloud_tiles.XXXXXX.
#
# For example, a 300 x 300 m forest in 25 m tiles, with ground between 1.2 and 3.5 m:
#   raydensity_tiles.sh forest.ply forest.zarr 400000,6100000 25 12,12 1,4.3 32 -- --ground forest_mesh.ply \
#     --stem_diameter 0.01
set -e
if [ $# -lt 6 ]; then
  awk 'NR >= 6 { if (!/^#/) exit; sub(/^# ?/, ""); print }' "$0"
  exit 1
fi
cloud=$1
store=$2
origin=$3
width=$4
tiles=$5
heights=$6
shift 6
jobs=$(nproc 2>/dev/null || echo 4)
if [ $# -gt 0 ] && [ "$1" != "--" ]; then
  jobs=$1
  shift
fi
if [ $# -gt 0 ] && [ "$1" == "--" ]; then
  shift
fi
case $heights in
  *,*) ;;
  *) echo "Error: zmin,zmax expected, not $heights" >&2; exit 1 ;;
esac

# the halo, computed as raydensity does, from its options or their defaults
voxel=0.05
pool_radius=0.05
halo=
args=("$@")
for ((k = 0; k < ${#args[@]}; k++)); do
  case ${args[k]} in
    --voxel | -v) voxel=${args[k + 1]} ;;
    --pool_radius | -q) pool_radius=${args[k + 1]} ;;
    --halo | -a) halo=${args[k + 1]} ;;
  esac
done
min_halo=$(awk -v p="$pool_radius" -v w="$voxel" 'BEGIN { printf "%.10g", p + 2 * w }')
if [ -z "$halo" ]; then
  halo=$min_halo
elif awk -v h="$halo" -v m="$min_halo" -v w="$voxel" 'BEGIN { exit !(h < m - 1e-6 * w) }'; then
  echo "Error: --halo $halo is less than --pool_radius + 2 --voxel ($min_halo)" >&2
  exit 1
fi
# raysplit's halo must be at least raydensity's. It is made to cover raydensity's grid, which extends the tile by the
# halo rounded up to whole voxels plus one voxel (for the neighbour prior), with the tile edges (where split rays end)
# off the voxel faces, so that the split's clipped ray ends are outside the grid
split_halo=$(awk -v h="$halo" -v w="$voxel" \
  'BEGIN { n = int(h / w); if (n * w < h - 1e-9 * w) n++; printf "%.10g", (n + 1.5) * w }')
echo "halo $halo m (split with $split_halo m)"
nx=${tiles%,*}
ny=${tiles#*,}
stub=${cloud%.ply}
workdir=$(mktemp -d "${stub}_tiles.XXXXXX")

# 1. split into tiles, keeping ray start points, within the domain only
# (a link, not a copy, names the tile files tile_i_j.ply without duplicating a large cloud)
cloud_dir=$(cd "$(dirname "$cloud")" && pwd)
ln -s "$cloud_dir/$(basename "$cloud")" "$workdir/tile.ply"
(cd "$workdir" && raysplit tile.ply tiles "$origin" "$width" "$split_halo" "$heights" "$tiles")
rm "$workdir/tile.ply"

# 2. estimate each tile in parallel, writing its chunk of the store
mkdir -p "$(dirname "$store")"
store_dir=$(cd "$(dirname "$store")" && pwd)
export RD_STORE=$store_dir/$(basename "$store")
export RD_DIR=$workdir RD_ORIGIN=$origin RD_WIDTH=$width RD_TILES=$tiles
# the options, quoted so that they can be re-read as an array by each job
RD_OPTIONS=$(printf '%q ' "$@")
export RD_OPTIONS
run_tile() {
  local i=$1 j=$2
  local -a options
  eval "options=($RD_OPTIONS)"
  # (raydensity uses the same halo: --halo if given, otherwise the same default)
  if ! raydensity "$RD_DIR/tile_${i}_${j}.ply" "$RD_STORE" --origin "$RD_ORIGIN" --tile_width "$RD_WIDTH" \
    --tiles "$RD_TILES" --tile "$i,$j" "${options[@]}" > "$RD_DIR/tile_${i}_${j}.log" 2>&1; then
    echo "tile $i,$j failed: see $RD_DIR/tile_${i}_${j}.log" >&2
    return 1
  fi
}
export -f run_tile
num_tiles=0
for ((i = 0; i < nx; i++)); do
  for ((j = 0; j < ny; j++)); do
    if [ -f "$workdir/tile_${i}_${j}.ply" ]; then
      num_tiles=$((num_tiles + 1))
    fi
  done
done
echo "estimating $num_tiles tiles (of $((nx * ny)) in the domain; tiles without rays are skipped)"
status=0
for ((i = 0; i < nx; i++)); do
  for ((j = 0; j < ny; j++)); do
    if [ -f "$workdir/tile_${i}_${j}.ply" ]; then
      echo "$i $j"
    fi
  done
done | xargs -P "$jobs" -L 1 bash -c 'run_tile "$0" "$1"' || status=$?
if [ $status -ne 0 ]; then
  echo "Error: some tiles failed (see above); the store $RD_STORE is incomplete" >&2
  exit 1
fi
echo "tiles written to $RD_STORE (tile clouds and logs in $workdir)"
