// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#include "raylib/raydensity.h"
#include <cstdlib>
#include <iostream>
#include <sstream>
#include "raylib/raylibconfig.h"
#include "raylib/rayparse.h"

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Estimate the volumetric density of vertical stems (or foliage) in one tile of a domain, writing the" << std::endl;
  std::cout << "tile's chunk of a Zarr store of cells (and optionally voxels). Tiles are independent, so can be run in parallel." << std::endl;
  std::cout << "usage:" << std::endl;
  std::cout << "raydensity cloud_3_7.ply forest.zarr --tile 3,7 --tiles 12,12 --origin 0,0 --tile_width 25" << std::endl;
  std::cout << "                     --voxel 0.05          - voxel width (m)" << std::endl;
  std::cout << "                     --cell 1,1,0.1        - output cell size (m). The tile width must be a multiple" << std::endl;
  std::cout << "                     --heights 0,0.6       - range of heights to output, above the ground mesh if given" << std::endl;
  std::cout << "                     --ground mesh.ply     - ground mesh (e.g. from rayextract terrain), for heights above ground" << std::endl;
  std::cout << "                     --halo 0.15           - rays within this distance of the tile are used (default" << std::endl;
  std::cout << "                                             pool_radius + 2 voxels). Use the same halo in raysplit tiles" << std::endl;
  std::cout << "                     --voxels              - also output the voxel densities" << std::endl;
  std::cout << "                     --density_model spherical - leaf area per volume (default vertical_cylinder: N d)" << std::endl;
  std::cout << "                     --estimator pooled    - also beer_lambert or crossing (see rayrender)" << std::endl;
  std::cout << "                     --stem_diameter 0.01  - mean stem diameter (m), for the occlusion correction" << std::endl;
  std::cout << "                     --footprint 0.01,0.003,0.05 - beam exit diameter (m), divergence (rad), threshold" << std::endl;
  std::cout << "                     --pool_radius 0.05, --min_gap 0, --min_sin_zenith 0.1, --max_sin_zenith 1," << std::endl;
  std::cout << "                     --no_overlap_correction, --flat_top - as for rayrender" << std::endl;
  std::cout << "The input should contain every ray crossing the tile and its halo, keeping ray start points, as given by:" << std::endl;
  std::cout << "raysplit cloud.ply tiles 0,0 25 0.15 -1,2" << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayDensity(int argc, char *argv[])
{
  ray::FileArgument cloud_file, store_file(false);
  ray::Vector2dArgument tile(0.0, 1e6), tiles(1.0, 1e6, Eigen::Vector2d(1, 1)), origin,
    heights(-1e6, 1e6, Eigen::Vector2d(0.0, 0.6));
  ray::DoubleArgument tile_width(0.001, 1e6, 25.0), voxel(0.0001, 100.0, 0.05), halo(0.0, 1e6),
    stem_diameter(0.0001, 10.0), min_sin(0.0, 1.0, 0.1), max_sin(0.0, 1.0, 1.0), pool_radius(0.0, 100.0, 0.05),
    min_gap(0.0, 100.0);
  ray::Vector3dArgument cell(0.0001, 1e6, Eigen::Vector3d(1.0, 1.0, 0.1)), footprint(0.0, 1.0);
  ray::FileArgument ground_file;
  ray::KeyChoice density_model({ "spherical", "vertical_cylinder" }),
    estimator({ "beer_lambert", "crossing", "pooled" });
  ray::OptionalKeyValueArgument tile_option("tile", 't', &tile), tiles_option("tiles", 'n', &tiles),
    origin_option("origin", 'o', &origin), width_option("tile_width", 'w', &tile_width),
    voxel_option("voxel", 'v', &voxel), cell_option("cell", 'c', &cell), heights_option("heights", 'h', &heights),
    ground_option("ground", 'g', &ground_file), halo_option("halo", 'a', &halo),
    model_option("density_model", 'd', &density_model), estimator_option("estimator", 'e', &estimator),
    stem_option("stem_diameter", 's', &stem_diameter), footprint_option("footprint", 'b', &footprint),
    min_sin_option("min_sin_zenith", 'z', &min_sin), max_sin_option("max_sin_zenith", 'y', &max_sin),
    pool_option("pool_radius", 'q', &pool_radius), gap_option("min_gap", 'm', &min_gap);
  ray::OptionalFlagArgument voxels_flag("voxels", 'x'), no_overlap("no_overlap_correction", 'k'),
    flat_top("flat_top", 'f');
  if (!ray::parseCommandLine(argc, argv, { &cloud_file, &store_file },
                             { &tile_option,      &tiles_option,   &origin_option,    &width_option,   &voxel_option,
                               &cell_option,      &heights_option, &ground_option,    &halo_option,    &model_option,
                               &estimator_option, &stem_option,    &footprint_option, &min_sin_option, &max_sin_option,
                               &pool_option,      &gap_option,     &voxels_flag,      &no_overlap,     &flat_top }))
  {
    usage();
  }
  ray::DensityTileParams params;
  params.origin = origin.value();
  params.tile_width = tile_width.value();
  params.tile =
    Eigen::Vector2i(static_cast<int>(std::round(tile.value()[0])), static_cast<int>(std::round(tile.value()[1])));
  params.num_tiles =
    Eigen::Vector2i(static_cast<int>(std::round(tiles.value()[0])), static_cast<int>(std::round(tiles.value()[1])));
  params.voxel_width = voxel.value();
  params.cell_size = cell.value();
  params.min_height = heights.value()[0];
  params.max_height = heights.value()[1];
  if (ground_option.isSet())
  {
    params.ground_mesh_file = ground_file.name();
  }
  if (halo_option.isSet())
  {
    params.halo = halo.value();
  }
  params.write_voxels = voxels_flag.isSet();
  ray::DensityOptions &options = params.options;
  options.distribution = ray::ElementDistribution::VerticalCylinder;
  if (model_option.isSet())
  {
    options.distribution = static_cast<ray::ElementDistribution>(density_model.selectedID());
  }
  options.path_estimator = ray::PathEstimator::PooledBeerLambert;
  if (estimator_option.isSet())
  {
    options.path_estimator = static_cast<ray::PathEstimator>(estimator.selectedID());
  }
  options.flat_top_compensation = flat_top.isSet();
  options.overlap_correction = !no_overlap.isSet();
  options.min_sin_zenith = min_sin.value();
  options.max_sin_zenith = max_sin.value();
  options.pool_radius = pool_radius.value();
  options.min_gap_ratio = min_gap.value();
  if (stem_option.isSet())
  {
    options.stem_diameter = stem_diameter.value();
  }
  if (footprint_option.isSet())
  {
    if (!stem_option.isSet())
    {
      std::cerr << "Error: --footprint requires --stem_diameter" << std::endl;
      usage();
    }
    options.footprint.exit_diameter = footprint.value()[0];
    options.footprint.divergence = footprint.value()[1];
    options.footprint.threshold = footprint.value()[2];
  }
  std::stringstream description;
  for (int i = 0; i < argc; i++)
  {
    description << (i > 0 ? " " : "") << argv[i];
  }
  params.description = description.str();
  return ray::estimateDensityTile(cloud_file.name(), params, store_file.name()) ? 0 : 1;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayDensity, argc, argv);
}
