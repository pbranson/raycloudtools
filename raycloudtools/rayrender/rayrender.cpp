// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "raylib/raycloud.h"
#include "raylib/raycuboid.h"
#include "raylib/raylibconfig.h"
#include "raylib/rayparse.h"
#include "raylib/rayrenderer.h"

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Render a ray cloud as an image, from a specified viewpoint" << std::endl;
  std::cout << "usage (--view / -v to view results):" << std::endl;
  std::cout << "rayrender raycloudfile.ply top ends        - render from the top (plan view) the end points" << std::endl;
  std::cout << "                           left            - facing negative x axis" << std::endl;
  std::cout << "                           right           - facing positive x axis" << std::endl;
  std::cout << "                           front           - facing negative y axis" << std::endl;
  std::cout << "                           back            - facing positive y axis" << std::endl;
  std::cout << "                               mean        - mean colour on axis" << std::endl;
  std::cout << "                               sum         - sum colours (globally scaled to colour range)" << std::endl;
  std::cout << "                               starts      - render the ray start points" << std::endl;
  std::cout << "                               rays        - render the full set of rays" << std::endl;
  std::cout << "                               height      - render the maximum heights in the view axis" << std::endl;
  std::cout << "                               density     - shade according to estimated density within pixel" << std::endl;
  std::cout << "                               density_rgb - r->g->b colour by estimated density" << std::endl;
  std::cout << "                     --resolution 512      - long axis resolution" << std::endl;
  std::cout << "                     --output name.png     - optional output file name. " << std::endl;
  std::cout << "                                             Supports .png, .tga, .hdr, .jpg, .bmp" << std::endl;
  std::cout << "                     --mark_origin         - place a 255,0,255 pixel at the coordinate origin. " << std::endl;
  std::cout << "                     --output_transform    - generate a yaml file containing the" << std::endl;
  std::cout << "                                             transform from the raycloud to" << std::endl;
  std::cout << "                                             pixels. Only compatible with top" << std::endl;
  std::cout << "                                             view." << std::endl;
  std::cout << "                     --georeference name.proj- projection file name, to output (geo)tif file. " << std::endl;
  std::cout << "                     --pixel_width 0.1     - optional pixel width in m, instead of resolution" << std::endl;
  std::cout << "                     --grid_width 100      - optionally bound to a grid cell width such that one cell centre is 0,0" << std::endl;
  std::cout << "                     --density_model vertical_cylinder - for density styles: model elements as vertical" << std::endl;
  std::cout << "                                             cylinders (e.g. pneumatophores), giving frontal area per" << std::endl;
  std::cout << "                                             volume N*d (1/m). Default: spherical (leaf area per volume)" << std::endl;
  std::cout << "                     --min_sin_zenith 0.1  - vertical_cylinder only: ignore rays steeper than this" << std::endl;
  std::cout << "                     --profile name.csv    - for density styles: output mean density per height layer" << std::endl;
  std::cout << "                     --cells name.csv      - for density styles: output density aggregated into cells" << std::endl;
  std::cout << "                     --cell_size 1,1,0.1   - size of the aggregated cells (m), tiling from the min bound" << std::endl;
  std::cout << "                     --max_sin_zenith 1    - vertical_cylinder only: ignore rays shallower than this" << std::endl;
  std::cout << "                     --stem_diameter 0.01  - vertical_cylinder only: mean true stem diameter (m)," << std::endl;
  std::cout << "                     --footprint 0.01,0.003,0.05 - with beam exit diameter (m), full divergence (rad)" << std::endl;
  std::cout << "                                             and detection threshold (fraction of beam power), to" << std::endl;
  std::cout << "                                             correct for the beam footprint widening thin stems" << std::endl;
  std::cout << "                     --estimator crossing  - for density styles: count hit rays' path to the voxel" << std::endl;
  std::cout << "                                             exit. Exact for few, discrete elements per voxel (small" << std::endl;
  std::cout << "                                             voxels). Default beer_lambert: path to hit (turbid medium)" << std::endl;
  std::cout << "                     --estimator pooled    - beer_lambert pooled over a neighbourhood of voxels, and" << std::endl;
  std::cout << "                                             corrected for the stems' (widened) cross-sections" << std::endl;
  std::cout << "                     --pool_radius 0.05    - half-width (m) of the neighbourhood for --estimator pooled" << std::endl;
  std::cout << "                                             and the stem occlusion correction" << std::endl;
  std::cout << "                     --min_gap 0           - minimum gap between stem surfaces (in stem diameters)" << std::endl;
  std::cout << "                                             assumed by the stem occlusion correction" << std::endl;
  std::cout << "                     --no_overlap_correction - with --stem_diameter: don't correct for stems hiding" << std::endl;
  std::cout << "                                             each other (footprint overlap and within-voxel occlusion)" << std::endl;
  std::cout << "                     --no_flat_top         - for density styles: disable the flat canopy top model" << std::endl;
  std::cout << "                                             (designed for grass/crops, it biases uneven stems)" << std::endl;
  std::cout << "                     --min_bound x,y,z     - optional bounds of the region to render, otherwise" << std::endl;
  std::cout << "                     --max_bound x,y,z       the bounds of the ray end points" << std::endl;
  std::cout << "Default output is raycloudfile.png" << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayRender(int argc, char *argv[])
{
  ray::KeyChoice viewpoint({ "top", "left", "right", "front", "back" });
  ray::KeyChoice style({ "ends", "mean", "sum", "starts", "rays", "height", "density", "density_rgb" });
  ray::DoubleArgument pixel_width(0.0001, 1000.0), grid_width(0.01, 1000000.0), min_sin_zenith(0.0, 1.0);
  ray::KeyChoice density_model({ "spherical", "vertical_cylinder" });
  ray::KeyChoice estimator({ "beer_lambert", "crossing", "pooled" });
  ray::DoubleArgument stem_diameter(0.0001, 10.0);
  ray::Vector3dArgument footprint(0.0, 1.0), min_bound(-1e10, 1e10), max_bound(-1e10, 1e10);
  ray::Vector3dArgument cell_size(0.0001, 1e6, Eigen::Vector3d(1.0, 1.0, 0.1));
  ray::DoubleArgument max_sin_zenith(0.0, 1.0), pool_radius(0.0, 100.0), min_gap(0.0, 100.0);
  ray::FileArgument cells_file(false);
  ray::IntArgument resolution(1, 20000, 512);
  ray::FileArgument cloud_file, image_file, transform_file, projection_file(false), profile_file(false);
  ray::OptionalFlagArgument mark_origin("mark_origin", 'm'), view_flag("view", 'v'), no_flat_top("no_flat_top", 'u'),
    no_overlap("no_overlap_correction", 'k');
  ray::OptionalKeyValueArgument resolution_option("resolution", 'r', &resolution);
  ray::OptionalKeyValueArgument pixel_width_option("pixel_width", 'p', &pixel_width);
  ray::OptionalKeyValueArgument grid_width_option("grid_width", 'g', &grid_width);
  ray::OptionalKeyValueArgument output_file_option("output", 'o', &image_file);
  ray::OptionalKeyValueArgument projection_file_option("georeference", 'g', &projection_file);
  ray::OptionalKeyValueArgument transform_file_option("output_transform", 't', &transform_file);
  ray::OptionalKeyValueArgument density_model_option("density_model", 'd', &density_model);
  ray::OptionalKeyValueArgument min_sin_zenith_option("min_sin_zenith", 'z', &min_sin_zenith);
  ray::OptionalKeyValueArgument profile_file_option("profile", 'f', &profile_file);
  ray::OptionalKeyValueArgument estimator_option("estimator", 'e', &estimator);
  ray::OptionalKeyValueArgument cells_option("cells", 'c', &cells_file);
  ray::OptionalKeyValueArgument cell_size_option("cell_size", 'l', &cell_size);
  ray::OptionalKeyValueArgument max_sin_zenith_option("max_sin_zenith", 'y', &max_sin_zenith);
  ray::OptionalKeyValueArgument pool_radius_option("pool_radius", 'q', &pool_radius);
  ray::OptionalKeyValueArgument min_gap_option("min_gap", 'a', &min_gap);
  ray::OptionalKeyValueArgument stem_diameter_option("stem_diameter", 's', &stem_diameter);
  ray::OptionalKeyValueArgument footprint_option("footprint", 'b', &footprint);
  ray::OptionalKeyValueArgument min_bound_option("min_bound", 'n', &min_bound);
  ray::OptionalKeyValueArgument max_bound_option("max_bound", 'x', &max_bound);
  if (!ray::parseCommandLine(
        argc, argv, { &cloud_file, &viewpoint, &style },
        { &resolution_option,    &pixel_width_option,    &output_file_option,  &view_flag,
          &mark_origin,          &transform_file_option, &grid_width_option,   &projection_file_option,
          &density_model_option, &min_sin_zenith_option, &profile_file_option, &stem_diameter_option,
          &footprint_option,     &min_bound_option,      &max_bound_option,    &no_flat_top,
          &estimator_option,     &cells_option,          &cell_size_option,    &max_sin_zenith_option,
          &no_overlap,           &pool_radius_option,    &min_gap_option }))
  {
    usage();
  }
  ray::DensityOptions density_options;
  density_options.flat_top_compensation = !no_flat_top.isSet();
  density_options.overlap_correction = !no_overlap.isSet();
  density_options.min_gap_ratio = min_gap.value();
  if (pool_radius_option.isSet())
  {
    density_options.pool_radius = pool_radius.value();
  }
  if (cells_option.isSet())
  {
    density_options.cells_file = cells_file.name();
    density_options.cell_size = cell_size.value();
  }
  if (max_sin_zenith_option.isSet())
  {
    density_options.max_sin_zenith = max_sin_zenith.value();
  }
  if (estimator_option.isSet())
  {
    density_options.path_estimator = static_cast<ray::PathEstimator>(estimator.selectedID());
  }
  if (density_model_option.isSet())
  {
    // the text and enums are in the same order
    density_options.distribution = static_cast<ray::ElementDistribution>(density_model.selectedID());
  }
  if (min_sin_zenith_option.isSet())
  {
    density_options.min_sin_zenith = min_sin_zenith.value();
  }
  if (profile_file_option.isSet())
  {
    density_options.profile_file = profile_file.name();
  }
  if (stem_diameter_option.isSet())
  {
    density_options.stem_diameter = stem_diameter.value();
  }
  if (footprint_option.isSet())
  {
    density_options.footprint.exit_diameter = footprint.value()[0];
    density_options.footprint.divergence = footprint.value()[1];
    density_options.footprint.threshold = footprint.value()[2];
    if (density_options.distribution != ray::ElementDistribution::VerticalCylinder || !stem_diameter_option.isSet())
    {
      std::cerr << "Error: --footprint requires --density_model vertical_cylinder and --stem_diameter" << std::endl;
      usage();
    }
  }
  if (!output_file_option.isSet())
  {
    image_file.name() = cloud_file.nameStub() + (projection_file_option.isSet() ? ".tif" : ".png");
  }
  // a projection file describes where the ray cloud is in the world, which allows
  // images to be output in geotiff (geolocalised tiff) format.
  if (projection_file_option.isSet())
  {
#if !RAYLIB_WITH_TIFF
    std::cerr << "Error: georeferencing requires the RAYLIB_WITH_TIFF build flag enabled. See README.md." << std::endl;
    usage();
#endif
    if (image_file.nameExt() != "tif")
    {
      std::cerr << "Error: projection files can only be used when outputting a .tif file" << std::endl;
      usage();
    }
    if (viewpoint.selectedKey() != "top")
    {
      std::cerr << "Error: can only geolocate a top-down render" << std::endl;
      usage();
    }
  }

  ray::Cloud::Info info;
  if (!ray::Cloud::getInfo(cloud_file.name(), info))
  {
    usage();
  }
  ray::Cuboid bounds = info.ends_bound;  // exclude the unbounded ray lengths (e.g. up into the sky)
  if (min_bound_option.isSet())
  {
    bounds.min_bound_ = min_bound.value();
  }
  if (max_bound_option.isSet())
  {
    bounds.max_bound_ = max_bound.value();
  }
  if ((bounds.max_bound_.array() <= bounds.min_bound_.array()).any())
  {
    std::cerr << "Error: max_bound must be greater than min_bound in each axis" << std::endl;
    usage();
  }
  if (grid_width_option.isSet())  // adjust min_bound to be well-aligned
  {
    Eigen::Vector3d mid = (bounds.min_bound_ + bounds.max_bound_) / 2.0;
    Eigen::Vector3d min_bound = bounds.min_bound_;
    Eigen::Vector3d max_bound = bounds.max_bound_;
    min_bound[0] = grid_width.value() * std::round(mid[0] / grid_width.value()) - 0.5 * grid_width.value();
    min_bound[1] = grid_width.value() * std::round(mid[1] / grid_width.value()) - 0.5 * grid_width.value();
    max_bound[0] = min_bound[0] + grid_width.value();
    max_bound[1] = min_bound[1] + grid_width.value();
    if (min_bound[0] > bounds.min_bound_[0] || min_bound[1] > bounds.min_bound_[1] ||
        max_bound[0] < bounds.max_bound_[0] || max_bound[1] < bounds.max_bound_[1])
    {
      std::cout << "Warning: cloud overlaps grid cell of width: " << grid_width.value() << " image bounds extended"
                << std::endl;
    }
    bounds.min_bound_ = ray::minVector(bounds.min_bound_, min_bound);
    bounds.max_bound_ = ray::maxVector(bounds.max_bound_, max_bound);
  }
  double pix_width = pixel_width.value();
  if (!pixel_width_option.isSet())
  {
    Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
    double length = std::max(extent[0], extent[1]);
    pix_width = length / resolution.value();
  }
  if (pix_width <= 0.0)
  {
    usage();
  }

  // quick casting allowed, taking care that the text and enums are in the same order
  const ray::ViewDirection view_dir = static_cast<ray::ViewDirection>(viewpoint.selectedID());
  const ray::RenderStyle render_style = static_cast<ray::RenderStyle>(style.selectedID());

  // an option to output the transformation from image to world frame
  if (transform_file_option.isSet() && (view_dir != ray::ViewDirection::Top))
  {
    std::cout << "--output_transform can only be used when view is top." << std::endl;
    usage();
  }

  if (!ray::renderCloud(cloud_file.name(), bounds, view_dir, render_style, pix_width, image_file.name(),
                        projection_file.name(), mark_origin.isSet(),
                        transform_file_option.isSet() ? &transform_file.name() : nullptr, density_options))
  {
    usage();
  }
  if (view_flag.isSet())
  {
    const std::string command = std::string(RAYLIB_IMAGETOOL) + " " + image_file.name();
    if (std::system(command.c_str()) != 0)
    {
      std::cerr << "failed to view file with command: " << command << std::endl;
    }
  }
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayRender, argc, argv);
}
