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
  std::cout << "                     --profile name.csv    - for density styles: output the density per height layer" << std::endl;
  std::cout << "                     --cells name.csv      - for density styles: output density aggregated into cells." << std::endl;
  std::cout << "                                             Both give the density pooled over the observed voxels" << std::endl;
  std::cout << "                                             (hits / path), the mean voxel density, and their ratio" << std::endl;
  std::cout << "                     --cell_size 1,1,0.1   - size of the aggregated cells (m), tiling from the min bound" << std::endl;
  std::cout << "                     --max_sin_zenith 1    - vertical_cylinder only: ignore rays shallower than this" << std::endl;
  std::cout << "                     --stem_diameter 0.01  - vertical_cylinder only: mean true stem diameter (m)," << std::endl;
  std::cout << "                     --footprint 0.01,0.003,0.05 - with beam exit diameter (m), full divergence (rad)" << std::endl;
  std::cout << "                                             and detection threshold (fraction of beam power), to" << std::endl;
  std::cout << "                                             correct for the beam footprint widening thin stems" << std::endl;
  std::cout << "                     --estimator crossing  - for density styles: count hit rays' path to the voxel" << std::endl;
  std::cout << "                                             exit (and through the stem). No occlusion correction, so" << std::endl;
  std::cout << "                                             biased low for dense stems with a beam footprint." << std::endl;
  std::cout << "                                             Default beer_lambert: path to hit (turbid medium)" << std::endl;
  std::cout << "                     --estimator pooled    - beer_lambert pooled over a neighbourhood of voxels, and" << std::endl;
  std::cout << "                                             corrected for the stems' (widened) cross-sections" << std::endl;
  std::cout << "                     --pool_radius 0.05    - nominal half-width (m) of the neighbourhood for --estimator" << std::endl;
  std::cout << "                                             pooled and the stem occlusion correction. The neighbourhood" << std::endl;
  std::cout << "                                             is (2k+1)^2 voxels for k = round(pool_radius/pixel_width)," << std::endl;
  std::cout << "                                             so +-(k+1/2) voxels: +-0.075 m for 0.05 m at 0.05 m voxels" << std::endl;
  std::cout << "                     --return_range outline - with --footprint: returns from partly intercepted stems are" << std::endl;
  std::cout << "                                             reported where the beam axis enters the widened stem" << std::endl;
  std::cout << "                                             outline. Default surface: the range of the stem surface" << std::endl;
  std::cout << "                     --min_gap 0           - minimum gap between stem surfaces (in stem diameters)" << std::endl;
  std::cout << "                                             assumed by the stem occlusion correction" << std::endl;
  std::cout << "                     --no_overlap_correction - with --stem_diameter: don't correct for stems hiding" << std::endl;
  std::cout << "                                             each other (footprint overlap and within-voxel occlusion)" << std::endl;
  std::cout << "                     --self_shadow_correction - with --estimator pooled and --stem_diameter: correct the" << std::endl;
  std::cout << "                                             neighbourhoods for the shadow each stem casts within its own" << std::endl;
  std::cout << "                                             neighbourhood (experimental; the mean density reads high without)" << std::endl;
  std::cout << "                     --self_shadow_second_order - with --self_shadow_correction: the second-order factor," << std::endl;
  std::cout << "                                             for small neighbourhoods (e.g. --pool_radius 0)" << std::endl;
  std::cout << "                     --occlusion_radius 0.05 - with --estimator pooled: invert the stem occlusion correction" << std::endl;
  std::cout << "                                             at the density of this (larger) neighbourhood, and apply its" << std::endl;
  std::cout << "                                             factor to the pooled one (default: the pooled neighbourhood)" << std::endl;
  std::cout << "                     --ratio_bias_correction - with --estimator pooled: correct each neighbourhood's" << std::endl;
  std::cout << "                                             hits/path for its finite number of rays (delta method)" << std::endl;
  std::cout << "                     --no_neighbour_priors - with --estimator pooled: don't blend voxels of under 10 rays" << std::endl;
  std::cout << "                                             with their neighbours (it biases small neighbourhoods low)" << std::endl;
  std::cout << "                                             Single-voxel neighbourhoods (experimental): --pool_radius 0" << std::endl;
  std::cout << "                                             --self_shadow_correction --self_shadow_second_order" << std::endl;
  std::cout << "                                             --occlusion_radius 0.05 --ratio_bias_correction" << std::endl;
  std::cout << "                                             --no_neighbour_priors" << std::endl;
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
  ray::DoubleArgument occlusion_radius(0.0, 100.0);
  ray::FileArgument cells_file(false);
  ray::IntArgument resolution(1, 20000, 512);
  ray::FileArgument cloud_file, image_file, transform_file, projection_file(false), profile_file(false);
  ray::OptionalFlagArgument mark_origin("mark_origin", 'm'), view_flag("view", 'v'), no_flat_top("no_flat_top", 'u'),
    no_overlap("no_overlap_correction", 'k'), self_shadow("self_shadow_correction", 'w'),
    second_order("self_shadow_second_order", 'S'), ratio_bias("ratio_bias_correction", 'R'),
    no_priors("no_neighbour_priors", 'N');
  ray::KeyChoice return_range({ "surface", "outline" });
  ray::OptionalKeyValueArgument return_range_option("return_range", 'j', &return_range);
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
  ray::OptionalKeyValueArgument occlusion_radius_option("occlusion_radius", 'O', &occlusion_radius);
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
          &no_overlap,           &pool_radius_option,    &min_gap_option,      &return_range_option,
          &self_shadow,          &second_order,          &occlusion_radius_option, &ratio_bias,
          &no_priors }))
  {
    usage();
  }
  // the density options only apply to the density styles, and some only to the vertical cylinder model
  const bool density_style = style.selectedKey() == "density" || style.selectedKey() == "density_rgb";
  const std::vector<std::pair<std::string, bool>> density_args = { { "density_model", density_model_option.isSet() },
                                                                   { "min_sin_zenith", min_sin_zenith_option.isSet() },
                                                                   { "max_sin_zenith", max_sin_zenith_option.isSet() },
                                                                   { "profile", profile_file_option.isSet() },
                                                                   { "cells", cells_option.isSet() },
                                                                   { "cell_size", cell_size_option.isSet() },
                                                                   { "stem_diameter", stem_diameter_option.isSet() },
                                                                   { "footprint", footprint_option.isSet() },
                                                                   { "estimator", estimator_option.isSet() },
                                                                   { "pool_radius", pool_radius_option.isSet() },
                                                                   { "min_gap", min_gap_option.isSet() },
                                                                   { "no_overlap_correction", no_overlap.isSet() },
                                                                   { "self_shadow_correction", self_shadow.isSet() },
                                                                   { "self_shadow_second_order", second_order.isSet() },
                                                                   { "occlusion_radius",
                                                                     occlusion_radius_option.isSet() },
                                                                   { "ratio_bias_correction", ratio_bias.isSet() },
                                                                   { "no_neighbour_priors", no_priors.isSet() },
                                                                   { "no_flat_top", no_flat_top.isSet() },
                                                                   { "return_range", return_range_option.isSet() } };
  const bool cylinders = density_model_option.isSet() && density_model.selectedKey() == "vertical_cylinder";
  for (const auto &arg : density_args)
  {
    if (arg.second && !density_style)
    {
      std::cerr << "Error: --" << arg.first << " only applies to the density and density_rgb styles" << std::endl;
      usage();
    }
    const bool cylinder_only = arg.first == "min_sin_zenith" || arg.first == "max_sin_zenith" ||
                               arg.first == "stem_diameter" || arg.first == "footprint" || arg.first == "min_gap" ||
                               arg.first == "no_overlap_correction" || arg.first == "return_range" ||
                               arg.first == "self_shadow_correction" || arg.first == "self_shadow_second_order" ||
                               arg.first == "occlusion_radius";
    if (arg.second && cylinder_only && !cylinders)
    {
      std::cerr << "Error: --" << arg.first << " requires --density_model vertical_cylinder" << std::endl;
      usage();
    }
  }
  if (cell_size_option.isSet() && !cells_option.isSet())
  {
    std::cerr << "Error: --cell_size requires --cells" << std::endl;
    usage();
  }
  if (return_range_option.isSet() && !footprint_option.isSet())
  {
    std::cerr << "Error: --return_range requires --footprint" << std::endl;
    usage();
  }
  ray::DensityOptions density_options;
  ray::DensityOutputs density_outputs;
  density_options.flat_top_compensation = !no_flat_top.isSet();
  density_options.overlap_correction = !no_overlap.isSet();
  density_options.self_shadow_correction = self_shadow.isSet();
  if (self_shadow.isSet() &&
      (!stem_diameter_option.isSet() || !estimator_option.isSet() || estimator.selectedKey() != "pooled"))
  {
    std::cerr << "Error: --self_shadow_correction requires --estimator pooled and --stem_diameter" << std::endl;
    usage();
  }
  if (second_order.isSet() && !self_shadow.isSet())
  {
    std::cerr << "Error: --self_shadow_second_order requires --self_shadow_correction" << std::endl;
    usage();
  }
  const bool pooled_estimator = estimator_option.isSet() && estimator.selectedKey() == "pooled";
  if ((occlusion_radius_option.isSet() || ratio_bias.isSet() || no_priors.isSet()) && !pooled_estimator)
  {
    std::cerr << "Error: --occlusion_radius, --ratio_bias_correction and --no_neighbour_priors require --estimator "
                 "pooled"
              << std::endl;
    usage();
  }
  density_options.self_shadow_second_order = second_order.isSet();
  density_options.occlusion_radius = occlusion_radius_option.isSet() ? occlusion_radius.value() : -1.0;
  density_options.ratio_bias_correction = ratio_bias.isSet();
  density_options.pooled_neighbour_priors = !no_priors.isSet();
  density_options.min_gap_ratio = min_gap.value();
  if (pool_radius_option.isSet())
  {
    density_options.pool_radius = pool_radius.value();
  }
  if (cells_option.isSet())
  {
    density_outputs.cells_file = cells_file.name();
    density_outputs.cell_size = cell_size.value();
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
    density_outputs.profile_file = profile_file.name();
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
    density_options.footprint.range_to_outline = return_range_option.isSet() && return_range.selectedKey() == "outline";
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
                        transform_file_option.isSet() ? &transform_file.name() : nullptr, density_options,
                        density_outputs))
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
