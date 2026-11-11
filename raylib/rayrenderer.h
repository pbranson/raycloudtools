// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYRENDERER_H
#define RAYLIB_RAYRENDERER_H

#include "raycuboid.h"
#include "raydensitygrid.h"
#include "raypose.h"
#include "rayutils.h"

namespace ray
{
/// Supported view directions on cloud data
enum class RAYLIB_EXPORT ViewDirection
{
  Top,
  Left,
  Right,
  Front,
  Back
};

/// Supported render styles on cloud data
enum class RAYLIB_EXPORT RenderStyle
{
  Ends,
  Mean,
  Sum,
  Starts,
  Rays,
  Height,
  Density,
  Density_rgb
};

/// Optional tables written by the density render styles, in addition to the image
struct RAYLIB_EXPORT DensityOutputs
{
  /// Optional .csv file to write the density of each horizontal layer of voxels (see aggregateDensityCells)
  std::string profile_file;
  /// Optional .csv file to write the density aggregated into cells of size @c cell_size (e.g. for hydrodynamic
  /// models), tiling from the minimum bound (see aggregateDensityCells)
  std::string cells_file;
  Eigen::Vector3d cell_size = Eigen::Vector3d(1.0, 1.0, 0.1);
  /// The cells (and layers of the profile) whose observed fraction is below this have nan density, pooled density
  /// and patchiness ratio (see DensityCell::mask). The observed fraction is always written. 0 (the default) masks none
  double min_observed_fraction = 0.0;
};

/// Render a ray cloud according to the supplied parameters
bool RAYLIB_EXPORT renderCloud(const std::string &cloud_file, const Cuboid &bounds, ViewDirection view_direction,
                               RenderStyle style, double pix_width, const std::string &image_file,
                               const std::string &projection_file, bool mark_origin,
                               const std::string *transform_file = nullptr,
                               const DensityOptions &density_options = DensityOptions(),
                               const DensityOutputs &density_outputs = DensityOutputs());

#if RAYLIB_WITH_TIFF
// save to geotif format using floating-point per-channel colour data. This function passes a projection file in order
// to geolocate the image
bool RAYLIB_EXPORT writeGeoTiffFloat(const std::string &filename, int x, int y, const float *data, double pixel_width,
                                     bool scalar, const std::string &projection_file, double origin_x, double origin_y);
#endif
}  // namespace ray
#endif  // RAYLIB_RAYRENDERER_H
