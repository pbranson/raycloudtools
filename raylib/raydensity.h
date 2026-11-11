// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#ifndef RAYLIB_RAYDENSITY_H
#define RAYLIB_RAYDENSITY_H

#include "raylib/raylibconfig.h"

#include "rayrenderer.h"
#include "rayutils.h"

#include <string>
#include <vector>

namespace ray
{
/// Minimal writer for a Zarr (version 2) store of uncompressed, little-endian float32 arrays in C order. Readable by
/// zarr-python, xarray (open_zarr) and GDAL. Files are written to a temporary name and then renamed, so concurrent
/// writers of the same (identical) metadata, or of different chunks, don't conflict.
class RAYLIB_EXPORT ZarrStore
{
public:
  explicit ZarrStore(const std::string &path)
    : path_(path)
  {}
  /// Writes the group metadata, with the attributes given as the members of a JSON object (without braces)
  bool writeGroup(const std::string &attributes) const;
  /// Writes an array's metadata. @c dimensions are the dimension names (for xarray), slowest varying first
  bool writeArray(const std::string &name, const std::vector<int> &shape, const std::vector<int> &chunks,
                  const std::vector<std::string> &dimensions, const std::string &attributes = "") const;
  /// Writes one chunk of an array, of the full chunk shape (C order)
  bool writeChunk(const std::string &name, const std::vector<int> &chunk_index, const std::vector<float> &data) const;
  /// Writes a 1D array as a single chunk (e.g. a coordinate)
  bool write1D(const std::string &name, const std::vector<float> &data, const std::string &attributes = "") const;

private:
  bool writeFile(const std::string &file_name, const char *data, size_t size) const;
  std::string path_;
};

/// Parameters for estimating the stem (or foliage) density in one square tile of a larger domain, and writing it into
/// a Zarr store holding the whole domain, as one chunk per tile.
struct RAYLIB_EXPORT DensityTileParams
{
  /// The domain origin (x, y). Tile (i, j) covers [origin + i tile_width, origin + (i+1) tile_width)
  Eigen::Vector2d origin = Eigen::Vector2d(0.0, 0.0);
  double tile_width = 25.0;
  /// This tile's index, from 0 to num_tiles - 1 in each axis
  Eigen::Vector2i tile = Eigen::Vector2i(0, 0);
  /// The number of tiles in the domain in x and y
  Eigen::Vector2i num_tiles = Eigen::Vector2i(1, 1);
  /// Rays within this distance outside the tile are included, so that the tile's voxels are estimated exactly as
  /// in a single run over the whole domain. Negative chooses the minimum: pool_radius plus two voxels
  double halo = -1.0;
  double voxel_width = 0.05;
  /// The output cells. The tile width must be a multiple of the cell width in x and y, and of the voxel width
  Eigen::Vector3d cell_size = Eigen::Vector3d(1.0, 1.0, 0.1);
  /// The range of heights to output, above the ground mesh if given, otherwise absolute. A multiple of the cell
  /// height
  double min_height = 0.0;
  double max_height = 0.6;
  /// Optional ground mesh (e.g. from rayextract terrain). Each voxel is assigned to a layer by the height of its
  /// centre above the ground beneath it
  std::string ground_mesh_file;
  /// Also output the voxel densities (as heights above ground), not just the cells
  bool write_voxels = false;
  DensityOptions options;
  /// Free text stored in the store's attributes, e.g. the command line
  std::string description;
};

/// Estimates the density in one tile from @c cloud_file (for example, from raysplit tiles, which keeps every ray
/// crossing the tile and its halo), and writes the tile's chunk of each array in the Zarr store at @c store_path.
/// The store's metadata (identical for every tile) is also written. Output arrays on dimensions (z, y, x) of cells:
///   density: the mean of the observed voxel densities in the cell (N d for vertical_cylinder, in m^-1)
///   single_voxel_density: total hits over total weighted path in the cell, as if one voxel the size of the cell
///   hits, weighted_path: the sums over the cell's observed voxels
///   observed_fraction: the fraction of the cell's voxels with enough rays for an estimate
///   solid_fraction: N pi d^2 / 4, when a stem diameter is given
/// and optionally voxel_density on dimensions (zv, yv, xv). Unobserved cells and voxels are NaN.
bool RAYLIB_EXPORT estimateDensityTile(const std::string &cloud_file, const DensityTileParams &params,
                                       const std::string &store_path);
}  // namespace ray

#endif  // RAYLIB_RAYDENSITY_H
