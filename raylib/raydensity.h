// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#ifndef RAYLIB_RAYDENSITY_H
#define RAYLIB_RAYDENSITY_H

#include "raylib/raylibconfig.h"

#include "raydensitygrid.h"
#include "rayutils.h"

#include <string>
#include <vector>

namespace ray
{
/// Minimal writer for a Zarr (version 2) store of uncompressed, little-endian float32 (or float64) arrays in C order.
/// Readable by zarr-python, xarray (open_zarr) and GDAL. Files are written to a temporary name (unique to the process
/// and call) and then renamed, so concurrent writers of the same (identical) metadata, or of different chunks, don't
/// conflict.
class RAYLIB_EXPORT ZarrStore
{
public:
  explicit ZarrStore(const std::string &path)
    : path_(path)
  {}
  /// Writes the group metadata, with the attributes given as the members of a JSON object (without braces)
  bool writeGroup(const std::string &attributes) const;
  /// Writes an array's metadata. @c dimensions are the dimension names (for xarray), slowest varying first. @c dtype
  /// is "<f4" (float32) or "<f8" (float64), and must match the data given to writeChunk
  bool writeArray(const std::string &name, const std::vector<int> &shape, const std::vector<int> &chunks,
                  const std::vector<std::string> &dimensions, const std::string &attributes = "",
                  const std::string &dtype = "<f4") const;
  /// Writes one chunk of a float32 array, of the full chunk shape (C order)
  bool writeChunk(const std::string &name, const std::vector<int> &chunk_index, const std::vector<float> &data) const;
  /// Writes one chunk of a float64 array, of the full chunk shape (C order)
  bool writeChunk(const std::string &name, const std::vector<int> &chunk_index, const std::vector<double> &data) const;
  /// Writes a 1D float32 array as a single chunk
  bool write1D(const std::string &name, const std::vector<float> &data, const std::string &attributes = "") const;
  /// Writes a 1D float64 array as a single chunk (e.g. a coordinate, which needs float64 for georeferenced values)
  bool write1D(const std::string &name, const std::vector<double> &data, const std::string &attributes = "") const;
  /// Writes a file that isn't part of the Zarr data model, at @c relative_path within the store (creating its
  /// directory). Zarr readers ignore directories without Zarr metadata
  bool writeText(const std::string &relative_path, const std::string &text) const;
  /// Reads a file at @c relative_path within the store. Returns false if it doesn't exist
  bool readText(const std::string &relative_path, std::string &text) const;

private:
  bool writeFile(const std::string &file_name, const char *data, size_t size) const;
  bool writeChunkData(const std::string &name, const std::vector<int> &chunk_index, const char *data,
                      size_t size) const;
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
  /// in a single run over the whole domain. It must be at least pool_radius plus two voxel widths (the pooled
  /// neighbourhood and the neighbour prior). Negative chooses this minimum. The voxel grid covers the halo rounded up
  /// to whole voxels, plus one voxel, so the input (e.g. from splitTiles) needs a halo at least 1.5 voxels larger,
  /// for the ray ends that it clips to lie outside the grid
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
  /// Free text describing this tile's run (e.g. the command line), written with the tile's index, input and halo to
  /// tiles/i_j.json in the store. It is not in the store's attributes, which are identical for every tile
  std::string description;
};

/// Estimates the density in one tile from @c cloud_file (for example, from raysplit tiles, which keeps every ray
/// crossing the tile and its halo), and writes the tile's chunk of each array in the Zarr store at @c store_path.
/// The store's metadata, identical for every tile, is written first. If the store already has metadata that differs
/// (it was written with different parameters), the tile fails without writing. Output arrays on dimensions (z, y, x)
/// of cells, with float64 coordinates x, y (cell centres) and z (cell centre height), from aggregateDensityCells:
///   density: pooled over the cell's observed voxels, total (corrected) hits over total weighted path (N d for
///     vertical_cylinder, in m^-1)
///   mean_voxel_density: the mean of the observed voxel densities in the cell
///   patchiness_ratio: mean_voxel_density / density, above 1 where the density varies within the cell
///   hits, weighted_path: the sums over the cell's observed voxels
///   observed_fraction: the fraction of the cell's voxels with enough rays for an estimate
///   solid_fraction: N pi d^2 / 4 for N d = density, when a stem diameter is given
/// and optionally voxel_density on dimensions (zv, yv, xv). Unobserved cells and voxels are NaN, except for the
/// observed fraction (0). A tile that no ray crosses is valid: it is written as unobserved, with a warning. Returns
/// false (with a message) on invalid parameters, a missing or unreadable input, or a failed write.
bool RAYLIB_EXPORT estimateDensityTile(const std::string &cloud_file, const DensityTileParams &params,
                                       const std::string &store_path);
}  // namespace ray

#endif  // RAYLIB_RAYDENSITY_H
