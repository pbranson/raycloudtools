// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#include "raydensity.h"

#include "raymesh.h"
#include "rayply.h"

#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

namespace ray
{
namespace
{
bool isLittleEndian()
{
  const uint16_t value = 1;
  uint8_t byte;
  std::memcpy(&byte, &value, 1);
  return byte == 1;
}

std::string jsonList(const std::vector<int> &values)
{
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < values.size(); i++)
  {
    ss << (i > 0 ? ", " : "") << values[i];
  }
  ss << "]";
  return ss.str();
}

std::string jsonString(const std::string &text)
{
  std::stringstream ss;
  ss << "\"";
  for (const char c : text)
  {
    if (c == '"' || c == '\\')
    {
      ss << '\\' << c;
    }
    else if (c == '\n')
    {
      ss << "\\n";
    }
    else if (c == '\t')
    {
      ss << "\\t";
    }
    else if (c == '\r')
    {
      ss << "\\r";
    }
    else if (static_cast<unsigned char>(c) < 0x20)
    {
      ss << "\\u00" << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(c) << std::dec;
    }
    else
    {
      ss << c;
    }
  }
  ss << "\"";
  return ss.str();
}

long processId()
{
#if defined(_WIN32)
  return static_cast<long>(_getpid());
#else
  return static_cast<long>(getpid());
#endif
}
}  // namespace

bool ZarrStore::writeFile(const std::string &file_name, const char *data, size_t size) const
{
  // write to a temporary file unique to this process and call, then rename (atomic on POSIX file systems), so that
  // concurrent writers don't conflict and readers never see a partial file
  static std::atomic<unsigned long> counter(0);
  std::stringstream tmp;
  tmp << file_name << ".tmp." << processId() << "." << counter++;
  const std::string tmp_name = tmp.str();
  bool ok;
  {
    std::ofstream ofs(tmp_name, std::ios::binary | std::ios::out);
    if (!ofs.is_open())
    {
      std::cerr << "Error: cannot write " << tmp_name << std::endl;
      return false;
    }
    ofs.write(data, static_cast<std::streamsize>(size));
    ofs.close();
    ok = !ofs.fail();
  }
  std::error_code error;
  if (!ok)
  {
    std::cerr << "Error: failed to write " << tmp_name << std::endl;
    std::filesystem::remove(tmp_name, error);
    return false;
  }
  std::filesystem::rename(tmp_name, file_name, error);
  if (error)
  {
    std::cerr << "Error: cannot rename " << tmp_name << " to " << file_name << ": " << error.message() << std::endl;
    std::filesystem::remove(tmp_name, error);
    return false;
  }
  return true;
}

bool ZarrStore::writeGroup(const std::string &attributes) const
{
  std::error_code error;
  std::filesystem::create_directories(path_, error);
  const std::string group = "{\n  \"zarr_format\": 2\n}\n";
  const std::string attrs = "{\n" + attributes + "\n}\n";
  return writeFile(path_ + "/.zgroup", group.data(), group.size()) &&
         writeFile(path_ + "/.zattrs", attrs.data(), attrs.size());
}

bool ZarrStore::writeArray(const std::string &name, const std::vector<int> &shape, const std::vector<int> &chunks,
                           const std::vector<std::string> &dimensions, const std::string &attributes,
                           const std::string &dtype) const
{
  if (dtype != "<f4" && dtype != "<f8")
  {
    std::cerr << "Error: unsupported Zarr dtype " << dtype << std::endl;
    return false;
  }
  const std::string dir = path_ + "/" + name;
  std::error_code error;
  std::filesystem::create_directories(dir, error);
  std::stringstream zarray;
  zarray << "{\n  \"chunks\": " << jsonList(chunks) << ",\n  \"compressor\": null,\n  \"dtype\": \"" << dtype
         << "\",\n  \"fill_value\": \"NaN\",\n  \"filters\": null,\n  \"order\": \"C\",\n  \"shape\": "
         << jsonList(shape) << ",\n  \"zarr_format\": 2\n}\n";
  std::stringstream zattrs;
  zattrs << "{\n  \"_ARRAY_DIMENSIONS\": [";
  for (size_t i = 0; i < dimensions.size(); i++)
  {
    zattrs << (i > 0 ? ", " : "") << jsonString(dimensions[i]);
  }
  zattrs << "]" << (attributes.empty() ? "" : ",\n" + attributes) << "\n}\n";
  const std::string a = zarray.str(), b = zattrs.str();
  return writeFile(dir + "/.zarray", a.data(), a.size()) && writeFile(dir + "/.zattrs", b.data(), b.size());
}

bool ZarrStore::writeChunkData(const std::string &name, const std::vector<int> &chunk_index, const char *data,
                               size_t size) const
{
  if (!isLittleEndian())
  {
    std::cerr << "Error: the Zarr writer requires a little-endian machine" << std::endl;
    return false;
  }
  std::stringstream key;
  for (size_t i = 0; i < chunk_index.size(); i++)
  {
    key << (i > 0 ? "." : "") << chunk_index[i];
  }
  return writeFile(path_ + "/" + name + "/" + key.str(), data, size);
}

bool ZarrStore::writeChunk(const std::string &name, const std::vector<int> &chunk_index,
                           const std::vector<float> &data) const
{
  return writeChunkData(name, chunk_index, reinterpret_cast<const char *>(data.data()), data.size() * sizeof(float));
}

bool ZarrStore::writeChunk(const std::string &name, const std::vector<int> &chunk_index,
                           const std::vector<double> &data) const
{
  return writeChunkData(name, chunk_index, reinterpret_cast<const char *>(data.data()), data.size() * sizeof(double));
}

bool ZarrStore::write1D(const std::string &name, const std::vector<float> &data, const std::string &attributes) const
{
  const int size = static_cast<int>(data.size());
  return writeArray(name, { size }, { size }, { name }, attributes, "<f4") && writeChunk(name, { 0 }, data);
}

bool ZarrStore::write1D(const std::string &name, const std::vector<double> &data, const std::string &attributes) const
{
  const int size = static_cast<int>(data.size());
  return writeArray(name, { size }, { size }, { name }, attributes, "<f8") && writeChunk(name, { 0 }, data);
}

bool ZarrStore::writeText(const std::string &relative_path, const std::string &text) const
{
  const std::filesystem::path file = std::filesystem::path(path_) / relative_path;
  std::error_code error;
  std::filesystem::create_directories(file.parent_path(), error);
  return writeFile(file.string(), text.data(), text.size());
}

bool ZarrStore::readText(const std::string &relative_path, std::string &text) const
{
  std::ifstream ifs(path_ + "/" + relative_path, std::ios::binary | std::ios::in);
  if (!ifs.is_open())
  {
    return false;
  }
  std::stringstream ss;
  ss << ifs.rdbuf();
  text = ss.str();
  return true;
}

namespace
{
/// whether x is (close to) an integer multiple of unit, and the multiple
bool multipleOf(double x, double unit, int &multiple)
{
  const double ratio = x / unit;
  multiple = static_cast<int>(std::round(ratio));
  return multiple > 0 && std::abs(ratio - multiple) < 1e-6 * std::max(1.0, ratio);
}

const char *estimatorName(PathEstimator estimator)
{
  switch (estimator)
  {
  case PathEstimator::BeerLambert:
    return "beer_lambert";
  case PathEstimator::Crossing:
    return "crossing";
  case PathEstimator::PooledBeerLambert:
    return "pooled";
  }
  return "unknown";
}

/// Whether @c file_name is a .ply file whose header declares no vertices (a ray cloud without rays)
bool isEmptyRayCloud(const std::string &file_name)
{
  std::ifstream ifs(file_name, std::ios::binary | std::ios::in);
  std::string line;
  if (!ifs.is_open() || !std::getline(ifs, line) || line.rfind("ply", 0) != 0)
  {
    return false;
  }
  bool empty = false;
  while (std::getline(ifs, line))
  {
    std::stringstream ss(line);
    std::string word, element;
    long long count = -1;
    ss >> word;
    if (word == "end_header")
    {
      return empty;
    }
    if (word == "element" && ss >> element >> count && element == "vertex")
    {
      empty = count == 0;
    }
  }
  return false;
}

/// Removes the triangles whose horizontal bounds don't overlap the box [box_min, box_max] in x and y, so that
/// Mesh::toHeightField only bins the triangles near the tile, rather than those of the whole domain
void cropMesh(Mesh &mesh, const Eigen::Vector2d &box_min, const Eigen::Vector2d &box_max)
{
  const std::vector<Eigen::Vector3d> &vertices = mesh.vertices();
  std::vector<Eigen::Vector3i> &triangles = mesh.indexList();
  size_t kept = 0;
  for (const Eigen::Vector3i &triangle : triangles)
  {
    Eigen::Vector2d tri_min = vertices[static_cast<size_t>(triangle[0])].head<2>(), tri_max = tri_min;
    for (int j = 1; j < 3; j++)
    {
      const Eigen::Vector2d v = vertices[static_cast<size_t>(triangle[j])].head<2>();
      tri_min = tri_min.cwiseMin(v);
      tri_max = tri_max.cwiseMax(v);
    }
    if ((tri_max.array() >= box_min.array()).all() && (tri_min.array() <= box_max.array()).all())
    {
      triangles[kept++] = triangle;
    }
  }
  triangles.resize(kept);
}
}  // namespace

bool estimateDensityTile(const std::string &cloud_file, const DensityTileParams &params, const std::string &store_path)
{
  const double w = params.voxel_width;
  int tile_voxels, cells_x, cells_y, cells_z, voxels_z;
  if (!multipleOf(params.tile_width, w, tile_voxels) || !multipleOf(params.tile_width, params.cell_size[0], cells_x) ||
      !multipleOf(params.tile_width, params.cell_size[1], cells_y) ||
      !multipleOf(params.max_height - params.min_height, params.cell_size[2], cells_z) ||
      !multipleOf(params.max_height - params.min_height, w, voxels_z))
  {
    std::cerr << "Error: the tile width must be a multiple of the voxel width and the cell width, and the height "
                 "range a multiple of the cell height and voxel width"
              << std::endl;
    return false;
  }
  if ((params.num_tiles.array() < 1).any() || (params.tile.array() < 0).any() ||
      (params.tile.array() >= params.num_tiles.array()).any())
  {
    std::cerr << "Error: tile " << params.tile.transpose() << " is outside the domain of "
              << params.num_tiles.transpose() << " tiles" << std::endl;
    return false;
  }
  // the halo must cover the pooled neighbourhood and the neighbour prior of the tile's edge voxels, for the tile to
  // match a single run over the whole domain
  const double min_halo = params.options.pool_radius + 2.0 * w;
  if (params.halo >= 0.0 && params.halo < min_halo - 1e-6 * w)
  {
    std::cerr << "Error: the halo " << params.halo << " m is less than the pool radius plus two voxel widths ("
              << min_halo << " m), so the tile would not match a single run over the domain" << std::endl;
    return false;
  }
  const double halo = params.halo >= 0.0 ? params.halo : min_halo;
  const int halo_voxels = static_cast<int>(std::ceil(halo / w - 1e-9));
  const int nx = tile_voxels + 2 * halo_voxels;
  const Eigen::Vector2d core_min(params.origin[0] + params.tile[0] * params.tile_width,
                                 params.origin[1] + params.tile[1] * params.tile_width);
  const Eigen::Vector2d grid_min = core_min - Eigen::Vector2d(halo_voxels * w, halo_voxels * w);
  const bool cylinders = params.options.distribution == ElementDistribution::VerticalCylinder;
  const double d = params.options.stem_diameter;
  const bool write_solid = cylinders && d > 0.0;

  // the store's metadata, identical for every tile of the domain (per-tile information goes in tiles/i_j.json)
  const ZarrStore store(store_path);
  const int shape_x = cells_x * params.num_tiles[0], shape_y = cells_y * params.num_tiles[1];
  std::stringstream attrs;
  attrs << std::setprecision(15);
  attrs << "  \"title\": \"raycloudtools raydensity\",\n"
        << "  \"quantity\": "
        << jsonString(cylinders ? "frontal area per volume N d (m^-1)" : "leaf area per volume (m^-1)")
        << ",\n  \"density_model\": " << jsonString(cylinders ? "vertical_cylinder" : "spherical")
        << ",\n  \"origin\": [" << params.origin[0] << ", " << params.origin[1] << "],\n"
        << "  \"tile_width\": " << params.tile_width << ",\n  \"num_tiles\": [" << params.num_tiles[0] << ", "
        << params.num_tiles[1] << "],\n  \"voxel_width\": " << w << ",\n  \"cell_size\": [" << params.cell_size[0]
        << ", " << params.cell_size[1] << ", " << params.cell_size[2] << "],\n  \"height_range\": ["
        << params.min_height << ", " << params.max_height
        << "],\n  \"heights_above_ground\": " << (params.ground_mesh_file.empty() ? "false" : "true")
        << ",\n  \"voxels\": " << (params.write_voxels ? "true" : "false") << ",\n  \"stem_diameter\": " << d
        << ",\n  \"footprint\": [" << params.options.footprint.exit_diameter << ", "
        << params.options.footprint.divergence << ", " << params.options.footprint.threshold
        << "],\n  \"return_range\": " << (params.options.footprint.range_to_outline ? "\"outline\"" : "\"surface\"")
        << ",\n  \"estimator\": " << jsonString(estimatorName(params.options.path_estimator))
        << ",\n  \"pool_radius\": " << params.options.pool_radius
        << ",\n  \"min_gap_ratio\": " << params.options.min_gap_ratio
        << ",\n  \"min_sin_zenith\": " << params.options.min_sin_zenith
        << ",\n  \"max_sin_zenith\": " << params.options.max_sin_zenith
        << ",\n  \"occlusion_correction\": " << (params.options.overlap_correction ? "true" : "false")
        << ",\n  \"flat_top_compensation\": " << (params.options.flat_top_compensation ? "true" : "false");
  const std::string group_attrs = "{\n" + attrs.str() + "\n}\n";
  std::string existing;
  if (store.readText(".zattrs", existing) && existing != group_attrs)
  {
    std::cerr << "Error: the store " << store_path
              << " was written with different parameters (see its .zattrs). Use the same domain and options for "
                 "every tile, or a new store"
              << std::endl;
    return false;
  }
  const std::vector<int> shape = { cells_z, shape_y, shape_x }, chunks = { cells_z, cells_y, cells_x };
  const std::vector<std::string> dims = { "z", "y", "x" };
  const std::vector<int> chunk_index = { 0, params.tile[1], params.tile[0] };
  struct Output
  {
    std::string name;
    std::string description;
    std::vector<float> data;
  };
  // the order of the outputs is that of the chunk data below
  std::vector<Output> outputs = {
    { "density", "total (corrected) hits / total weighted path length over the cell's observed voxels", {} },
    { "mean_voxel_density", "mean of the observed voxel densities in the cell", {} },
    { "patchiness_ratio", "mean_voxel_density / density, above 1 where the density varies within the cell", {} },
    { "hits", "sum of (corrected) hits over the observed voxels", {} },
    { "weighted_path", "sum of the G-weighted path lengths (m) over the observed voxels", {} },
    { "observed_fraction", "fraction of the cell's voxels with enough rays for an estimate", {} }
  };
  if (write_solid)
  {
    outputs.push_back(
      { "solid_fraction", "N pi d^2 / 4 for N d = density, the fraction of the ground area occupied by stems", {} });
  }
  // coordinates in float64, as float32 can't resolve cells or voxels at georeferenced (e.g. UTM) coordinates
  std::vector<double> xs(static_cast<size_t>(shape_x)), ys(static_cast<size_t>(shape_y)),
    zs(static_cast<size_t>(cells_z));
  for (int i = 0; i < shape_x; i++)
  {
    xs[static_cast<size_t>(i)] = params.origin[0] + (i + 0.5) * params.cell_size[0];
  }
  for (int i = 0; i < shape_y; i++)
  {
    ys[static_cast<size_t>(i)] = params.origin[1] + (i + 0.5) * params.cell_size[1];
  }
  for (int i = 0; i < cells_z; i++)
  {
    zs[static_cast<size_t>(i)] = params.min_height + (i + 0.5) * params.cell_size[2];
  }
  bool ok = store.writeGroup(attrs.str()) && store.write1D("x", xs, "  \"units\": \"m\"") &&
            store.write1D("y", ys, "  \"units\": \"m\"") &&
            store.write1D("z", zs, "  \"units\": \"m\", \"long_name\": \"cell centre height\"");
  for (const auto &output : outputs)
  {
    ok =
      ok && store.writeArray(output.name, shape, chunks, dims, "  \"description\": " + jsonString(output.description));
  }
  const int vx = tile_voxels * params.num_tiles[0], vy = tile_voxels * params.num_tiles[1];
  if (ok && params.write_voxels)
  {
    std::vector<double> xv(static_cast<size_t>(vx)), yv(static_cast<size_t>(vy)), zv(static_cast<size_t>(voxels_z));
    for (int i = 0; i < vx; i++)
    {
      xv[static_cast<size_t>(i)] = params.origin[0] + (i + 0.5) * w;
    }
    for (int i = 0; i < vy; i++)
    {
      yv[static_cast<size_t>(i)] = params.origin[1] + (i + 0.5) * w;
    }
    for (int i = 0; i < voxels_z; i++)
    {
      zv[static_cast<size_t>(i)] = params.min_height + (i + 0.5) * w;
    }
    ok = store.write1D("xv", xv, "  \"units\": \"m\"") && store.write1D("yv", yv, "  \"units\": \"m\"") &&
         store.write1D("zv", zv, "  \"units\": \"m\", \"long_name\": \"voxel centre height\"") &&
         store.writeArray("voxel_density", { voxels_z, vy, vx }, { voxels_z, tile_voxels, tile_voxels },
                          { "zv", "yv", "xv" }, "  \"description\": \"voxel density, NaN where unobserved\"");
  }
  if (!ok)
  {
    return false;
  }

  // ground height at each voxel column centre
  std::vector<double> ground(static_cast<size_t>(nx * nx), 0.0);
  if (!params.ground_mesh_file.empty())
  {
    Mesh mesh;
    if (!readPlyMesh(params.ground_mesh_file, mesh))
    {
      return false;
    }
    const Eigen::Vector3d box_min(grid_min[0], grid_min[1], -1e6);
    const Eigen::Vector3d box_max(grid_min[0] + nx * w, grid_min[1] + nx * w, 1e6);
    cropMesh(mesh, box_min.head<2>(), box_max.head<2>());
    if (mesh.indexList().empty())
    {
      std::cerr << "Error: the ground mesh doesn't cover tile " << params.tile.transpose() << std::endl;
      return false;
    }
    Eigen::ArrayXXd field;
    mesh.toHeightField(field, box_min, box_max, w);
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < field.rows(); i++)
    {
      for (int j = 0; j < field.cols(); j++)
      {
        if (field(i, j) > std::numeric_limits<double>::lowest())
        {
          sum += field(i, j);
          count++;
        }
      }
    }
    if (count == 0)
    {
      std::cerr << "Error: the ground mesh doesn't cover tile " << params.tile.transpose() << std::endl;
      return false;
    }
    for (int x = 0; x < nx; x++)
    {
      for (int y = 0; y < nx; y++)
      {
        const double h =
          field(std::min(x, static_cast<int>(field.rows()) - 1), std::min(y, static_cast<int>(field.cols()) - 1));
        ground[static_cast<size_t>(x + nx * y)] = h > std::numeric_limits<double>::lowest() ? h : sum / count;
      }
    }
  }
  double ground_min = std::numeric_limits<double>::max(), ground_max = std::numeric_limits<double>::lowest();
  for (const double g : ground)
  {
    ground_min = std::min(ground_min, g);
    ground_max = std::max(ground_max, g);
  }
  // the vertical extent of the voxel grid, aligned to the global voxel lattice, with a voxel's margin
  const double z0 = std::floor((ground_min + params.min_height) / w) * w - w;
  const int nz = static_cast<int>(std::ceil((ground_max + params.max_height - z0) / w - 1e-9)) + 1;

  const Cuboid bounds(Eigen::Vector3d(grid_min[0], grid_min[1], z0),
                      Eigen::Vector3d(grid_min[0] + nx * w, grid_min[1] + nx * w, z0 + nz * w));

  // aggregate the tile's core voxels into cells, by height above ground (of the voxel centres)
  CellLayout layout;
  layout.min_bound = Eigen::Vector3d(core_min[0], core_min[1], params.min_height);
  layout.cell_size = params.cell_size;
  layout.num_cells = Eigen::Vector3i(cells_x, cells_y, cells_z);
  const size_t num_voxel_out = static_cast<size_t>(tile_voxels * tile_voxels * voxels_z);
  const float nan = std::numeric_limits<float>::quiet_NaN();
  std::vector<float> voxel_out(params.write_voxels ? num_voxel_out : 0, nan);
  // every cell has voxels (the grid spans the height range above the lowest and highest ground)
  std::vector<DensityCell> cells(static_cast<size_t>(cells_x * cells_y * cells_z));
  for (DensityCell &cell : cells)
  {
    cell.total = 1;
  }
  // a tile that no ray crosses is valid (e.g. a tile of the domain beyond the scan), and is written as unobserved.
  // The ray cloud reader fails on a file without rays, so that is detected first
  const bool no_rays = isEmptyRayCloud(cloud_file);
  size_t used_rays = 0;
  if (!no_rays)
  {
    int shift = 0;
    std::unique_ptr<DensityGrid> grid_ptr =
      estimateDensityGrid(cloud_file, bounds, w, Eigen::Vector3i(nx, nx, nz), params.options, shift);
    if (!grid_ptr)
    {
      return false;
    }
    const DensityGrid &grid = *grid_ptr;
    used_rays = grid.numUsedRays();
    auto ground_height = [&](int x, int y) { return ground[static_cast<size_t>(x + nx * y)]; };
    auto write_voxel = [&](const Eigen::Vector3i &voxel, const Eigen::Vector3d &centre, double density) {
      const int vz = std::min(voxels_z - 1, static_cast<int>(std::floor(centre[2] / w)));
      voxel_out[static_cast<size_t>((voxel[0] - halo_voxels) +
                                    tile_voxels * ((voxel[1] - halo_voxels) + tile_voxels * vz))] =
        static_cast<float>(density);
    };
    cells = aggregateDensityCells(
      grid, shift, Eigen::Vector3i(halo_voxels, halo_voxels, 0),
      Eigen::Vector3i(halo_voxels + tile_voxels, halo_voxels + tile_voxels, nz), layout, ground_height,
      params.write_voxels ? std::function<void(const Eigen::Vector3i &, const Eigen::Vector3d &, double)>(write_voxel) :
                            nullptr);
  }
  size_t num_observed = 0;
  for (const DensityCell &cell : cells)
  {
    num_observed += cell.observed > 0 ? 1 : 0;
  }
  if (num_observed == 0)
  {
    std::cout << "Warning: "
              << (no_rays ? "the input has no rays" :
                            (used_rays == 0 ? "no (usable) rays of the input cross the tile and its halo" :
                                              "no voxel of the tile has enough rays"))
              << ", so tile " << params.tile.transpose() << " is written as unobserved (NaN, observed_fraction 0)"
              << std::endl;
  }

  // chunk data in (z, y, x) C order, which is the order of the cells. Unobserved cells are NaN, except for their
  // observed fraction, which is 0
  for (auto &output : outputs)
  {
    output.data.assign(cells.size(), nan);
  }
  for (size_t c = 0; c < cells.size(); c++)
  {
    const DensityCell &cell = cells[c];
    outputs[5].data[c] = static_cast<float>(cell.observedFraction());
    if (cell.observed == 0)
    {
      continue;
    }
    outputs[0].data[c] = static_cast<float>(cell.density());
    outputs[1].data[c] = static_cast<float>(cell.meanVoxelDensity());
    outputs[2].data[c] = static_cast<float>(cell.patchinessRatio());
    outputs[3].data[c] = static_cast<float>(cell.hits);
    outputs[4].data[c] = static_cast<float>(cell.weighted_path);
    if (write_solid)
    {
      outputs[6].data[c] = static_cast<float>(0.25 * kPi * d * cell.density());
    }
  }

  // this tile's chunks, then its description
  for (const auto &output : outputs)
  {
    ok = ok && store.writeChunk(output.name, chunk_index, output.data);
  }
  if (params.write_voxels)
  {
    ok = ok && store.writeChunk("voxel_density", chunk_index, voxel_out);
  }
  std::stringstream tile_name, tile_info;
  tile_name << "tiles/" << params.tile[0] << "_" << params.tile[1] << ".json";
  tile_info << std::setprecision(15) << "{\n  \"tile\": [" << params.tile[0] << ", " << params.tile[1]
            << "],\n  \"input\": " << jsonString(cloud_file) << ",\n  \"halo\": " << halo
            << ",\n  \"description\": " << jsonString(params.description) << "\n}\n";
  ok = ok && store.writeText(tile_name.str(), tile_info.str());
  if (ok)
  {
    std::cout << "wrote tile " << params.tile.transpose() << " to " << store_path << std::endl;
  }
  return ok;
}
}  // namespace ray
