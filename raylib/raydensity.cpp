// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#include "raydensity.h"

#include "raymesh.h"
#include "rayply.h"

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
    else if (static_cast<unsigned char>(c) >= 0x20)
    {
      ss << c;
    }
  }
  ss << "\"";
  return ss.str();
}
}  // namespace

bool ZarrStore::writeFile(const std::string &file_name, const char *data, size_t size) const
{
  // write to a uniquely named temporary file, then rename, so that concurrent writers don't conflict
  static int counter = 0;
  std::stringstream tmp;
  tmp << file_name << ".tmp" << reinterpret_cast<uintptr_t>(this) << "_" << counter++ << "_"
      << std::hash<std::string>()(file_name);
  {
    std::ofstream ofs(tmp.str(), std::ios::binary | std::ios::out);
    if (!ofs.is_open())
    {
      std::cerr << "Error: cannot write " << tmp.str() << std::endl;
      return false;
    }
    ofs.write(data, static_cast<std::streamsize>(size));
    if (!ofs.good())
    {
      std::cerr << "Error: failed to write " << tmp.str() << std::endl;
      return false;
    }
  }
  std::error_code error;
  std::filesystem::rename(tmp.str(), file_name, error);
  if (error)
  {
    std::cerr << "Error: cannot rename " << tmp.str() << " to " << file_name << ": " << error.message() << std::endl;
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
                           const std::vector<std::string> &dimensions, const std::string &attributes) const
{
  const std::string dir = path_ + "/" + name;
  std::error_code error;
  std::filesystem::create_directories(dir, error);
  std::stringstream zarray;
  zarray << "{\n  \"chunks\": " << jsonList(chunks) << ",\n  \"compressor\": null,\n  \"dtype\": \"<f4\",\n"
         << "  \"fill_value\": \"NaN\",\n  \"filters\": null,\n  \"order\": \"C\",\n  \"shape\": " << jsonList(shape)
         << ",\n  \"zarr_format\": 2\n}\n";
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

bool ZarrStore::writeChunk(const std::string &name, const std::vector<int> &chunk_index,
                           const std::vector<float> &data) const
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
  return writeFile(path_ + "/" + name + "/" + key.str(), reinterpret_cast<const char *>(data.data()),
                   data.size() * sizeof(float));
}

bool ZarrStore::write1D(const std::string &name, const std::vector<float> &data, const std::string &attributes) const
{
  const int size = static_cast<int>(data.size());
  return writeArray(name, { size }, { size }, { name }, attributes) && writeChunk(name, { 0 }, data);
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
  if ((params.tile.array() < 0).any() || (params.tile.array() >= params.num_tiles.array()).any())
  {
    std::cerr << "Error: tile " << params.tile.transpose() << " is outside the domain of "
              << params.num_tiles.transpose() << " tiles" << std::endl;
    return false;
  }
  const double halo = params.halo >= 0.0 ? params.halo : params.options.pool_radius + 2.0 * w;
  const int halo_voxels = static_cast<int>(std::ceil(halo / w - 1e-9));
  const int nx = tile_voxels + 2 * halo_voxels;
  const Eigen::Vector2d core_min(params.origin[0] + params.tile[0] * params.tile_width,
                                 params.origin[1] + params.tile[1] * params.tile_width);
  const Eigen::Vector2d grid_min = core_min - Eigen::Vector2d(halo_voxels * w, halo_voxels * w);

  // ground height at each voxel column centre
  std::vector<double> ground(static_cast<size_t>(nx * nx), 0.0);
  if (!params.ground_mesh_file.empty())
  {
    Mesh mesh;
    if (!readPlyMesh(params.ground_mesh_file, mesh))
    {
      return false;
    }
    Eigen::ArrayXXd field;
    const Eigen::Vector3d box_min(grid_min[0], grid_min[1], -1e6);
    const Eigen::Vector3d box_max(grid_min[0] + nx * w, grid_min[1] + nx * w, 1e6);
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
  int shift = 0;
  std::unique_ptr<DensityGrid> grid_ptr =
    estimateDensityGrid(cloud_file, bounds, w, Eigen::Vector3i(nx, nx, nz), params.options, shift);
  if (!grid_ptr)
  {
    return false;
  }
  const DensityGrid &grid = *grid_ptr;
  const Eigen::Vector3i full_dims = grid.dimensions();

  // aggregate the tile's core voxels into cells, by height above ground
  const size_t num_cells = static_cast<size_t>(cells_x * cells_y * cells_z);
  std::vector<double> sum(num_cells, 0.0), hits(num_cells, 0.0), path(num_cells, 0.0);
  std::vector<int> observed(num_cells, 0), total(num_cells, 0);
  const size_t num_voxel_out = static_cast<size_t>(tile_voxels * tile_voxels * voxels_z);
  std::vector<float> voxel_out(params.write_voxels ? num_voxel_out : 0, std::numeric_limits<float>::quiet_NaN());
  for (int x = halo_voxels; x < halo_voxels + tile_voxels; x++)
  {
    for (int y = halo_voxels; y < halo_voxels + tile_voxels; y++)
    {
      const double g = ground[static_cast<size_t>(x + nx * y)];
      const int cx = static_cast<int>(((x - halo_voxels) + 0.5) * w / params.cell_size[0]);
      const int cy = static_cast<int>(((y - halo_voxels) + 0.5) * w / params.cell_size[1]);
      for (int z = 0; z < nz; z++)
      {
        const double height = z0 + (z + 0.5) * w - g - params.min_height;  // of the voxel centre, from min_height
        if (height < 0.0 || height >= params.max_height - params.min_height)
        {
          continue;
        }
        const int cz = std::min(cells_z - 1, static_cast<int>(height / params.cell_size[2]));
        const size_t c =
          static_cast<size_t>(std::min(cx, cells_x - 1) + cells_x * (std::min(cy, cells_y - 1) + cells_y * cz));
        total[c]++;
        const DensityGrid::Voxel &voxel = grid.voxels()[static_cast<size_t>(grid.getIndex(Eigen::Vector3i(x, y, z)))];
        if (voxel.numRays() <= DensityGrid::min_voxel_hits)
        {
          continue;
        }
        const double density = voxelAverageDensity(grid, voxel, x, y, z, full_dims, shift);
        observed[c]++;
        sum[c] += density;
        hits[c] += voxel.numHits();
        path[c] += voxel.pathLength();
        if (params.write_voxels)
        {
          const int vz = std::min(voxels_z - 1, static_cast<int>(height / w));
          voxel_out[static_cast<size_t>((x - halo_voxels) + tile_voxels * ((y - halo_voxels) + tile_voxels * vz))] =
            static_cast<float>(density);
        }
      }
    }
  }

  // chunk data in (z, y, x) C order
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const bool cylinders = params.options.distribution == ElementDistribution::VerticalCylinder;
  const double d = params.options.stem_diameter;
  std::vector<float> density(num_cells, nan), single(num_cells, nan), hits_out(num_cells), path_out(num_cells),
    observed_fraction(num_cells, 0.0f), solid(num_cells, nan);
  for (int z = 0; z < cells_z; z++)
  {
    for (int y = 0; y < cells_y; y++)
    {
      for (int x = 0; x < cells_x; x++)
      {
        const size_t c = static_cast<size_t>(x + cells_x * (y + cells_y * z));  // also the C order index for (z,y,x)
        hits_out[c] = static_cast<float>(hits[c]);
        path_out[c] = static_cast<float>(path[c]);
        if (total[c] > 0)
        {
          observed_fraction[c] = static_cast<float>(observed[c]) / static_cast<float>(total[c]);
        }
        if (observed[c] > 0)
        {
          density[c] = static_cast<float>(sum[c] / observed[c]);
          if (cylinders && d > 0.0)
          {
            solid[c] = static_cast<float>(0.25 * kPi * d * density[c]);
          }
        }
        if (path[c] > 0.0)
        {
          single[c] = static_cast<float>(hits[c] / path[c]);
        }
      }
    }
  }

  // write the metadata (identical for every tile), then this tile's chunks
  const ZarrStore store(store_path);
  const int shape_x = cells_x * params.num_tiles[0], shape_y = cells_y * params.num_tiles[1];
  std::stringstream attrs;
  attrs << std::setprecision(10);
  attrs << "  \"title\": \"raycloudtools raydensity\",\n"
        << "  \"quantity\": "
        << jsonString(cylinders ? "frontal area per volume N d (m^-1)" : "leaf area per volume (m^-1)")
        << ",\n  \"origin\": [" << params.origin[0] << ", " << params.origin[1] << "],\n"
        << "  \"tile_width\": " << params.tile_width << ",\n  \"num_tiles\": [" << params.num_tiles[0] << ", "
        << params.num_tiles[1] << "],\n  \"voxel_width\": " << w << ",\n  \"cell_size\": [" << params.cell_size[0]
        << ", " << params.cell_size[1] << ", " << params.cell_size[2] << "],\n  \"height_range\": ["
        << params.min_height << ", " << params.max_height
        << "],\n  \"heights_above_ground\": " << (params.ground_mesh_file.empty() ? "false" : "true")
        << ",\n  \"stem_diameter\": " << d << ",\n  \"footprint\": [" << params.options.footprint.exit_diameter << ", "
        << params.options.footprint.divergence << ", " << params.options.footprint.threshold
        << "],\n  \"estimator\": " << static_cast<int>(params.options.path_estimator)
        << ",\n  \"pool_radius\": " << params.options.pool_radius
        << ",\n  \"min_gap_ratio\": " << params.options.min_gap_ratio
        << ",\n  \"occlusion_correction\": " << (params.options.overlap_correction ? "true" : "false")
        << ",\n  \"description\": " << jsonString(params.description);
  const std::vector<int> shape = { cells_z, shape_y, shape_x }, chunks = { cells_z, cells_y, cells_x };
  const std::vector<std::string> dims = { "z", "y", "x" };
  const std::vector<int> chunk_index = { 0, params.tile[1], params.tile[0] };
  struct Output
  {
    std::string name;
    const std::vector<float> *data;
    std::string description;
  };
  std::vector<Output> outputs = {
    { "density", &density, "mean of the observed voxel densities in the cell" },
    { "single_voxel_density", &single, "total hits / total weighted path length, as if the cell were one voxel" },
    { "hits", &hits_out, "sum of (corrected) hits over the observed voxels" },
    { "weighted_path", &path_out, "sum of the G-weighted path lengths (m) over the observed voxels" },
    { "observed_fraction", &observed_fraction, "fraction of the cell's voxels with enough rays for an estimate" }
  };
  if (cylinders && d > 0.0)
  {
    outputs.push_back({ "solid_fraction", &solid, "N pi d^2 / 4, the fraction of the ground area occupied by stems" });
  }
  bool ok = store.writeGroup(attrs.str());
  std::vector<float> xs(static_cast<size_t>(shape_x)), ys(static_cast<size_t>(shape_y)),
    zs(static_cast<size_t>(cells_z));
  for (int i = 0; i < shape_x; i++)
  {
    xs[static_cast<size_t>(i)] = static_cast<float>(params.origin[0] + (i + 0.5) * params.cell_size[0]);
  }
  for (int i = 0; i < shape_y; i++)
  {
    ys[static_cast<size_t>(i)] = static_cast<float>(params.origin[1] + (i + 0.5) * params.cell_size[1]);
  }
  for (int i = 0; i < cells_z; i++)
  {
    zs[static_cast<size_t>(i)] = static_cast<float>(params.min_height + (i + 0.5) * params.cell_size[2]);
  }
  ok = ok && store.write1D("x", xs, "  \"units\": \"m\"") && store.write1D("y", ys, "  \"units\": \"m\"") &&
       store.write1D("z", zs, "  \"units\": \"m\", \"long_name\": \"cell centre height\"");
  for (const auto &output : outputs)
  {
    ok =
      ok && store.writeArray(output.name, shape, chunks, dims, "  \"description\": " + jsonString(output.description));
    ok = ok && store.writeChunk(output.name, chunk_index, *output.data);
  }
  if (params.write_voxels)
  {
    const int vx = tile_voxels * params.num_tiles[0], vy = tile_voxels * params.num_tiles[1];
    std::vector<float> xv(static_cast<size_t>(vx)), yv(static_cast<size_t>(vy)), zv(static_cast<size_t>(voxels_z));
    for (int i = 0; i < vx; i++)
    {
      xv[static_cast<size_t>(i)] = static_cast<float>(params.origin[0] + (i + 0.5) * w);
    }
    for (int i = 0; i < vy; i++)
    {
      yv[static_cast<size_t>(i)] = static_cast<float>(params.origin[1] + (i + 0.5) * w);
    }
    for (int i = 0; i < voxels_z; i++)
    {
      zv[static_cast<size_t>(i)] = static_cast<float>(params.min_height + (i + 0.5) * w);
    }
    ok = ok && store.write1D("xv", xv) && store.write1D("yv", yv) && store.write1D("zv", zv);
    ok = ok && store.writeArray("voxel_density", { voxels_z, vy, vx }, { voxels_z, tile_voxels, tile_voxels },
                                { "zv", "yv", "xv" }, "  \"description\": \"voxel density, NaN where unobserved\"");
    ok = ok && store.writeChunk("voxel_density", chunk_index, voxel_out);
  }
  if (ok)
  {
    std::cout << "wrote tile " << params.tile.transpose() << " to " << store_path << std::endl;
  }
  return ok;
}
}  // namespace ray
