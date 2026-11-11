// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raydensitygrid.h"

#include "raycloud.h"
#include "rayunused.h"

#define DENSITY_MIN_RAYS 10  // larger is more accurate but more blurred. 0 for no adaptive blending

namespace ray
{
bool DensityGrid::calculatePeaks(const std::string &file_name)
{
  auto calc_peaks = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                        std::vector<RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); ++i)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      if (!bounds_.clipRay(start, end, 1e-10))
      {
        continue;  // ray is outside of bounds
      }
      // as in calculateDensities, a ray whose end point was clipped by the bounds has no return inside the grid
      const double clip_tolerance = 1e-6;
      const bool end_clipped = (end - ends[i]).squaredNorm() > clip_tolerance * clip_tolerance;
      if (colours[i].alpha > 0 && !end_clipped)
      {
        const Eigen::Vector3d vox_end = (end - bounds_.min_bound_) / voxel_width_;
        const Eigen::Vector3i target =
          Eigen::Vector3d(std::floor(vox_end[0]), std::floor(vox_end[1]), std::floor(vox_end[2])).cast<int>();
        int peak_id = target[0] + target[1] * voxel_dims_[0];
        peaks_[peak_id] = std::max(peaks_[peak_id], vox_end[2]);
      }
    }
  };
  return Cloud::read(file_name, calc_peaks);
}

/// Calculate the density within each voxel of the grid: by default the surface area per cubic metre, assuming an
/// unbiased distribution of surface angles. See @c ElementDistribution for the alternatives.
bool DensityGrid::calculateDensities(const std::string &file_name)
{
  auto calculate = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                       std::vector<RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); ++i)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      if (!bounds_.clipRay(start, end, 1e-10))
      {
        continue;  // ray is outside of bounds
      }
      // a ray whose end point was clipped by the bounds has passed through the grid without hitting anything inside
      // it, so it is treated as a miss. Otherwise its hit would be wrongly counted in the boundary voxel.
      const double clip_tolerance = 1e-6;
      const bool end_clipped = (end - ends[i]).squaredNorm() > clip_tolerance * clip_tolerance;
      range_offset_ = (start - starts[i]).norm();
      bounded_ = colours[i].alpha > 0 && !end_clipped;
      const Eigen::Vector3d vox_start = (start - bounds_.min_bound_) / voxel_width_;
      const Eigen::Vector3d vox_end = (end - bounds_.min_bound_) / voxel_width_;
      source_ = vox_start;
      dir_ = (vox_end - vox_start).normalized();
      ray_g_ = gFunction(distribution_, dir_);
      if (ray_g_ < min_g_ || ray_g_ > max_g_ || ray_g_ <= 0.0)
      {
        num_ignored_rays_++;  // selecting rays on direction alone doesn't bias the estimate
        continue;
      }
      num_used_rays_++;
      hit_voxel_ = Eigen::Vector3d(std::floor(vox_end[0]), std::floor(vox_end[1]), std::floor(vox_end[2])).cast<int>();
      hit_length_ = (vox_end - vox_start).norm();
      Eigen::Vector3d walk_end = vox_end;
      if (bounded_ && path_estimator_ == PathEstimator::Crossing && stem_diameter_ > 0.0)
      {
        // extend the ray through the mean chord of the stem's widened outline beyond the return (horizontally). This
        // is pi w / 4 from the outline, or pi (w^2 + d^2) / (8 w) from the physical surface or closest approach
        const double hit_range = range_offset_ + (end - start).norm();
        const double width = stem_width_.enabled() ? stem_width_.width(hit_range) : stem_diameter_;
        const double chord =
          range_to_outline_ ?
            0.25 * kPi * width :
            (width > 0.0 ? 0.125 * kPi * (width * width + stem_diameter_ * stem_diameter_) / width : 0.0);
        const double extension = std::min(1.0, chord / ray_g_);
        Eigen::Vector3d ext_start = start;
        Eigen::Vector3d ext_end = end + (end - start).normalized() * extension;
        if (bounds_.clipRay(ext_start, ext_end, 1e-10))
        {
          walk_end = (ext_end - bounds_.min_bound_) / voxel_width_;
        }
      }
      walkGrid(vox_start, walk_end, *this);
    }
  };
  return Cloud::read(file_name, calculate);
}

void DensityGrid::enableOverlapCorrection()
{
  if (distribution_ != ElementDistribution::VerticalCylinder || stem_diameter_ <= 0.0)
  {
    return;
  }
  overlap_path_.assign(voxels_.size(), 0.0f);
}

void DensityGrid::applyOverlapCorrection()
{
  const bool pooled = path_estimator_ == PathEstimator::PooledBeerLambert;
  if (overlap_path_.empty() && !pooled)
  {
    return;
  }
  // the neighbourhood is the (2k+1) x (2k+1) voxels centred on each voxel, horizontally, to estimate the local
  // density with less noise than a single voxel
  const int radius = std::max(0, static_cast<int>(std::round(pool_radius_ / voxel_width_)));
  std::cout << "pooling over neighbourhoods of " << 2 * radius + 1 << " x " << 2 * radius + 1 << " voxels (+-"
            << (radius + 0.5) * voxel_width_ << " m) for a pool radius of " << pool_radius_ << " m" << std::endl;
  const bool correct = !overlap_path_.empty() && path_estimator_ != PathEstimator::Crossing;
  StemResponseTable table(min_gap_ratio_, range_to_outline_);
  size_t num_corrected = 0, num_saturated = 0;
  double sum_factor = 0.0, sum_k_path = 0.0, sum_path = 0.0;
  // per-layer summed-area tables of hits, path and rho-weighted path, so each neighbourhood sum is O(1). The
  // neighbourhoods are horizontal, so each layer can be updated in place once its tables are built
  const int nx = voxel_dims_[0], ny = voxel_dims_[1];
  const size_t table_size = static_cast<size_t>(nx + 1) * static_cast<size_t>(ny + 1);
  std::vector<double> sum_hits(table_size), sum_paths(table_size), sum_k_paths(table_size);
  auto at = [&](int x, int y) { return static_cast<size_t>(x) + static_cast<size_t>(nx + 1) * static_cast<size_t>(y); };
  auto box = [&](const std::vector<double> &sat, int x0, int y0, int x1, int y1) {
    // sum over x0 <= x < x1, y0 <= y < y1
    return sat[at(x1, y1)] - sat[at(x0, y1)] - sat[at(x1, y0)] + sat[at(x0, y0)];
  };
  for (int z = 0; z < voxel_dims_[2]; z++)
  {
    for (int y = 0; y < ny; y++)
    {
      double row_hits = 0.0, row_path = 0.0, row_k_path = 0.0;
      for (int x = 0; x < nx; x++)
      {
        const int ind = getIndex(Eigen::Vector3i(x, y, z));
        row_hits += voxels_[ind].numHits();
        row_path += voxels_[ind].pathLength();
        row_k_path += overlap_path_.empty() ? 0.0 : overlap_path_[static_cast<size_t>(ind)];
        sum_hits[at(x + 1, y + 1)] = sum_hits[at(x + 1, y)] + row_hits;
        sum_paths[at(x + 1, y + 1)] = sum_paths[at(x + 1, y)] + row_path;
        sum_k_paths[at(x + 1, y + 1)] = sum_k_paths[at(x + 1, y)] + row_k_path;
      }
    }
    const double path_tolerance = 1e-9 * sum_paths[at(nx, ny)] + 1e-12;
    const double hits_tolerance = 1e-9 * sum_hits[at(nx, ny)] + 1e-6;
    for (int y = 0; y < ny; y++)
    {
      for (int x = 0; x < nx; x++)
      {
        const int ind = getIndex(Eigen::Vector3i(x, y, z));
        if (pooled ? voxels_[ind].numRays() <= 0.0f : voxels_[ind].numHits() <= 0.0f)
        {
          continue;
        }
        const int x0 = std::max(0, x - radius), x1 = std::min(nx, x + radius + 1);
        const int y0 = std::max(0, y - radius), y1 = std::min(ny, y + radius + 1);
        // the differences of the sums have rounding errors relative to the layer's totals, so ignore residuals
        double hits = box(sum_hits, x0, y0, x1, y1);
        const double path = box(sum_paths, x0, y0, x1, y1);
        const double k_path = box(sum_k_paths, x0, y0, x1, y1);
        if (path <= path_tolerance)
        {
          continue;
        }
        hits = hits > hits_tolerance ? hits : 0.0;
        const double measured = hits / path;  // local measured density
        double factor = 1.0;
        if (correct && k_path > 0.0 && measured > 0.0)
        {
          // measured = a R(pi d a / 4, rho) for true density a = N d
          const double rho = k_path / path;  // the (G rho L)-weighted mean ratio of effective width to diameter
          bool saturated = false;
          const double a =
            invertStemResponse(measured, rho, stem_diameter_, min_gap_ratio_, max_solid_fraction_, table, saturated);
          num_saturated += saturated ? 1 : 0;
          factor = a / measured;
          num_corrected++;
          sum_factor += factor;
          sum_k_path += overlap_path_[static_cast<size_t>(ind)];
          sum_path += voxels_[ind].pathLength();
        }
        // for the pooled estimator, the voxel takes the neighbourhood's density, keeping its own path length
        voxels_[ind].numHits() =
          static_cast<float>(factor * (pooled ? measured * voxels_[ind].pathLength() : voxels_[ind].numHits()));
      }
    }
  }
  pooled_densities_ = pooled;
  if (num_corrected > 0 && sum_path > 0.0)
  {
    std::cout << "stem occlusion correction: mean factor " << sum_factor / static_cast<double>(num_corrected)
              << " over " << num_corrected << " voxels, with mean effective width / diameter of "
              << sum_k_path / sum_path << ", returns reported at the "
              << (range_to_outline_ ? "widened outline" : "physical stem surface") << std::endl;
    if (num_saturated > 0)
    {
      const double scale = (1.0 + min_gap_ratio_) * (1.0 + min_gap_ratio_);
      std::cout << "Warning: " << num_saturated << " voxels are beyond the maximum solid fraction ("
                << max_solid_fraction_ / scale << ") for " << stem_diameter_ << " m stems with gaps of "
                << min_gap_ratio_
                << " diameters, of N*d = " << max_solid_fraction_ * 4.0 / (kPi * stem_diameter_ * scale)
                << " m^-1. Their densities are set to this limit. Check the stem diameter" << std::endl;
    }
  }
}

// When the cloud has a sharp change in density at the top (e.g. grass or wheat field) between air and the crop, then
// this function adjusts the density estimation based on this two-phase density, rather than assuming the top voxel is
// uniform density
void DensityGrid::flatTopCompensation()
{
  for (int x = 0; x < voxel_dims_[0]; x++)
  {
    for (int y = 0; y < voxel_dims_[1]; y++)
    {
      double p = peaks_[x + voxel_dims_[0] * y];
      if (p < -10000.0)
        continue;
      int z = (int)std::floor(p);
      if (z < 0 || z >= voxel_dims_[2])
        continue;
      int ind = getIndex(Eigen::Vector3i(x, y, z));

      float r = voxels_[ind].numRays();
      if (r > 2.0f)
        voxels_[ind].numHits() *=
          (r - 2.0f) / (r - 1.0f);  // unbiased estimator reduces density estimation slightly in a manner that is
                                    // equivalent to choosing a grass height slightly above the peak point height
    }
  }
}

// This is a form of windowed average over the Moore neighbourhood (3x3x3) window.
void DensityGrid::addNeighbourPriors()
{
#if DENSITY_MIN_RAYS > 0
  const int X = 1;
  const int Y = voxel_dims_[0];
  const int Z = voxel_dims_[0] * voxel_dims_[1];
  DensityGrid::Voxel neighbours;
  double num_hit_points = 0.0;
  double num_hit_points_unsatisfied = 0.0;

  // This simple 3x3x3 convolution needs to be a bit sneaky to avoid having to double the memory cost.
  // well, not that sneaky, we just shift the output -1,-1,-1 for each cell
  for (int x = 1; x < voxel_dims_[0] - 1; x++)
  {
    for (int y = 1; y < voxel_dims_[1] - 1; y++)
    {
      for (int z = 1; z < voxel_dims_[2] - 1; z++)
      {
        const int ind = getIndex(Eigen::Vector3i(x, y, z));
        if (voxels_[ind].numHits() > 0)
          num_hit_points++;
        float needed = DENSITY_MIN_RAYS - voxels_[ind].numRays();
        const DensityGrid::Voxel corner_vox = voxels_[ind - X - Y - Z];
        voxels_[ind - X - Y - Z] = voxels_[ind];  // move centre up to corner
        DensityGrid::Voxel &voxel = voxels_[ind - X - Y - Z];
        if (needed < 0.0)
          continue;
        neighbours = voxels_[ind - X];
        neighbours += voxels_[ind + X];
        neighbours += voxels_[ind - Y];
        neighbours += voxels_[ind + Y];
        neighbours += voxels_[ind - Z];
        neighbours += voxels_[ind + Z];
        if (neighbours.numRays() >= needed)
        {
          voxel += neighbours * (needed / neighbours.numRays());  // add minimal amount to reach DENSITY_MIN_RAYS
          continue;
        }
        voxel += neighbours;
        needed -= neighbours.numRays();

        neighbours = voxels_[ind - X - Y];
        neighbours += voxels_[ind - X + Y];
        neighbours += voxels_[ind + X - Y];
        neighbours += voxels_[ind + X + Y];

        neighbours += voxels_[ind - X - Z];
        neighbours += voxels_[ind - X + Z];
        neighbours += voxels_[ind + X - Z];
        neighbours += voxels_[ind + X + Z];

        neighbours += voxels_[ind - Y - Z];
        neighbours += voxels_[ind - Y + Z];
        neighbours += voxels_[ind + Y - Z];
        neighbours += voxels_[ind + Y + Z];
        if (neighbours.numRays() >= needed)
        {
          voxel += neighbours * (needed / neighbours.numRays());  // add minimal amount to reach DENSITY_MIN_RAYS
          continue;
        }
        voxel += neighbours;
        needed -= neighbours.numRays();

        neighbours = corner_vox;
        neighbours += voxels_[ind - X - Y + Z];
        neighbours += voxels_[ind - X + Y - Z];
        neighbours += voxels_[ind + X - Y - Z];
        neighbours += voxels_[ind - X + Y + Z];
        neighbours += voxels_[ind + X - Y + Z];
        neighbours += voxels_[ind + X + Y - Z];
        neighbours += voxels_[ind + X + Y + Z];
        if (neighbours.numRays() >= needed)
        {
          voxel += neighbours * (needed / neighbours.numRays());  // add minimal amount to reach DENSITY_MIN_RAYS
          continue;
        }
        voxel += neighbours;
        if (voxels_[ind].numHits() > 0)
          num_hit_points_unsatisfied++;
      }
    }
  }
  const double percentage = 100.0 * num_hit_points_unsatisfied / num_hit_points;
  std::cout << "Density calculation: " << percentage << "% of voxels had insufficient (<" << DENSITY_MIN_RAYS
            << ") rays within them" << std::endl;
  if (percentage > 50.0)
  {
    std::cout << "This is high. Consider using a larger pixel size, or a denser cloud, or reducing DENSITY_MIN_RAYS, "
                 "for consistent results"
              << std::endl;
  }
  else if (percentage < 1.0)
  {
    std::cout << "This is low enough that you could get more fidelity from using a smaller pixel size" << std::endl;
    std::cout << "or more accuracy by increasing DENSITY_MIN_RAYS" << std::endl;
  }
#endif
}

namespace
{
/// The fraction of voxel x, y, z (shifted indices) that the density estimate applies to: with flat top compensation,
/// the top voxel of each column is estimated for its part below the column's highest point, so its density is scaled
/// by this fraction to give the voxel average, as for the rendered image. Otherwise 1
double flatTopFraction(const DensityGrid &grid, int x, int y, int z, const Eigen::Vector3i &full_dims, int shift)
{
#if defined FLAT_TOP_COMPENSATION
  const double p = grid.peaks_[static_cast<size_t>((x + shift) + full_dims[0] * (y + shift))];
  const int top_z = p < -10000.0 ? -1 : static_cast<int>(std::floor(p));
  if (z == top_z - shift)
  {
    return p - static_cast<double>(top_z);
  }
#else
  RAYLIB_UNUSED(grid);
  RAYLIB_UNUSED(x);
  RAYLIB_UNUSED(y);
  RAYLIB_UNUSED(z);
  RAYLIB_UNUSED(full_dims);
  RAYLIB_UNUSED(shift);
#endif
  return 1.0;
}
}  // namespace

/// The density of a voxel at valid (shifted) index x,y,z, as an average over the voxel
double voxelAverageDensity(const DensityGrid &grid, const DensityGrid::Voxel &voxel, int x, int y, int z,
                           const Eigen::Vector3i &full_dims, int shift)
{
  return grid.voxelDensity(voxel) * flatTopFraction(grid, x, y, z, full_dims, shift);
}

std::unique_ptr<DensityGrid> estimateDensityGrid(const std::string &cloud_file, const Cuboid &bounds,
                                                 double voxel_width, const Eigen::Vector3i &dims_in,
                                                 const DensityOptions &options, int &shift)
{
  Eigen::Vector3i dims = dims_in;
  Cuboid grid_bounds = bounds;
  shift = 0;
#if DENSITY_MIN_RAYS > 0
  shift = 1;
  dims += 2 * Eigen::Vector3i(shift, shift, shift);  // so that we have extra space to convolve
  grid_bounds.min_bound_ -= Eigen::Vector3d(voxel_width, voxel_width, voxel_width);
#endif
  const bool cylinders = options.distribution == ElementDistribution::VerticalCylinder;
  auto grid_ptr = std::make_unique<DensityGrid>(grid_bounds, voxel_width, dims, options.distribution,
                                                cylinders ? options.min_sin_zenith : 0.0);
  DensityGrid &grid = *grid_ptr;
  grid.setPathEstimator(options.path_estimator);
  grid.setPoolRadius(options.pool_radius);
  grid.setMinGapRatio(options.min_gap_ratio);
  if (cylinders)
  {
    grid.setMaxG(options.max_sin_zenith);
  }
  if (cylinders && options.stem_diameter > 0.0)
  {
    grid.setStemDiameter(options.stem_diameter);
  }
  if (cylinders && options.footprint.enabled())
  {
    if (options.stem_diameter <= 0.0)
    {
      std::cerr << "Error: the beam footprint correction requires the stem diameter" << std::endl;
      return nullptr;
    }
    grid.setFootprint(options.footprint, options.stem_diameter);
    std::cout << "correcting for beam footprint, effective stem width of " << options.stem_diameter << " m stems is "
              << effectiveStemWidth(options.stem_diameter, options.footprint.diameter(1.0), options.footprint.threshold)
              << " m at 1 m range, "
              << effectiveStemWidth(options.stem_diameter, options.footprint.diameter(5.0), options.footprint.threshold)
              << " m at 5 m range" << std::endl;
  }

#if defined FLAT_TOP_COMPENSATION
  if (options.flat_top_compensation)  // without the peaks, no compensation is applied
  {
    if (!grid.calculatePeaks(cloud_file))
    {
      return nullptr;
    }
  }
#endif
  if (cylinders && options.overlap_correction && options.stem_diameter > 0.0)
  {
    grid.enableOverlapCorrection();
  }
  if (!grid.calculateDensities(cloud_file))
  {
    return nullptr;
  }
  if (grid.numIgnoredRays() > 0)
  {
    std::cout << grid.numIgnoredRays() << " rays ignored outside of the sin zenith angle range "
              << options.min_sin_zenith << " to " << options.max_sin_zenith << std::endl;
  }
#if defined FLAT_TOP_COMPENSATION
  if (options.flat_top_compensation)
  {
    grid.flatTopCompensation();
  }
#endif
  grid.applyOverlapCorrection();
#if DENSITY_MIN_RAYS > 0
  grid.addNeighbourPriors();
#endif
  return grid_ptr;
}

std::vector<DensityCell> aggregateDensityCells(
  const DensityGrid &grid, int shift, const Eigen::Vector3i &voxel_min, const Eigen::Vector3i &voxel_max,
  const CellLayout &layout, const std::function<double(int, int)> &ground_height,
  const std::function<void(const Eigen::Vector3i &, const Eigen::Vector3d &, double)> &observed_voxel)
{
  std::vector<DensityCell> cells(static_cast<size_t>(layout.num_cells[0]) * static_cast<size_t>(layout.num_cells[1]) *
                                 static_cast<size_t>(layout.num_cells[2]));
  const Eigen::Vector3i full_dims = grid.dimensions();
  const double w = grid.voxelWidth();
  // the grid's bounds include the padding of shift voxels, so (shifted) voxel 0, 0, 0 has its centre at
  // bounds().min_bound_ + w (shift + 1/2) in each axis
  const Eigen::Vector3d first_centre =
    grid.bounds().min_bound_ + w * Eigen::Vector3d(shift + 0.5, shift + 0.5, shift + 0.5) - layout.min_bound;
  for (int y = voxel_min[1]; y < voxel_max[1]; y++)
  {
    for (int x = voxel_min[0]; x < voxel_max[0]; x++)
    {
      Eigen::Vector3d centre = first_centre + w * Eigen::Vector3d(x, y, 0.0);
      const int cx = static_cast<int>(std::floor(centre[0] / layout.cell_size[0]));
      const int cy = static_cast<int>(std::floor(centre[1] / layout.cell_size[1]));
      if (cx < 0 || cy < 0 || cx >= layout.num_cells[0] || cy >= layout.num_cells[1])
      {
        continue;
      }
      const double ground = ground_height ? ground_height(x, y) : 0.0;
      for (int z = voxel_min[2]; z < voxel_max[2]; z++)
      {
        centre[2] = first_centre[2] + w * static_cast<double>(z) - ground;
        const int cz = static_cast<int>(std::floor(centre[2] / layout.cell_size[2]));
        if (cz < 0 || cz >= layout.num_cells[2])
        {
          continue;
        }
        DensityCell &cell = cells[layout.index(cx, cy, cz)];
        cell.total++;
        const DensityGrid::Voxel &voxel = grid.voxels()[static_cast<size_t>(grid.getIndex(Eigen::Vector3i(x, y, z)))];
        if (voxel.numRays() <= DensityGrid::min_voxel_hits)
        {
          continue;
        }
        const double fraction = flatTopFraction(grid, x, y, z, full_dims, shift);
        const double density = grid.voxelDensity(voxel) * fraction;
        cell.observed++;
        cell.density_sum += density;
        cell.hits += voxel.numHits() * fraction;
        cell.weighted_path += voxel.pathLength();
        if (observed_voxel)
        {
          observed_voxel(Eigen::Vector3i(x, y, z), centre, density);
        }
      }
    }
  }
  return cells;
}
}  // namespace ray
