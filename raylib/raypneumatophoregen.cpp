// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#include "raypneumatophoregen.h"

#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ray
{
namespace
{
/// Spatial hash of stem positions, to reject candidate stems that would overlap existing ones (a hard-core process)
class StemHash
{
public:
  StemHash(double field_min, double field_width, double diameter)
    : min_(field_min)
    , cell_(std::max(diameter, 1e-4))
    , diameter2_(diameter * diameter)
  {
    cells_ = std::max(1, static_cast<int>(std::ceil(field_width / cell_)));
    buckets_.resize(static_cast<size_t>(cells_) * static_cast<size_t>(cells_));
  }
  bool overlaps(const Eigen::Vector2d &pos) const
  {
    const int bx = index(pos[0]), by = index(pos[1]);
    for (int x = std::max(0, bx - 1); x <= std::min(cells_ - 1, bx + 1); x++)
    {
      for (int y = std::max(0, by - 1); y <= std::min(cells_ - 1, by + 1); y++)
      {
        for (const auto &c : buckets_[static_cast<size_t>(x) + static_cast<size_t>(cells_) * static_cast<size_t>(y)])
        {
          if ((pos - c).squaredNorm() < diameter2_)
          {
            return true;
          }
        }
      }
    }
    return false;
  }
  void add(const Eigen::Vector2d &pos)
  {
    buckets_[static_cast<size_t>(index(pos[0])) + static_cast<size_t>(cells_) * static_cast<size_t>(index(pos[1]))]
      .push_back(pos);
  }

private:
  int index(double x) const { return std::min(cells_ - 1, std::max(0, static_cast<int>((x - min_) / cell_))); }
  double min_, cell_, diameter2_;
  int cells_;
  std::vector<std::vector<Eigen::Vector2d>> buckets_;
};
}  // namespace

void PneumatophoreGen::generateStems(std::mt19937 &gen)
{
  const PneumatophoreParams &p = params_;
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 1.0);
  const double half_field = 0.5 * p.plot_width + p.margin;
  const double field_area = 4.0 * half_field * half_field;
  auto random_height = [&]() {
    double h;
    do
    {
      h = p.mean_height + p.height_sd * normal(gen);
    } while (h < 0.02);
    return h;
  };
  auto random_position = [&]() {
    return Eigen::Vector2d(-half_field + 2.0 * half_field * uniform(gen),
                           -half_field + 2.0 * half_field * uniform(gen));
  };
  // stems as a fraction of ground area, at the (local) density stems_per_m2
  const double solid_fraction = p.stems_per_m2 * kPi * p.diameter * p.diameter / 4.0;
  if (p.hard_core && solid_fraction > max_solid_fraction)
  {
    std::stringstream ss;
    ss << "stems of diameter " << p.diameter << " m cannot be placed without overlap at " << p.stems_per_m2
       << " per m^2: the solid fraction (N pi d^2 / 4) of " << solid_fraction << " exceeds " << max_solid_fraction
       << ", near the random packing limit of 0.547. That is, at most "
       << max_solid_fraction * 4.0 / (kPi * p.diameter * p.diameter) << " stems per m^2";
    throw std::runtime_error(ss.str());
  }
  StemHash hash(-half_field, 2.0 * half_field, p.diameter);
  // accepts a stem unless it overlaps an existing one (when hard_core)
  auto try_add = [&](const Eigen::Vector2d &pos) {
    if (p.hard_core && hash.overlaps(pos))
    {
      return false;
    }
    hash.add(pos);
    stems_.emplace_back(pos[0], pos[1], random_height());
    return true;
  };
  // the number of stems to place in an area: exact (rounded) when hard_core, otherwise Poisson distributed
  auto stem_count = [&](double expected) {
    if (p.hard_core)
    {
      return static_cast<int>(std::round(expected));
    }
    std::poisson_distribution<int> count(expected);
    return count(gen);
  };
  // random sequential addition of a given number of stems, to positions sampled by the given acceptance function
  auto place = [&](int target, const std::function<bool(const Eigen::Vector2d &)> &accept) {
    const long long max_attempts = 1000LL * std::max(target, 1) + 100000LL;
    int placed = 0;
    for (long long attempt = 0; placed < target && attempt < max_attempts; attempt++)
    {
      const Eigen::Vector2d pos = random_position();
      if (accept(pos) && try_add(pos))
      {
        placed++;
      }
    }
    if (placed < target)
    {
      std::cerr << "Warning: only " << placed << " of " << target << " stems could be placed without overlap"
                << std::endl;
    }
  };

  stems_.clear();
  if (p.patchy)
  {
    // a Boolean model of discs, with centre density giving the required ground cover
    const double half_patch_field = half_field + p.patch_radius;
    const double cover = std::min(p.patch_cover, 0.999);
    const double centre_density = -std::log(1.0 - cover) / (kPi * p.patch_radius * p.patch_radius);
    std::poisson_distribution<int> patch_count(centre_density * 4.0 * half_patch_field * half_patch_field);
    const int num_patches = patch_count(gen);
    // bucket the patch centres for a fast inside test
    const double cell = std::max(p.patch_radius, 0.05);
    const int cells = static_cast<int>(std::ceil(2.0 * half_patch_field / cell));
    std::vector<std::vector<Eigen::Vector2d>> patch_buckets(static_cast<size_t>(cells * cells));
    for (int i = 0; i < num_patches; i++)
    {
      const Eigen::Vector2d c(-half_patch_field + 2.0 * half_patch_field * uniform(gen),
                              -half_patch_field + 2.0 * half_patch_field * uniform(gen));
      const int bx = std::min(cells - 1, static_cast<int>((c[0] + half_patch_field) / cell));
      const int by = std::min(cells - 1, static_cast<int>((c[1] + half_patch_field) / cell));
      patch_buckets[static_cast<size_t>(bx + cells * by)].push_back(c);
    }
    auto in_patch = [&](const Eigen::Vector2d &pos) {
      const int bx = static_cast<int>((pos[0] + half_patch_field) / cell);
      const int by = static_cast<int>((pos[1] + half_patch_field) / cell);
      for (int x = std::max(0, bx - 1); x <= std::min(cells - 1, bx + 1); x++)
      {
        for (int y = std::max(0, by - 1); y <= std::min(cells - 1, by + 1); y++)
        {
          for (const auto &c : patch_buckets[static_cast<size_t>(x + cells * y)])
          {
            if ((pos - c).squaredNorm() < p.patch_radius * p.patch_radius)
            {
              return true;
            }
          }
        }
      }
      return false;
    };
    // measure the patch area within the field, so that the density within patches is as specified
    const double sample_width = std::min(0.02, 0.1 * p.patch_radius);
    const int samples = static_cast<int>(std::ceil(2.0 * half_field / sample_width));
    long long inside = 0;
    for (int i = 0; i < samples; i++)
    {
      for (int j = 0; j < samples; j++)
      {
        const Eigen::Vector2d pos(-half_field + sample_width * (i + 0.5), -half_field + sample_width * (j + 0.5));
        inside += in_patch(pos) ? 1 : 0;
      }
    }
    const double patch_area = field_area * static_cast<double>(inside) / (static_cast<double>(samples) * samples);
    place(stem_count(p.stems_per_m2 * patch_area), in_patch);
    if (p.matrix_ratio > 0.0)
    {
      place(stem_count(p.matrix_ratio * p.stems_per_m2 * (field_area - patch_area)),
            [&](const Eigen::Vector2d &pos) { return !in_patch(pos); });
    }
    return;
  }
  if (!p.clustered)
  {
    // complete spatial randomness, or (when hard_core) random sequential addition of non-overlapping stems
    place(stem_count(p.stems_per_m2 * field_area), [](const Eigen::Vector2d &) { return true; });
    return;
  }
  // clustered: stems lie along cable roots radiating from trees. Trees are generated beyond the field so that roots
  // from outside reach into it, keeping the mean density uniform
  const double half_tree_field = half_field + p.root_length;
  const double tree_area = 4.0 * half_tree_field * half_tree_field;
  std::poisson_distribution<int> tree_count(p.trees_per_m2 * tree_area);
  const int num_trees = tree_count(gen);
  const double stems_per_metre_of_root = p.stems_per_m2 / (p.trees_per_m2 * p.roots_per_tree * p.root_length);
  int num_failed = 0;
  for (int t = 0; t < num_trees; t++)
  {
    const Eigen::Vector2d tree(-half_tree_field + 2.0 * half_tree_field * uniform(gen),
                               -half_tree_field + 2.0 * half_tree_field * uniform(gen));
    const double angle_offset = 2.0 * kPi * uniform(gen);
    for (int r = 0; r < p.roots_per_tree; r++)
    {
      // evenly spread root directions, with some jitter
      const double angle =
        angle_offset + 2.0 * kPi * (static_cast<double>(r) + 0.3 * (uniform(gen) - 0.5)) / p.roots_per_tree;
      const Eigen::Vector2d along(std::cos(angle), std::sin(angle));
      const Eigen::Vector2d across(-along[1], along[0]);
      const int num = stem_count(stems_per_metre_of_root * p.root_length);
      for (int i = 0; i < num; i++)
      {
        bool placed = false;
        for (int attempt = 0; attempt < 100 && !placed; attempt++)
        {
          const Eigen::Vector2d pos =
            tree + along * (p.root_length * uniform(gen)) + across * (p.lateral_spread * normal(gen));
          if (std::abs(pos[0]) >= half_field || std::abs(pos[1]) >= half_field)
          {
            placed = true;  // outside of the field, so not needed
          }
          else
          {
            placed = try_add(pos);
          }
        }
        num_failed += placed ? 0 : 1;
      }
    }
  }
  if (num_failed > 0)
  {
    std::cerr << "Warning: " << num_failed << " stems could not be placed along cable roots without overlap"
              << std::endl;
  }
}

void PneumatophoreGen::generate(const PneumatophoreParams &params, unsigned int seed)
{
  params_ = params;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  generateStems(gen);
  stem_width_ = StemWidthModel(params_.diameter, params_.footprint);

  // build the bucket grid. Each stem is added to every bucket that its maximum effective radius overlaps
  const double half_field = 0.5 * params_.plot_width + params_.margin;
  field_min_ = -half_field;
  num_buckets_ = static_cast<int>(std::ceil(2.0 * half_field / bucket_width_));
  buckets_.assign(static_cast<size_t>(num_buckets_ * num_buckets_), std::vector<int>());
  const double max_radius = 0.5 * (stem_width_.enabled() ? stem_width_.maxWidth(params_.max_range) : params_.diameter);
  max_height_ = 0.0;
  for (size_t i = 0; i < stems_.size(); i++)
  {
    const Eigen::Vector3d &stem = stems_[i];
    max_height_ = std::max(max_height_, stem[2]);
    const int x0 = std::max(0, static_cast<int>((stem[0] - max_radius - field_min_) / bucket_width_));
    const int x1 = std::min(num_buckets_ - 1, static_cast<int>((stem[0] + max_radius - field_min_) / bucket_width_));
    const int y0 = std::max(0, static_cast<int>((stem[1] - max_radius - field_min_) / bucket_width_));
    const int y1 = std::min(num_buckets_ - 1, static_cast<int>((stem[1] + max_radius - field_min_) / bucket_width_));
    for (int x = x0; x <= x1; x++)
    {
      for (int y = y0; y <= y1; y++)
      {
        buckets_[static_cast<size_t>(x + num_buckets_ * y)].push_back(static_cast<int>(i));
      }
    }
  }

  // a lawnmower walking path within the plot
  const double half_plot = 0.5 * params_.plot_width;
  const int num_lines = std::max(1, static_cast<int>(params_.plot_width / params_.walk_spacing));
  const double path_length = params_.plot_width * num_lines;
  cloud_.clear();
  const RGBA stem_colour(40, 160, 40, 255), ground_colour(120, 80, 40, 255), none_colour(128, 128, 128, 0);
  int num_stem_hits = 0, num_ground_hits = 0;
  for (int i = 0; i < params_.num_rays; i++)
  {
    const double s = path_length * (static_cast<double>(i) + 0.5) / static_cast<double>(params_.num_rays);
    const int line = std::min(num_lines - 1, static_cast<int>(s / params_.plot_width));
    double along = s - params_.plot_width * line;
    if (line % 2 == 1)
    {
      along = params_.plot_width - along;
    }
    const double y = -half_plot + params_.walk_spacing * (static_cast<double>(line) + 0.5);
    const Eigen::Vector3d start(-half_plot + along, y,
                                params_.scanner_height + 0.05 * std::sin(10.0 * s));  // some walking bob

    // a rotating hand-held scanner covers most directions. Directions are uniform over the band of zenith angles
    const double cos_min = std::cos(params_.min_zenith_deg * kPi / 180.0);
    const double cos_max = std::cos(params_.max_zenith_deg * kPi / 180.0);
    const double dz = -(cos_max + (cos_min - cos_max) * uniform(gen));
    const double phi = 2.0 * kPi * uniform(gen);
    const double horizontal = std::sqrt(std::max(0.0, 1.0 - dz * dz));
    const Eigen::Vector3d dir(horizontal * std::cos(phi), horizontal * std::sin(phi), dz);

    double length;
    bool hit_stem, hit_ground;
    castRay(start, dir, length, hit_stem, hit_ground);
    // stems green, ground brown, and rays with no return (at maximum range) grey
    RGBA colour = hit_stem ? stem_colour : (hit_ground ? ground_colour : none_colour);
    if (!hit_stem && (!hit_ground || params_.ground_unbounded))
    {
      colour.alpha = 0;
    }
    num_stem_hits += hit_stem ? 1 : 0;
    num_ground_hits += hit_ground ? 1 : 0;
    cloud_.addRay(start, start + dir * length, 1e-4 * static_cast<double>(i), colour);
  }
  std::cout << "generated " << stems_.size() << " stems and " << params_.num_rays << " rays, of which " << num_stem_hits
            << " hit stems and " << num_ground_hits << " hit the ground" << std::endl;
}

void PneumatophoreGen::castRay(const Eigen::Vector3d &start, const Eigen::Vector3d &dir, double &length, bool &hit_stem,
                               bool &hit_ground) const
{
  hit_stem = hit_ground = false;
  const double inf = std::numeric_limits<double>::max();
  const double t_ground = dir[2] < 0.0 ? start[2] / -dir[2] : inf;
  double t_end = std::min(params_.max_range, t_ground);
  length = t_end;
  hit_ground = t_ground <= params_.max_range;

  // the part of the ray below the tallest stem
  double t_lo = 0.0;
  if (start[2] > max_height_)
  {
    if (dir[2] >= 0.0)
    {
      return;
    }
    t_lo = (start[2] - max_height_) / -dir[2];
  }
  if (t_lo >= t_end)
  {
    return;
  }
  const Eigen::Vector2d s2(start[0], start[1]);
  const Eigen::Vector2d d2(dir[0], dir[1]);
  const double h2 = d2.squaredNorm();
  const double h = std::sqrt(h2);

  double best_t = inf;
  auto test_bucket = [&](int bx, int by) {
    for (const int index : buckets_[static_cast<size_t>(bx + num_buckets_ * by)])
    {
      const Eigen::Vector3d &stem = stems_[static_cast<size_t>(index)];
      const Eigen::Vector2d rel = Eigen::Vector2d(stem[0], stem[1]) - s2;
      const double along = h2 > 1e-12 ? rel.dot(d2) / h2 : 0.0;  // range of closest approach
      const double offset2 = std::max(0.0, rel.squaredNorm() - (h2 > 1e-12 ? along * along * h2 : 0.0));
      const double radius = 0.5 * (stem_width_.enabled() ? stem_width_.width(std::max(0.0, along)) : params_.diameter);
      if (offset2 >= radius * radius)
      {
        continue;
      }
      const double half_chord = h > 1e-6 ? std::sqrt(radius * radius - offset2) / h : inf;
      const double t_in = std::max(0.0, along - half_chord);
      const double t_out = along + half_chord;
      const double z_in = start[2] + t_in * dir[2];
      double t_hit = inf;
      if (z_in >= 0.0 && z_in <= stem[2])
      {
        t_hit = t_in;  // hits the side of the stem
      }
      else if (z_in > stem[2] && dir[2] < 0.0)
      {
        const double t_top = (start[2] - stem[2]) / -dir[2];
        if (t_top <= t_out)
        {
          t_hit = t_top;  // hits the top of the stem
        }
      }
      if (t_hit >= t_lo - 1e-9 && t_hit < best_t)
      {
        best_t = t_hit;
      }
    }
  };

  // walk the bucket grid along the ray, from t_lo to t_end, stopping once a hit is found before the bucket exit
  const Eigen::Vector2d p0 = s2 + d2 * t_lo;
  int bx = static_cast<int>(std::floor((p0[0] - field_min_) / bucket_width_));
  int by = static_cast<int>(std::floor((p0[1] - field_min_) / bucket_width_));
  const int step_x = d2[0] > 0.0 ? 1 : -1;
  const int step_y = d2[1] > 0.0 ? 1 : -1;
  const double next_x = field_min_ + bucket_width_ * (bx + (step_x > 0 ? 1 : 0));
  const double next_y = field_min_ + bucket_width_ * (by + (step_y > 0 ? 1 : 0));
  double t_max_x = std::abs(d2[0]) > 1e-12 ? t_lo + (next_x - p0[0]) / d2[0] : inf;
  double t_max_y = std::abs(d2[1]) > 1e-12 ? t_lo + (next_y - p0[1]) / d2[1] : inf;
  const double t_delta_x = std::abs(d2[0]) > 1e-12 ? bucket_width_ / std::abs(d2[0]) : inf;
  const double t_delta_y = std::abs(d2[1]) > 1e-12 ? bucket_width_ / std::abs(d2[1]) : inf;
  while (bx >= 0 && by >= 0 && bx < num_buckets_ && by < num_buckets_)
  {
    test_bucket(bx, by);
    const double t_exit = std::min(t_max_x, t_max_y);
    if (best_t <= t_exit || t_exit >= t_end)
    {
      break;
    }
    if (t_max_x < t_max_y)
    {
      bx += step_x;
      t_max_x += t_delta_x;
    }
    else
    {
      by += step_y;
      t_max_y += t_delta_y;
    }
  }
  if (best_t < t_end)
  {
    length = best_t;
    hit_stem = true;
    hit_ground = false;
  }
}

bool PneumatophoreGen::writeTruthProfile(const std::string &file_name, double layer_width) const
{
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing" << std::endl;
    return false;
  }
  const double half_plot = 0.5 * params_.plot_width;
  const double area = params_.plot_width * params_.plot_width;
  const int num_layers = static_cast<int>(std::ceil(max_height_ / layer_width)) + 1;
  std::vector<double> frontal_area(static_cast<size_t>(num_layers), 0.0);
  int num_stems = 0;
  for (const auto &stem : stems_)
  {
    if (std::abs(stem[0]) >= half_plot || std::abs(stem[1]) >= half_plot)
    {
      continue;
    }
    num_stems++;
    for (int i = 0; i < num_layers; i++)
    {
      const double overlap = std::max(0.0, std::min(stem[2], layer_width * (i + 1)) - layer_width * i);
      frontal_area[static_cast<size_t>(i)] += params_.diameter * overlap;
    }
  }
  ofs << "z_min,z_max,true_frontal_area_per_volume_Nd" << std::endl;
  for (int i = 0; i < num_layers; i++)
  {
    ofs << layer_width * i << "," << layer_width * (i + 1) << ","
        << frontal_area[static_cast<size_t>(i)] / (area * layer_width) << std::endl;
  }
  std::cout << "truth profile from " << num_stems << " stems in the " << params_.plot_width << " m plot ("
            << num_stems / area << " stems per m^2) written to " << file_name << std::endl;
  return true;
}
bool PneumatophoreGen::writeStems(const std::string &file_name) const
{
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing" << std::endl;
    return false;
  }
  ofs << "x,y,height,diameter" << std::endl;
  for (const auto &stem : stems_)
  {
    ofs << stem[0] << "," << stem[1] << "," << stem[2] << "," << params_.diameter << std::endl;
  }
  return true;
}

bool PneumatophoreGen::writeTruthCells(const std::string &file_name, const Eigen::Vector3d &cell_size) const
{
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing" << std::endl;
    return false;
  }
  const double half_plot = 0.5 * params_.plot_width;
  const int nx = static_cast<int>(std::floor(params_.plot_width / cell_size[0] + 1e-9));
  const int ny = static_cast<int>(std::floor(params_.plot_width / cell_size[1] + 1e-9));
  const int nz = static_cast<int>(std::ceil(max_height_ / cell_size[2]));
  std::vector<double> frontal_area(static_cast<size_t>(nx * ny * nz), 0.0);
  for (const auto &stem : stems_)
  {
    const int x = static_cast<int>(std::floor((stem[0] + half_plot) / cell_size[0]));
    const int y = static_cast<int>(std::floor((stem[1] + half_plot) / cell_size[1]));
    if (x < 0 || y < 0 || x >= nx || y >= ny)
    {
      continue;
    }
    for (int z = 0; z < nz; z++)
    {
      const double overlap =
        std::max(0.0, std::min(stem[2], cell_size[2] * (z + 1)) - cell_size[2] * static_cast<double>(z));
      frontal_area[static_cast<size_t>(x + nx * (y + ny * z))] += params_.diameter * overlap;
    }
  }
  const double volume = cell_size[0] * cell_size[1] * cell_size[2];
  ofs << "x_min,y_min,z_min,x_max,y_max,z_max,true_frontal_area_per_volume_Nd" << std::endl;
  for (int z = 0; z < nz; z++)
  {
    for (int y = 0; y < ny; y++)
    {
      for (int x = 0; x < nx; x++)
      {
        const Eigen::Vector3d min_bound(-half_plot + cell_size[0] * x, -half_plot + cell_size[1] * y, cell_size[2] * z);
        const Eigen::Vector3d max_bound = min_bound + cell_size;
        ofs << min_bound[0] << "," << min_bound[1] << "," << min_bound[2] << "," << max_bound[0] << "," << max_bound[1]
            << "," << max_bound[2] << "," << frontal_area[static_cast<size_t>(x + nx * (y + ny * z))] / volume
            << std::endl;
      }
    }
  }
  return true;
}
}  // namespace ray
