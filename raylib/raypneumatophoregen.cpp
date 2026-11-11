// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#include "raypneumatophoregen.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace ray
{
namespace
{
/// Random numbers drawn from std::mt19937, whose output sequence is specified by the C++ standard. The std::
/// distributions are implementation-defined, so would give different scenes with different standard libraries
class SceneRandom
{
public:
  explicit SceneRandom(std::mt19937 &gen)
    : gen_(gen)
  {}
  /// uniform in the open interval (0, 1)
  double uniform() { return (static_cast<double>(gen_()) + 0.5) / 4294967296.0; }
  /// standard normal, by the Box-Muller transform. Values are generated in pairs
  double normal()
  {
    if (has_spare_)
    {
      has_spare_ = false;
      return spare_;
    }
    const double radius = std::sqrt(-2.0 * std::log(uniform()));
    const double angle = 2.0 * kPi * uniform();
    spare_ = radius * std::sin(angle);
    has_spare_ = true;
    return radius * std::cos(angle);
  }
  /// Poisson distributed count of the given mean, by inversion. Larger means are split into parts of mean at most 64,
  /// so that exp(-mean) can't underflow. The sum of independent Poisson counts is Poisson, with the total mean
  int poisson(double mean)
  {
    int count = 0;
    while (mean > 0.0)
    {
      const double part = std::min(mean, 64.0);
      mean -= part;
      const double u = uniform();
      double probability = std::exp(-part);
      double cumulative = probability;
      int k = 0;
      // probability > 0 ends the search if rounding leaves the total probability just below u
      while (u > cumulative && probability > 0.0)
      {
        k++;
        probability *= part / static_cast<double>(k);
        cumulative += probability;
      }
      count += k;
    }
    return count;
  }

private:
  std::mt19937 &gen_;
  double spare_ = 0.0;
  bool has_spare_ = false;
};

/// A uniform 2D grid of buckets of items (such as positions or indices) over the square [min, min + width)^2, for
/// finding the items near a position. Positions beyond the grid fall in its edge buckets
template <class T>
class BucketGrid2D
{
public:
  BucketGrid2D(double min, double width, double bucket_width)
    : min_(min)
    , bucket_width_(bucket_width)
    , num_(std::max(1, static_cast<int>(std::ceil(width / bucket_width))))
    , buckets_(static_cast<size_t>(num_) * static_cast<size_t>(num_))
  {}
  /// the bucket index of a coordinate, along either axis
  int index(double x) const
  {
    return std::min(num_ - 1, std::max(0, static_cast<int>(std::floor((x - min_) / bucket_width_))));
  }
  /// add the item to the bucket containing pos
  void add(const Eigen::Vector2d &pos, const T &item) { bucket(index(pos[0]), index(pos[1])).push_back(item); }
  /// add the item to every bucket that overlaps the square of half-width @c radius centred on pos
  void addWithin(const Eigen::Vector2d &pos, double radius, const T &item)
  {
    for (int x = index(pos[0] - radius); x <= index(pos[0] + radius); x++)
    {
      for (int y = index(pos[1] - radius); y <= index(pos[1] + radius); y++)
      {
        bucket(x, y).push_back(item);
      }
    }
  }
  /// whether @c predicate is true for any item in the 3 x 3 buckets around pos, which include every item added
  /// within one bucket width of pos
  template <class Predicate>
  bool anyNear(const Eigen::Vector2d &pos, const Predicate &predicate) const
  {
    const int bx = index(pos[0]), by = index(pos[1]);
    for (int x = std::max(0, bx - 1); x <= std::min(num_ - 1, bx + 1); x++)
    {
      for (int y = std::max(0, by - 1); y <= std::min(num_ - 1, by + 1); y++)
      {
        for (const auto &item : bucket(x, y))
        {
          if (predicate(item))
          {
            return true;
          }
        }
      }
    }
    return false;
  }
  /// the number of buckets along each axis
  inline int size() const { return num_; }
  /// the buckets, indexed by x + size() * y
  inline std::vector<std::vector<T>> &buckets() { return buckets_; }

private:
  inline size_t id(int x, int y) const
  {
    return static_cast<size_t>(x) + static_cast<size_t>(num_) * static_cast<size_t>(y);
  }
  inline std::vector<T> &bucket(int x, int y) { return buckets_[id(x, y)]; }
  inline const std::vector<T> &bucket(int x, int y) const { return buckets_[id(x, y)]; }
  double min_, bucket_width_;
  int num_;
  std::vector<std::vector<T>> buckets_;
};
}  // namespace

void PneumatophoreGen::generateStems(std::mt19937 &gen)
{
  const PneumatophoreParams &p = params_;
  SceneRandom rng(gen);
  const double half_field = 0.5 * p.plot_width + p.margin;
  const double field_area = 4.0 * half_field * half_field;
  auto random_height = [&]() {
    double h;
    do
    {
      h = p.mean_height + p.height_sd * rng.normal();
    } while (h < 0.02);
    return h;
  };
  auto random_position = [&]() {
    return Eigen::Vector2d(-half_field + 2.0 * half_field * rng.uniform(),
                           -half_field + 2.0 * half_field * rng.uniform());
  };
  // stems as a fraction of ground area, at the (local) density stems_per_m2
  const double solid_fraction = p.stems_per_m2 * kPi * p.diameter * p.diameter / 4.0;
  // stem centres are at least this far apart
  const double separation = p.diameter * (1.0 + std::max(0.0, p.min_gap_ratio));
  const double separation_scale = (separation / p.diameter) * (separation / p.diameter);
  if (p.hard_core && solid_fraction * separation_scale > max_solid_fraction)
  {
    std::stringstream ss;
    ss << "stems of diameter " << p.diameter << " m cannot be placed with gaps of at least " << p.min_gap_ratio
       << " diameters at " << p.stems_per_m2 << " per m^2: the solid fraction (N pi d^2 / 4) of " << solid_fraction
       << " exceeds " << max_solid_fraction / separation_scale
       << " (for the exclusion discs, near the random packing limit of 0.547). That is, at most "
       << max_solid_fraction * 4.0 / (kPi * separation * separation) << " stems per m^2";
    throw std::runtime_error(ss.str());
  }
  // the placed stems' positions, to reject candidate stems that would overlap them (a hard-core process)
  BucketGrid2D<Eigen::Vector2d> placed_stems(-half_field, 2.0 * half_field, std::max(separation, 1e-4));
  // sparse trunks, not overlapping each other
  trunks_.clear();
  const int num_trunks = static_cast<int>(std::round(p.trunks_per_m2 * field_area));
  for (int attempt = 0; static_cast<int>(trunks_.size()) < num_trunks && attempt < 1000 * num_trunks; attempt++)
  {
    const Eigen::Vector2d pos = random_position();
    double diameter;
    do
    {
      diameter = p.trunk_diameter + p.trunk_diameter_sd * rng.normal();
    } while (diameter < 0.05);
    bool clear = true;
    for (const auto &t : trunks_)
    {
      clear = clear && (pos - t.head<2>()).norm() > 0.5 * (diameter + t[2]) + 0.1;
    }
    if (clear)
    {
      trunks_.emplace_back(pos[0], pos[1], diameter);
    }
  }
  // accepts a stem unless it overlaps an existing one (when hard_core), or a trunk
  auto try_add = [&](const Eigen::Vector2d &pos) {
    if (p.hard_core && placed_stems.anyNear(pos, [&](const Eigen::Vector2d &stem) {
          return (pos - stem).squaredNorm() < separation * separation;
        }))
    {
      return false;
    }
    for (const auto &t : trunks_)
    {
      if ((pos - t.head<2>()).norm() < 0.5 * (t[2] + p.diameter))
      {
        return false;
      }
    }
    placed_stems.add(pos, pos);
    stems_.emplace_back(pos[0], pos[1], random_height());
    return true;
  };
  // the number of stems to place in an area: exact (rounded) when hard_core, otherwise Poisson distributed
  auto stem_count = [&](double expected) {
    if (p.hard_core)
    {
      return static_cast<int>(std::round(expected));
    }
    return rng.poisson(expected);
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
    const int num_patches = rng.poisson(centre_density * 4.0 * half_patch_field * half_patch_field);
    // bucket the patch centres for a fast inside test
    BucketGrid2D<Eigen::Vector2d> patches(-half_patch_field, 2.0 * half_patch_field, std::max(p.patch_radius, 0.05));
    for (int i = 0; i < num_patches; i++)
    {
      const Eigen::Vector2d centre(-half_patch_field + 2.0 * half_patch_field * rng.uniform(),
                                   -half_patch_field + 2.0 * half_patch_field * rng.uniform());
      patches.add(centre, centre);
    }
    auto in_patch = [&](const Eigen::Vector2d &pos) {
      return patches.anyNear(pos, [&](const Eigen::Vector2d &centre) {
        return (pos - centre).squaredNorm() < p.patch_radius * p.patch_radius;
      });
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
  const int num_trees = rng.poisson(p.trees_per_m2 * tree_area);
  const double stems_per_metre_of_root = p.stems_per_m2 / (p.trees_per_m2 * p.roots_per_tree * p.root_length);
  int num_failed = 0;
  for (int t = 0; t < num_trees; t++)
  {
    const Eigen::Vector2d tree(-half_tree_field + 2.0 * half_tree_field * rng.uniform(),
                               -half_tree_field + 2.0 * half_tree_field * rng.uniform());
    const double angle_offset = 2.0 * kPi * rng.uniform();
    for (int r = 0; r < p.roots_per_tree; r++)
    {
      // evenly spread root directions, with some jitter
      const double angle =
        angle_offset + 2.0 * kPi * (static_cast<double>(r) + 0.3 * (rng.uniform() - 0.5)) / p.roots_per_tree;
      const Eigen::Vector2d along(std::cos(angle), std::sin(angle));
      const Eigen::Vector2d across(-along[1], along[0]);
      const int num = stem_count(stems_per_metre_of_root * p.root_length);
      for (int i = 0; i < num; i++)
      {
        bool placed = false;
        for (int attempt = 0; attempt < 100 && !placed; attempt++)
        {
          const Eigen::Vector2d pos =
            tree + along * (p.root_length * rng.uniform()) + across * (p.lateral_spread * rng.normal());
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
  generateStems(gen);
  SceneRandom rng(gen);
  stem_width_ = StemWidthModel(params_.diameter, params_.footprint);

  // build the bucket grid. Each stem is added to every bucket that its maximum effective radius overlaps
  const double half_field = 0.5 * params_.plot_width + params_.margin;
  field_min_ = -half_field;
  BucketGrid2D<int> stem_buckets(field_min_, 2.0 * half_field, bucket_width_);
  const double max_radius = 0.5 * (stem_width_.enabled() ? stem_width_.maxWidth(params_.max_range) : params_.diameter);
  max_height_ = 0.0;
  for (size_t i = 0; i < stems_.size(); i++)
  {
    const Eigen::Vector3d &stem = stems_[i];
    max_height_ = std::max(max_height_, stem[2]);
    stem_buckets.addWithin(stem.head<2>(), max_radius, static_cast<int>(i));
  }
  num_buckets_ = stem_buckets.size();
  buckets_ = std::move(stem_buckets.buckets());

  // a lawnmower walking path within the plot
  const double half_plot = 0.5 * params_.plot_width;
  const int num_lines = std::max(1, static_cast<int>(params_.plot_width / params_.walk_spacing));
  const double path_length = params_.plot_width * num_lines;
  cloud_.clear();
  const RGBA stem_colour(40, 160, 40, 255), ground_colour(120, 80, 40, 255), none_colour(128, 128, 128, 0),
    trunk_colour(150, 110, 90, 255);
  int num_stem_hits = 0, num_ground_hits = 0, num_omitted = 0;
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
    // the operator walks around trunks, so omit scanner positions within 0.2 m of a trunk
    bool near_trunk = false;
    for (const auto &trunk : trunks_)
    {
      near_trunk = near_trunk || (start.head<2>() - trunk.head<2>()).norm() < 0.5 * trunk[2] + 0.2;
    }
    if (near_trunk)
    {
      num_omitted++;
      continue;
    }

    // a rotating hand-held scanner covers most directions. Directions are uniform over the band of zenith angles
    const double cos_min = std::cos(params_.min_zenith_deg * kPi / 180.0);
    const double cos_max = std::cos(params_.max_zenith_deg * kPi / 180.0);
    const double dz = -(cos_max + (cos_min - cos_max) * rng.uniform());
    const double phi = 2.0 * kPi * rng.uniform();
    const double horizontal = std::sqrt(std::max(0.0, 1.0 - dz * dz));
    const Eigen::Vector3d dir(horizontal * std::cos(phi), horizontal * std::sin(phi), dz);

    double length;
    bool hit_stem, hit_ground, hit_trunk;
    castRay(start, dir, length, hit_stem, hit_ground, hit_trunk);
    // stems green, trunks brown-grey, ground brown, and rays with no return (at maximum range) grey
    RGBA colour = hit_stem ? stem_colour : (hit_trunk ? trunk_colour : (hit_ground ? ground_colour : none_colour));
    if (!hit_stem && !hit_trunk && (!hit_ground || params_.ground_unbounded))
    {
      colour.alpha = 0;
    }
    num_stem_hits += hit_stem ? 1 : 0;
    num_ground_hits += hit_ground ? 1 : 0;
    cloud_.addRay(start, start + dir * length, 1e-4 * static_cast<double>(i), colour);
  }
  if (num_omitted > 0)
  {
    std::cout << num_omitted << " rays omitted, from scanner positions within 0.2 m of a trunk" << std::endl;
  }
  std::cout << "generated " << stems_.size() << " stems and " << params_.num_rays - num_omitted << " rays, of which "
            << num_stem_hits << " hit stems and " << num_ground_hits << " hit the ground" << std::endl;
}

void PneumatophoreGen::castRay(const Eigen::Vector3d &start, const Eigen::Vector3d &dir, double &length, bool &hit_stem,
                               bool &hit_ground, bool &hit_trunk) const
{
  hit_stem = hit_ground = hit_trunk = false;
  const double inf = std::numeric_limits<double>::max();
  const double t_ground = dir[2] < 0.0 ? start[2] / -dir[2] : inf;
  double t_end = std::min(params_.max_range, t_ground);
  length = t_end;
  hit_ground = t_ground <= params_.max_range;

  // trunks (few, so tested directly), with the same footprint widening as the stems
  if (!trunks_.empty())
  {
    const Eigen::Vector2d s2(start[0], start[1]), d2(dir[0], dir[1]);
    const double h2 = d2.squaredNorm();
    const double max_widening = params_.footprint.enabled() ? params_.footprint.diameter(params_.max_range) : 0.0;
    for (const auto &trunk : trunks_)
    {
      const Eigen::Vector2d rel = trunk.head<2>() - s2;
      const double along = h2 > 1e-12 ? rel.dot(d2) / h2 : 0.0;
      const double offset2 = std::max(0.0, rel.squaredNorm() - along * along * h2);
      const double bound = 0.5 * trunk[2] + max_widening;
      if (offset2 >= bound * bound || h2 <= 1e-12)
      {
        continue;
      }
      const double width =
        params_.footprint.enabled() ?
          effectiveStemWidth(trunk[2], params_.footprint.diameter(std::max(0.0, along)), params_.footprint.threshold) :
          trunk[2];
      const double radius = 0.5 * width;
      if (offset2 >= radius * radius)
      {
        continue;
      }
      // the trunk is intercepted. Its return is reported at the physical surface (or the closest approach of the
      // beam axis, when it misses the trunk), or at the widened outline (BeamFootprint::range_to_outline).
      // The trunk must be ahead of the scanner (which is never inside a trunk)
      const double report_radius = params_.footprint.range_to_outline ? radius : std::min(radius, 0.5 * trunk[2]);
      const double t_in = along - std::sqrt(std::max(0.0, report_radius * report_radius - offset2)) / std::sqrt(h2);
      const double z_in = start[2] + t_in * dir[2];
      if (t_in > 0.0 && t_in < t_end && z_in >= 0.0 && z_in <= params_.trunk_height)
      {
        t_end = t_in;
        length = t_in;
        hit_trunk = true;
        hit_ground = false;
      }
    }
  }

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
        // hits the side of the stem. The return is reported where the beam axis enters the widened outline, or
        // (by default) at the physical stem surface, or the closest approach when the axis misses the stem
        t_hit = t_in;
        if (!params_.footprint.range_to_outline && h > 1e-6)
        {
          const double stem_radius = 0.5 * params_.diameter;
          t_hit = std::max(t_in, along - std::sqrt(std::max(0.0, stem_radius * stem_radius - offset2)) / h);
          if (dir[2] > 0.0)
          {
            t_hit = std::min(t_hit, (stem[2] - start[2]) / dir[2]);  // a rising ray's return is not above the stem
          }
        }
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
    hit_trunk = false;
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
bool PneumatophoreGen::writeTrunks(const std::string &file_name) const
{
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing" << std::endl;
    return false;
  }
  ofs << "x,y,diameter" << std::endl;
  for (const auto &trunk : trunks_)
  {
    ofs << trunk[0] << "," << trunk[1] << "," << trunk[2] << std::endl;
  }
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
  // trunks' frontal area, in the cell holding each trunk's centre
  std::vector<double> trunk_area(static_cast<size_t>(nx * ny * nz), 0.0);
  for (const auto &trunk : trunks_)
  {
    const int x = static_cast<int>(std::floor((trunk[0] + half_plot) / cell_size[0]));
    const int y = static_cast<int>(std::floor((trunk[1] + half_plot) / cell_size[1]));
    if (x < 0 || y < 0 || x >= nx || y >= ny)
    {
      continue;
    }
    for (int z = 0; z < nz; z++)
    {
      trunk_area[static_cast<size_t>(x + nx * (y + ny * z))] += trunk[2] * cell_size[2];
    }
  }
  ofs << "x_min,y_min,z_min,x_max,y_max,z_max,true_frontal_area_per_volume_Nd,true_trunk_frontal_area_per_volume"
      << std::endl;
  for (int z = 0; z < nz; z++)
  {
    for (int y = 0; y < ny; y++)
    {
      for (int x = 0; x < nx; x++)
      {
        const Eigen::Vector3d min_bound(-half_plot + cell_size[0] * x, -half_plot + cell_size[1] * y, cell_size[2] * z);
        const Eigen::Vector3d max_bound = min_bound + cell_size;
        ofs << min_bound[0] << "," << min_bound[1] << "," << min_bound[2] << "," << max_bound[0] << "," << max_bound[1]
            << "," << max_bound[2] << "," << frontal_area[static_cast<size_t>(x + nx * (y + ny * z))] / volume << ","
            << trunk_area[static_cast<size_t>(x + nx * (y + ny * z))] / volume << std::endl;
      }
    }
  }
  return true;
}
}  // namespace ray
