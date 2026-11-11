// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#include "raystemmodel.h"

#include "raystemresponsetable.h"

#include <random>

namespace ray
{
double stripInterceptFraction(double width, double beam_diameter, double offset)
{
  if (beam_diameter <= 0.0)
  {
    return std::abs(offset) < 0.5 * width ? 1.0 : 0.0;
  }
  // the 1/e^2 diameter is 4 standard deviations. The power across the strip is the 1D marginal of the 2D Gaussian
  const double sigma = 0.25 * beam_diameter;
  const double scale = 1.0 / (std::sqrt(2.0) * sigma);
  return 0.5 * (std::erf((offset + 0.5 * width) * scale) - std::erf((offset - 0.5 * width) * scale));
}

double effectiveStemWidth(double diameter, double beam_diameter, double threshold)
{
  if (beam_diameter <= 0.0)
  {
    return threshold <= 1.0 ? diameter : 0.0;
  }
  if (stripInterceptFraction(diameter, beam_diameter, 0.0) < threshold)
  {
    return 0.0;  // never intercepts enough of the beam to be detected
  }
  // the intercepted fraction decreases monotonically with offset, so bisect for where it equals the threshold
  double low = 0.0;
  double high = 0.5 * diameter + 4.0 * beam_diameter;
  for (int i = 0; i < 60; i++)
  {
    const double mid = 0.5 * (low + high);
    if (stripInterceptFraction(diameter, beam_diameter, mid) >= threshold)
    {
      low = mid;
    }
    else
    {
      high = mid;
    }
  }
  return 2.0 * low;
}

double meanExitDistance(double width_x, double width_y)
{
  if (width_x <= 0.0 || width_y <= 0.0)
  {
    return 0.0;
  }
  // 4 F / (2 pi a b), where F is the integral over 0 <= u <= a, 0 <= v <= b of (a - u)(b - v) / sqrt(u^2 + v^2),
  // evaluated in polar coordinates: the radial integral is analytic, the angular one by Simpson's rule
  const double a = width_x, b = width_y;
  const double corner = std::atan2(b, a);
  auto radial = [&](double theta) {
    const double c = std::cos(theta), s = std::sin(theta);
    const double r = theta < corner ? a / c : b / s;
    return a * b * r - 0.5 * (a * s + b * c) * r * r + c * s * r * r * r / 3.0;
  };
  // integrate each side of the corner angle separately, as the integrand has a kink there
  auto simpson = [&](double t0, double t1) {
    const int n = 512;  // even
    const double h = (t1 - t0) / n;
    double sum = radial(t0) + radial(t1);
    for (int i = 1; i < n; i++)
    {
      sum += (i % 2 == 1 ? 4.0 : 2.0) * radial(t0 + h * i);
    }
    return sum * h / 3.0;
  };
  const double integral = simpson(0.0, corner) + simpson(corner, 0.5 * kPi);
  return 4.0 * integral / (2.0 * kPi * a * b);
}

StemWidthModel::StemWidthModel(double diameter, const BeamFootprint &footprint)
  : diameter_(diameter)
  , enabled_(diameter > 0.0 && footprint.enabled())
{
  if (!enabled_)
  {
    return;
  }
  auto scale_at = [&](double range) {
    return effectiveStemWidth(diameter_, footprint.diameter(range), footprint.threshold) / diameter_;
  };
  const int num = static_cast<int>(std::round(table_max_range_ / table_step_)) + 1;
  table_.resize(static_cast<size_t>(num));
  for (int i = 0; i < num; i++)
  {
    table_[static_cast<size_t>(i)] = scale_at(table_step_ * static_cast<double>(i));
  }
  table_max_index_ = static_cast<double>(num - 1);
  // beyond the regular table, the width varies slowly with range, so use logarithmic intervals
  far_table_.resize(static_cast<size_t>(far_steps_) + 1);
  for (int i = 0; i <= far_steps_; i++)
  {
    far_table_[static_cast<size_t>(i)] =
      scale_at(table_max_range_ * std::pow(far_range_ / table_max_range_, static_cast<double>(i) / far_steps_));
  }
}

double StemWidthModel::farScale(double range) const
{
  const double x = far_steps_ * std::log(std::max(range, table_max_range_) / table_max_range_) /
                   std::log(far_range_ / table_max_range_);
  if (!(x < static_cast<double>(far_steps_)))
  {
    return far_table_.back();
  }
  const size_t i = static_cast<size_t>(x);
  const double blend = x - static_cast<double>(i);
  return far_table_[i] * (1.0 - blend) + far_table_[i + 1] * blend;
}

double StemWidthModel::maxWidth(double max_range) const
{
  double max_width = diameter_;
  for (double range = 0.0; range <= max_range; range += table_step_)
  {
    max_width = std::max(max_width, width(range));
  }
  return max_width;
}

namespace
{
/// Deterministic uniform random number in (0, 1), independent of the standard library's distributions
inline double uniform01(std::mt19937 &gen)
{
  return (static_cast<double>(gen()) + 0.5) / 4294967296.0;
}

/// One realisation of the free path simulation for freePathResponse, accumulating hits and free path length.
/// Returns false when the free space is too small to find the lines' start points (saturated)
bool freePathRealisation(double solid_fraction, double width_ratio, double min_gap_ratio, bool range_to_outline,
                         unsigned int seed, double &hits, double &path, double &density_out)
{
  // non-overlapping stems of unit diameter, placed by random sequential addition in a square
  const int num_stems = 1500;
  const int num_lines = 10000;
  // start points are found by rejection, so give up when under 1% of the candidate points are free
  const long long max_attempts = 100LL * num_lines;
  const double diameter = 1.0;
  const double side = std::sqrt(num_stems * kPi * diameter * diameter / (4.0 * solid_fraction));
  density_out = num_stems / (side * side);
  std::mt19937 gen(seed);
  const double separation = diameter * (1.0 + min_gap_ratio);  // minimum distance between stem centres
  const double radius = 0.5 * width_ratio * diameter;          // of the widened stems
  const double stem_radius = 0.5 * diameter;
  const double cell = std::max(separation, radius);  // so stems within these distances are in adjacent cells
  const int cells = std::max(1, static_cast<int>(std::ceil(side / cell)));
  std::vector<std::vector<Eigen::Vector2d>> buckets(static_cast<size_t>(cells * cells));
  auto bucket = [&](double x) { return std::min(cells - 1, std::max(0, static_cast<int>(x / cell))); };
  // whether any stem centre is closer than dist to pos (dist no more than the cell width)
  auto near = [&](const Eigen::Vector2d &pos, double dist) {
    const int bx = bucket(pos[0]), by = bucket(pos[1]);
    for (int x = std::max(0, bx - 1); x <= std::min(cells - 1, bx + 1); x++)
    {
      for (int y = std::max(0, by - 1); y <= std::min(cells - 1, by + 1); y++)
      {
        for (const auto &c : buckets[static_cast<size_t>(x + cells * y)])
        {
          if ((pos - c).squaredNorm() < dist * dist)
          {
            return true;
          }
        }
      }
    }
    return false;
  };
  // The range at which a line reports a return from a stem, given that the line passes within the widened radius
  // of its centre. b is the range of closest approach and offset2 the squared offset there
  auto reported_range = [&](double b, double offset2) {
    const double r = range_to_outline ? radius : stem_radius;
    return b - std::sqrt(std::max(0.0, r * r - offset2));
  };
  // with returns at the physical range, a line can start within a widened outline (outside the physical stem), but
  // only before its reported range, as a ray beyond that would have already returned
  auto past_return = [&](const Eigen::Vector2d &pos, const Eigen::Vector2d &dir) {
    const int bx = bucket(pos[0]), by = bucket(pos[1]);
    for (int x = std::max(0, bx - 1); x <= std::min(cells - 1, bx + 1); x++)
    {
      for (int y = std::max(0, by - 1); y <= std::min(cells - 1, by + 1); y++)
      {
        for (const auto &c : buckets[static_cast<size_t>(x + cells * y)])
        {
          const Eigen::Vector2d rel = c - pos;
          if (rel.squaredNorm() < radius * radius)
          {
            const double b = rel.dot(dir);
            if (reported_range(b, rel.squaredNorm() - b * b) <= 0.0)
            {
              return true;
            }
          }
        }
      }
    }
    return false;
  };
  int placed = 0;
  for (long long attempt = 0; placed < num_stems && attempt < 100000000LL; attempt++)
  {
    const Eigen::Vector2d pos(side * uniform01(gen), side * uniform01(gen));
    if (!near(pos, separation))
    {
      buckets[static_cast<size_t>(bucket(pos[0]) + cells * bucket(pos[1]))].push_back(pos);
      placed++;
    }
  }
  // lines from random points of the free space, until their first reported return or the edge of the square
  long long attempts = 0;
  for (int l = 0; l < num_lines; l++)
  {
    Eigen::Vector2d start, dir;
    while (true)
    {
      if (++attempts > max_attempts)
      {
        return false;
      }
      start = Eigen::Vector2d(side * (0.1 + 0.8 * uniform01(gen)), side * (0.1 + 0.8 * uniform01(gen)));
      if (near(start, range_to_outline ? radius : stem_radius))
      {
        continue;
      }
      const double angle = 2.0 * kPi * uniform01(gen);
      dir = Eigen::Vector2d(std::cos(angle), std::sin(angle));
      if (range_to_outline || width_ratio <= 1.0 || !past_return(start, dir))
      {
        break;
      }
    }
    double t_max = std::numeric_limits<double>::max();
    for (int i = 0; i < 2; i++)
    {
      if (dir[i] > 1e-12)
      {
        t_max = std::min(t_max, (side - start[i]) / dir[i]);
      }
      else if (dir[i] < -1e-12)
      {
        t_max = std::min(t_max, -start[i] / dir[i]);
      }
    }
    // walk the cells along the line. Any stem whose widened disc reaches a cell has its centre in that cell or an
    // adjacent one, and its reported range is no less than its outline entry, so once a reported return lies within
    // the current cell, no later cell can have an earlier one
    double best = t_max;
    int cx = bucket(start[0]), cy = bucket(start[1]);
    const int step_x = dir[0] > 0.0 ? 1 : -1, step_y = dir[1] > 0.0 ? 1 : -1;
    double next_x = std::abs(dir[0]) > 1e-12 ? ((cx + (step_x > 0 ? 1 : 0)) * cell - start[0]) / dir[0] :
                                               std::numeric_limits<double>::max();
    double next_y = std::abs(dir[1]) > 1e-12 ? ((cy + (step_y > 0 ? 1 : 0)) * cell - start[1]) / dir[1] :
                                               std::numeric_limits<double>::max();
    const double delta_x = std::abs(dir[0]) > 1e-12 ? cell / std::abs(dir[0]) : std::numeric_limits<double>::max();
    const double delta_y = std::abs(dir[1]) > 1e-12 ? cell / std::abs(dir[1]) : std::numeric_limits<double>::max();
    while (true)
    {
      for (int x = std::max(0, cx - 1); x <= std::min(cells - 1, cx + 1); x++)
      {
        for (int y = std::max(0, cy - 1); y <= std::min(cells - 1, cy + 1); y++)
        {
          for (const auto &c : buckets[static_cast<size_t>(x + cells * y)])
          {
            const Eigen::Vector2d rel = c - start;
            const double b = rel.dot(dir);
            const double offset2 = rel.squaredNorm() - b * b;
            if (offset2 <= radius * radius)
            {
              const double t = reported_range(b, offset2);
              if (t > 0.0 && t < best)
              {
                best = t;
              }
            }
          }
        }
      }
      const double t_exit = std::min(next_x, next_y);
      if (best <= t_exit || t_exit >= t_max)
      {
        break;
      }
      if (next_x < next_y)
      {
        cx += step_x;
        next_x += delta_x;
      }
      else
      {
        cy += step_y;
        next_y += delta_y;
      }
      if (cx < 0 || cy < 0 || cx >= cells || cy >= cells)
      {
        break;
      }
    }
    path += best;
    hits += best < t_max ? 1.0 : 0.0;
  }
  return true;
}
}  // namespace

double freePathResponse(double solid_fraction, double width_ratio, double min_gap_ratio, bool range_to_outline,
                        unsigned int seed)
{
  if (solid_fraction <= 0.0)
  {
    return 1.0;
  }
  min_gap_ratio = std::max(0.0, min_gap_ratio);
  const double separation_scale = (1.0 + min_gap_ratio) * (1.0 + min_gap_ratio);
  solid_fraction = std::min(solid_fraction, 0.4 / separation_scale);
  width_ratio = std::max(width_ratio, 1.0);
  // pool several independent realisations, as one realisation's response varies by about 2-3% between seeds
  const int num_realisations = 64;
  double hits = 0.0, path = 0.0, density = 0.0;
  for (int i = 0; i < num_realisations; i++)
  {
    if (!freePathRealisation(solid_fraction, width_ratio, min_gap_ratio, range_to_outline,
                             seed + 7919u * static_cast<unsigned int>(i), hits, path, density))
    {
      return std::numeric_limits<double>::infinity();  // effectively no free space
    }
  }
  return hits / (path * density * width_ratio);
}

StemResponseTable::StemResponseTable(double min_gap_ratio, bool range_to_outline)
  : min_gap_ratio_(std::max(0.0, min_gap_ratio))
  , range_to_outline_(range_to_outline)
{
  using namespace stem_response;
  table_.assign(static_cast<size_t>(kNumFractions * kNumRatios), std::numeric_limits<double>::quiet_NaN());
  if (min_gap_ratio_ == 0.0)
  {
    const int convention = range_to_outline_ ? 1 : 0;
    for (int f = 0; f < kNumFractions; f++)
    {
      for (int r = 0; r < kNumRatios; r++)
      {
        table_[static_cast<size_t>(f * kNumRatios + r)] = kGapZero[convention][f][r];
      }
    }
  }
}

double StemResponseTable::value(size_t fraction_index, size_t ratio_index)
{
  using namespace stem_response;
  double &entry = table_[fraction_index * static_cast<size_t>(kNumRatios) + ratio_index];
  if (std::isnan(entry))
  {
    const double separation_scale = (1.0 + min_gap_ratio_) * (1.0 + min_gap_ratio_);
    entry = kFractions[fraction_index] <= 0.0 ?
              1.0 :
              freePathResponse(kFractions[fraction_index] / separation_scale,
                               1.0 + kRatioStep * static_cast<double>(ratio_index), min_gap_ratio_, range_to_outline_);
  }
  return entry;
}

double StemResponseTable::response(double solid_fraction, double width_ratio)
{
  using namespace stem_response;
  // bilinear interpolation over a grid of solid fractions and width ratios
  const double max_ratio = 1.0 + kRatioStep * (kNumRatios - 1);
  // the table is over the solid fraction of the exclusion discs (diameter d (1 + gap)), which have the packing limit
  const double separation_scale = (1.0 + min_gap_ratio_) * (1.0 + min_gap_ratio_);
  solid_fraction = std::max(0.0, std::min(solid_fraction * separation_scale, kFractions[kNumFractions - 1]));
  const double r = (std::max(1.0, std::min(width_ratio, max_ratio)) - 1.0) / kRatioStep;
  const size_t r0 = static_cast<size_t>(std::min(static_cast<int>(r), kNumRatios - 2));
  size_t f0 = 0;
  while (f0 + 2 < static_cast<size_t>(kNumFractions) && kFractions[f0 + 1] < solid_fraction)
  {
    f0++;
  }
  const double fb = (solid_fraction - kFractions[f0]) / (kFractions[f0 + 1] - kFractions[f0]);
  const double rb = r - static_cast<double>(r0);
  // skip zero weights, so that a saturated (infinite) neighbour is only used where it contributes
  double sum = 0.0;
  for (size_t i = 0; i < 2; i++)
  {
    for (size_t j = 0; j < 2; j++)
    {
      const double weight = (i == 0 ? 1.0 - fb : fb) * (j == 0 ? 1.0 - rb : rb);
      if (weight > 0.0)
      {
        sum += weight * value(f0 + i, r0 + j);
      }
    }
  }
  return sum;
}

double invertStemResponse(double measured, double width_ratio, double stem_diameter, double min_gap_ratio,
                          double max_solid_fraction, StemResponseTable &table, bool &saturated)
{
  const double separation_scale = (1.0 + min_gap_ratio) * (1.0 + min_gap_ratio);
  const double max_density = max_solid_fraction * 4.0 / (kPi * stem_diameter * separation_scale);
  auto response = [&](double a) { return table.response(0.25 * kPi * stem_diameter * a, width_ratio); };
  return invertDensityResponse(measured, max_density, response, saturated);
}
}  // namespace ray
