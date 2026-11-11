// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#ifndef RAYLIB_RAYSTEMMODEL_H
#define RAYLIB_RAYSTEMMODEL_H

#include "raylib/raylibconfig.h"

#include "rayutils.h"

#include <vector>

/// Models of thin vertical stems (e.g. mangrove pneumatophores) as seen by a lidar beam: the beam footprint and the
/// stems' effective (detectable) width, and the response of the free path density estimate to stems that don't
/// overlap (the occlusion correction). Used by DensityGrid (raydensitygrid.h) and PneumatophoreGen.
namespace ray
{
/// Gaussian model of a lidar beam's footprint, used to estimate the detectable width of thin elements
struct RAYLIB_EXPORT BeamFootprint
{
  /// beam diameter at the scanner (m), defined at the 1/e^2 intensity level
  double exit_diameter = 0.0;
  /// full-angle beam divergence (radians), defined at the 1/e^2 intensity level
  double divergence = 0.0;
  /// minimum fraction of the beam's power that an element must intercept to trigger a return. In practice this is
  /// a calibration parameter: it depends on the scanner's detection threshold, the element reflectance and range.
  double threshold = 0.05;
  /// Where along the beam axis a return from a partly intercepted stem is reported. false (default): at the range
  /// of the physical stem surface (the nearest illuminated part of the stem), as a range measurement would give.
  /// true: where the beam axis enters the footprint-widened outline of the stem (an idealisation used before).
  /// The correct choice for a given scanner can be checked from the spread of returns around isolated stems.
  bool range_to_outline = false;
  inline bool enabled() const { return exit_diameter > 0.0 || divergence > 0.0; }
  /// the 1/e^2 beam diameter at the given range
  inline double diameter(double range) const { return exit_diameter + divergence * range; }
};

/// Fraction of a Gaussian beam's power (1/e^2 diameter @c beam_diameter) that is intercepted by an infinitely long
/// strip of width @c width whose centreline is offset by @c offset from the beam centre. A vertical cylinder of
/// diameter d is seen as such a strip of width d from any ray direction.
double RAYLIB_EXPORT stripInterceptFraction(double width, double beam_diameter, double offset);

/// The effective (detectable) width of a long cylinder of diameter @c diameter: the width of the band of beam-centre
/// offsets for which the intercepted power fraction reaches @c threshold. It is larger than @c diameter for small
/// thresholds, and is zero when the cylinder is too thin to ever trigger a return.
double RAYLIB_EXPORT effectiveStemWidth(double diameter, double beam_diameter, double threshold);

/// The response of the Beer-Lambert (free path) estimate to non-overlapping vertical stems: the hit rate per unit
/// free path divided by N d_eff, for stems at solid fraction N pi d^2 / 4, seen by the beam with effective width
/// width_ratio * d. A stem is intercepted when the beam axis passes within d_eff / 2 of its centre, and the return
/// is reported at the range given by @c range_to_outline (see BeamFootprint::range_to_outline):
/// - true: where the beam axis enters the widened outline (diameter d_eff). The free space, where paths can start,
///   is outside the widened outlines.
/// - false: at the physical stem surface (or at the closest approach of the beam axis, when the axis misses the
///   stem). The free space is outside the physical stems, except where a ray would already have returned: inside a
///   widened outline beyond its reported range (the far side of the outline, in the ray direction). So paths can
///   start within a widened outline, before its return. Starting paths only outside the widened outlines instead
///   would bias the response low (below 1), as no path would then end in the short distance between the outline
///   and the reported range.
/// For overlapping (Poisson) stems the response is exactly 1 under both conventions. Non-overlapping stems leave
/// less free space, so the response is at least 1 / (1 - solid fraction), and grows with width_ratio. The two
/// conventions give responses within the simulation noise of each other (about 1%).
/// Paths start at random points of the free space, as rays enter the stem layer at random points. Computed
/// deterministically by a 2D Monte-Carlo simulation (fixed seed), pooled over 64 independent realisations (each of
/// 1500 stems and 10000 lines), giving a standard error of under 1%.
/// min_gap_ratio is the minimum gap between stem surfaces in stem diameters, so that stem centres are at least
/// d (1 + min_gap_ratio) apart.
/// Returns infinity (saturated) when there is effectively no free space (under 1% of the area), e.g. for widened
/// outlines that cover the plane.
/// @c seed sets the simulation's random number sequence (the correction uses the default).
double RAYLIB_EXPORT freePathResponse(double solid_fraction, double width_ratio, double min_gap_ratio = 0.0,
                                      bool range_to_outline = false, unsigned int seed = 12345);

/// The mean distance from a uniformly random point of a rectangle of sides @c width_x and @c width_y to the
/// rectangle's boundary, along a uniformly random direction in the plane: (1 / (2 pi A)) times the integral of
/// 1 / |p - q| over pairs of points p, q of the rectangle (area A). 0.4732 times the side for a square. This is the
/// mean length of the shadow that an element casts within the rectangle, behind it, in a uniform field of horizontal
/// ray directions (see DensityOptions::self_shadow_correction).
double RAYLIB_EXPORT meanExitDistance(double width_x, double width_y);

/// The second-order self-shadow factor of the pooled estimator's neighbourhoods (see
/// DensityOptions::self_shadow_second_order): the mean over neighbourhoods of their measured density (hits / free
/// path) divided by the pooled measured density, for stems of effective width @c effective_width and diameter
/// @c diameter in neighbourhoods of area @c area and mean exit distance @c exit_distance (meanExitDistance). With
/// x = exit_distance * effective_width / area (the first-order factor is 1 + x):
///   1 / factor = 1 - x - k_f pi w^2 / (16 area) + 0.76 x^3,
/// where w and k_f are the diameter and 0.55 for returns at the physical surface, and the effective width and 1.2 for
/// returns at the widened outline (@c range_to_outline). 1 - x is derived: a stem's shadow within its own
/// neighbourhood removes x of the neighbourhood's free path on average (the variation of the shadow's length with the
/// stem's position adds under 0.4% for x <= 0.3). The other terms, for the spread of the reported returns along the
/// ray and for the width of the shadow (and returns falling in neighbouring voxels), are fitted to a 2D Monte-Carlo
/// simulation of tall stems at low density (N d_eff W <= 0.3 for neighbourhoods of width W). For returns at the
/// surface the fit is within about 1.5% of the factor (the simulation's noise is about 1%) for d_eff / W up to 1 with
/// diameter / effective width 0.3 to 0.65, and for d_eff = d up to d / W = 0.45. For returns at the outline it rests
/// on two simulations at d_eff / W = 0.45.
/// Denser neighbourhoods have a larger factor, as their free path is concentrated near their edges. The factor adds
/// 0.35 x tau for the optical depth tau = N d_eff exit_distance (@c optical_depth), fitted to the simulation (for
/// d_eff / W = 0.2 to 0.45 and N d_eff W up to 1, where it rises by up to 4%) and to 3D scans (600 and 1000 stems/m^2
/// in single 0.05 m voxels). Beyond N d_eff W of about 1 the simulated factor grows much faster (by 50% at 2 for
/// d_eff / W = 0.45), and in 3D scans that term also depends on how the rays enter the neighbourhood (from above as
/// well as from the sides), so the factor is too small there
double RAYLIB_EXPORT selfShadowFactor(double exit_distance, double area, double effective_width, double diameter,
                                      double optical_depth, bool range_to_outline);

/// Lookup of freePathResponse, bilinearly interpolated over a grid of solid fractions and width ratios (1 to 5).
/// For no minimum gap, the grid values are precomputed (raystemresponsetable.h, generated by
/// scripts/generate_stem_response_table.cpp). Other values are computed when first needed.
class RAYLIB_EXPORT StemResponseTable
{
public:
  explicit StemResponseTable(double min_gap_ratio = 0.0, bool range_to_outline = false);
  /// The interpolated response. Infinity when the interpolation uses a saturated grid value
  double response(double solid_fraction, double width_ratio);

private:
  double value(size_t fraction_index, size_t ratio_index);
  std::vector<double> table_;  // grid values, NaN where not yet computed
  double min_gap_ratio_;
  bool range_to_outline_;
};

/// Solves measured = a R(a) for the true density a in [0, max_density], where a R(a) increases with a.
/// R (@c response) may be below 1, so the solution may exceed the measured density: the solution is bracketed upwards
/// from the measured density, doubling up to max_density. Starting from the measured density only evaluates the
/// response up to the density needed (it can be slow to compute). @c saturated is set when the measured density is
/// beyond the response at max_density, and max_density is returned.
template <class Response>
double invertDensityResponse(double measured, double max_density, const Response &response, bool &saturated)
{
  saturated = false;
  if (measured <= 0.0)
  {
    return 0.0;
  }
  double low = 0.0;
  double high = std::min(measured, max_density);
  while (high * response(high) < measured)
  {
    if (high >= max_density)
    {
      saturated = true;
      return max_density;
    }
    low = high;
    high = std::min(2.0 * high, max_density);
  }
  for (int it = 0; it < 40; it++)
  {
    const double mid = 0.5 * (low + high);
    (mid * response(mid) < measured ? low : high) = mid;
  }
  return 0.5 * (low + high);
}

/// Inverts measured = a R(pi d a / 4, rho) for the true density a = N d of stems of diameter @c stem_diameter, with
/// the response R from @c table, and the solid fraction limited to @c max_solid_fraction / (1 + min_gap_ratio)^2.
/// See invertDensityResponse.
double RAYLIB_EXPORT invertStemResponse(double measured, double width_ratio, double stem_diameter, double min_gap_ratio,
                                        double max_solid_fraction, StemResponseTable &table, bool &saturated);

/// Tabulated effective width of a stem of given diameter as a function of range, for fast lookup
class RAYLIB_EXPORT StemWidthModel
{
public:
  StemWidthModel() {}
  StemWidthModel(double diameter, const BeamFootprint &footprint);
  /// effective (detectable) width of the stem at the given range from the scanner
  inline double width(double range) const { return enabled_ ? scale(range) * diameter_ : diameter_; }
  /// effective width divided by true diameter
  inline double scale(double range) const
  {
    if (!enabled_)
    {
      return 1.0;
    }
    const double x = std::max(0.0, range) * inv_table_step_;
    if (x < table_max_index_)
    {
      const size_t i = static_cast<size_t>(x);
      const double blend = x - static_cast<double>(i);
      return table_[i] * (1.0 - blend) + table_[i + 1] * blend;
    }
    return farScale(range);
  }
  /// the maximum effective width over ranges 0 to @c max_range
  double maxWidth(double max_range) const;
  inline bool enabled() const { return enabled_; }

private:
  /// the scale beyond table_max_range_, interpolated over logarithmic range intervals
  double farScale(double range) const;
  double diameter_ = 0.0;
  bool enabled_ = false;
  std::vector<double> table_;      // effective width / diameter at regular range intervals up to table_max_range_
  std::vector<double> far_table_;  // effective width / diameter at logarithmic range intervals beyond, to far_range_
  double table_max_index_ = 0.0;   // the interpolation index limit of table_
  static constexpr double table_step_ = 0.01;
  static constexpr double inv_table_step_ = 1.0 / table_step_;
  static constexpr double table_max_range_ = 100.0;
  static constexpr double far_range_ = 100000.0;  // ranges beyond this take the value at far_range_
  static constexpr int far_steps_ = 1000;
};
}  // namespace ray

#endif  // RAYLIB_RAYSTEMMODEL_H
