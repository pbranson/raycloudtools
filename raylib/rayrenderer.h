// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYRENDERER_H
#define RAYLIB_RAYRENDERER_H

#include "raycuboid.h"
#include "raypose.h"
#include "rayutils.h"

#include <map>
#include <memory>
#define FLAT_TOP_COMPENSATION  // use a 2-medium model for the top of the cloud, so only the cloud's density is
                               // estimated, not the air above it

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

/// The assumed orientation distribution of the intercepting elements (leaves, stems, roots) within each voxel.
/// This selects the G function: the ratio of the elements' mean projected area along a ray direction to the
/// reference area that the density is reported in.
enum class RAYLIB_EXPORT ElementDistribution
{
  /// Uniform spherical (isotropic) distribution of flat leaf angles. G = 0.5 for every ray direction, and density
  /// is reported as one-sided leaf area per unit volume (m^2/m^3). This is the original raycloudtools behaviour.
  Spherical,
  /// Vertical cylinders (e.g. stems, mangrove pneumatophores) of diameter d at N stems per m^2 of ground. The
  /// interception rate of a ray at zenith angle theta is N d sin(theta), so G = sin(theta) and density is reported
  /// as the frontal area per unit volume a = N d (m^-1), as used in vegetation drag models. The integral over
  /// height (the rendered pixel value) is the frontal area index N d h (m^2/m^2).
  /// Stem tops (end caps) are not modelled.
  VerticalCylinder
};

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

/// How the path length of a ray that hits within a voxel is counted
enum class RAYLIB_EXPORT PathEstimator
{
  /// The path up to the hit point. This is the maximum likelihood estimate for a turbid medium (many small elements
  /// randomly placed within each voxel, Beer-Lambert attenuation). It is biased high when voxels contain only a few
  /// discrete elements, by about 1/(1 - p/2) where p is the chance of a ray crossing a voxel being intercepted.
  BeerLambert,
  /// The path to where the ray would have exited the voxel. density = hits / crossing length is exact for any
  /// spatial arrangement of elements (including clumped), provided elements rarely occlude each other within a
  /// single voxel. So it suits voxels that are small relative to the mean free path, e.g. pneumatophores.
  Crossing,
  /// The path up to the hit point (as BeerLambert), but each voxel takes the density of its horizontal
  /// neighbourhood (see DensityOptions::pool_radius), pooled as total hits over total path length. Pooling avoids the
  /// bias of averaging noisy per-voxel ratios, and the Beer-Lambert estimate needs no correction for stems occluding
  /// each other within a voxel, as it is exact for randomly placed stems. The neighbourhood limits it to resolving
  /// patchiness at scales smaller than the neighbourhood.
  PooledBeerLambert
};

/// Options controlling the voxel density estimate
struct RAYLIB_EXPORT DensityOptions
{
  PathEstimator path_estimator = PathEstimator::BeerLambert;
  ElementDistribution distribution = ElementDistribution::Spherical;
  /// For VerticalCylinder only: rays whose sin(zenith angle) is below this are ignored entirely, since near-vertical
  /// rays carry almost no information on vertical elements, and their few (end cap) hits would give outliers.
  /// Ignoring rays based on their direction alone does not bias the estimate.
  double min_sin_zenith = 0.1;
  /// For VerticalCylinder only: rays whose sin(zenith angle) is above this are ignored, to analyse a band of angles
  double max_sin_zenith = 1.0;
  /// Optional .csv file to write the mean density per horizontal layer of voxels
  std::string profile_file;
  /// Optional .csv file to write the density aggregated into cells of size @c cell_size (e.g. for hydrodynamic
  /// models), tiling from the minimum bound
  std::string cells_file;
  Eigen::Vector3d cell_size = Eigen::Vector3d(1.0, 1.0, 0.1);
  /// For VerticalCylinder only: when the footprint is enabled and stem_diameter > 0, each ray is weighted by the
  /// ratio of the stem's effective (detectable) width at that range to its true diameter. This corrects the
  /// apparent widening of thin stems by the beam footprint, so the density remains N d for the true diameter d.
  /// footprint.range_to_outline selects the return range convention assumed by the occlusion correction and the
  /// crossing estimator (see BeamFootprint::range_to_outline and freePathResponse).
  BeamFootprint footprint;
  double stem_diameter = 0.0;
  /// For VerticalCylinder with a stem_diameter, and the beer_lambert or pooled estimators: correct the free-path
  /// estimate for the space occupied by the stems, and for the footprint widening. The free path excludes the space
  /// that rays can't reach before their return (the physical stems, or the widened outlines, depending on the
  /// return range convention), so the measured density is a R(phi, rho) for true density a = N d, solid fraction
  /// phi = N pi d^2 / 4 and widening rho = d_eff / d (see freePathResponse). This is inverted using the density of
  /// the neighbourhood of pool_radius horizontally (see pool_radius). The solid fraction is limited to 0.4 (below
  /// the random packing limit of 0.547 for non-overlapping discs), so N d is at most 1.6 / (pi d).
  bool overlap_correction = true;
  /// Nominal half-width (m) of the horizontal neighbourhood used by the pooled estimator and the occlusion
  /// correction. The neighbourhood is the (2k+1) x (2k+1) voxels centred on each voxel, for
  /// k = round(pool_radius / voxel width), so its actual half-width is (k + 1/2) voxel widths: +-0.075 m for 0.05 m
  /// at 0.05 m voxels, and +-0.125 m for 0.1 m. Smaller resolves finer patchiness, but with fewer hits per
  /// neighbourhood, which biases the estimate high in sparse areas (about +6% at 30 stems/m^2 for 0.05 m, +3% for 0.1
  /// m). Don't adapt it to the number of hits, as that is a stopping rule that biases the estimate high.
  double pool_radius = 0.05;
  /// The minimum gap between stem surfaces assumed by the occlusion correction, in stem diameters (stem centres are
  /// at least d (1 + min_gap_ratio) apart). More regular spacing leaves less free space between the widened stems.
  double min_gap_ratio = 0.0;
  /// Model the top of the cloud as a sharp boundary between air and a uniform medium (e.g. a grass or crop
  /// canopy), so the air above the highest point in each column is excluded. Not suited to elements of varying
  /// height such as pneumatophores, where it biases the upper layers.
  bool flat_top_compensation = true;
};

/// The G function: the mean projected area of the elements in the (unit) ray direction @c dir, per unit reference
/// area. The interception rate (hits per metre of ray) is G * density.
inline double gFunction(ElementDistribution distribution, const Eigen::Vector3d &dir)
{
  switch (distribution)
  {
  case ElementDistribution::VerticalCylinder:
    return std::sqrt(dir[0] * dir[0] + dir[1] * dir[1]);  // sin(zenith angle)
  case ElementDistribution::Spherical:
  default:
    return 0.5;
  }
}

/// Render a ray cloud according to the supplied parameters
bool RAYLIB_EXPORT renderCloud(const std::string &cloud_file, const Cuboid &bounds, ViewDirection view_direction,
                               RenderStyle style, double pix_width, const std::string &image_file,
                               const std::string &projection_file, bool mark_origin,
                               const std::string *transform_file = nullptr,
                               const DensityOptions &density_options = DensityOptions());

#if RAYLIB_WITH_TIFF
// save to geotif format using floating-point per-channel colour data. This function passes a projection file in order
// to geolocate the image
bool RAYLIB_EXPORT writeGeoTiffFloat(const std::string &filename, int x, int y, const float *data, double pixel_width,
                                     bool scalar, const std::string &projection_file, double origin_x, double origin_y);
#endif

/// This is used for estimating the per-voxel density of a ray cloud
/// By default density represents the surface area per volume, assuming an unbiased distribution of surface angles
/// It is most effective as a measure of leaf area per volume on vegetation, and is described in:
/// Lowe, Thomas, et al. "Canopy Density Estimation in Perennial Horticulture Crops Using 3D Spinning LiDAR SLAM."
/// arXiv preprint arXiv:2007.15652 (2020).
/// Other element distributions are supported via @c ElementDistribution. Each ray's path length is weighted by the
/// G function of its direction, so the per-voxel maximum likelihood estimate is density = hits / sum_i(G_i L_i),
/// which reduces to the original 2 hits / sum_i(L_i) for the spherical distribution.
struct RAYLIB_EXPORT DensityGrid
{
  static const int min_voxel_hits = 2;

  DensityGrid(const Cuboid &grid_bounds, double vox_width, const Eigen::Vector3i &dims,
              ElementDistribution distribution = ElementDistribution::Spherical, double min_g = 0.0)
    : bounds_(grid_bounds)
    , voxel_width_(vox_width)
    , voxel_dims_(dims)
    , distribution_(distribution)
    , min_g_(min_g)
    , ray_g_(0.5)
  {
    // 3D dimensions are stored in a vector of 32-bit signed ints. It's fairly
    // easy to overflow when calculating voxel counts or 1D indices. Make sure
    // this won't happen.
    const std::size_t voxel_count =
      static_cast<std::size_t>(dims[0]) * static_cast<std::size_t>(dims[1]) * static_cast<std::size_t>(dims[2]);
    if (voxel_count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
      std::cerr << "                        dims[0]: " << dims[0] << '\n';
      std::cerr << "                        dims[1]: " << dims[1] << '\n';
      std::cerr << "                        dims[2]: " << dims[2] << '\n';
      std::cerr << "    dims[0] * dims[1] * dims[2]: " << voxel_count << '\n';
      std::cerr << "std::numeric_limits<int>::max(): " << std::numeric_limits<int>::max() << '\n';
      throw std::runtime_error("Dimensions given to ray::DensityGrid are too large. "
                               "Their product must be less than std::numeric_limits<int>::max().");
    }

    // Allocate.
    voxels_.resize(dims[0] * dims[1] * dims[2]);
    peaks_.resize(dims[0] * dims[1], std::numeric_limits<double>::lowest());
  }

  /// This specific voxel class represents a density
  class Voxel
  {
  public:
    Voxel() { num_hits_ = num_rays_ = path_length_ = 0.0; }
    /// return the density that the voxel represents
    inline double density() const;
    inline double numerator() const;
    inline double denominator() const;
    /// the densities can be summed element-wise
    inline void operator+=(const Voxel &other);
    /// the densities can be multiplied by a scalar, element-wise
    inline Voxel operator*(float scale) const;
    /// Add a ray which enters the voxel and hits within it. @c length is the G-weighted ray length within the voxel
    inline void addHitRay(float length);
    /// Add a ray which enters and exits the voxel. @c length is the G-weighted ray path length within the voxel
    inline void addMissRay(float length);
    inline const float &numHits() const { return num_hits_; }
    inline float &numHits() { return num_hits_; }
    inline const float &numRays() const { return num_rays_; }
    /// The sum over rays of the G-weighted path length (G_i L_i) within the voxel
    inline const float &pathLength() const { return path_length_; }
    inline float &pathLength() { return path_length_; }

  private:
    float num_hits_;
    float num_rays_;
    float path_length_;  // G-weighted
  };

  /// The highest end point in each column, for flatTopCompensation. Returns false if the file can't be read
  bool calculatePeaks(const std::string &file_name);
  /// This streams in a ray cloud file, and fills in the voxel density information. Returns false if the file can't be
  /// read
  bool calculateDensities(const std::string &file_name);
  void flatTopCompensation();
  /// To void low-ray-count voxels giving unstable density estimates, we fuse with neighbour information
  /// up to a specified minimum number of rays. Specified in DENSITY_MIN_RAYS
  void addNeighbourPriors();
  /// Note, for performance, this index function does not check that the specified indices are in valid bounds.
  /// It is up to the calling function to assure this condition
  inline int getIndex(const Eigen::Vector3i &inds) const;
  inline int getIndexFromPos(const Eigen::Vector3d &pos) const;
  /// Return the vector of density voxels
  inline const std::vector<Voxel> &voxels() const { return voxels_; }
  inline Eigen::Vector3i dimensions() const { return voxel_dims_; }
  inline Cuboid bounds() { return bounds_; }
  inline double voxelWidth() const { return voxel_width_; }
  inline ElementDistribution distribution() const { return distribution_; }
  inline void setPathEstimator(PathEstimator estimator) { path_estimator_ = estimator; }
  /// Rays with G above this value are ignored
  inline void setMaxG(double max_g) { max_g_ = max_g; }
  /// Half-width (m) of the horizontal neighbourhood for pooling and the occlusion correction
  inline void setPoolRadius(double radius) { pool_radius_ = radius; }
  /// Minimum gap between stem surfaces (in stem diameters) assumed by the occlusion correction
  inline void setMinGapRatio(double gap) { min_gap_ratio_ = gap; }
  /// Accumulate the effective width ratio of each ray segment, for applyOverlapCorrection. Requires the stem diameter
  void enableOverlapCorrection();
  /// Pool the voxel densities over neighbourhoods (PooledBeerLambert), and correct them for the space occupied by the
  /// stems (see DensityOptions::overlap_correction). Call after calculateDensities and before addNeighbourPriors.
  /// The neighbourhood is the (2k+1) x (2k+1) voxels centred on each voxel, for k = round(pool_radius / voxel width),
  /// so it spans +-(k + 1/2) voxel widths, e.g. +-0.075 m for the default 0.05 m and 0.05 m voxels.
  void applyOverlapCorrection();
  /// The density of a voxel of this grid. This is voxel.density(), except for the pooled estimator, whose voxels take
  /// the pooled neighbourhood ratio of hits to path length, without the small-sample factor (n-1)/n for the voxel's
  /// own number of rays n.
  inline double voxelDensity(const Voxel &voxel) const;
  /// Correct for the beam footprint widening thin vertical stems, using their (mean) true diameter.
  /// Only applies to the VerticalCylinder distribution
  inline void setFootprint(const BeamFootprint &footprint, double stem_diameter)
  {
    if (distribution_ == ElementDistribution::VerticalCylinder && footprint.enabled() && stem_diameter > 0.0)
    {
      stem_width_ = StemWidthModel(stem_diameter, footprint);
    }
    range_to_outline_ = footprint.range_to_outline;
    setStemDiameter(stem_diameter);
  }
  /// The (mean) true stem diameter. For VerticalCylinder with the Crossing estimator, this extends each hit ray by
  /// the mean chord of the stem's widened outline (effective width w) that remains beyond the reported return,
  /// horizontally, counted as unobstructed path: pi w / 4 when returns are reported at the outline, or
  /// pi (w^2 + d^2) / (8 w) when they are reported at the physical surface (or closest approach) of the stem. This
  /// removes a positive bias for voxels not much larger than the stems: without it, a stem straddling a voxel
  /// boundary blocks rays that would only have passed through its own interior in the neighbouring voxel.
  inline void setStemDiameter(double stem_diameter)
  {
    if (distribution_ == ElementDistribution::VerticalCylinder)
    {
      stem_diameter_ = stem_diameter;
    }
  }
  /// Number of rays ignored because their G value was below min_g
  inline size_t numIgnoredRays() const { return num_ignored_rays_; }
  // used in walking grid only
  inline bool operator()(const Eigen::Vector3i &p, const Eigen::Vector3i &target, double in_length, double out_length,
                         double max_length);

  std::vector<double> peaks_;  // highest points
private:
  Cuboid bounds_;
  std::vector<Voxel> voxels_;
  double voxel_width_;
  Eigen::Vector3i voxel_dims_;
  ElementDistribution distribution_;
  double min_g_;                                       // rays with G below this are ignored
  double max_g_ = std::numeric_limits<double>::max();  // rays with G above this are ignored
  size_t num_ignored_rays_ = 0;
  bool bounded_;
  double ray_g_;                     // G function value of the ray currently being walked
  double range_offset_ = 0.0;        // distance from the ray's original start to its (clipped) start in the grid
  StemWidthModel stem_width_;        // optional beam footprint correction
  double stem_diameter_ = 0.0;       // optional mean stem diameter
  std::vector<float> overlap_path_;  // per voxel, sum over segments of rho * (G rho L), rho = effective width / d
  bool range_to_outline_ = false;    // the return range convention (BeamFootprint::range_to_outline)
  bool pooled_densities_ = false;    // whether applyOverlapCorrection has pooled the voxels
  double pool_radius_ = 0.05;
  double min_gap_ratio_ = 0.0;
  static constexpr double max_solid_fraction_ = 0.4;
  Eigen::Vector3i hit_voxel_;  // voxel of the end point of the ray currently being walked
  double hit_length_ = 0.0;    // length (in voxels) to the end point of the ray currently being walked
  PathEstimator path_estimator_ = PathEstimator::BeerLambert;

  // used when waling the grid
  Eigen::Vector3d source_;
  Eigen::Vector3d dir_;
};

// inline functions
double DensityGrid::Voxel::numerator() const
{
  return (num_rays_ - 1.0) * num_hits_;
}
double DensityGrid::Voxel::denominator() const
{
  const double eps = 1e-10;  // avoid division by 0
  return eps + num_rays_ * path_length_;
}
double DensityGrid::Voxel::density() const
{
  if (num_rays_ <= min_voxel_hits)
  {
    return 0.0;
  }
  const double eps = 1e-10;  // avoid division by 0
  // below -1.0 should be -2.0 when min length is estimated (e.g. when initially air)
  return (num_rays_ - 1.0) * num_hits_ / (eps + num_rays_ * path_length_);
}
void DensityGrid::Voxel::operator+=(const DensityGrid::Voxel &other)
{
  num_hits_ += other.num_hits_;
  num_rays_ += other.num_rays_;
  path_length_ += other.path_length_;
}
double DensityGrid::voxelDensity(const Voxel &voxel) const
{
  if (!pooled_densities_)
  {
    return voxel.density();
  }
  if (voxel.numRays() <= min_voxel_hits)
  {
    return 0.0;
  }
  return voxel.numHits() / (1e-10 + voxel.pathLength());
}
DensityGrid::Voxel DensityGrid::Voxel::operator*(float scale) const
{
  DensityGrid::Voxel voxel;
  voxel.num_hits_ = num_hits_ * scale;
  voxel.num_rays_ = num_rays_ * scale;
  voxel.path_length_ = path_length_ * scale;
  return voxel;
}
void DensityGrid::Voxel::addHitRay(float length)
{
  path_length_ += length;
  num_hits_++;
  num_rays_++;
}
void DensityGrid::Voxel::addMissRay(float length)
{
  path_length_ += length;
  num_rays_++;
}
int DensityGrid::getIndex(const Eigen::Vector3i &inds) const
{
  return inds[0] + inds[1] * voxel_dims_[0] + inds[2] * voxel_dims_[0] * voxel_dims_[1];
}
int DensityGrid::getIndexFromPos(const Eigen::Vector3d &pos) const
{
  Eigen::Vector3d gridspace = (pos - bounds_.min_bound_) / voxel_width_;
  return getIndex(gridspace.cast<int>());
}
inline bool DensityGrid::operator()(const Eigen::Vector3i &p, const Eigen::Vector3i &, double in_length,
                                    double out_length, double max_length)
{
  int index = getIndex(p);
  double end_length = std::min(out_length, max_length);

#if defined FLAT_TOP_COMPENSATION
  int peak_id = p[0] + p[1] * voxel_dims_[0];
  double peak = peaks_[peak_id];
  double in_height = source_[2] + dir_[2] * in_length;
  double end_height = source_[2] + dir_[2] * end_length;
  if (dir_[2] < 0.0 && in_height > peak && end_height <= peak)  // currently only supported on downwards rays
  {
    double t = (in_height - peak) / (in_height - end_height);
    in_length +=
      (end_length - in_length) *
      std::max(0.0,
               std::min(t, 0.99));  // 0.99 because 1 can cause very tiny distances, which can give density outliers
  }
#endif

  // the ratio of the stems' effective width to their diameter at the range of this segment of the ray
  const double rho =
    stem_width_.enabled() ? stem_width_.scale(range_offset_ + 0.5 * (in_length + end_length) * voxel_width_) : 1.0;
  const double g = ray_g_ * rho;
  // the hit voxel is that of the ray's end point. It is stored, as the walked ray may be extended beyond its hit
  double weighted_length;
  if (p == hit_voxel_ && bounded_)
  {
    const double length_in_voxel =
      (path_estimator_ == PathEstimator::Crossing ? out_length : std::min(out_length, hit_length_)) - in_length;
    weighted_length = g * length_in_voxel * voxel_width_;
    voxels_[index].addHitRay(static_cast<float>(weighted_length));
  }
  else
  {
    // end_length, so that rays which end within the voxel without hitting (e.g. zero-alpha ground returns) only
    // count the path length up to their end point
    weighted_length = g * (end_length - in_length) * voxel_width_;
    voxels_[index].addMissRay(static_cast<float>(weighted_length));
  }
  if (!overlap_path_.empty())
  {
    // accumulate rho weighted by G rho L (weighted_length), which is the weighting of the expected hits
    overlap_path_[static_cast<size_t>(index)] += static_cast<float>(rho * weighted_length);
  }
  return false;
}
/// Runs the density estimation pipeline (rays, occlusion correction and neighbour priors) on a grid of
/// @c dims voxels of width @c voxel_width from bounds.min_bound_. The returned grid is padded by @c shift voxels each
/// side, and its contents are shifted by -shift, so that after this call voxel (x, y, z) for 0 <= x < dims[0] etc.
/// covers bounds.min_bound_ + voxel_width * [x, x+1] etc. Returns nullptr on error, including when the cloud file
/// can't be read.
RAYLIB_EXPORT std::unique_ptr<DensityGrid> estimateDensityGrid(const std::string &cloud_file, const Cuboid &bounds,
                                                               double voxel_width, const Eigen::Vector3i &dims,
                                                               const DensityOptions &options, int &shift);

/// The density of the voxel at (shifted) index x, y, z of a grid from estimateDensityGrid, as an average over the
/// voxel (scaling the top voxel of each column when flat top compensation was used). full_dims is grid.dimensions().
/// It uses grid.voxelDensity, so pooled voxels have no small-sample factor
double RAYLIB_EXPORT voxelAverageDensity(const DensityGrid &grid, const DensityGrid::Voxel &voxel, int x, int y, int z,
                                         const Eigen::Vector3i &full_dims, int shift);
}  // namespace ray
#endif  // RAYLIB_RAYRENDERER_H
