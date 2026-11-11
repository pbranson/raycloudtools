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
/// width_ratio * d. It exceeds 1 because the free path excludes the stems' (widened) cross-sections, which leave less
/// free space than if they overlapped (for which it is exactly 1). Paths start at random free points, as rays enter
/// the stem layer at random points. Computed deterministically by a 2D Monte-Carlo simulation (fixed seed), pooled
/// over 64 independent realisations (each of 1500 stems and 10000 lines), giving a standard error of under 1%.
/// min_gap_ratio is the minimum gap between stem surfaces in stem diameters, so that stem centres are at least
/// d (1 + min_gap_ratio) apart.
/// @c seed sets the simulation's random number sequence (the correction uses the default).
double RAYLIB_EXPORT freePathResponse(double solid_fraction, double width_ratio, double min_gap_ratio = 0.0,
                                      unsigned int seed = 12345);

/// Lookup of freePathResponse, interpolated over a grid of values that are computed when first needed
class RAYLIB_EXPORT StemResponseTable
{
public:
  explicit StemResponseTable(double min_gap_ratio = 0.0)
    : min_gap_ratio_(min_gap_ratio)
  {}
  double response(double solid_fraction, double width_ratio);

private:
  std::map<std::pair<int, int>, double> table_;
  double min_gap_ratio_;
};

/// Tabulated effective width of a stem of given diameter as a function of range, for fast lookup
class RAYLIB_EXPORT StemWidthModel
{
public:
  StemWidthModel() {}
  StemWidthModel(double diameter, const BeamFootprint &footprint);
  /// effective (detectable) width of the stem at the given range from the scanner
  double width(double range) const;
  /// effective width divided by true diameter
  inline double scale(double range) const { return width(range) / diameter_; }
  /// the maximum effective width over ranges 0 to @c max_range
  double maxWidth(double max_range) const;
  inline bool enabled() const { return enabled_; }

private:
  double diameter_ = 0.0;
  BeamFootprint footprint_;
  bool enabled_ = false;
  std::vector<double> table_;  // effective width at regular range intervals
  static constexpr double table_step_ = 0.01;
  static constexpr double table_max_range_ = 100.0;
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
  /// neighbourhood of pool_radius (default +-0.05 m), pooled as total hits over total path length. Pooling avoids the
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
  BeamFootprint footprint;
  double stem_diameter = 0.0;
  /// For VerticalCylinder with a stem_diameter, and the beer_lambert or pooled estimators: correct the free-path
  /// estimate for the space occupied by the stems, as widened by the footprint. The free path excludes the stems'
  /// widened cross-sections, so the measured density is a R(phi, rho) for true density a = N d, solid fraction
  /// phi = N pi d^2 / 4 and widening rho = d_eff / d (see freePathResponse). This is inverted using the density of
  /// a neighbourhood of pool_radius (default +-0.05 m) horizontally. The solid fraction is limited to 0.4 (below the
  /// random packing limit of 0.547 for non-overlapping discs), so N d is at most 1.6 / (pi d).
  bool overlap_correction = true;
  /// Half-width (m) of the horizontal neighbourhood used by the pooled estimator and the occlusion correction.
  /// Smaller resolves finer patchiness, but with fewer hits per neighbourhood, which biases the estimate high in
  /// sparse areas (about +6% at 30 stems/m^2 for 0.05 m, +3% for 0.1 m). Don't adapt it to the number of hits, as
  /// that is a stopping rule that biases the estimate high.
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

  void calculatePeaks(const std::string &file_name);
  /// This streams in a ray cloud file, and fills in the voxel density information
  void calculateDensities(const std::string &file_name);
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
  inline void setMinGapRatio(double gap)
  {
    min_gap_ratio_ = gap;
    response_ = StemResponseTable(gap);
  }
  /// Accumulate the effective width ratio of each ray segment, for applyOverlapCorrection. Requires the stem diameter
  void enableOverlapCorrection();
  /// Pool the voxel densities over neighbourhoods (PooledBeerLambert), and correct them for the space occupied by the
  /// stems (see DensityOptions::overlap_correction). Call after calculateDensities and before addNeighbourPriors.
  void applyOverlapCorrection();
  /// Correct for the beam footprint widening thin vertical stems, using their (mean) true diameter.
  /// Only applies to the VerticalCylinder distribution
  inline void setFootprint(const BeamFootprint &footprint, double stem_diameter)
  {
    if (distribution_ == ElementDistribution::VerticalCylinder && footprint.enabled() && stem_diameter > 0.0)
    {
      stem_width_ = StemWidthModel(stem_diameter, footprint);
    }
    setStemDiameter(stem_diameter);
  }
  /// The (mean) true stem diameter. For VerticalCylinder with the Crossing estimator, this extends each hit ray by
  /// the stem's mean internal chord (pi/4 of its effective width, horizontally), counted as unobstructed path. This
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
  std::vector<float> overlap_path_;  // per voxel, sum of the path-weighted overlap coefficient k of each segment
  StemResponseTable response_;       // free path response to the stems, for the occlusion correction
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

  double g = ray_g_;
  if (stem_width_.enabled())
  {
    // the effective width of the stems at the range of this segment of the ray
    const double range = range_offset_ + 0.5 * (in_length + end_length) * voxel_width_;
    g *= stem_width_.scale(range);
  }
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
    // accumulate the path-weighted ratio of effective stem width to diameter
    const double rho =
      stem_width_.enabled() ? stem_width_.scale(range_offset_ + 0.5 * (in_length + end_length) * voxel_width_) : 1.0;
    overlap_path_[static_cast<size_t>(index)] += static_cast<float>(rho * weighted_length);
  }
  return false;
}
/// Runs the density estimation pipeline (rays, occlusion correction and neighbour priors) on a grid of
/// @c dims voxels of width @c voxel_width from bounds.min_bound_. The returned grid is padded by @c shift voxels each
/// side, and its contents are shifted by -shift, so that after this call voxel (x, y, z) for 0 <= x < dims[0] etc.
/// covers bounds.min_bound_ + voxel_width * [x, x+1] etc. Returns nullptr on error.
RAYLIB_EXPORT std::unique_ptr<DensityGrid> estimateDensityGrid(const std::string &cloud_file, const Cuboid &bounds,
                                                               double voxel_width, const Eigen::Vector3i &dims,
                                                               const DensityOptions &options, int &shift);

/// The density of the voxel at (shifted) index x, y, z of a grid from estimateDensityGrid, as an average over the
/// voxel (scaling the top voxel of each column when flat top compensation was used). full_dims is grid.dimensions()
double RAYLIB_EXPORT voxelAverageDensity(const DensityGrid &grid, const DensityGrid::Voxel &voxel, int x, int y, int z,
                                         const Eigen::Vector3i &full_dims, int shift);
}  // namespace ray
#endif  // RAYLIB_RAYRENDERER_H
