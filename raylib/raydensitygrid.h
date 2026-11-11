// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYDENSITYGRID_H
#define RAYLIB_RAYDENSITYGRID_H

#include "raylib/raylibconfig.h"

#include "raycuboid.h"
#include "raystemmodel.h"
#include "rayutils.h"

#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#define FLAT_TOP_COMPENSATION  // use a 2-medium model for the top of the cloud, so only the cloud's density is
                               // estimated, not the air above it

namespace ray
{
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
  /// at 0.05 m voxels, and +-0.125 m for 0.1 m. Smaller resolves finer patchiness, but each stem's shadow within its
  /// own neighbourhood biases the mean of the voxel densities high by about 0.47 d_eff / W for a neighbourhood of width
  /// W (+7% for 0.05 m, +4% for 0.1 m, with a 2.3 cm effective stem width; see self_shadow_correction), whatever the
  /// number of rays. Don't adapt it to the number of hits, as that is a stopping rule that biases the estimate high.
  double pool_radius = 0.05;
  /// The minimum gap between stem surfaces assumed by the occlusion correction, in stem diameters (stem centres are
  /// at least d (1 + min_gap_ratio) apart). More regular spacing leaves less free space between the widened stems.
  double min_gap_ratio = 0.0;
  /// For the pooled estimator with VerticalCylinder and a stem_diameter: correct each neighbourhood's density for
  /// the stems' self-shadowing (experimental, off by default). A stem removes the free path behind it within its own
  /// neighbourhood but not its own hits, so neighbourhoods holding stems have a higher hits / path than their share,
  /// and the mean of the voxel densities (DensityCell::density) is biased high by about s d_eff / W^2, where W^2 is
  /// the neighbourhood's area, s its meanExitDistance (0.473 W for a square) and d_eff the (G rho L weighted) mean
  /// effective stem width (about +7% for 0.05 m voxels, the default pool radius and a 2.3 cm effective width). With
  /// this option the neighbourhood's measured density is divided by 1 + s d_eff / W^2 before the occlusion correction.
  /// The pooled cell density (DensityCell::pooledDensity), which is a path-weighted mean without this bias, is
  /// unchanged. Derived for tall stems in a uniform field of horizontal ray directions; it overcorrects slightly
  /// where stems end (their shadows are shorter) and in very regular stands
  bool self_shadow_correction = false;
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
  /// Whether addNeighbourPriors also fills voxels that no ray reached (true, the original behaviour). When false,
  /// those voxels stay unobserved, so voxels beyond the scanned volume (e.g. below the ground) aren't reported with
  /// their neighbours' densities
  inline void setFillUnobservedPriors(bool fill) { fill_unobserved_priors_ = fill; }
  /// Note, for performance, this index function does not check that the specified indices are in valid bounds.
  /// It is up to the calling function to assure this condition
  inline int getIndex(const Eigen::Vector3i &inds) const;
  inline int getIndexFromPos(const Eigen::Vector3d &pos) const;
  /// Return the vector of density voxels
  inline const std::vector<Voxel> &voxels() const { return voxels_; }
  inline Eigen::Vector3i dimensions() const { return voxel_dims_; }
  inline Cuboid bounds() const { return bounds_; }
  inline double voxelWidth() const { return voxel_width_; }
  inline ElementDistribution distribution() const { return distribution_; }
  inline void setPathEstimator(PathEstimator estimator) { path_estimator_ = estimator; }
  /// Rays with G above this value are ignored
  inline void setMaxG(double max_g) { max_g_ = max_g; }
  /// Half-width (m) of the horizontal neighbourhood for pooling and the occlusion correction
  inline void setPoolRadius(double radius) { pool_radius_ = radius; }
  /// Minimum gap between stem surfaces (in stem diameters) assumed by the occlusion correction
  inline void setMinGapRatio(double gap) { min_gap_ratio_ = gap; }
  /// Correct the pooled neighbourhood densities for the stems' self-shadowing (DensityOptions::self_shadow_correction)
  inline void setSelfShadowCorrection(bool correct) { self_shadow_correction_ = correct; }
  /// The factor that converts a voxel's (corrected) hits to those without the self-shadow correction, for the pooled
  /// cell density. 1 without the correction. @c index is the voxel's index in voxels()
  inline double uncorrectedHitsScale(size_t index) const
  {
    return index < uncorrected_scale_.size() ? static_cast<double>(uncorrected_scale_[index]) : 1.0;
  }
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
  /// Number of rays that crossed the grid and were used by calculateDensities
  inline size_t numUsedRays() const { return num_used_rays_; }
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
  size_t num_used_rays_ = 0;
  bool bounded_;
  double ray_g_;                     // G function value of the ray currently being walked
  double range_offset_ = 0.0;        // distance from the ray's original start to its (clipped) start in the grid
  StemWidthModel stem_width_;        // optional beam footprint correction
  double stem_diameter_ = 0.0;       // optional mean stem diameter
  std::vector<float> overlap_path_;  // per voxel, sum over segments of rho * (G rho L), rho = effective width / d
  bool range_to_outline_ = false;    // the return range convention (BeamFootprint::range_to_outline)
  bool fill_unobserved_priors_ = true;
  bool pooled_densities_ = false;  // whether applyOverlapCorrection has pooled the voxels
  double pool_radius_ = 0.05;
  double min_gap_ratio_ = 0.0;
  bool self_shadow_correction_ = false;
  std::vector<float> uncorrected_scale_;  // per voxel, uncorrected / corrected hits, with the self-shadow correction
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
    // A ray that passes exactly through a voxel edge or corner, or ends exactly on a voxel face, is also walked
    // through a voxel with (close to) zero length. Which of the adjacent voxels that is depends on rounding (e.g. of a
    // ray end clipped onto a voxel face by raysplit), so counting it as one of the voxel's rays would make the ray
    // counts differ by one between a tiled and a single run, though the hits and path lengths are the same. The ray
    // count selects the observed voxels (min_voxel_hits) and the neighbour priors (DENSITY_MIN_RAYS), and gives the
    // small-sample factor (n-1)/n (which changed two tile edge voxels of the pooled estimator by up to 0.44% before
    // it dropped the factor). It adds no path, so it isn't counted
    if (end_length <= in_length + 1e-9)
    {
      return false;
    }
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
std::unique_ptr<DensityGrid> RAYLIB_EXPORT estimateDensityGrid(const std::string &cloud_file, const Cuboid &bounds,
                                                               double voxel_width, const Eigen::Vector3i &dims,
                                                               const DensityOptions &options, int &shift);

/// The density of the voxel at (shifted) index x, y, z of a grid from estimateDensityGrid, as an average over the
/// voxel (scaling the top voxel of each column when flat top compensation was used). full_dims is grid.dimensions().
/// It uses grid.voxelDensity, so pooled voxels have no small-sample factor
double RAYLIB_EXPORT voxelAverageDensity(const DensityGrid &grid, const DensityGrid::Voxel &voxel, int x, int y, int z,
                                         const Eigen::Vector3i &full_dims, int shift);

/// The sums over the observed voxels of one cell (or layer) of voxels, from aggregateDensityCells. A voxel is observed
/// when it has more than DensityGrid::min_voxel_hits rays (after the neighbour priors).
struct RAYLIB_EXPORT DensityCell
{
  /// Sum of the voxels' (corrected) hits. With flat top compensation, the hits of each column's top voxel are scaled
  /// by the fraction of the voxel below the column's highest point, as in voxelAverageDensity
  double hits = 0.0;
  /// Sum of the voxels' G-weighted path lengths (m)
  double weighted_path = 0.0;
  /// Sum of the voxel densities (voxelAverageDensity)
  double density_sum = 0.0;
  int observed = 0;  // number of observed voxels
  int total = 0;     // number of voxels in the cell
  /// The cell's density: the mean of its observed voxel densities (each voxel weighted equally, so by volume). NaN
  /// when no voxel is observed. This is the main output: it represents patchy cells correctly, but each stem's shadow
  /// within its own neighbourhood biases it high, by about 0.47 d_eff / W for neighbourhoods of width W, unless
  /// DensityOptions::self_shadow_correction is set
  inline double density() const
  {
    return observed > 0 ? density_sum / static_cast<double>(observed) : std::numeric_limits<double>::quiet_NaN();
  }
  /// The cell's density pooled over its observed voxels: hits / weighted_path. This is the weighted_path-weighted
  /// mean of the (neighbourhood-pooled and corrected) voxel densities. It is accurate for uniform stems, but much
  /// too low where stems are patchy within the cell, as rays penetrate dense patches less, so they carry less
  /// weight. NaN when no voxel is observed
  inline double pooledDensity() const
  {
    return observed > 0 && weighted_path > 0.0 ? hits / weighted_path : std::numeric_limits<double>::quiet_NaN();
  }
  /// density / pooledDensity, which exceeds 1 where the density varies within the cell (about 2 for 0.3 m patches
  /// in 1 m cells). NaN when either is NaN, or the pooled density is 0
  inline double patchinessRatio() const
  {
    const double pooled = pooledDensity();
    return pooled > 0.0 ? density() / pooled : std::numeric_limits<double>::quiet_NaN();
  }
  /// Mean hits per observed voxel. NaN when no voxel is observed
  inline double meanHitsPerVoxel() const
  {
    return observed > 0 ? hits / static_cast<double>(observed) : std::numeric_limits<double>::quiet_NaN();
  }
  /// The fraction of the cell's voxels that are observed. NaN for a cell without voxels
  inline double observedFraction() const
  {
    return total > 0 ? static_cast<double>(observed) / static_cast<double>(total) :
                       std::numeric_limits<double>::quiet_NaN();
  }
};

/// A regular grid of cells, for aggregateDensityCells
struct RAYLIB_EXPORT CellLayout
{
  /// The minimum corner of cell (0, 0, 0). Its z is a height above the ground, when aggregating with ground heights
  Eigen::Vector3d min_bound = Eigen::Vector3d(0.0, 0.0, 0.0);
  Eigen::Vector3d cell_size = Eigen::Vector3d(1.0, 1.0, 0.1);
  Eigen::Vector3i num_cells = Eigen::Vector3i(1, 1, 1);
  /// The index of cell (x, y, z) in the result of aggregateDensityCells, x fastest (also the C order of (z, y, x))
  inline size_t index(int x, int y, int z) const
  {
    return static_cast<size_t>(x) +
           static_cast<size_t>(num_cells[0]) * (static_cast<size_t>(y) + static_cast<size_t>(num_cells[1]) * z);
  }
};

/// Aggregates the voxels of @c grid (from estimateDensityGrid, with its @c shift) with (shifted) indices in
/// [voxel_min, voxel_max) into the cells of @c layout. Each voxel is in the cell that contains its centre, where the
/// centre's z is its height above @c ground_height(x, y) (the ground height beneath voxel column x, y, in shifted
/// indices) when given, and absolute otherwise. Voxels outside every cell are ignored. This is the common aggregation
/// for the cells and profile of renderCloud and for estimateDensityTile.
/// @c observed_voxel, when given, is called for each observed voxel in a cell, with its (shifted) index, its centre
/// relative to layout.min_bound (z above the ground when ground heights are given), and its voxelAverageDensity.
std::vector<DensityCell> RAYLIB_EXPORT aggregateDensityCells(
  const DensityGrid &grid, int shift, const Eigen::Vector3i &voxel_min, const Eigen::Vector3i &voxel_max,
  const CellLayout &layout, const std::function<double(int, int)> &ground_height = nullptr,
  const std::function<void(const Eigen::Vector3i &, const Eigen::Vector3d &, double)> &observed_voxel = nullptr);
}  // namespace ray
#endif  // RAYLIB_RAYDENSITYGRID_H
