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

/// Options controlling the voxel density estimate
struct RAYLIB_EXPORT DensityOptions
{
  ElementDistribution distribution = ElementDistribution::Spherical;
  /// For VerticalCylinder only: rays whose sin(zenith angle) is below this are ignored entirely, since near-vertical
  /// rays carry almost no information on vertical elements, and their few (end cap) hits would give outliers.
  /// Ignoring rays based on their direction alone does not bias the estimate.
  double min_sin_zenith = 0.1;
  /// Optional .csv file to write the mean density per horizontal layer of voxels
  std::string profile_file;
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
  inline Eigen::Vector3i dimensions() { return voxel_dims_; }
  inline Cuboid bounds() { return bounds_; }
  inline double voxelWidth() const { return voxel_width_; }
  inline ElementDistribution distribution() const { return distribution_; }
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
  double min_g_;  // rays with G below this are ignored
  size_t num_ignored_rays_ = 0;
  bool bounded_;
  double ray_g_;  // G function value of the ray currently being walked

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
inline bool DensityGrid::operator()(const Eigen::Vector3i &p, const Eigen::Vector3i &target, double in_length,
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

  if (p == target && bounded_)
  {
    double length_in_voxel = end_length - in_length;
    voxels_[index].addHitRay(static_cast<float>(ray_g_ * length_in_voxel * voxel_width_));
  }
  else
  {
    // end_length, so that rays which end within the voxel without hitting (e.g. zero-alpha ground returns) only
    // count the path length up to their end point
    voxels_[index].addMissRay(static_cast<float>(ray_g_ * (end_length - in_length) * voxel_width_));
  }
  return false;
}
}  // namespace ray
#endif  // RAYLIB_RAYRENDERER_H
