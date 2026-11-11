// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#ifndef RAYLIB_RAYPNEUMATOPHOREGEN_H
#define RAYLIB_RAYPNEUMATOPHOREGEN_H

#include "raylib/raylibconfig.h"

#include "raycloud.h"
#include "raystemmodel.h"
#include "rayutils.h"

#include <random>

namespace ray
{
/// Parameters of a synthetic field of mangrove pneumatophores (vertical cylinders on flat ground at z = 0), scanned
/// by a hand-held lidar walked along a lawnmower path. Used to validate density estimation.
struct RAYLIB_EXPORT PneumatophoreParams
{
  // stems
  double stems_per_m2 = 100.0;  // mean number of stems per square metre of ground
  double diameter = 0.01;       // stem diameter (m)
  double mean_height = 0.25;    // mean stem height (m)
  double height_sd = 0.08;      // standard deviation of stem height (m), truncated to be positive
  /// The distribution of stem heights. Normal: mean_height and height_sd, truncated below 0.02 m (a height_sd of 0
  /// gives stems of one height, and a height above the scanner gives stems without tops in view). Uniform: between
  /// min_height and max_height
  enum class HeightDistribution
  {
    Normal,
    Uniform
  };
  HeightDistribution height_distribution = HeightDistribution::Normal;
  double min_height = 0.1;  // for the uniform height distribution (m)
  double max_height = 0.4;  // for the uniform height distribution (m)

  // Stems are placed without overlapping each other (a hard-core process, by random sequential addition), with the
  // exact number of stems for the density. Otherwise stem positions are independent (Poisson), so stems may overlap
  // and the number of stems is Poisson distributed. Hard-core stems limit the solid fraction N pi d^2 / 4 to the
  // random packing limit of about 0.547
  bool hard_core = true;
  /// When hard_core: the minimum gap between neighbouring stems' surfaces, as a ratio of the stem diameter, so stem
  /// centres are at least d (1 + min_gap_ratio) apart. This limits the solid fraction to
  /// PneumatophoreGen::max_solid_fraction / (1 + min_gap_ratio)^2
  double min_gap_ratio = 0.0;

  // sparse tree trunks: tall vertical cylinders at random positions (not overlapping each other), which the
  // pneumatophores don't overlap. Trunk returns are coloured brown-grey
  double trunks_per_m2 = 0.0;
  double trunk_diameter = 0.2;      // mean trunk diameter (m)
  double trunk_diameter_sd = 0.05;  // standard deviation of trunk diameter (m), truncated at 0.05 m
  double trunk_height = 10.0;

  // spatial arrangement. Random places stems uniformly at random over the ground. Clustered places stems along
  // cable roots radiating from trees, with the same mean density
  bool clustered = false;
  double trees_per_m2 = 0.04;    // tree density, when clustered
  int roots_per_tree = 12;       // cable roots per tree
  double root_length = 6.0;      // length of each cable root (m)
  double lateral_spread = 0.03;  // standard deviation of stem offsets from the cable root line (m)
  // Patchy places stems at stems_per_m2 within circular patches, and at matrix_ratio * stems_per_m2 elsewhere, as in
  // degraded areas. Patch centres are random, so the patches overlap
  bool patchy = false;
  double patch_radius = 0.5;  // radius of each patch (m)
  double patch_cover = 0.3;   // fraction of the ground covered by patches
  double matrix_ratio = 0.0;  // stem density between patches, relative to within them

  // scene and scan
  double plot_width = 6.0;       // the scanner walks within this square plot, centred on 0,0
  double margin = 8.0;           // stems are generated this far beyond the plot, to avoid edge effects
  double scanner_height = 1.2;   // height of the hand-held scanner (m)
  double walk_spacing = 1.0;     // spacing of the lawnmower walking lines (m)
  int num_rays = 2000000;        // number of rays
  double max_range = 30.0;       // range of rays with no return (m)
  bool ground_unbounded = true;  // ground returns are given zero alpha (as if perfectly classified and removed)
  BeamFootprint footprint;       // beam model. Disabled by default (an infinitely thin beam)
  // range of ray zenith angles (degrees from straight down). Rays are uniformly distributed over this band of
  // directions. The default omits upward rays, which see only the sky or canopy here
  double min_zenith_deg = 0.0;
  double max_zenith_deg = 101.5;
  // measurement errors, drawn from a separate random sequence so that the scene and the ideal rays are unchanged
  /// standard deviation (m) of Gaussian noise added to the range of every return (stem, trunk or ground), along its ray
  double range_noise = 0.0;
  /// standard deviation (m), per axis, of a slowly varying trajectory (pose) error: the recorded start and end of each
  /// ray are both shifted by an offset interpolated linearly between independent Gaussian values every
  /// pose_noise_interval metres along the walk (a translation error, as from SLAM drift; rotation is not modelled)
  double pose_noise = 0.0;
  double pose_noise_interval = 1.0;
};

/// Generates a synthetic pneumatophore field and its ray cloud. Stem returns are coloured green, ground returns
/// brown and rays with no return grey. The scene depends only on the parameters and seed, not on the platform
/// or standard library
class RAYLIB_EXPORT PneumatophoreGen
{
public:
  /// The highest solid fraction N pi d^2 / 4 at which non-overlapping stems are generated, just below the random
  /// sequential addition (jamming) limit of 0.547 for discs
  static constexpr double max_solid_fraction = 0.5;

  /// Generate the stems and the ray cloud
  void generate(const PneumatophoreParams &params, unsigned int seed);

  /// The stems, as x, y, height
  inline const std::vector<Eigen::Vector3d> &stems() const { return stems_; }
  /// The trunks, as x, y, diameter
  inline const std::vector<Eigen::Vector3d> &trunks() const { return trunks_; }
  /// The trunks as a .csv file of x, y, diameter
  bool writeTrunks(const std::string &file_name) const;
  inline const Cloud &cloud() const { return cloud_; }

  /// The true frontal area per unit volume N d, averaged over each height layer within the plot, as a .csv file
  bool writeTruthProfile(const std::string &file_name, double layer_width) const;

  /// The stems as a .csv file of x, y, height, with the stem diameter
  bool writeStems(const std::string &file_name) const;

  /// The true N d averaged over each cell of size @c cell_size, tiling the plot from (-w/2, -w/2, 0) for plot width w.
  /// Written as a .csv file
  bool writeTruthCells(const std::string &file_name, const Eigen::Vector3d &cell_size) const;

  /// Cast a single ray, returning its length and whether it hit a stem, or the ground
  void castRay(const Eigen::Vector3d &start, const Eigen::Vector3d &dir, double &length, bool &hit_stem,
               bool &hit_ground, bool &hit_trunk) const;

private:
  void generateStems(std::mt19937 &gen);
  PneumatophoreParams params_;
  std::vector<Eigen::Vector3d> stems_;
  std::vector<Eigen::Vector3d> trunks_;
  Cloud cloud_;
  // bucket grid of stem indices, for accelerating ray intersection
  std::vector<std::vector<int>> buckets_;
  double bucket_width_ = 0.25;
  double field_min_ = 0.0;
  int num_buckets_ = 0;
  double max_height_ = 0.0;
  StemWidthModel stem_width_;
};
}  // namespace ray

#endif  // RAYLIB_RAYPNEUMATOPHOREGEN_H
