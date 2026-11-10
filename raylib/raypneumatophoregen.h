// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
#ifndef RAYLIB_RAYPNEUMATOPHOREGEN_H
#define RAYLIB_RAYPNEUMATOPHOREGEN_H

#include "raylib/raylibconfig.h"

#include "raycloud.h"
#include "rayrenderer.h"
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

  // spatial arrangement. Random is a Poisson (complete spatial randomness) process. Clustered places stems along
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
};

/// Generates a synthetic pneumatophore field and its ray cloud. Stems are coloured green, ground returns brown
class RAYLIB_EXPORT PneumatophoreGen
{
public:
  /// Generate the stems and the ray cloud
  void generate(const PneumatophoreParams &params, unsigned int seed);

  /// The stems, as x, y, height
  inline const std::vector<Eigen::Vector3d> &stems() const { return stems_; }
  inline const Cloud &cloud() const { return cloud_; }

  /// The true frontal area per unit volume N d, averaged over each height layer within the plot, as a .csv file
  bool writeTruthProfile(const std::string &file_name, double layer_width) const;

  /// The true N d averaged over each cell of size @c cell_size, tiling the plot from (-w/2, -w/2, 0) for plot width w.
  /// Written as a .csv file
  bool writeTruthCells(const std::string &file_name, const Eigen::Vector3d &cell_size) const;

  /// Cast a single ray, returning its length and whether it hit a stem, or the ground
  void castRay(const Eigen::Vector3d &start, const Eigen::Vector3d &dir, double &length, bool &hit_stem,
               bool &hit_ground) const;

private:
  void generateStems(std::mt19937 &gen);
  PneumatophoreParams params_;
  std::vector<Eigen::Vector3d> stems_;
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
