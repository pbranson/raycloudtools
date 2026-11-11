// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include <cstdlib>
#include <vector>

#include "gtest/gtest.h"

#include "raylib/raylibconfig.h"

#if RAYLIB_WITH_TBB
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/parallel_for.h>
#endif  // RAYLIB_WITH_TBB

#include "raylib/raycloud.h"
#include "raylib/raydensity.h"
#include "raylib/rayforeststructure.h"
#include "raylib/raylaz.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/raypneumatophoregen.h"
#include "raylib/rayrenderer.h"
#include "raylib/raysplitter.h"
#include "raylib/raythreads.h"

#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <sstream>

#include "environment.h"

/// Raycloud testing framework. In each test, the statistics of the resulting clouds are compared to the statistics
/// of the cloud when it was confirmed to be operating correctly.
namespace raytest
{
/// Issues the specified system command, including the required prefix on non-windows systems.
int command(const std::string &system_command)
{
  return raytest::Environment::Command(system_command);
}

/// Issues the command to copy a file, which is a platform dependent system command.
int copy(const std::string &copy_command)
{
#ifdef _WIN32
  return system("copy " + copy_command);
#else
  return system(("cp " + copy_command).c_str());
#endif  // _WIN32
}

/// Compare the statistical (1st and 2nd order) moments of the two ray clouds. This almost surely
/// detects differing clouds, and always equal clouds, given a tolerance @c eps.
std::size_t compareMoments(const Eigen::ArrayXd &m1, const std::vector<double> &m2, double eps = 0.1)
{
  size_t i;
  for (i = 0; i < m2.size(); i++)
  {
    EXPECT_NEAR(m1[i], m2[i], eps);
    if (std::abs(m1[i] - m2[i]) > eps)
    {
      return i;
    }
  }
  return i;
}

/// Compare the statistical (1st and 2nd order) moments of the two ray clouds. This almost surely
/// detects differing clouds, and always equal clouds, given a tolerance @c eps.
std::size_t compareMomentsPercentageError(const Eigen::ArrayXd &m1, const std::vector<double> &m2,
                                          double percentage = 5.0)
{
  size_t i;
  for (i = 0; i < m2.size(); i++)
  {
    double eps = std::abs(m2[i]) * percentage * 0.01;
    EXPECT_NEAR(m1[i], m2[i], eps);
    if (std::abs(m1[i] - m2[i]) > eps)
    {
      return i;
    }
  }
  return i;
}

/// Creates two copies of the same room with a rotational difference, then aligns the first onto the second
TEST(Basic, RayAlign)
{
  EXPECT_EQ(command("raycreate room 1"), 0);
  EXPECT_EQ(copy("room.ply room2.ply"), 0);
  EXPECT_EQ(command("rayrotate room2.ply 0,0,35"), 0);
  EXPECT_EQ(command("rayalign room.ply room2.ply"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("room_aligned.ply"));
  EXPECT_EQ(
    compareMoments(cloud.getMoments(),
                   { -0.0618268, -0.077552, 0.0531072, 7.58334e-08, 7.97642e-08, 1.93877e-08, -0.180532, -0.219257,
                     0.0654452,  2.47241,   2.08183,   1.28226,     17.539,      10.1994,     0.304682,  0.761892,
                     0.429502,   0.987362,  0.318932,  0.225742,    0.389901,    0.111705 }),
    cloud.getMoments().size());
}
/// Colours a room according to the normal direction of the surfaces, comparing to the expected results
TEST(Basic, RayColour)
{
  EXPECT_EQ(command("raycreate room 1"), 0);
  EXPECT_EQ(command("raycolour room.ply normal"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("room_coloured.ply"));
  EXPECT_EQ(compareMoments(cloud.getMoments(), { -0.108066, -0.0410134, 0.052168, 7.05134e-08, 8.45038e-08, 1.93877e-08,
                                                 -0.276144, -0.0760758, 0.065631, 2.42455,     2.13738,     1.28226,
                                                 17.539,    10.1994,    0.497919, 0.496369,    0.490293,    0.987362,
                                                 0.248361,  0.203648,   0.385192, 0.111705 }),
            cloud.getMoments().size());
}

/// Creates two rooms, with different transformations, then combines them, and compares to the expected result.
TEST(Basic, RayCombine)
{
  EXPECT_EQ(command("raycreate room 1"), 0);
  EXPECT_EQ(copy("room.ply room2.ply"), 0);
  EXPECT_EQ(command("raytranslate room2.ply 0,0,1"), 0);
  EXPECT_EQ(command("rayrotate room2.ply 0,0,35"), 0);
  EXPECT_EQ(command("raycombine min room.ply room2.ply 1 rays"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("room_combined.ply"));
  EXPECT_EQ(compareMoments(cloud.getMoments(),
                           { -0.0867714, -0.0679941, 0.546619, 0.0215326, 0.0272819, 0.499969, -0.305657, -0.186353,
                             0.582642,   2.95777,    2.47531,  1.63323,   17.4967,   10.1789,  0.305355,  0.763356,
                             0.427376,   0.979005,   0.318409, 0.225661,  0.389366,  0.143369 }),
            cloud.getMoments().size());
}

/// Creates a building with random seed 1, and compares to the expected results
TEST(Basic, RayCreate)
{
  EXPECT_EQ(command("raycreate building 1"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("building.ply"));
  EXPECT_EQ(
    compareMoments(cloud.getMoments(), { -3.168,   16.5472,  7.04175,  3.28539,  18.0871,  2.96654,  -3.19551, 16.564,
                                         7.30826,  4.14359,  18.2463,  3.34431,  935.715,  540.236,  0.499997, 0.500408,
                                         0.429416, 0.998894, 0.372134, 0.372389, 0.390499, 0.0332398 }),
    cloud.getMoments().size());
}

/// Creates a forest and decimates it, comparing to the expected result
TEST(Basic, RayDecimate)
{
  EXPECT_EQ(command("raycreate forest 1"), 0);
  EXPECT_EQ(command("raydecimate forest.ply 10 cm"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("forest_decimated.ply"));
  // Below does not compare the time values (or the colour values, which are based on time here)
  // because spatial decimation does not constraint which time it picks points from.
  EXPECT_EQ(compareMoments(cloud.getMoments(),
                           { -0.222571, 1.08156, 1.67264,  6.00755,  5.78731,  0.508713, -0.202668, 1.09517,
                             2.6238,    6.0285,  5.85715,  3.22093,  69.0574,  35.2775,  0.48969,   0.498403,
                             0.443549,  1,       0.379062, 0.366963, 0.389535, 0 }),
            cloud.getMoments().size());
}

/// Creates a room, and calls denoise using a fixed distance threshols, and compares to expected result
TEST(Basic, RayDenoise)
{
  EXPECT_EQ(command("raycreate room 1"), 0);
  EXPECT_EQ(command("raydenoise room.ply 3 cm"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("room_denoised.ply"));
  EXPECT_EQ(compareMoments(cloud.getMoments(), { -0.108066, -0.0410134, 0.052168, 8.67026e-08, 8.81787e-08, 2.24394e-08,
                                                 -0.464107, -0.113806,  0.161496, 2.82122,     2.34281,     1.35279,
                                                 17.81,     10.2005,    0.297047, 0.758802,    0.440232,    0.975166,
                                                 0.317215,  0.226682,   0.390971, 0.155618 }),
            cloud.getMoments().size());
}

#if RAYLIB_WITH_LAS
/// Create a forest as .ply, load it as .ply, save it as .laz, load it as
/// .laz, make sure the clouds loaded as .ply and .laz match.
TEST(Basic, RayLaz)
{
  // Write .ply.
  EXPECT_EQ(command("raycreate forest 1"), 0);

  // Read .ply.
  ray::Cloud ply_cloud;
  EXPECT_TRUE(ply_cloud.load("forest.ply"));

  // Check.
  const auto size = ply_cloud.ends.size();
  EXPECT_EQ(ply_cloud.ends.size(), size);
  EXPECT_EQ(ply_cloud.times.size(), size);
  EXPECT_EQ(ply_cloud.colours.size(), size);

  // Write .laz.
  EXPECT_TRUE(ray::writeLas("forest.laz", ply_cloud.ends, ply_cloud.times, ply_cloud.colours));

  // Read .laz.
  ray::Cloud laz_cloud;
  const double max_intensity = 255.0;
  EXPECT_TRUE(ray::readLas("forest.laz", laz_cloud.ends, laz_cloud.times, laz_cloud.colours, max_intensity));

  // Check.
  EXPECT_EQ(laz_cloud.ends.size(), size);
  EXPECT_EQ(laz_cloud.times.size(), size);
  EXPECT_EQ(laz_cloud.colours.size(), size);

  // Compare.
  for (size_t index = 0; index < size; ++index)
  {
    EXPECT_NEAR(ply_cloud.ends[index].x(), laz_cloud.ends[index].x(), 1e-4);
    EXPECT_NEAR(ply_cloud.ends[index].y(), laz_cloud.ends[index].y(), 1e-4);
    EXPECT_NEAR(ply_cloud.ends[index].z(), laz_cloud.ends[index].z(), 1e-4);
    EXPECT_EQ(ply_cloud.colours[index].red, laz_cloud.colours[index].red);
    EXPECT_EQ(ply_cloud.colours[index].green, laz_cloud.colours[index].green);
    EXPECT_EQ(ply_cloud.colours[index].blue, laz_cloud.colours[index].blue);
    EXPECT_EQ(ply_cloud.colours[index].alpha, laz_cloud.colours[index].alpha);
    EXPECT_EQ(ply_cloud.times[index], laz_cloud.times[index]);
  }
}
#endif  // RAYLIB_WITH_LAS

/// Creates two rooms, the second is decimated and transformed, then rayrestore is called to apply this transformation
/// to the first (high resolution) room
TEST(Basic, RayRestore)
{
  EXPECT_EQ(command("raycreate room 1"), 0);
  EXPECT_EQ(copy("room.ply room2.ply"), 0);
  EXPECT_EQ(command("raydecimate room2.ply 10 cm"), 0);
  EXPECT_EQ(command("raytranslate room2_decimated.ply 1,2,3"), 0);
  EXPECT_EQ(command("rayrotate room2_decimated.ply 0,0,-50"), 0);
  EXPECT_EQ(command("rayrestore room2_decimated.ply 10 cm room.ply"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("room_restored.ply"));
  EXPECT_EQ(compareMoments(cloud.getMoments(),
                           { 2.07399,  0.575952, 3.05217,  7.85442e-08, 7.70963e-08, 1.93877e-08, 1.9391,   0.682169,
                             3.06563,  2.10068,  2.45642,  1.28226,     17.539,      10.1994,     0.304682, 0.761892,
                             0.429502, 0.987362, 0.318932, 0.225742,    0.389901,    0.111705 }),
            cloud.getMoments().size());
}

/// Creates a forest and rotates it in all three axes, comparing to the expected result
TEST(Basic, RayRotate)
{
  EXPECT_EQ(command("raycreate forest 1"), 0);
  EXPECT_EQ(command("rayrotate forest.ply 10,20,30"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("forest.ply"));
  EXPECT_EQ(
    compareMoments(cloud.getMoments(), { -0.254879, 0.846076, 2.02322,  5.62648,  5.68622,  2.56306, 0.266873, 0.772016,
                                         3.2888,    5.54595,  5.69021,  4.28394,  62.683,   36.1903, 0.514327, 0.504407,
                                         0.413534,  1,        0.372377, 0.365965, 0.391709, 0 }),
    cloud.getMoments().size());
}

/// Creates a room and smooths this ray cloud, comparing to the expected result
TEST(Basic, RaySmooth)
{
  EXPECT_EQ(command("raycreate room 1"), 0);
  EXPECT_EQ(command("raysmooth room.ply"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("room_smooth.ply"));
  EXPECT_EQ(
    compareMoments(cloud.getMoments(),
                   { -0.108066, -0.0410134, 0.052168, 7.05134e-08, 8.45038e-08, 1.93877e-08, -0.27615, -0.0761079,
                     0.0656267, 2.42413,    2.13691,  1.28163,     17.539,      10.1994,     0.304682, 0.761892,
                     0.429502,  0.987362,   0.318932, 0.225742,    0.389901,    0.111705 }),
    cloud.getMoments().size());
}

/// Creates a room, then splits it around a plane, comparing agaisnt the expected result
TEST(Basic, RaySplit)
{
  EXPECT_EQ(command("raycreate room 1"), 0);
  EXPECT_EQ(command("raysplit room.ply plane 0,0.1,1.5"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("room_outside.ply"));
  EXPECT_EQ(
    compareMoments(cloud.getMoments(), { -0.467731, 1.05075,  1.43662,  2.20441,  1.60162, 0.106775, -0.77974, 1.03139,
                                         1.57353,   3.67521,  2.64766,  0.485084, 17.3995, 10.279,   0.311066, 0.759795,
                                         0.425206,  0.951355, 0.321609, 0.226785, 0.39073, 0.215125 }),
    cloud.getMoments().size());
}

/// Call functions defined in raythreads and then, if TBB is enabled, try to use
/// TBB and check if it actually works.
TEST(Basic, RayThreads)
{
  const int original_available_threads = ray::Threads::availableThreads();
  const int original_recommended_thread_count = ray::Threads::recommendedThreadCount();
  for (int thread_count = -2; thread_count <= 2; ++thread_count)
  {
    if (thread_count >= -1)
    {
      EXPECT_NO_THROW(ray::Threads::init(thread_count));
    }
    const int available_threads = ray::Threads::availableThreads();
    const int recommended_thread_count = ray::Threads::recommendedThreadCount();
    EXPECT_EQ(available_threads, original_available_threads);
    EXPECT_EQ(recommended_thread_count, original_recommended_thread_count);
#if RAYLIB_WITH_TBB
    // Check that the global limit was applied.
    const int tbb_thread_count =
      static_cast<int>(oneapi::tbb::global_control::active_value(oneapi::tbb::global_control::max_allowed_parallelism));
    switch (thread_count)
    {
    case -2:
      EXPECT_EQ(tbb_thread_count, original_available_threads);
      break;
    case -1:
      EXPECT_EQ(tbb_thread_count, original_available_threads);
      break;
    case 0:
      EXPECT_EQ(tbb_thread_count, original_recommended_thread_count);
      break;
    case 1:
      EXPECT_EQ(tbb_thread_count, 1);
      break;
    default:
      EXPECT_EQ(tbb_thread_count, std::min(thread_count, original_available_threads));
      break;
    }
    // Check that the threads can actually be used.
    std::atomic<size_t> active_threads = 0;
    std::atomic<size_t> active_threads_max = 0;
    oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, 1024),
      [&active_threads, &active_threads_max](const oneapi::tbb::blocked_range<size_t> &) {
        // Spin a bit while other threads start.
        ++active_threads;
        for (size_t count = 0; count < 9999; ++count)
        {
          std::sqrt(rand());
        }
        // Track max number of active threads.
        while (true)
        {
          auto old_active_threads_max = active_threads_max.load();
          const auto new_active_threads_max = std::max({
            active_threads.load(),
            old_active_threads_max,
          });
          if (active_threads_max.compare_exchange_weak(old_active_threads_max, new_active_threads_max))
          {
            break;
          }
        }
        // Spin a bit while other threads track max number of active threads.
        for (size_t count = 0; count < 9999; ++count)
        {
          std::sqrt(rand());
        }
        --active_threads;
      },
      oneapi::tbb::static_partitioner());
    // This isn't particularly reliable so just check it's in the expected
    // range.
    EXPECT_GE(active_threads_max, 1);
    EXPECT_LE(active_threads_max, tbb_thread_count);
#endif  // RAYLIB_WITH_TBB
  }
}

/// Creates a room and runs raytransients, comparing the identified transients ray cloud to the expected results
TEST(Basic, RayTransients)
{
  EXPECT_EQ(command("raycreate room 2"), 0);
  EXPECT_EQ(command("raytransients min room.ply 1 rays"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("room_transient.ply"));
  EXPECT_EQ(compareMoments(cloud.getMoments(),
                           { -1.05406,  -0.240721, -0.0629182, 5.05649e-08, 3.32941e-08, 2.54759e-08, 0.268724,
                             -0.136746, -0.596782, 1.04798,    0.921776,    0.527205,    32.1452,     6.7491,
                             0.205871,  0.395641,  0.884296,   1,           0.225501,    0.296487,    0.153923,
                             0 }),
            cloud.getMoments().size());
}

/// Creates a forest and translates it in all three axes, comparing to the expected result
TEST(Basic, RayTranslate)
{
  EXPECT_EQ(command("raycreate forest 1"), 0);
  EXPECT_EQ(command("raytranslate forest.ply 10,20,30"), 0);
  ray::Cloud cloud;
  EXPECT_TRUE(cloud.load("forest.ply"));
  EXPECT_EQ(
    compareMoments(cloud.getMoments(), { 9.66298,  21.3454, 31.7177,  6.0926,   5.75511,  0.56438, 9.69155,  21.3605,
                                         33.0883,  6.10555, 5.82564,  3.20507,  62.683,   36.1903, 0.514327, 0.504407,
                                         0.413534, 1,       0.372377, 0.365965, 0.391709, 0 }),
    cloud.getMoments().size());
}

/// Generates a ray cloud of rays at a fixed zenith angle through a random (Poisson) field of vertical cylinders
/// of diameter @c diameter and @c stems_per_m2 per square metre. Returns the cloud, with each ray ending at its
/// first cylinder intersection (bounded) or at its maximum length (unbounded)
ray::Cloud generateCylinderFieldCloud(double zenith_deg, double stems_per_m2, double diameter, int num_rays,
                                      unsigned seed)
{
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  const double field_min = -2.0, field_max = 6.0;  // cylinders extend beyond where rays can reach
  const double cell = 0.5;                         // bucket grid of cylinder centres, to accelerate intersection
  const int cells = static_cast<int>((field_max - field_min) / cell);
  std::vector<std::vector<Eigen::Vector2d>> buckets(static_cast<size_t>(cells * cells));
  const int num_stems = static_cast<int>(stems_per_m2 * (field_max - field_min) * (field_max - field_min));
  for (int i = 0; i < num_stems; i++)
  {
    Eigen::Vector2d c(field_min + (field_max - field_min) * uniform(gen),
                      field_min + (field_max - field_min) * uniform(gen));
    int bx = std::min(cells - 1, static_cast<int>((c[0] - field_min) / cell));
    int by = std::min(cells - 1, static_cast<int>((c[1] - field_min) / cell));
    buckets[static_cast<size_t>(bx + cells * by)].push_back(c);
  }
  const double radius = 0.5 * diameter;
  const double max_length = 1.5;
  const double theta = zenith_deg * ray::kPi / 180.0;
  ray::Cloud cloud;
  int ray_count = 0;
  while (ray_count < num_rays)
  {
    const Eigen::Vector3d start(0.5 + 3.0 * uniform(gen), 0.5 + 3.0 * uniform(gen), 1.5 + uniform(gen));
    const double phi = 2.0 * ray::kPi * uniform(gen);
    const Eigen::Vector3d dir(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), -std::cos(theta));
    const Eigen::Vector2d s2(start[0], start[1]);
    const Eigen::Vector2d d2(dir[0], dir[1]);
    const double h2 = d2.squaredNorm();
    const Eigen::Vector2d e2 = s2 + d2 * max_length;
    const int x0 = static_cast<int>((std::min(s2[0], e2[0]) - radius - field_min) / cell);
    const int x1 = static_cast<int>((std::max(s2[0], e2[0]) + radius - field_min) / cell);
    const int y0 = static_cast<int>((std::min(s2[1], e2[1]) - radius - field_min) / cell);
    const int y1 = static_cast<int>((std::max(s2[1], e2[1]) + radius - field_min) / cell);
    double hit_length = max_length;
    bool inside = false;
    for (int x = x0; x <= x1; x++)
    {
      for (int y = y0; y <= y1; y++)
      {
        for (const auto &c : buckets[static_cast<size_t>(x + cells * y)])
        {
          const Eigen::Vector2d to = c - s2;
          if (to.squaredNorm() < radius * radius)
          {
            inside = true;
          }
          if (h2 <= 0.0)
          {
            continue;
          }
          // solve |s2 + t d2 - c| = radius for the smallest t > 0
          const double b = to.dot(d2);
          const double disc = b * b - h2 * (to.squaredNorm() - radius * radius);
          if (disc < 0.0)
          {
            continue;
          }
          const double t = (b - std::sqrt(disc)) / h2;
          if (t > 0.0 && t < hit_length)
          {
            hit_length = t;
          }
        }
      }
    }
    if (inside)
    {
      continue;  // a scanner can't start inside a stem
    }
    const bool hit = hit_length < max_length;
    cloud.addRay(start, start + dir * hit_length, static_cast<double>(ray_count),
                 ray::RGBA(127, 127, 127, hit ? 255 : 0));
    ray_count++;
  }
  return cloud;
}

/// Returns the mean voxel density over the voxels that have sufficient rays
double meanVoxelDensity(const std::string &cloud_file, ray::ElementDistribution distribution)
{
  const double vox_width = 0.5;
  ray::Cuboid bounds(Eigen::Vector3d(0.5, 0.5, 0.5), Eigen::Vector3d(3.5, 3.5, 2.5));
  const Eigen::Vector3i dims(6, 6, 4);
  ray::DensityGrid grid(bounds, vox_width, dims, distribution);
  grid.calculateDensities(cloud_file);
  double sum = 0.0;
  int count = 0;
  for (const auto &voxel : grid.voxels())
  {
    if (voxel.numRays() > 100)
    {
      sum += voxel.density();
      count++;
    }
  }
  EXPECT_GT(count, 20);
  return count > 0 ? sum / static_cast<double>(count) : 0.0;
}

/// Monte-Carlo verification of the vertical cylinder G function: the frontal area per volume N*d of a random field
/// of vertical cylinders should be recovered regardless of the ray zenith angle. The spherical (isotropic leaf)
/// model, by contrast, gives a result proportional to sin(zenith).
TEST(Basic, RayDensityVerticalCylinders)
{
  const double stems_per_m2 = 100.0;
  const double diameter = 0.01;
  const double frontal_area_density = stems_per_m2 * diameter;  // a = N d = 1 m^-1
  const double zeniths[] = { 15.0, 45.0, 80.0, 90.0 };
  unsigned seed = 1;
  for (const double zenith : zeniths)
  {
    ray::Cloud cloud = generateCylinderFieldCloud(zenith, stems_per_m2, diameter, 100000, seed++);
    const std::string file_name = "density_cylinders.ply";
    cloud.save(file_name);
    const double cylinder_estimate = meanVoxelDensity(file_name, ray::ElementDistribution::VerticalCylinder);
    const double spherical_estimate = meanVoxelDensity(file_name, ray::ElementDistribution::Spherical);
    const double sin_zenith = std::sin(zenith * ray::kPi / 180.0);
    std::cout << "zenith " << zenith << " deg: true N*d = " << frontal_area_density
              << ", vertical cylinder estimate = " << cylinder_estimate
              << ", spherical estimate = " << spherical_estimate << std::endl;
    EXPECT_NEAR(cylinder_estimate, frontal_area_density, 0.05 * frontal_area_density);
    // spherical model returns 2 * G * a with G = sin(zenith)
    EXPECT_NEAR(spherical_estimate, 2.0 * sin_zenith * frontal_area_density, 0.1 * frontal_area_density);
  }
}

/// Checks the Gaussian beam footprint model against a Monte-Carlo integration of the beam power over a strip
TEST(Basic, RayDensityFootprintWidth)
{
  std::mt19937 gen(5);
  std::normal_distribution<double> normal(0.0, 1.0);
  const double diameter = 0.01, beam_diameter = 0.025;
  const double sigma = beam_diameter / 4.0;  // 1/e^2 diameter is 4 standard deviations
  for (const double offset : { 0.0, 0.005, 0.012 })
  {
    int inside = 0;
    const int num = 200000;
    for (int i = 0; i < num; i++)
    {
      const double x = sigma * normal(gen);  // the strip is vertical, so only the horizontal coordinate matters
      inside += std::abs(x - offset) < 0.5 * diameter ? 1 : 0;
    }
    EXPECT_NEAR(ray::stripInterceptFraction(diameter, beam_diameter, offset), static_cast<double>(inside) / num, 0.005);
  }
  // the effective width is where the intercepted fraction equals the threshold
  const double threshold = 0.05;
  const double width = ray::effectiveStemWidth(diameter, beam_diameter, threshold);
  EXPECT_GT(width, diameter);
  EXPECT_NEAR(ray::stripInterceptFraction(diameter, beam_diameter, 0.5 * width), threshold, 1e-6);
  // a negligible beam gives the true diameter, and too thin a stem is never detected
  EXPECT_NEAR(ray::effectiveStemWidth(diameter, 1e-6, 0.5), diameter, 1e-6);
  EXPECT_EQ(ray::effectiveStemWidth(diameter, 0.1, 0.5), 0.0);
  // the tabulated model matches the direct calculation
  ray::BeamFootprint footprint;
  footprint.exit_diameter = 0.01;
  footprint.divergence = 0.003;
  footprint.threshold = threshold;
  ray::StemWidthModel model(diameter, footprint);
  EXPECT_NEAR(model.width(5.003), ray::effectiveStemWidth(diameter, footprint.diameter(5.003), threshold), 1e-5);
}

/// The free path response to non-overlapping stems is deterministic. For a negligible footprint it is at least
/// 1 / (1 - solid fraction), from the stems' area being excluded from the free path. It is a little more, as the
/// paths start at random free points, so sample the gaps between stems in proportion to their length, and these
/// gaps are more regular than the exponential gaps between overlapping (Poisson) stems. Without widening, the two
/// return range conventions are the same
TEST(Basic, RayDensityFreePathResponse)
{
  const double response = ray::freePathResponse(0.1, 1.0, 0.0, true);
  EXPECT_EQ(response, ray::freePathResponse(0.1, 1.0, 0.0, true));
  EXPECT_EQ(response, ray::freePathResponse(0.1, 1.0, 0.0, false));
  EXPECT_GT(response, 1.0 / (1.0 - 0.1));
  EXPECT_LT(response, 1.0 / (1.0 - 0.1) + 0.05);
  EXPECT_EQ(ray::freePathResponse(0.0, 2.0), 1.0);
  // with returns at the widened outline, widening the stems leaves less free space, so increases the response
  const double outline = ray::freePathResponse(0.1, 2.25, 0.0, true);
  EXPECT_GT(outline, response + 0.03);
  // with returns at the physical surface, paths can also start within the widened outlines (before the return), and
  // end further on. The response is similar
  const double surface = ray::freePathResponse(0.1, 2.25, 0.0, false);
  std::cout << "free path response at solid fraction 0.1, width ratio 2.25: " << outline << " (outline), " << surface
            << " (surface)" << std::endl;
  EXPECT_GT(surface, response + 0.03);
  EXPECT_NEAR(surface, outline, 0.03);
  // the precomputed table for no minimum gap matches the simulation (within its rounding across platforms)
  ray::StemResponseTable outline_table(0.0, true), surface_table(0.0, false);
  EXPECT_NEAR(outline_table.response(0.1, 2.25), outline, 2e-3);
  EXPECT_NEAR(surface_table.response(0.1, 2.25), surface, 2e-3);
  EXPECT_NEAR(surface_table.response(0.25, 1.5), ray::freePathResponse(0.25, 1.5, 0.0, false), 2e-3);
  // widened outlines that cover the plane leave no free space, and are saturated
  EXPECT_TRUE(std::isinf(ray::freePathResponse(0.4, 5.0, 0.0, true)));
  EXPECT_TRUE(std::isinf(outline_table.response(0.4, 5.0)));
  // a non-zero gap is computed on demand, so is the simulation at the grid nodes
  ray::StemResponseTable gap_table(0.5, false);
  EXPECT_NEAR(gap_table.response(0.04 / 2.25, 2.0), ray::freePathResponse(0.04 / 2.25, 2.0, 0.5, false), 1e-9);
}

/// The inversion of the stem response doesn't assume the response is at least 1: a sparse voxel with a response
/// below 1 has a true density above its measured density, and is not saturated
TEST(Basic, RayDensityResponseInversion)
{
  bool saturated = true;
  const double max_density = 50.0;
  // a response below 1, e.g. from simulation noise, or a model in which part of the stems' outline is free path
  auto low_response = [](double a) { return 0.9 + 0.001 * a; };
  const double measured = 3.0;
  double a = ray::invertDensityResponse(measured, max_density, low_response, saturated);
  EXPECT_FALSE(saturated);
  EXPECT_GT(a, measured);
  EXPECT_NEAR(a * low_response(a), measured, 1e-9);
  // a response that only reaches the measured density beyond twice it is still bracketed
  auto very_low_response = [](double) { return 0.3; };
  a = ray::invertDensityResponse(measured, max_density, very_low_response, saturated);
  EXPECT_FALSE(saturated);
  EXPECT_NEAR(a, measured / 0.3, 1e-9);
  // beyond the response at the maximum density, it is saturated
  a = ray::invertDensityResponse(100.0, max_density, low_response, saturated);
  EXPECT_TRUE(saturated);
  EXPECT_EQ(a, max_density);

  // the stem response table inversion, for 1 cm stems
  const double diameter = 0.01, width_ratio = 3.0;
  ray::StemResponseTable table(0.0, false);
  a = ray::invertStemResponse(measured, width_ratio, diameter, 0.0, 0.4, table, saturated);
  EXPECT_FALSE(saturated);
  EXPECT_LT(a, measured);
  EXPECT_NEAR(a * table.response(0.25 * ray::kPi * diameter * a, width_ratio), measured, 1e-6);
  ray::invertStemResponse(1000.0, width_ratio, diameter, 0.0, 0.4, table, saturated);
  EXPECT_TRUE(saturated);
}

/// A missing input cloud is an error, rather than an empty density grid
TEST(Basic, RayDensityMissingFile)
{
  const ray::Cuboid bounds(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1, 1, 1));
  ray::DensityOptions options;
  int shift = 0;
  EXPECT_EQ(ray::estimateDensityGrid("missing_cloud.ply", bounds, 0.1, Eigen::Vector3i(10, 10, 10), options, shift),
            nullptr);
  EXPECT_FALSE(ray::renderCloud("missing_cloud.ply", bounds, ray::ViewDirection::Top, ray::RenderStyle::Density, 0.1,
                                "missing_cloud.png", "", false, nullptr, options));
  ray::DensityTileParams params;
  params.tile_width = 1.0;
  EXPECT_FALSE(ray::estimateDensityTile("missing_cloud.ply", params, "missing_cloud.zarr"));
}

/// The generator reports returns from stems at the range of the physical stem surface (or the closest approach of
/// the beam axis, when it misses the stem) by default, or at the widened outline with range_to_outline
TEST(Basic, RayDensityReturnRange)
{
  ray::PneumatophoreParams params;
  params.num_rays = 20000;
  params.stems_per_m2 = 30.0;
  params.plot_width = 2.0;
  params.margin = 1.5;
  params.footprint.exit_diameter = 0.01;
  params.footprint.divergence = 0.003;
  params.footprint.threshold = 0.05;
  const ray::StemWidthModel width_model(params.diameter, params.footprint);
  for (const bool outline : { false, true })
  {
    params.footprint.range_to_outline = outline;
    ray::PneumatophoreGen gen;
    gen.generate(params, 11);
    const ray::Cloud &cloud = gen.cloud();
    int side_hits = 0, on_surface = 0, within = 0, on_outline = 0;
    for (size_t i = 0; i < cloud.rayCount(); i++)
    {
      if (cloud.colours[i].green != 160 || cloud.colours[i].alpha == 0)
      {
        continue;  // not a stem return
      }
      // the horizontal distance from the return to the nearest stem's axis
      const Eigen::Vector3d &end = cloud.ends[i];
      double min_dist = std::numeric_limits<double>::max();
      Eigen::Vector3d nearest(0, 0, 0);
      for (const auto &stem : gen.stems())
      {
        const double dist = (stem.head<2>() - end.head<2>()).norm();
        if (dist < min_dist)
        {
          min_dist = dist;
          nearest = stem;
        }
      }
      if (end[2] > nearest[2] - 1e-6)
      {
        continue;  // a return from the top of the stem
      }
      side_hits++;
      const double radius = 0.5 * params.diameter;
      // the effective width is at the range of the beam axis' closest approach to the stem
      const Eigen::Vector3d dir = (end - cloud.starts[i]).normalized();
      const Eigen::Vector2d d2 = dir.head<2>();
      const double along = (nearest.head<2>() - cloud.starts[i].head<2>()).dot(d2) / d2.squaredNorm();
      const double outline_radius = 0.5 * width_model.width(along);
      on_surface += std::abs(min_dist - radius) < 1e-6 ? 1 : 0;
      within += min_dist > radius - 1e-6 && min_dist < outline_radius + 1e-6 ? 1 : 0;
      on_outline += std::abs(min_dist - outline_radius) < 1e-6 ? 1 : 0;
    }
    std::cout << (outline ? "outline" : "surface") << " convention: " << side_hits << " side returns, " << on_surface
              << " on the stem surface, " << on_outline << " on the widened outline" << std::endl;
    ASSERT_GT(side_hits, 200);
    if (outline)
    {
      EXPECT_GT(on_outline, 0.95 * side_hits);
    }
    else
    {
      // the beam axis hits the physical stem for about d / d_eff of the returns, and otherwise the return is at the
      // closest approach, between the surface and the outline
      EXPECT_GE(within, 0.99 * side_hits);
      EXPECT_GT(on_surface, 0.25 * side_hits);
      EXPECT_LT(on_outline, 0.05 * side_hits);
    }
  }
}

/// Generated stems don't overlap, and have the exact density requested, up to the packing limit
TEST(Basic, RayDensityStemPlacement)
{
  ray::PneumatophoreParams params;
  params.num_rays = 1;
  params.stems_per_m2 = 2000.0;  // solid fraction 0.157 for 1 cm stems
  params.plot_width = 2.0;
  params.margin = 0.5;
  ray::PneumatophoreGen gen;
  gen.generate(params, 3);
  const double field_width = params.plot_width + 2.0 * params.margin;
  EXPECT_EQ(static_cast<double>(gen.stems().size()), std::round(params.stems_per_m2 * field_width * field_width));
  double min_separation2 = std::numeric_limits<double>::max();
  for (size_t i = 0; i < gen.stems().size(); i++)
  {
    for (size_t j = i + 1; j < gen.stems().size(); j++)
    {
      min_separation2 = std::min(min_separation2, (gen.stems()[i] - gen.stems()[j]).head<2>().squaredNorm());
    }
  }
  EXPECT_GE(std::sqrt(min_separation2), params.diameter);
  params.stems_per_m2 = 7000.0;  // solid fraction 0.55, beyond the packing limit
  EXPECT_THROW(gen.generate(params, 3), std::runtime_error);

  // with a minimum gap of half a diameter, stem centres are at least 1.5 diameters apart
  params.stems_per_m2 = 1000.0;
  params.min_gap_ratio = 0.5;
  gen.generate(params, 4);
  EXPECT_EQ(static_cast<double>(gen.stems().size()), std::round(params.stems_per_m2 * field_width * field_width));
  min_separation2 = std::numeric_limits<double>::max();
  for (size_t i = 0; i < gen.stems().size(); i++)
  {
    for (size_t j = i + 1; j < gen.stems().size(); j++)
    {
      min_separation2 = std::min(min_separation2, (gen.stems()[i] - gen.stems()[j]).head<2>().squaredNorm());
    }
  }
  EXPECT_GE(std::sqrt(min_separation2), 1.5 * params.diameter);
  params.stems_per_m2 = 3000.0;  // 0.24 solid fraction, but 0.53 for the exclusion discs of 1.5 diameters
  EXPECT_THROW(gen.generate(params, 3), std::runtime_error);

  // trunks: stems don't overlap them, and rays only hit trunks ahead of the scanner (which is above the stems)
  params.stems_per_m2 = 300.0;
  params.min_gap_ratio = 0.0;
  params.trunks_per_m2 = 0.5;
  params.num_rays = 20000;
  gen.generate(params, 5);
  ASSERT_FALSE(gen.trunks().empty());
  for (const auto &stem : gen.stems())
  {
    for (const auto &trunk : gen.trunks())
    {
      EXPECT_GE((stem.head<2>() - trunk.head<2>()).norm(), 0.5 * (trunk[2] + params.diameter));
    }
  }
  int trunk_hits = 0;
  const ray::Cloud &cloud = gen.cloud();
  for (size_t i = 0; i < cloud.rayCount(); i++)
  {
    if (cloud.colours[i].red == 150 && cloud.colours[i].alpha > 0)  // trunk returns
    {
      trunk_hits++;
      EXPECT_GT((cloud.ends[i] - cloud.starts[i]).norm(), 0.15);
    }
  }
  EXPECT_GT(trunk_hits, 0);
}

/// Reads a chunk of float32 values from a Zarr store
std::vector<float> readZarrChunk(const std::string &file_name)
{
  std::ifstream ifs(file_name, std::ios::binary);
  std::vector<float> data;
  float value;
  while (ifs.read(reinterpret_cast<char *>(&value), sizeof(float)))
  {
    data.push_back(value);
  }
  return data;
}

/// A tile without rays is valid, and is written as unobserved
TEST(Basic, RayDensityEmptyTile)
{
  ray::Cloud cloud;
  cloud.addRay(Eigen::Vector3d(10, 10, 1), Eigen::Vector3d(10, 11, 0), 0.0, ray::RGBA(255, 255, 255, 255));
  cloud.save("far_tile.ply");  // a ray far from the tile
  ray::Cloud empty;
  empty.save("empty_tile.ply");
  for (const std::string name : { "far_tile", "empty_tile" })
  {
    ray::DensityTileParams params;
    params.tile_width = 1.0;
    params.options.distribution = ray::ElementDistribution::VerticalCylinder;
    params.options.path_estimator = ray::PathEstimator::PooledBeerLambert;
    EXPECT_TRUE(ray::estimateDensityTile(name + ".ply", params, name + ".zarr"));
    const std::vector<float> density = readZarrChunk(name + ".zarr/density/0.0.0");
    const std::vector<float> observed = readZarrChunk(name + ".zarr/observed_fraction/0.0.0");
    ASSERT_EQ(density.size(), 6u);  // 1 x 1 m cells of 0.1 m from 0 to 0.6 m
    ASSERT_EQ(observed.size(), 6u);
    for (size_t i = 0; i < density.size(); i++)
    {
      EXPECT_TRUE(std::isnan(density[i]));
      EXPECT_EQ(observed[i], 0.0f);
    }
  }
}

/// Density estimated over 2 x 2 tiles (split by raysplit tiles, keeping ray starts) matches a single run over the
/// whole domain, as the halo gives each tile's voxels exactly the rays that a single run would. This holds for a split
/// halo equal to raydensity's halo (the pool radius plus two voxels, when the split's clipped ray ends lie on voxel
/// faces of the grid), and for the larger split halo of scripts/raydensity_tiles.sh (1.5 voxels more than the halo
/// rounded up to whole voxels)
TEST(Basic, RayDensityTiles)
{
  ray::PneumatophoreParams scene;
  scene.plot_width = 2.0;
  scene.margin = 3.0;
  scene.num_rays = 300000;
  scene.patchy = true;
  scene.patch_radius = 0.3;
  scene.patch_cover = 0.5;
  scene.footprint.exit_diameter = 0.01;
  scene.footprint.divergence = 0.003;
  scene.footprint.threshold = 0.05;
  ray::PneumatophoreGen gen;
  gen.generate(scene, 9);
  gen.cloud().save("tiles_test.ply");

  ray::DensityTileParams params;
  params.origin = Eigen::Vector2d(-1.0, -1.0);
  params.cell_size = Eigen::Vector3d(0.5, 0.5, 0.1);
  params.max_height = 0.4;
  params.options.distribution = ray::ElementDistribution::VerticalCylinder;
  params.options.path_estimator = ray::PathEstimator::PooledBeerLambert;
  params.options.flat_top_compensation = false;
  params.options.stem_diameter = scene.diameter;
  params.options.footprint = scene.footprint;
  const double w = params.voxel_width;
  const double halo = params.options.pool_radius + 2.0 * w;  // raydensity's (default) halo
  // the split halo of scripts/raydensity_tiles.sh
  const double script_halo = (std::ceil(halo / w - 1e-9) + 1.5) * w;

  // a single tile over the whole domain
  params.tile_width = 2.0;
  EXPECT_TRUE(ray::estimateDensityTile("tiles_test.ply", params, "tiles_single.zarr"));
  const std::vector<float> single = readZarrChunk("tiles_single.zarr/density/0.0.0");
  const std::vector<float> single_hits = readZarrChunk("tiles_single.zarr/hits/0.0.0");
  // the single chunk is (z, y, x) = (4, 4, 4) cells, each split chunk is (4, 2, 2)
  ASSERT_EQ(single.size(), 64u);
  ASSERT_EQ(single_hits.size(), 64u);

  params.tile_width = 1.0;
  params.num_tiles = Eigen::Vector2i(2, 2);
  for (const double split_halo : { halo, script_halo })
  {
    std::stringstream stub;
    stub << "tiles_test_halo" << static_cast<int>(std::round(1000.0 * split_halo));
    const std::string store = stub.str() + ".zarr";
    EXPECT_TRUE(
      ray::splitTiles("tiles_test.ply", stub.str(), params.origin, 1.0, split_halo, -1.0, 2.0, Eigen::Vector2i(2, 2)));
    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < 2; j++)
      {
        params.tile = Eigen::Vector2i(i, j);
        std::stringstream name;
        name << stub.str() << "_" << i << "_" << j << ".ply";
        EXPECT_TRUE(ray::estimateDensityTile(name.str(), params, store));
      }
    }
    int compared = 0;
    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < 2; j++)
      {
        std::stringstream chunk;
        chunk << "/0." << j << "." << i;
        const std::vector<float> tile = readZarrChunk(store + "/density" + chunk.str());
        const std::vector<float> hits = readZarrChunk(store + "/hits" + chunk.str());
        ASSERT_EQ(tile.size(), 16u);
        ASSERT_EQ(hits.size(), 16u);
        for (int z = 0; z < 4; z++)
        {
          for (int y = 0; y < 2; y++)
          {
            for (int x = 0; x < 2; x++)
            {
              const size_t t = static_cast<size_t>(x + 2 * (y + 2 * z));
              const size_t s = static_cast<size_t>((x + 2 * i) + 4 * ((y + 2 * j) + 4 * z));
              const float a = tile[t], b = single[s];
              EXPECT_EQ(std::isnan(a), std::isnan(b)) << "split halo " << split_halo;
              if (!std::isnan(a) && !std::isnan(b))
              {
                EXPECT_NEAR(a, b, 1e-4 * std::max(1.0f, std::abs(b))) << "split halo " << split_halo;
                EXPECT_NEAR(hits[t], single_hits[s], 1e-4 * std::max(1.0f, std::abs(single_hits[s])))
                  << "split halo " << split_halo;
                compared++;
              }
            }
          }
        }
      }
    }
    EXPECT_GT(compared, 40) << "split halo " << split_halo;
  }
}

/// The ratio of the total hits to the total G-weighted path length over all voxels of a grid
double pooledGridDensity(const ray::DensityGrid &grid)
{
  double hits = 0.0, path = 0.0;
  for (const auto &voxel : grid.voxels())
  {
    hits += voxel.numHits();
    path += voxel.pathLength();
  }
  return path > 0.0 ? hits / path : 0.0;
}

/// The crossing estimator counts each hit ray's path to where it would have exited its voxel (plus the stem interior
/// beyond the return), so in voxels much shorter than the mean free path it recovers N d of a field of vertical
/// cylinders. Its bias is about -a G L / 2 for a chord L, from rays occluded within the voxel, whereas Beer-Lambert
/// (path to the hit) is unbiased for randomly placed stems
TEST(Basic, RayDensityCrossing)
{
  // many thin stems, so that the field's density over the scanned region varies little between realisations
  const double stems_per_m2 = 300.0, diameter = 0.001;
  const double frontal_area_density = stems_per_m2 * diameter;  // a = N d = 0.3 m^-1
  const ray::Cuboid bounds(Eigen::Vector3d(0.5, 0.5, 0.5), Eigen::Vector3d(3.5, 3.5, 2.5));
  const double vox_width = 0.1;
  const Eigen::Vector3i dims(30, 30, 20);
  unsigned seed = 21;
  for (const double zenith : { 45.0, 90.0 })
  {
    ray::Cloud cloud = generateCylinderFieldCloud(zenith, stems_per_m2, diameter, 100000, seed++);
    cloud.save("density_crossing.ply");
    double estimates[2];
    for (int k = 0; k < 2; k++)
    {
      ray::DensityGrid grid(bounds, vox_width, dims, ray::ElementDistribution::VerticalCylinder);
      grid.setPathEstimator(k == 0 ? ray::PathEstimator::Crossing : ray::PathEstimator::BeerLambert);
      grid.setStemDiameter(diameter);
      ASSERT_TRUE(grid.calculateDensities("density_crossing.ply"));
      estimates[k] = pooledGridDensity(grid);
    }
    std::cout << "zenith " << zenith << " deg: true N*d = " << frontal_area_density
              << ", crossing estimate = " << estimates[0] << ", beer-lambert estimate = " << estimates[1] << std::endl;
    // within 6%, which is mostly the variation of the (Poisson) field's density over the scanned region (about 2%)
    EXPECT_NEAR(estimates[0], frontal_area_density, 0.06 * frontal_area_density);
    EXPECT_NEAR(estimates[1], frontal_area_density, 0.06 * frontal_area_density);
    // on the same rays, hit rays count a longer path, by about a G L / 2 for 0.1 m voxels (0.5% to 1.5%)
    EXPECT_LT(estimates[0], estimates[1]);
    EXPECT_GT(estimates[0], 0.97 * estimates[1]);
  }
}

/// Reads a whole file as text, or an empty string if it doesn't exist
std::string readTextFile(const std::string &file_name)
{
  std::ifstream ifs(file_name, std::ios::binary);
  std::stringstream ss;
  ss << ifs.rdbuf();
  return ss.str();
}

/// Reads a chunk of float64 values from a Zarr store
std::vector<double> readZarrChunkDouble(const std::string &file_name)
{
  std::ifstream ifs(file_name, std::ios::binary);
  std::vector<double> data;
  double value;
  while (ifs.read(reinterpret_cast<char *>(&value), sizeof(double)))
  {
    data.push_back(value);
  }
  return data;
}

/// Writes a planar ground mesh z = slope . (x, y), as a grid of squares of width @c spacing over [min, max]^2
void writePlaneMesh(const std::string &file_name, const Eigen::Vector2d &slope, double min, double max, double spacing)
{
  ray::Mesh mesh;
  const int n = static_cast<int>(std::round((max - min) / spacing));
  for (int j = 0; j <= n; j++)
  {
    for (int i = 0; i <= n; i++)
    {
      const double x = min + i * spacing, y = min + j * spacing;
      mesh.vertices().push_back(Eigen::Vector3d(x, y, slope[0] * x + slope[1] * y));
    }
  }
  for (int j = 0; j < n; j++)
  {
    for (int i = 0; i < n; i++)
    {
      const int v = i + (n + 1) * j;
      mesh.indexList().push_back(Eigen::Vector3i(v, v + 1, v + n + 2));
      mesh.indexList().push_back(Eigen::Vector3i(v, v + n + 2, v + n + 1));
    }
  }
  EXPECT_TRUE(ray::writePlyMesh(file_name, mesh));
}

/// Layers are heights above the ground mesh: a sloped scene (the flat scene sheared vertically by the ground plane,
/// which keeps each ray's horizontal path, so its hits and its G-weighted path length) estimated with the ground mesh
/// matches the flat scene, up to the registration of the voxels to the ground. The same holds with a ground mesh much
/// larger than the tile, which is cropped to the tile
TEST(Basic, RayDensityGroundMesh)
{
  ray::PneumatophoreParams scene;
  scene.plot_width = 2.0;
  scene.margin = 1.5;
  scene.num_rays = 200000;
  ray::PneumatophoreGen gen;
  gen.generate(scene, 13);
  gen.cloud().save("ground_flat.ply");
  const Eigen::Vector2d slope(0.1, -0.05);
  ray::Cloud sloped = gen.cloud();
  for (size_t i = 0; i < sloped.rayCount(); i++)
  {
    sloped.starts[i][2] += slope.dot(sloped.starts[i].head<2>());
    sloped.ends[i][2] += slope.dot(sloped.ends[i].head<2>());
  }
  sloped.save("ground_sloped.ply");
  writePlaneMesh("ground_small_mesh.ply", slope, -5.0, 5.0, 10.0);
  writePlaneMesh("ground_large_mesh.ply", slope, -200.0, 200.0, 2.0);  // 80000 triangles, mostly far from the tile

  ray::DensityTileParams params;
  params.origin = Eigen::Vector2d(-1.0, -1.0);
  params.tile_width = 2.0;
  params.max_height = 0.3;
  params.options.distribution = ray::ElementDistribution::VerticalCylinder;
  params.options.path_estimator = ray::PathEstimator::PooledBeerLambert;
  params.options.flat_top_compensation = false;
  params.options.stem_diameter = scene.diameter;
  EXPECT_TRUE(ray::estimateDensityTile("ground_flat.ply", params, "ground_flat.zarr"));
  const std::vector<float> flat = readZarrChunk("ground_flat.zarr/density/0.0.0");
  ASSERT_EQ(flat.size(), 12u);  // 2 x 2 cells of 1 m, 3 layers
  // without the ground mesh, the sloped scene's layers are absolute heights, so differ
  EXPECT_TRUE(ray::estimateDensityTile("ground_sloped.ply", params, "ground_absolute.zarr"));
  const std::vector<float> absolute = readZarrChunk("ground_absolute.zarr/density/0.0.0");
  ASSERT_EQ(absolute.size(), 12u);
  for (const std::string mesh : { "ground_small_mesh", "ground_large_mesh" })
  {
    params.ground_mesh_file = mesh + ".ply";
    EXPECT_TRUE(ray::estimateDensityTile("ground_sloped.ply", params, mesh + ".zarr"));
    const std::vector<float> relative = readZarrChunk(mesh + ".zarr/density/0.0.0");
    ASSERT_EQ(relative.size(), 12u);
    double max_error = 0.0, max_absolute_error = 0.0;
    for (size_t c = 0; c < flat.size(); c++)
    {
      ASSERT_FALSE(std::isnan(flat[c]));
      ASSERT_FALSE(std::isnan(relative[c]));
      max_error = std::max(max_error, std::abs(relative[c] / flat[c] - 1.0));
      max_absolute_error = std::max(max_absolute_error, std::abs(absolute[c] / flat[c] - 1.0));
      EXPECT_NEAR(relative[c], flat[c], 0.05 * flat[c]) << mesh << " cell " << c;
    }
    std::cout << mesh << ": largest relative difference from the flat scene " << max_error
              << " with heights above ground, " << max_absolute_error << " with absolute heights" << std::endl;
    EXPECT_GT(max_absolute_error, 3.0 * max_error);
  }
  EXPECT_NE(readTextFile("ground_large_mesh.zarr/.zattrs").find("\"heights_above_ground\": true"), std::string::npos);
}

/// The occlusion correction for a minimum gap between stems. A gap makes the stems more regular, leaving less free
/// space between the widened stems, so increases the free path response. A scene generated with a minimum gap is
/// estimated more accurately when the correction assumes the same gap
TEST(Basic, RayDensityMinGap)
{
  for (const double rho : { 1.0, 2.0 })
  {
    const double no_gap = ray::freePathResponse(0.1, rho, 0.0, false);
    const double gap = ray::freePathResponse(0.1, rho, 1.0, false);
    std::cout << "free path response at solid fraction 0.1, width ratio " << rho << ": " << no_gap << " (no gap), "
              << gap << " (gap of one diameter)" << std::endl;
    EXPECT_GT(gap, no_gap + 0.01);
  }

  ray::PneumatophoreParams scene;
  scene.plot_width = 1.0;
  scene.margin = 1.5;
  scene.num_rays = 300000;
  scene.stems_per_m2 = 1000.0;
  scene.min_gap_ratio = 1.0;
  scene.footprint.exit_diameter = 0.01;
  scene.footprint.divergence = 0.003;
  scene.footprint.threshold = 0.05;
  ray::PneumatophoreGen gen;
  gen.generate(scene, 17);
  gen.cloud().save("min_gap.ply");
  const double half_plot = 0.5 * scene.plot_width;
  double frontal_area = 0.0;
  for (const auto &stem : gen.stems())
  {
    if (std::abs(stem[0]) < half_plot && std::abs(stem[1]) < half_plot)
    {
      frontal_area += scene.diameter * std::min(stem[2], 0.1);
    }
  }
  const double truth = frontal_area / (scene.plot_width * scene.plot_width * 0.1);

  ray::DensityTileParams params;
  params.origin = Eigen::Vector2d(-half_plot, -half_plot);
  params.tile_width = scene.plot_width;
  params.max_height = 0.1;
  params.options.distribution = ray::ElementDistribution::VerticalCylinder;
  params.options.path_estimator = ray::PathEstimator::PooledBeerLambert;
  params.options.flat_top_compensation = false;
  params.options.stem_diameter = scene.diameter;
  params.options.footprint = scene.footprint;
  double estimates[2];
  for (int k = 0; k < 2; k++)
  {
    params.options.min_gap_ratio = k == 0 ? 0.0 : scene.min_gap_ratio;
    const std::string store = k == 0 ? "min_gap_0.zarr" : "min_gap_1.zarr";
    EXPECT_TRUE(ray::estimateDensityTile("min_gap.ply", params, store));
    const std::vector<float> density = readZarrChunk(store + "/density/0.0.0");
    ASSERT_EQ(density.size(), 1u);
    estimates[k] = density[0];
  }
  std::cout << "true N*d near ground at a minimum gap of one diameter: " << truth
            << ", estimated assuming no gap: " << estimates[0] << ", assuming the gap: " << estimates[1] << std::endl;
  EXPECT_LT(std::abs(estimates[1] - truth), std::abs(estimates[0] - truth));
  EXPECT_LT(estimates[1], estimates[0]);
}

/// The Zarr store of a tile at georeferenced (e.g. UTM) coordinates: valid metadata, float32 data and float64
/// coordinates that resolve the cell centres exactly, NaN for unobserved cells (here below the ground), the tile's
/// description, and a second tile with different parameters is refused
TEST(Basic, RayDensityZarr)
{
  ray::PneumatophoreParams scene;
  scene.plot_width = 2.0;
  scene.margin = 1.0;
  scene.num_rays = 50000;
  ray::PneumatophoreGen gen;
  gen.generate(scene, 19);
  const Eigen::Vector3d offset(400000.0, 6100000.0, 0.0);
  ray::Cloud cloud = gen.cloud();
  for (size_t i = 0; i < cloud.rayCount(); i++)
  {
    cloud.starts[i] += offset;
    cloud.ends[i] += offset;
  }
  cloud.save("zarr_utm.ply");

  ray::DensityTileParams params;
  params.origin = Eigen::Vector2d(offset[0] - 1.0, offset[1] - 1.0);
  params.tile_width = 1.0;
  params.num_tiles = Eigen::Vector2i(2, 2);
  params.tile = Eigen::Vector2i(1, 0);
  params.cell_size = Eigen::Vector3d(0.5, 0.5, 0.1);
  // three layers below the ground, which no ray reaches. The voxels just below the ground (within one voxel) take
  // the neighbour prior of the voxels above it (DensityGrid::addNeighbourPriors), so are observed, and the layer just
  // below the ground is partly observed
  params.min_height = -0.3;
  params.max_height = 0.3;
  params.options.distribution = ray::ElementDistribution::VerticalCylinder;
  params.options.path_estimator = ray::PathEstimator::PooledBeerLambert;
  params.options.flat_top_compensation = false;
  params.options.stem_diameter = scene.diameter;
  params.description = "test tile";
  const std::string store = "zarr_utm.zarr";
  ASSERT_TRUE(ray::estimateDensityTile("zarr_utm.ply", params, store));

  EXPECT_EQ(readTextFile(store + "/.zgroup"), "{\n  \"zarr_format\": 2\n}\n");
  const std::string attrs = readTextFile(store + "/.zattrs");
  ASSERT_GT(attrs.size(), 2u);
  EXPECT_EQ(attrs.front(), '{');
  EXPECT_EQ(attrs.substr(attrs.size() - 2), "}\n");
  EXPECT_NE(attrs.find("\"origin\": [399999, 6099999]"), std::string::npos) << attrs;
  EXPECT_NE(attrs.find("\"density_model\": \"vertical_cylinder\""), std::string::npos);
  EXPECT_NE(attrs.find("\"estimator\": \"pooled\""), std::string::npos);
  for (const std::string name : { "density", "mean_voxel_density", "patchiness_ratio", "hits", "weighted_path",
                                  "observed_fraction", "solid_fraction" })
  {
    const std::string zarray = readTextFile(store + "/" + name + "/.zarray");
    EXPECT_NE(zarray.find("\"dtype\": \"<f4\""), std::string::npos) << name;
    EXPECT_NE(zarray.find("\"shape\": [6, 4, 4]"), std::string::npos) << name << ": " << zarray;
    EXPECT_NE(zarray.find("\"chunks\": [6, 2, 2]"), std::string::npos) << name;
    EXPECT_NE(zarray.find("\"fill_value\": \"NaN\""), std::string::npos) << name;
    EXPECT_NE(zarray.find("\"order\": \"C\""), std::string::npos) << name;
    EXPECT_NE(readTextFile(store + "/" + name + "/.zattrs").find("\"_ARRAY_DIMENSIONS\": [\"z\", \"y\", \"x\"]"),
              std::string::npos)
      << name;
  }
  // coordinates: float64 cell centres, exact at these values
  for (const std::string name : { "x", "y", "z" })
  {
    EXPECT_NE(readTextFile(store + "/" + name + "/.zarray").find("\"dtype\": \"<f8\""), std::string::npos) << name;
  }
  const std::vector<double> xs = readZarrChunkDouble(store + "/x/0");
  const std::vector<double> ys = readZarrChunkDouble(store + "/y/0");
  const std::vector<double> zs = readZarrChunkDouble(store + "/z/0");
  ASSERT_EQ(xs.size(), 4u);
  ASSERT_EQ(ys.size(), 4u);
  ASSERT_EQ(zs.size(), 6u);
  for (int i = 0; i < 4; i++)
  {
    EXPECT_EQ(xs[static_cast<size_t>(i)], 399999.25 + 0.5 * i);
    EXPECT_EQ(ys[static_cast<size_t>(i)], 6099999.25 + 0.5 * i);
  }
  for (int i = 0; i < 6; i++)
  {
    EXPECT_NEAR(zs[static_cast<size_t>(i)], -0.25 + 0.1 * i, 1e-12);
  }
  // this tile's chunk (z, y, x) = (6, 2, 2). The two lowest layers are unobserved: NaN, with observed fraction 0
  const std::vector<float> density = readZarrChunk(store + "/density/0.0.1");
  const std::vector<float> mean = readZarrChunk(store + "/mean_voxel_density/0.0.1");
  const std::vector<float> observed = readZarrChunk(store + "/observed_fraction/0.0.1");
  ASSERT_EQ(density.size(), 24u);
  ASSERT_EQ(mean.size(), 24u);
  ASSERT_EQ(observed.size(), 24u);
  for (size_t c = 0; c < density.size(); c++)
  {
    if (c >= 8 && c < 12)
    {
      continue;  // the layer just below the ground
    }
    const bool below_ground = c < 8;
    EXPECT_EQ(std::isnan(density[c]), below_ground) << c;
    EXPECT_EQ(std::isnan(mean[c]), below_ground) << c;
    EXPECT_EQ(observed[c] == 0.0f, below_ground) << c;
    if (!below_ground)
    {
      EXPECT_GT(density[c], 0.0f);
    }
  }
  // only this tile's chunk is written
  EXPECT_FALSE(std::ifstream(store + "/density/0.0.0").good());
  const std::string tile_info = readTextFile(store + "/tiles/1_0.json");
  EXPECT_NE(tile_info.find("\"tile\": [1, 0]"), std::string::npos) << tile_info;
  EXPECT_NE(tile_info.find("\"input\": \"zarr_utm.ply\""), std::string::npos);
  EXPECT_NE(tile_info.find("\"description\": \"test tile\""), std::string::npos);

  // another tile of the same domain and parameters is accepted, but not one with different parameters
  params.tile = Eigen::Vector2i(0, 0);
  EXPECT_TRUE(ray::estimateDensityTile("zarr_utm.ply", params, store));
  EXPECT_TRUE(std::ifstream(store + "/tiles/0_0.json").good());
  params.tile = Eigen::Vector2i(0, 1);
  params.options.pool_radius = 0.1;
  EXPECT_FALSE(ray::estimateDensityTile("zarr_utm.ply", params, store));
  EXPECT_FALSE(std::ifstream(store + "/density/0.1.0").good());
  EXPECT_FALSE(std::ifstream(store + "/tiles/0_1.json").good());
  EXPECT_EQ(readTextFile(store + "/.zattrs"), attrs);
}

/// Splits a line of comma-separated values
std::vector<std::string> splitCsvLine(const std::string &line)
{
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, ','))
  {
    fields.push_back(field);
  }
  return fields;
}

/// The density command line tools: rayrender density cells and profile, the option checks of rayrender, raydensity
/// and raycreate, and raysplit tiles within a domain
TEST(Basic, RayDensityTools)
{
  EXPECT_EQ(command("raycreate pneumatophores 1 --rays 100000 --stems_per_m2 50"), 0);
  // cells from below the ground (unobserved, so nan, except within a voxel of the ground, from the neighbour prior)
  // to above the stems, with the crossing estimator
  EXPECT_EQ(
    command("rayrender pneumatophores.ply top density --density_model vertical_cylinder --estimator crossing "
            "--stem_diameter 0.01 --no_flat_top --pixel_width 0.05 --min_bound -1,-1,-0.3 --max_bound 1,1,0.6 "
            "--cells tools_cells.csv --cell_size 1,1,0.1 --profile tools_profile.csv --output tools_density.png"),
    0);
  {
    std::ifstream ifs("tools_cells.csv");
    std::string line;
    ASSERT_TRUE(std::getline(ifs, line));
    EXPECT_EQ(line, "x_min,y_min,z_min,x_max,y_max,z_max,observed_fraction,hits,weighted_path,"
                    "density_frontal_area_per_volume_Nd,mean_voxel_density_frontal_area_per_volume_Nd,"
                    "patchiness_ratio");
    int rows = 0, unobserved = 0, observed = 0;
    while (std::getline(ifs, line))
    {
      const std::vector<std::string> fields = splitCsvLine(line);
      ASSERT_EQ(fields.size(), 12u) << line;
      rows++;
      const double z_min = std::stod(fields[2]);
      if (z_min < -0.2 + 1e-6)  // well below the ground
      {
        unobserved++;
        EXPECT_EQ(std::stod(fields[6]), 0.0) << line;
        for (size_t i = 7; i < 12; i++)
        {
          EXPECT_EQ(fields[i], "nan") << line;
        }
      }
      else if (z_min > -1e-6 && z_min < 0.15)  // near-ground layers, within the stems
      {
        observed++;
        EXPECT_GT(std::stod(fields[9]), 0.0) << line;
        EXPECT_GT(std::stod(fields[10]), 0.0) << line;
      }
    }
    EXPECT_EQ(rows, 2 * 2 * 9);
    EXPECT_EQ(unobserved, 2 * 2 * 2);
    EXPECT_EQ(observed, 2 * 2 * 2);
  }
  {
    std::ifstream ifs("tools_profile.csv");
    std::string line;
    ASSERT_TRUE(std::getline(ifs, line));
    EXPECT_EQ(line, "z_min,z_max,observed_voxels,total_voxels,mean_hits_per_voxel,density_frontal_area_per_volume_Nd,"
                    "mean_voxel_density_frontal_area_per_volume_Nd,patchiness_ratio");
    int rows = 0, unobserved = 0;
    while (std::getline(ifs, line))
    {
      const std::vector<std::string> fields = splitCsvLine(line);
      ASSERT_EQ(fields.size(), 8u) << line;
      if (std::stod(fields[1]) < -0.05 + 1e-6)  // more than a voxel below the ground
      {
        EXPECT_EQ(fields[2], "0") << line;
        EXPECT_EQ(fields[5], "nan") << line;
        EXPECT_EQ(fields[6], "nan") << line;
        unobserved++;
      }
      rows++;
    }
    EXPECT_EQ(rows, 18);  // 0.05 m layers over 0.9 m
    EXPECT_EQ(unobserved, 5);
  }

  // density options require a density style, and raydensity requires the tile
  EXPECT_NE(command("rayrender pneumatophores.ply top ends --cells tools_bad.csv"), 0);
  EXPECT_NE(command("rayrender pneumatophores.ply top ends --estimator crossing"), 0);
  EXPECT_NE(command("rayrender pneumatophores.ply top ends --profile tools_bad.csv"), 0);
  EXPECT_NE(command("raydensity pneumatophores.ply tools_bad.zarr --origin 0,0 --tile_width 1 --tiles 1,1"), 0);
  EXPECT_FALSE(std::ifstream("tools_bad.csv").good());
  EXPECT_FALSE(std::ifstream("tools_bad.zarr/.zgroup").good());
  // pneumatophore options only apply to pneumatophores
  EXPECT_NE(command("raycreate room 1 --stems_per_m2 100"), 0);
  EXPECT_NE(command("raycreate forest 1 --rays 1000"), 0);

  // raysplit tiles within a domain of 2 x 2 tiles writes only those tiles, though the scan extends beyond them
  copy("pneumatophores.ply tools_split.ply");
  EXPECT_EQ(command("raysplit tools_split.ply tiles -1,-1 1 0.225 -1,2 2,2"), 0);
  std::set<std::string> tiles;
  for (const auto &entry : std::filesystem::directory_iterator("."))
  {
    const std::string name = entry.path().filename().string();
    if (name.rfind("tools_split_", 0) == 0)
    {
      tiles.insert(name);
    }
  }
  EXPECT_EQ(tiles, std::set<std::string>(
                     { "tools_split_0_0.ply", "tools_split_0_1.ply", "tools_split_1_0.ply", "tools_split_1_1.ply" }));
  EXPECT_EQ(command("raydensity tools_split_1_0.ply tools_split.zarr --origin -1,-1 --tile_width 1 --tiles 2,2 "
                    "--tile 1,0 --stem_diameter 0.01"),
            0);
  EXPECT_TRUE(std::ifstream("tools_split.zarr/density/0.0.1").good());
  EXPECT_TRUE(std::ifstream("tools_split.zarr/tiles/1_0.json").good());
}

/// End-to-end check on a synthetic hand-held scan of a random pneumatophore field with a finite beam footprint:
/// the near-ground frontal area per volume profile should match the truth when corrected for the footprint, and be
/// substantially overestimated when not
TEST(Basic, RayDensityPneumatophores)
{
  ray::PneumatophoreParams params;
  params.num_rays = 2000000;
  params.footprint.exit_diameter = 0.01;
  params.footprint.divergence = 0.003;
  params.footprint.threshold = 0.05;
  ray::PneumatophoreGen gen;
  gen.generate(params, 7);
  gen.cloud().save("pneumatophores_test.ply");
  // true N d in the lowest 0.1 m, from the stems within the plot
  const double half_plot = 0.5 * params.plot_width;
  double frontal_area = 0.0;
  for (const auto &stem : gen.stems())
  {
    if (std::abs(stem[0]) < half_plot && std::abs(stem[1]) < half_plot)
    {
      frontal_area += params.diameter * std::min(stem[2], 0.1);
    }
  }
  const double truth = frontal_area / (params.plot_width * params.plot_width * 0.1);

  const ray::Cuboid bounds(Eigen::Vector3d(-half_plot, -half_plot, 0.0), Eigen::Vector3d(half_plot, half_plot, 0.6));
  auto lowest_layer = [&](bool correct_footprint) {
    ray::DensityOptions options;
    options.distribution = ray::ElementDistribution::VerticalCylinder;
    options.path_estimator = ray::PathEstimator::PooledBeerLambert;
    options.flat_top_compensation = false;
    // a uniform (random) stand at 100 stems/m^2, which suits the larger pool (the default 0.05 m is about 7% high
    // here, from few hits per neighbourhood; see docs/pneumatophore_density_method.md)
    options.pool_radius = 0.1;
    ray::DensityOutputs outputs;
    outputs.profile_file = "pneumatophores_test_profile.csv";
    if (correct_footprint)
    {
      options.footprint = params.footprint;
      options.stem_diameter = params.diameter;
    }
    EXPECT_TRUE(ray::renderCloud("pneumatophores_test.ply", bounds, ray::ViewDirection::Top, ray::RenderStyle::Density,
                                 0.05, "pneumatophores_test.png", "", false, nullptr, options, outputs));
    std::ifstream ifs(outputs.profile_file);
    std::string line;
    std::getline(ifs, line);  // header
    EXPECT_EQ(line, "z_min,z_max,observed_voxels,total_voxels,mean_hits_per_voxel,density_frontal_area_per_volume_Nd,"
                    "mean_voxel_density_frontal_area_per_volume_Nd,patchiness_ratio");
    double pooled_sum = 0.0, mean_sum = 0.0;
    for (int i = 0; i < 2; i++)  // the two 0.05 m layers making up the lowest 0.1 m
    {
      std::getline(ifs, line);
      std::stringstream ss(line);
      std::string field;
      for (int j = 0; j < 6; j++)
      {
        std::getline(ss, field, ',');
      }
      pooled_sum += std::stod(field);  // pooled over the layer's voxels
      std::getline(ss, field, ',');
      mean_sum += std::stod(field);  // mean of the voxel densities
    }
    std::cout << (correct_footprint ? "corrected" : "uncorrected") << " mean of the voxel densities: " << 0.5 * mean_sum
              << std::endl;
    return 0.5 * pooled_sum;
  };
  const double corrected = lowest_layer(true);
  const double uncorrected = lowest_layer(false);
  std::cout << "true N*d near ground: " << truth << ", footprint corrected (pooled) estimate: " << corrected
            << ", uncorrected: " << uncorrected << std::endl;
  EXPECT_NEAR(corrected, truth, 0.05 * truth);
  EXPECT_GT(uncorrected, 1.5 * truth);
}

#if RAYLIB_WITH_QHULL
/// Creates a terrain ray cloud, then wraps it from below, comparing the mesh to the expected results
TEST(Basic, RayWrap)
{
  EXPECT_EQ(command("raycreate terrain 1"), 0);
  EXPECT_EQ(command("raywrap terrain.ply upwards 1.0"), 0);
  ray::Mesh mesh;
  EXPECT_TRUE(ray::readPlyMesh("terrain_mesh.ply", mesh));
  EXPECT_EQ(compareMoments(mesh.getMoments(), { 0.0386662, -1.52168, -0.139079, 3.30621, 3.35391, 0.705937 }),
            mesh.getMoments().size());
}

/// Tests extraction of terrain and extraction of trees
TEST(Basic, RayExtract)
{
  EXPECT_EQ(command("raycreate forest 2"), 0);
  EXPECT_EQ(command("rayextract terrain forest.ply"), 0);
  ray::Mesh mesh;
  EXPECT_TRUE(ray::readPlyMesh("forest_mesh.ply", mesh));
  EXPECT_EQ(compareMoments(mesh.getMoments(), { -0.00147491, -0.00191917, -0.0617946, 5.77995, 5.80266, 0.0426993 }),
            mesh.getMoments().size());
  EXPECT_EQ(command("rayextract trees forest.ply forest_mesh.ply"), 0);

  ray::ForestStructure forest;
  EXPECT_TRUE(forest.load("forest_trees.txt"));
  EXPECT_EQ(compareMomentsPercentageError(forest.getMoments(),
                                          { 20, 22.2172, 1058.42, 1.39416, 0.112794, 2.2863, 21918, 0, 105.4 }),
            forest.getMoments().size());

  EXPECT_EQ(command("rayextract forest forest.ply --ground forest_mesh.ply"), 0);

  ray::ForestStructure forest2;
  EXPECT_TRUE(forest2.load("forest_forest.txt"));
  EXPECT_EQ(compareMomentsPercentageError(forest2.getMoments(),
                                          { 11, 8.43829, 586.427, 1.40054, 0.200644, 0, 16697, 8.49419, 3.09917 }),
            forest2.getMoments().size());

  EXPECT_EQ(command("rayextract trunks forest.ply"), 0);

  ray::ForestStructure forest3;
  EXPECT_TRUE(forest3.load("forest_trunks.txt"));
  EXPECT_EQ(
    compareMomentsPercentageError(forest3.getMoments(), { 21, 20.0797, 1124.61, 1.60427, 0.135159, 0, 0, 0, 0 }),
    forest3.getMoments().size());
}
#endif  // RAYLIB_WITH_QHULL
}  // namespace raytest
