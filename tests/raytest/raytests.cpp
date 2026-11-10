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
#include "raylib/rayforeststructure.h"
#include "raylib/raylaz.h"
#include "raylib/raymesh.h"
#include "raylib/rayply.h"
#include "raylib/raypneumatophoregen.h"
#include "raylib/rayrenderer.h"
#include "raylib/raythreads.h"

#include <fstream>
#include <random>
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
    options.path_estimator = ray::PathEstimator::Crossing;
    options.flat_top_compensation = false;
    options.min_sin_zenith = 0.9;
    options.profile_file = "pneumatophores_test_profile.csv";
    if (correct_footprint)
    {
      options.footprint = params.footprint;
      options.stem_diameter = params.diameter;
    }
    EXPECT_TRUE(ray::renderCloud("pneumatophores_test.ply", bounds, ray::ViewDirection::Top, ray::RenderStyle::Density,
                                 0.05, "pneumatophores_test.png", "", false, nullptr, options));
    std::ifstream ifs(options.profile_file);
    std::string line;
    std::getline(ifs, line);  // header
    double layer_sum = 0.0;
    for (int i = 0; i < 2; i++)  // the two 0.05 m layers making up the lowest 0.1 m
    {
      std::getline(ifs, line);
      std::stringstream ss(line);
      std::string field;
      for (int j = 0; j < 6; j++)
      {
        std::getline(ss, field, ',');
      }
      layer_sum += std::stod(field);  // mean of the voxel densities
    }
    return 0.5 * layer_sum;
  };
  const double corrected = lowest_layer(true);
  const double uncorrected = lowest_layer(false);
  std::cout << "true N*d near ground: " << truth << ", footprint corrected estimate: " << corrected
            << ", uncorrected: " << uncorrected << std::endl;
  // the crossing estimator is biased low by stems occluding each other within a voxel, by about half the chance
  // that a ray crossing a voxel is intercepted (here N d_eff * 0.05 m ~ 0.12, so about -6%)
  EXPECT_NEAR(corrected, truth, 0.07 * truth);
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
