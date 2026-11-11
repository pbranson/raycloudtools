// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <optional>

#include "raylib/raycloud.h"
#include "raylib/raycloudwriter.h"
#include "raylib/raymesh.h"
#include "raylib/rayparse.h"
#include "raylib/rayply.h"

#define DETAILED_GROUND_SUBTRACTION  // ground lookup per-point. Otherwise it quantises to a height grid

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Translate a raycloud" << std::endl;
  std::cout << "usage (--view / -v to view results):" << std::endl;
  std::cout << "raytranslate raycloud 0,0,1 - translation (x,y,z) in metres" << std::endl;
  std::cout << "                      0,0,1,24.3 - optional 4th component translates time" << std::endl;
  std::cout << "                      subtract ground_mesh.ply  - translate vertically to remove ground_mesh heights" << std::endl;
  std::cout << "                      add ground_mesh.ply       - translate vertically to add ground_mesh heights" << std::endl;
  std::cout << "                     --shear - shift each ray start by the ground height beneath the start, and each" << std::endl;
  std::cout << "                               end by the ground height beneath the end. This shears the rays as if" << std::endl;
  std::cout << "                               the scene were physically flattened, so the whole ray path keeps its" << std::endl;
  std::cout << "                               height above ground. Recommended before ray-based density or free-space" << std::endl;
  std::cout << "                               analysis on sloping ground. Starts with no ground mesh beneath them use" << std::endl;
  std::cout << "                               the end's height instead (and are counted)." << std::endl;
  std::cout << "                               Without --shear (default) each ray is moved rigidly by the height beneath" << std::endl;
  std::cout << "                               its end, which preserves ray directions and lengths, but misplaces points" << std::endl;
  std::cout << "                               along the ray by (ground slope x horizontal distance from the end)." << std::endl;
  std::cout << "                     --split_untranslated - write any points that could not be" << std::endl;
  std::cout << "                                            translated to a new file (otherwise" << std::endl;
  std::cout << "                                            leave as they were in output raycloud)" << std::endl;
  // clang-format on
  exit(exit_code);
}

int rayTranslate(int argc, char *argv[])
{
  ray::FileArgument cloud_file, ground_file;
  ray::TextArgument subtract("subtract"), add("add");
  ray::OptionalFlagArgument view_flag("view", 'v');
  ray::Vector3dArgument translation3;
  ray::Vector4dArgument translation4;
  ray::OptionalFlagArgument split_untranslated("split_untranslated", 'u');
  ray::OptionalFlagArgument shear("shear", 's');

  bool vec3_format = ray::parseCommandLine(argc, argv, { &cloud_file, &translation3 }, { &view_flag });
  bool vec4_format = ray::parseCommandLine(argc, argv, { &cloud_file, &translation4 }, { &view_flag });
  bool ground_subtract_format = ray::parseCommandLine(argc, argv, { &cloud_file, &subtract, &ground_file },
                                                      { &view_flag, &split_untranslated, &shear });
  bool ground_add_format = ray::parseCommandLine(argc, argv, { &cloud_file, &add, &ground_file },
                                                 { &view_flag, &split_untranslated, &shear });
  if (!vec3_format && !vec4_format && !ground_subtract_format && !ground_add_format)
    usage();

  Eigen::Vector3d translation(0, 0, 0);
  double time_delta = 0.0;
  ray::Mesh ground_mesh;
  ray::Cloud::Info info;
  Eigen::Vector3d min_bound;
  Eigen::Vector3d max_bound;
#if defined DETAILED_GROUND_SUBTRACTION
  std::vector<ray::Triangle> triangles;
  ray::Grid<ray::Triangle *> grid;
#else
  Eigen::ArrayXXd ground_heights;
#endif
  double voxel_width = 0.25;
  if (vec3_format)
    translation = translation3.value();
  else if (vec4_format)
  {
    translation = translation4.value().head<3>();
    time_delta = translation4.value()[3];
  }
  else
  {
    if (!ray::readPlyMesh(ground_file.name(), ground_mesh))
    {
      usage();
    }
    if (!ray::Cloud::getInfo(cloud_file.name(), info))
    {
      usage();
    }
    min_bound = info.ends_bound.min_bound_;
    max_bound = info.ends_bound.max_bound_;
    // the toHeightfield function only casts rays from max_bound[2] down to min_bound[2], so set these according to
    // ground_file's extents, not cloud_file's extents.
    min_bound[2] = std::numeric_limits<double>::max();
    max_bound[2] = std::numeric_limits<double>::lowest();
    Eigen::Vector3d mesh_min(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), 0.0);
    Eigen::Vector3d mesh_max(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), 0.0);
    for (auto &vert : ground_mesh.vertices())
    {
      min_bound[2] = std::min(min_bound[2], vert[2]);
      max_bound[2] = std::max(max_bound[2], vert[2]);
      mesh_min = ray::minVector(mesh_min, vert);
      mesh_max = ray::maxVector(mesh_max, vert);
    }
    min_bound[2] -= 0.1;
    max_bound[2] += 0.1;
    if (shear.isSet())
    {
      // Ray starts can lie outside the lateral bounds of the ray ends (e.g. a hand-held scanner walked around the
      // edge of the scan), so extend the lookup grid to cover the starts, but only where the ground mesh is, so that
      // distant starts (e.g. airborne) don't create a huge grid. Starts outside the grid fall back to the end height.
      for (int i = 0; i < 2; i++)
      {
        const double lo = std::max(info.starts_bound.min_bound_[i], mesh_min[i]);
        const double hi = std::min(info.starts_bound.max_bound_[i], mesh_max[i]);
        if (lo <= hi)
        {
          min_bound[i] = std::min(min_bound[i], lo);
          max_bound[i] = std::max(max_bound[i], hi);
        }
      }
    }
#if defined DETAILED_GROUND_SUBTRACTION
    Eigen::Vector3d box_min = min_bound;
    Eigen::Vector3d box_max = max_bound;
    box_max[2] = box_min[2] + 0.5 * voxel_width;  // ensure that the grid is only 1 voxel high
    triangles.resize(ground_mesh.indexList().size());
    for (int i = 0; i < (int)ground_mesh.indexList().size(); i++)
    {
      ray::Triangle &tri = triangles[i];
      for (int j = 0; j < 3; j++)
      {
        tri.corners[j] = ground_mesh.vertices()[ground_mesh.indexList()[i][j]];
      }
      tri.tested = false;
      tri.normal = (tri.corners[1] - tri.corners[0]).cross(tri.corners[2] - tri.corners[0]);
    }

    // put the triangles into a grid
    grid.init(box_min, box_max, voxel_width);
    for (auto &tri : triangles)
    {
      Eigen::Vector3d tri_min =
        (ray::minVector(tri.corners[0], ray::minVector(tri.corners[1], tri.corners[2])) - box_min) / voxel_width;
      Eigen::Vector3d tri_max =
        (ray::maxVector(tri.corners[0], ray::maxVector(tri.corners[1], tri.corners[2])) - box_min) / voxel_width;
      // only insert into cells within the grid. Out of range (e.g. negative) indices are not valid grid cells.
      const int x_min = std::max(0, (int)std::floor(tri_min[0]));
      const int y_min = std::max(0, (int)std::floor(tri_min[1]));
      const int x_max = std::min(grid.dims[0] - 1, (int)std::floor(tri_max[0]));
      const int y_max = std::min(grid.dims[1] - 1, (int)std::floor(tri_max[1]));
      for (int x = x_min; x <= x_max; x++)
      {
        for (int y = y_min; y <= y_max; y++)
        {
          grid.insert(x, y, 0, &tri);
        }
      }
    }
#else
    ground_mesh.toHeightField(ground_heights, min_bound, max_bound, voxel_width);
#endif
  }

  // Stats to print at end.
  int num_missed_triangles = 0, num_totally_missed = 0, num_start_fallbacks = 0;

  // Look up the ground height directly beneath @c pos. Returns no value if there is no ground mesh there.
  // @c approximated is set when no triangle lies directly beneath @c pos, so the height is estimated from the planes
  // of the triangles in the same grid cell.
  auto groundHeight = [&](const Eigen::Vector3d &pos, bool &approximated) -> std::optional<double> {
    approximated = false;
    if (pos[0] < min_bound[0] || pos[1] < min_bound[1] || pos[0] > max_bound[0] || pos[1] > max_bound[1])
    {
      return std::nullopt;
    }
#if defined DETAILED_GROUND_SUBTRACTION
    // clamp so that points on the upper bound use the last cell
    const int x = std::min((int)std::floor((pos[0] - min_bound[0]) / voxel_width), grid.dims[0] - 1);
    const int y = std::min((int)std::floor((pos[1] - min_bound[1]) / voxel_width), grid.dims[1] - 1);
    if (x < 0 || y < 0)
    {
      return std::nullopt;
    }
    // intersect a vertical ray directly through the point (not the cell centre), so the height is exact per point
    Eigen::Vector3d pos_top(pos[0], pos[1], max_bound[2]);
    Eigen::Vector3d pos_base(pos[0], pos[1], min_bound[2]);
    auto &tris = grid.cell(x, y, 0).data;
    // search the triangles in this cell 'bucket'
    double mean_depth = 0.0;
    double mean_count = 0.0;
    for (auto &tri : tris)
    {
      double depth = 0.0;
      if (tri->intersectsRay(pos_top, pos_base, depth))
      {
        return pos_top[2] + (pos_base[2] - pos_top[2]) * depth;
      }
      if (depth != 0.0)
      {
        mean_depth += depth;
        mean_count++;
      }
    }
    if (mean_count > 0)
    {
      mean_depth /= mean_count;
      approximated = true;
      return pos_top[2] + (pos_base[2] - pos_top[2]) * mean_depth;
    }
    return std::nullopt;
#else
    const int x = std::min((int)std::floor((pos[0] - min_bound[0]) / voxel_width), (int)ground_heights.rows() - 1);
    const int y = std::min((int)std::floor((pos[1] - min_bound[1]) / voxel_width), (int)ground_heights.cols() - 1);
    if (x < 0 || y < 0 || ground_heights(x, y) == std::numeric_limits<double>::lowest())
    {
      return std::nullopt;
    }
    return ground_heights(x, y);
#endif
  };

  // Normal output of translated rays.
  ray::Cloud translated_chunk;
  ray::CloudWriter translated_cloud_writer;
  const std::string temp_name = cloud_file.nameStub() + "~.ply";  // tilde is a common suffix for temporary files
  if (!translated_cloud_writer.begin(temp_name))
  {
    usage();
  }

  // Alternative output for rays that could not be translated if
  // --split_untranslated arg is provided.
  ray::Cloud untranslated_chunk;
  ray::CloudWriter untranslated_cloud_writer;
  const std::string untranslated_name = cloud_file.nameStub() + "_untranslated.ply";
  if (split_untranslated.isSet())
  {
    if (!untranslated_cloud_writer.begin(untranslated_name))
    {
      usage();
    }
  }

  // Translate a ray then push it into the appropriate chunk.
  auto translate = [&](Eigen::Vector3d &start, Eigen::Vector3d &end, double &time, ray::RGBA &colour) {
    if (ground_subtract_format || ground_add_format)
    {
      bool approximated = false;
      const std::optional<double> end_height = groundHeight(end, approximated);
      if (approximated)
      {
        num_missed_triangles++;
      }
      if (end_height.has_value())
      {
        double start_height = end_height.value();
        if (shear.isSet())
        {
          bool start_approximated = false;
          const std::optional<double> height = groundHeight(start, start_approximated);
          if (height.has_value())
          {
            start_height = height.value();
          }
          else
          {
            num_start_fallbacks++;
          }
        }
        const double sign = ground_subtract_format ? -1.0 : 1.0;
        start[2] += sign * start_height;
        end[2] += sign * end_height.value();
      }
      else
      {
        num_totally_missed++;
        if (split_untranslated.isSet())
        {
          // Push into chunk to be written to untranslated output.
          untranslated_chunk.starts.push_back(start);
          untranslated_chunk.ends.push_back(end);
          untranslated_chunk.times.push_back(time);
          untranslated_chunk.colours.push_back(colour);
          // Return to avoid also pushing into the normal output.
          return;
        }
      }
    }
    else
    {
      start += translation;
      end += translation;
      time += time_delta;
    }
    // Push into chunk to be written to normal output.
    translated_chunk.starts.push_back(start);
    translated_chunk.ends.push_back(end);
    translated_chunk.times.push_back(time);
    translated_chunk.colours.push_back(colour);
  };
  // Translate all rays in a chunk then write the output chunks.
  auto per_chunk = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends,
                       std::vector<double> &times, std::vector<ray::RGBA> &colours) {
    for (size_t index = 0; index < ends.size(); ++index)
    {
      translate(starts[index], ends[index], times[index], colours[index]);
    }
    translated_cloud_writer.writeChunk(translated_chunk);
    translated_chunk.clear();
    if (split_untranslated.isSet())
    {
      untranslated_cloud_writer.writeChunk(untranslated_chunk);
      untranslated_chunk.clear();
    }
  };
  // Process all chunks.
  if (!ray::readPly(cloud_file.name(), true, per_chunk, 0))
    usage();

  // Finish writing.
  translated_cloud_writer.end();
  if (split_untranslated.isSet())
  {
    untranslated_cloud_writer.end();
  }

  // Log stats.
  if (num_missed_triangles > 2)
  {
    std::cout << num_missed_triangles
              << " cloud points not overlapping triangles, if this is large the ground file may not be a lateral "
                 "coverage of the ray cloud."
              << std::endl;
  }
  if (num_start_fallbacks > 0)
  {
    std::cout << num_start_fallbacks
              << " ray starts have no ground mesh beneath them, so were shifted by the ground height beneath their end "
                 "point instead."
              << std::endl;
  }
  if (num_totally_missed > 0)
  {
    std::cout << "Warning: " << num_totally_missed
              << " cloud points have no laterally overlapping triangles, so the ground_file is not a full coverage."
              << std::endl;
    if (split_untranslated.isSet())
    {
      std::cout << "Wrote " << num_totally_missed << " rays that had no laterally overlapping triangles to "
                << untranslated_name << std::endl;
    }
    else
    {
      std::cout << "Left " << num_totally_missed << " rays that had no laterally overlapping triangles as they were in "
                << cloud_file.name() << std::endl;
    }
  }
  std::rename(temp_name.c_str(), cloud_file.name().c_str());
  if (view_flag.isSet())
    ray::viewFile(cloud_file.name());
  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayTranslate, argc, argv);
}
