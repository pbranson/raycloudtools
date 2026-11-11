// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raylib/raybuildinggen.h"
#include "raylib/raycloud.h"
#include "raylib/rayforestgen.h"
#include "raylib/rayparse.h"
#include "raylib/raypneumatophoregen.h"
#include "raylib/rayroomgen.h"
#include "raylib/rayterraingen.h"
#include "raylib/raytreegen.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

void usage(int exit_code = 1)
{
  // clang-format off
  std::cout << "Generates simple example ray clouds" << std::endl;
  std::cout << "usage (--view / -v to view results):" << std::endl;
  std::cout << "raycreate room 3 - generates a room using the seed 3. Also:" << std::endl;
  std::cout << "          building" << std::endl;
  std::cout << "          tree" << std::endl;
  std::cout << "          forest" << std::endl;
  std::cout << "          terrain" << std::endl;
  std::cout << "          field" << std::endl;
  std::cout << "          pneumatophores - vertical stems on flat ground, scanned by a hand-held lidar. Also outputs" << std::endl;
  std::cout << "                           pneumatophores_truth.csv, the true frontal area per volume N*d profile," << std::endl;
  std::cout << "                           and pneumatophores_stems.csv, the stem positions and heights" << std::endl;
  std::cout << "             --stems_per_m2 100     - mean stem density" << std::endl;
  std::cout << "             --diameter 0.01        - stem diameter (m)" << std::endl;
  std::cout << "             --clustered            - stems along cable roots radiating from trees, else random" << std::endl;
  std::cout << "             --cluster_spread 0.03  - standard deviation of stem distance from the cable roots (m)" << std::endl;
  std::cout << "             --patchy 0.5,0.3,0     - stems only in random circular patches (e.g. degraded areas), given" << std::endl;
  std::cout << "                                      patch radius (m), ground cover fraction and the stem density" << std::endl;
  std::cout << "                                      between patches, relative to within them" << std::endl;
  std::cout << "             --walk_spacing 1       - spacing of the lawnmower walking lines (m)" << std::endl;
  std::cout << "             --zenith 0,101.5       - range of ray zenith angles (degrees from straight down)" << std::endl;
  std::cout << "             --truth_cells 1,1,0.1  - also output pneumatophores_truth_cells.csv for cells of this size" << std::endl;
  std::cout << "             --footprint 0.01,0.003,0.05 - beam exit diameter (m), divergence (rad), detection threshold" << std::endl;
  std::cout << "             --return_range surface - with --footprint, returns at range to the stem surface (default), or with" << std::endl;
  std::cout << "                                      outline, where the beam axis meets the footprint-widened stem outline" << std::endl;
  std::cout << "             --rays 2000000         - number of rays" << std::endl;
  std::cout << "             --height 0.25,0.08     - mean and standard deviation of the (normal) stem heights (m). A" << std::endl;
  std::cout << "                                      standard deviation of 0 gives stems of one height" << std::endl;
  std::cout << "             --uniform_height 0.1,0.4 - stem heights uniform between these (m), instead of normal" << std::endl;
  std::cout << "             --scanner_height 1.2   - height of the hand-held scanner (m)" << std::endl;
  std::cout << "             --plot_width 6         - width of the square plot that the scanner walks (m)" << std::endl;
  std::cout << "             --range_noise 0.01     - standard deviation of the range noise of returns (m)" << std::endl;
  std::cout << "             --pose_noise 0.02,1    - standard deviation (m) of a smooth trajectory error, per axis, and" << std::endl;
  std::cout << "                                      the walk distance (m) between its independent values" << std::endl;
  std::cout << "             --ground_returns       - keep ground returns as bounded rays (default: zero alpha)" << std::endl;
  std::cout << "             --trunks 0.1,0.2       - sparse trunks per m^2 and mean diameter (m), with a trunks.csv file" << std::endl;
  std::cout << "             --min_gap 0.5          - minimum gap between stem surfaces, in stem diameters (default 0)" << std::endl;
  std::cout << "             --overlapping          - independent (Poisson) stem positions, which may overlap, and a Poisson" << std::endl;
  std::cout << "                                      number of stems. Default: exact number of non-overlapping stems" << std::endl;
  std::cout << "             These options apply only to pneumatophores" << std::endl;
  std::cout << std::endl;
  std::cout << "          forest trees.txt - generate from a comma-separated list of x,y,z,radius trees" << std::endl;
  std::cout << "          terrain mesh.ply - generate from a ground mesh" << std::endl;
  // clang-format on
  exit(exit_code);
}

double random(double x, double y){ return ray::random(x, y); }

int rayCreate(int argc, char *argv[])
{
  ray::KeyChoice cloud_type({ "room", "building", "tree", "forest", "terrain", "field" });
  ray::TextArgument pneumatophores_type("pneumatophores");
  ray::IntArgument seed(1, 1000000);
  ray::OptionalFlagArgument view_flag("view", 'v');
  ray::FileArgument input_file;
  // pneumatophore options
  ray::DoubleArgument stems_per_m2(0.01, 100000.0, 100.0), diameter(0.0001, 1.0, 0.01), spread(0.0, 100.0, 0.03);
  ray::IntArgument num_rays(1, 1000000000, 2000000);
  ray::Vector3dArgument footprint(0.0, 1.0), patchy(0.0, 1000.0), truth_cells(0.0001, 1000.0);
  ray::DoubleArgument walk_spacing(0.01, 1000.0, 1.0), min_gap(0.0, 100.0);
  ray::Vector2dArgument zenith(0.0, 180.0), trunks(0.0, 100.0);
  ray::Vector2dArgument height(0.0, 100.0), uniform_height(0.0, 100.0), pose_noise(0.0, 100.0);
  ray::DoubleArgument scanner_height(0.01, 1000.0, 1.2), plot_width(0.1, 10000.0, 6.0), range_noise(0.0, 10.0);
  ray::KeyChoice return_range({ "surface", "outline" });
  ray::OptionalKeyValueArgument stems_option("stems_per_m2", 's', &stems_per_m2);
  ray::OptionalKeyValueArgument diameter_option("diameter", 'd', &diameter);
  ray::OptionalKeyValueArgument rays_option("rays", 'r', &num_rays);
  ray::OptionalKeyValueArgument footprint_option("footprint", 'b', &footprint);
  ray::OptionalKeyValueArgument spread_option("cluster_spread", 'p', &spread);
  ray::OptionalKeyValueArgument patchy_option("patchy", 'a', &patchy);
  ray::OptionalKeyValueArgument walk_option("walk_spacing", 'w', &walk_spacing);
  ray::OptionalKeyValueArgument min_gap_option("min_gap", 'm', &min_gap);
  ray::OptionalKeyValueArgument trunks_option("trunks", 'k', &trunks);
  ray::OptionalKeyValueArgument zenith_option("zenith", 'z', &zenith);
  ray::OptionalKeyValueArgument truth_cells_option("truth_cells", 't', &truth_cells);
  ray::OptionalKeyValueArgument return_range_option("return_range", 'e', &return_range);
  ray::OptionalKeyValueArgument height_option("height", 'h', &height);
  ray::OptionalKeyValueArgument uniform_height_option("uniform_height", 'u', &uniform_height);
  ray::OptionalKeyValueArgument scanner_height_option("scanner_height", 'n', &scanner_height);
  ray::OptionalKeyValueArgument plot_width_option("plot_width", 'l', &plot_width);
  ray::OptionalKeyValueArgument range_noise_option("range_noise", 'x', &range_noise);
  ray::OptionalKeyValueArgument pose_noise_option("pose_noise", 'j', &pose_noise);
  ray::OptionalFlagArgument clustered_flag("clustered", 'c'), ground_flag("ground_returns", 'g'),
    overlapping_flag("overlapping", 'o');
  bool from_seed = ray::parseCommandLine(argc, argv, { &cloud_type, &seed }, { &view_flag });
  bool from_file = ray::parseCommandLine(argc, argv, { &cloud_type, &input_file }, { &view_flag });
  // the pneumatophore options are accepted only for pneumatophores, which are generated from a seed
  const bool pneumatophores = ray::parseCommandLine(
    argc, argv, { &pneumatophores_type, &seed },
    { &view_flag, &stems_option, &diameter_option, &rays_option, &footprint_option, &clustered_flag, &ground_flag,
      &spread_option, &patchy_option, &walk_option, &zenith_option, &truth_cells_option, &overlapping_flag,
      &min_gap_option, &trunks_option, &return_range_option, &height_option, &uniform_height_option,
      &scanner_height_option, &plot_width_option, &range_noise_option, &pose_noise_option });
  if (!from_seed && !from_file && !pneumatophores)
    usage();

  if (from_seed || pneumatophores)
  {
    ray::srand(seed.value());
  }

  ray::Cloud cloud;
  const double time_delta = 0.001;  // between rays
  std::string type = pneumatophores ? "pneumatophores" : cloud_type.selectedKey();
  if (type == "room")
  {
    // create room
    ray::RoomGen room_gen;
    room_gen.generate();
    cloud.starts = room_gen.rayStarts();
    cloud.ends = room_gen.rayEnds();
    cloud.times.resize(cloud.starts.size());
    double time = 0.0;
    for (size_t i = 0; i < cloud.starts.size(); i++)
    {
      cloud.times[i] = time;
      if (i == cloud.starts.size() / 2)
        time += 0.5;
      time += time_delta;
    }
    colourByTime(cloud.times, cloud.colours);
    const std::vector<bool> &bounded = room_gen.rayBounded();
    for (int i = 0; i < (int)cloud.colours.size(); i++) cloud.colours[i].alpha = bounded[i] ? 255 : 0;
  }
  else if (type == "field") // like a field of wheat. The rays follow a Beer/Lambert law in their depth
  {
    double half_width = random(5.0, 20.0);
    double crop_height = random(0.2, 1.5);
    double height_gradient = 0.0;
    double density_gradient = 0.0;
    double crop_density = random(0.04, 0.5);
    double ray_density = random(400.0, 1600.0); // rays per square metre
    std::cout << "Creating field of width " << 2.0*half_width << " m, height: " << crop_height << " m, crop density " << crop_density << " and " << ray_density << " vertical rays per m^2" << std::endl;
    int num_rays = (int)(ray_density * 4.0 * half_width*half_width);
    if (random(0.0,1.0) < 0.5 && seed.value() != 1)
    {
      height_gradient = random(-crop_height/half_width, crop_height/half_width);
      std::cout << "height gradient in y: " << height_gradient << std::endl;
    }
    if (random(0.0, 1.0) < 0.5)
    {
      density_gradient = random(-crop_density/half_width, crop_density/half_width);  
      std::cout << "density gradient in x: " << density_gradient << std::endl;
    }

    for (int i = 0; i<num_rays; i++)
    {
      Eigen::Vector3d pos(random(-half_width, half_width), random(-half_width, half_width), 0.0);
      double lambda = crop_density + pos[0]*density_gradient;
      double height = crop_height + pos[1]*height_gradient;
      Eigen::Vector3d start(pos[0], pos[1], 2.0*height);

      double prob_ground = 1.0 - std::exp(-lambda * height);
      double rnd = random(0.0, 1.0);
      if (rnd < prob_ground)
      {
        double h = std::log(1.0 - rnd)/-lambda;
        pos[2] = height - h;
        if (pos[2] < -0.01)
          std::cout << "error: " << pos[2] << ", rnd: " << rnd << ", h: " << h << std::endl;
      }
      cloud.addRay(start, pos, 0.001*(double)i, ray::RGBA(255,255,255,255));
    }
    colourByTime(cloud.times, cloud.colours);
  }
  else if (type == "building")
  {
    // create building...
    ray::BuildingGen building_gen;
    building_gen.generate();
    cloud.starts = building_gen.rayStarts();
    cloud.ends = building_gen.rayEnds();
    cloud.times.resize(cloud.ends.size());
    double time = 0.0;
    for (size_t i = 0; i < cloud.starts.size(); i++)
    {
      cloud.times[i] = time;
      time += time_delta;
    }
    colourByTime(cloud.times, cloud.colours);
    const std::vector<bool> &bounded = building_gen.rayBounded();
    for (int i = 0; i < (int)cloud.colours.size(); i++) cloud.colours[i].alpha = bounded[i] ? 255 : 0;
  }
  else if (type == "tree" || type == "forest")
  {
    const double density = 500.0;                   // density of points on the branches of the trees
    const double tree_ground_extent = 2.0;          // extent (from centre) is the half-width
    const double ground_noise_extent = 0.025;       // vertical noise in the ground
    const double forest_ground_multiplier = 5.0;    // how much larger the forest ground is
    const double ground_ray_deviation = 0.1;        // lateral deviation in start of ray, from end point
    const double ground_ray_vertical_height = 1.5;  // height above ground for ray start

    ray::fillBranchAngleLookup();
    Eigen::Vector3d box_min(-tree_ground_extent, -tree_ground_extent, -ground_noise_extent);
    Eigen::Vector3d box_max(tree_ground_extent, tree_ground_extent, ground_noise_extent);
    double time = 0.0;
    if (type == "tree")  // create a single tree
    {
      ray::TreeGen tree_gen;
      ray::TreeParams params;
      params.random_factor = 0.25;
      tree_gen.segments().resize(1);
      tree_gen.segments()[0].tip = Eigen::Vector3d(0, 0, 0);
      tree_gen.segments()[0].radius = 0.1;
      tree_gen.make(params);
      tree_gen.generateRays(density);
      cloud.starts = tree_gen.rayStarts();
      cloud.ends = tree_gen.rayEnds();
      cloud.times.resize(cloud.ends.size());
      for (size_t i = 0; i < cloud.starts.size(); i++)
      {
        cloud.times[i] = time;
        time += time_delta;
      }
      colourByTime(cloud.times, cloud.colours);
    }
    else if (type == "forest")  // create multiple trees on a plane
    {
      ray::ForestParams params;
      params.random_factor = 0.25;
      ray::ForestGen forest_gen;
      if (from_file)  // load a forest from an _trees.txt file
      {
        if (!forest_gen.makeFromFile(input_file.name(), params))
        {
          usage();
        }
      }
      else
      {
        forest_gen.make(params);
      }
      forest_gen.generateRays(density);
      for (auto &tree : forest_gen.trees())
      {
        const std::vector<Eigen::Vector3d> &ray_starts = tree.rayStarts();
        const std::vector<Eigen::Vector3d> &ray_ends = tree.rayEnds();
        cloud.starts.insert(cloud.starts.end(), ray_starts.begin(), ray_starts.end());
        cloud.ends.insert(cloud.ends.end(), ray_ends.begin(), ray_ends.end());
      }
      cloud.times.resize(cloud.starts.size());
      for (size_t i = 0; i < cloud.starts.size(); i++)
      {
        cloud.times[i] = time;
        time += time_delta;
      }
      box_min *= forest_ground_multiplier;  // for a forest, we need a larger ground
      box_max *= forest_ground_multiplier;
    }
    if (!from_file)
    {
      int num = int(0.25 * density * (box_max[0] - box_min[0]) * (box_max[1] - box_min[1]));
      for (int i = 0; i < num; i++)
      {
        Eigen::Vector3d pos(ray::random(box_min[0], box_max[0]), ray::random(box_min[1], box_max[1]),
                            ray::random(box_min[2], box_max[2]));
        cloud.ends.push_back(pos);
        cloud.starts.push_back(pos + Eigen::Vector3d(ray::random(-ground_ray_deviation, ground_ray_deviation),
                                                     ray::random(-ground_ray_deviation, ground_ray_deviation),
                                                     ground_ray_vertical_height));
        cloud.times.push_back(time);
        time += time_delta;
      }
    }
    colourByTime(cloud.times, cloud.colours);
  }
  else if (type == "terrain")
  {
    ray::TerrainGen terrain;
    if (from_file)  // generate ray cloud terrain from a .ply mesh file
    {
      terrain.generateFromFile(input_file.name());
    }
    else
    {
      terrain.generate();
    }
    cloud.starts = terrain.rayStarts();
    cloud.ends = terrain.rayEnds();
    cloud.times.resize(cloud.starts.size());
    double time = 0.0;
    for (size_t i = 0; i < cloud.starts.size(); i++)
    {
      cloud.times[i] = time;
      time += time_delta;
    }
    colourByTime(cloud.times, cloud.colours);
  }
  else if (type == "pneumatophores")
  {
    ray::PneumatophoreParams params;
    params.stems_per_m2 = stems_per_m2.value();
    params.diameter = diameter.value();
    params.num_rays = num_rays.value();
    params.clustered = clustered_flag.isSet();
    params.lateral_spread = spread.value();
    params.walk_spacing = walk_spacing.value();
    params.hard_core = !overlapping_flag.isSet();
    params.min_gap_ratio = min_gap.value();
    if (trunks_option.isSet())
    {
      params.trunks_per_m2 = trunks.value()[0];
      params.trunk_diameter = trunks.value()[1];
    }
    if (patchy_option.isSet())
    {
      params.patchy = true;
      params.patch_radius = patchy.value()[0];
      params.patch_cover = patchy.value()[1];
      params.matrix_ratio = patchy.value()[2];
    }
    if (zenith_option.isSet())
    {
      params.min_zenith_deg = zenith.value()[0];
      params.max_zenith_deg = zenith.value()[1];
    }
    params.ground_unbounded = !ground_flag.isSet();
    if (height_option.isSet() && uniform_height_option.isSet())
    {
      std::cerr << "Error: --height and --uniform_height are alternatives" << std::endl;
      return 1;
    }
    if (height_option.isSet())
    {
      params.mean_height = height.value()[0];
      params.height_sd = height.value()[1];
      if (params.height_sd == 0.0 && params.mean_height < 0.02)
      {
        std::cerr << "Error: stems must be at least 0.02 m tall" << std::endl;
        return 1;
      }
    }
    if (uniform_height_option.isSet())
    {
      params.height_distribution = ray::PneumatophoreParams::HeightDistribution::Uniform;
      params.min_height = uniform_height.value()[0];
      params.max_height = uniform_height.value()[1];
      if (params.max_height < params.min_height)
      {
        std::cerr << "Error: --uniform_height needs the minimum height first" << std::endl;
        return 1;
      }
    }
    params.scanner_height = scanner_height.value();
    params.plot_width = plot_width.value();
    params.range_noise = range_noise.value();
    if (pose_noise_option.isSet())
    {
      params.pose_noise = pose_noise.value()[0];
      params.pose_noise_interval = std::max(0.01, pose_noise.value()[1]);
    }
    if (footprint_option.isSet())
    {
      params.footprint.exit_diameter = footprint.value()[0];
      params.footprint.divergence = footprint.value()[1];
      params.footprint.threshold = footprint.value()[2];
    }
    if (return_range_option.isSet())
    {
      if (!footprint_option.isSet())
      {
        std::cerr << "Error: --return_range applies only with --footprint" << std::endl;
        return 1;
      }
      params.footprint.range_to_outline = return_range.selectedKey() == "outline";
    }
    ray::PneumatophoreGen gen;
    try
    {
      gen.generate(params, static_cast<unsigned int>(seed.value()));
    }
    catch (const std::runtime_error &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;
    }
    cloud = gen.cloud();
    if (!gen.writeTruthProfile(type + "_truth.csv", 0.01))
    {
      return 1;
    }
    if (!gen.writeStems(type + "_stems.csv") || !gen.writeTrunks(type + "_trunks.csv"))
    {
      return 1;
    }
    if (truth_cells_option.isSet() && !gen.writeTruthCells(type + "_truth_cells.csv", truth_cells.value()))
    {
      return 1;
    }
  }
  else
    usage();
  cloud.save(type + ".ply");
  if (view_flag.isSet())
    ray::viewFile(type + ".ply");

  return 0;
}

int main(int argc, char *argv[])
{
  return ray::runWithMemoryCheck(rayCreate, argc, argv);
}