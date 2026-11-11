// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "rayrenderer.h"
#include "imagewrite.h"
#include "raycloud.h"
#include "raylib/raylibconfig.h"
#include "rayparse.h"
#if RAYLIB_WITH_TIFF   // build option to support outputting to geotif (.tif) format
#include "geotiffio.h" /* for GeoTIFF */
#include "xtiffio.h"   /* for TIFF */
#endif
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include "rayunused.h"

namespace ray
{
#if RAYLIB_WITH_TIFF

namespace
{
struct Proj4Data
{
  std::unordered_map<std::string, std::string> values;
  std::unordered_set<std::string> flags;
};

std::string trim(const std::string &s)
{
  std::string::size_type b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b])))
  {
    ++b;
  }
  std::string::size_type e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1])))
  {
    --e;
  }
  return s.substr(b, e - b);
}

std::string toLower(std::string s)
{
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

bool tryParseInt(const std::string &s, int &value)
{
  char *end = nullptr;
  const long v = std::strtol(s.c_str(), &end, 10);
  if (end == s.c_str() || *end != '\0')
    return false;
  value = static_cast<int>(v);
  return true;
}

bool tryParseDouble(const std::string &s, double &value)
{
  char *end = nullptr;
  const double v = std::strtod(s.c_str(), &end);
  if (end == s.c_str() || *end != '\0')
    return false;
  value = v;
  return true;
}

bool parseProj4File(const std::string &projection_file, Proj4Data &proj4)
{
  std::ifstream ifs(projection_file.c_str(), std::ios::in);
  if (ifs.fail())
    return false;

  std::string token;
  while (ifs >> token)
  {
    const std::string::size_type comment_pos = token.find('#');
    if (comment_pos == 0)
    {
      std::string ignored;
      std::getline(ifs, ignored);
      continue;
    }

    std::string clean = comment_pos == std::string::npos ? token : token.substr(0, comment_pos);
    clean = trim(clean);
    if (clean.empty())
      continue;

    if (clean.front() == '+')
      clean.erase(clean.begin());

    const std::string::size_type eq_pos = clean.find('=');
    if (eq_pos == std::string::npos)
    {
      proj4.flags.insert(toLower(clean));
      continue;
    }

    std::string key = toLower(trim(clean.substr(0, eq_pos)));
    std::string value = trim(clean.substr(eq_pos + 1));
    if (!value.empty() && value.front() == '"' && value.back() == '"' && value.size() >= 2)
      value = value.substr(1, value.size() - 2);
    if (!key.empty())
      proj4.values[key] = value;
  }

  return !proj4.values.empty() || !proj4.flags.empty();
}

int parseProjCoordTrans(const std::string &proj_name)
{
  const std::string p = toLower(proj_name);
  if (p == "utm" || p == "tmerc")
    return CT_TransverseMercator;
  if (p == "ortho")
    return CT_Orthographic;

  int projcoord = KvUserDefined;
  if (tryParseInt(p, projcoord))
    return projcoord;
  return KvUserDefined;
}

int parseDatumCode(const std::string &datum)
{
  const std::string d = toLower(datum);
  if (d == "wgs84")
    return Datum_WGS84;
  if (d == "gda2020")
    return 1168;
  if (d == "nad83")
    return 6269;
  if (d == "nad27")
    return 6267;

  int code = KvUserDefined;
  if (tryParseInt(d, code))
    return code;
  return KvUserDefined;
}

int parseEllpsCode(const std::string &ellps)
{
  const std::string e = toLower(ellps);
  if (e == "wgs84")
    return 7030;
  if (e == "grs80")
    return 7019;
  if (e == "clrk66")
    return 7008;

  int code = KvUserDefined;
  if (tryParseInt(e, code))
    return code;
  return KvUserDefined;
}

int parseGeographicTypeCode(const std::string &datum)
{
  const std::string d = toLower(datum);
  if (d == "wgs84")
    return GCS_WGS_84;
  if (d == "gda2020")
    return 7844;
  if (d == "nad83")
    return 4269;
  if (d == "nad27")
    return 4267;
  return KvUserDefined;
}

std::string defaultEllpsForDatum(const std::string &datum)
{
  const std::string d = toLower(datum);
  if (d == "wgs84")
    return "wgs84";
  if (d == "gda2020" || d == "nad83")
    return "grs80";
  if (d == "nad27")
    return "clrk66";
  return datum;
}

}  // namespace

// save to geotif format using floating-point per-channel colour data. This function passes a projection file in order
// to geolocate the image
bool writeGeoTiffFloat(const std::string &filename, int x, int y, const float *data, double pixel_width, bool scalar,
                       const std::string &projection_file, double origin_x, double origin_y)
{
  /* Open TIFF descriptor to write GeoTIFF tags */
  TIFF *tif = XTIFFOpen(filename.c_str(), "w");
  if (!tif)
  {
    return false;
  }

  /* Open GTIF Key parser */
  GTIF *gtif = GTIFNew(tif);
  if (!gtif)
  {
    XTIFFClose(tif);
    return false;
  }

  auto fail = [&](const std::string &msg) {
    std::cerr << msg << std::endl;
    GTIFFree(gtif);
    XTIFFClose(tif);
    return false;
  };

  if (x < 0 || y < 0)
  {
    return fail("Bad image size: " + std::to_string(x) + ", " + std::to_string(y));
  }
  const uint32_t w = (uint32_t)x;
  const uint32_t h = (uint32_t)y;
  const int channels = scalar ? 2 : 4;

  /* Set up standard TIFF file */
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
  /* set other TIFF tags and write out image ... */
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, scalar ? PHOTOMETRIC_MINISBLACK : PHOTOMETRIC_RGB);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
  TIFFSetField(tif, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, channels);
  const uint16_t ex_samp[] = { EXTRASAMPLE_ASSOCALPHA };
  TIFFSetField(tif, TIFFTAG_EXTRASAMPLES, 1, &ex_samp);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

  // now go line by line to write out the image data
  for (uint32_t row = 0; row < h; row++)
  {
    std::vector<float> pdst(channels * w);

    // moving the data from the dib to a row structure that
    // can be used by the tiff library
    for (uint32_t col = 0; col < w; col++)
    {
      const uint32_t index = 3 * ((h - 1 - row) * w + col);
      const float shade = (data[index + 0] + data[index + 1] + data[index + 2]) / 3.0f;
      if (scalar)
      {
        pdst[2 * col] = shade;
        pdst[2 * col + 1] = shade == 0.0 ? 0.0 : 255.0;
      }
      else
      {
        pdst[4 * col + 0] = data[index + 0];
        pdst[4 * col + 1] = data[index + 1];
        pdst[4 * col + 2] = data[index + 2];
        pdst[4 * col + 3] = shade == 0.0 ? 0.0 : 255.0;
      }
    }

    // now actually write the row data
    if (TIFFWriteScanline(tif, &pdst[0], row, 0) != 1)
    {
      return fail("Failed writing TIFF scanline " + std::to_string(row));
    }
  }

  // read in the projection parameters
  if (!projection_file.empty())
  {
    Proj4Data proj4;
    if (!parseProj4File(projection_file, proj4))
      return fail("Cannot parse projection file: " + projection_file);

    const auto getValue = [&](const std::string &k) -> std::string {
      const auto it = proj4.values.find(k);
      return it == proj4.values.end() ? std::string() : it->second;
    };

    const std::string proj_value = getValue("proj");
    const std::string datum_value = getValue("datum");
    std::string ellps_value = getValue("ellps");
    if (ellps_value.empty())
      ellps_value = defaultEllpsForDatum(datum_value);
    const std::string units_value = getValue("units");

    bool is_south = proj4.flags.count("south") > 0;
    int zone_value = KvUserDefined;
    if (!getValue("zone").empty())
    {
      if (!tryParseInt(getValue("zone"), zone_value))
        return fail("Invalid +zone in projection file: " + getValue("zone"));
    }

    double coord_lat = 0.0;
    if (!getValue("lat_0").empty() && !tryParseDouble(getValue("lat_0"), coord_lat))
      return fail("Invalid +lat_0 in projection file: " + getValue("lat_0"));

    double coord_long = 0.0;
    if (!getValue("lon_0").empty() && !tryParseDouble(getValue("lon_0"), coord_long))
      return fail("Invalid +lon_0 in projection file: " + getValue("lon_0"));

    const std::string x0_value = getValue("x_0");
    const std::string y0_value = getValue("y_0");
    const bool has_x0 = !x0_value.empty();
    const bool has_y0 = !y0_value.empty();
    double false_easting = 0.0;
    double false_northing = 0.0;
    if (has_x0 && !tryParseDouble(x0_value, false_easting))
      return fail("Invalid +x_0 in projection file: " + x0_value);
    if (has_y0 && !tryParseDouble(y0_value, false_northing))
      return fail("Invalid +y_0 in projection file: " + y0_value);

    int projected_cs_type = KvUserDefined;
    const std::string init_value = toLower(getValue("init"));
    const std::string epsg_value = getValue("epsg");
    if (!epsg_value.empty())
    {
      if (!tryParseInt(epsg_value, projected_cs_type))
        return fail("Invalid +epsg in projection file: " + epsg_value);
    }
    else if (!init_value.empty())
    {
      const std::string prefix = "epsg:";
      if (init_value.find(prefix) == 0)
      {
        if (!tryParseInt(init_value.substr(prefix.size()), projected_cs_type))
          return fail("Invalid +init EPSG code in projection file: " + init_value);
      }
    }

    if (projected_cs_type == KvUserDefined && toLower(proj_value) == "utm" && zone_value >= 1 && zone_value <= 60)
    {
      if (toLower(datum_value) == "wgs84")
        projected_cs_type = is_south ? (32700 + zone_value) : (32600 + zone_value);
      else if (toLower(datum_value) == "gda2020" && is_south)
        projected_cs_type = 7800 + zone_value;  // MGA2020 zones, e.g. zone 55 -> EPSG:7855
    }

    std::cout << "proj: " << proj_value << ", false_easting: " << false_easting
              << ", false_northing: " << false_northing << ", using_south_flag: " << (is_south ? "true" : "false")
              << ", ellps: " << ellps_value << ", datum: " << datum_value << ", coord_long: " << coord_long
              << ", coord_lat: " << coord_lat << " zone: " << zone_value << std::endl;

    const double scales[3] = { pixel_width, pixel_width, pixel_width };
    TIFFSetField(tif, TIFFTAG_GEOPIXELSCALE, 3, scales);  // set the width of a pixel

    // Set GeoTIFF information
    // We are only supporting a limited set of projection types, so here we assume standard settings
    GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelTypeProjected);
    GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);

    GTIFKeySet(gtif, ProjLinearUnitsGeoKey, TYPE_SHORT, 1, Linear_Meter);
    GTIFKeySet(gtif, VerticalUnitsGeoKey, TYPE_SHORT, 1, Linear_Meter);

    GTIFKeySet(gtif, ProjectionGeoKey, TYPE_SHORT, 1, KvUserDefined);
    if (!proj_value.empty())
    {
      const int projcoord = parseProjCoordTrans(proj_value);
      if (projcoord != KvUserDefined)
        GTIFKeySet(gtif, ProjCoordTransGeoKey, TYPE_SHORT, 1, projcoord);
      else
        std::cout << "warning: unknown projection type: " << proj_value << std::endl;
    }
    if (has_x0)
      GTIFKeySet(gtif, ProjFalseEastingGeoKey, TYPE_DOUBLE, 1, false_easting);

    if (has_y0)
      GTIFKeySet(gtif, ProjFalseNorthingGeoKey, TYPE_DOUBLE, 1, false_northing);
    else if (is_south)
    {
      GTIFKeySet(gtif, ProjFalseNorthingGeoKey, TYPE_DOUBLE, 1, 10000000.0);  // false northing
      // Usually paired with the Latitude of Natural Origin (Equator)
      GTIFKeySet(gtif, ProjNatOriginLatGeoKey, TYPE_DOUBLE, 1, 0.0);
    }

    // describe the coordinates of the image corners
    const double tiepoints[6] = { 0, 0, 0, origin_x, origin_y, 0 };
    TIFFSetField(tif, TIFFTAG_GEOTIEPOINTS, 6, tiepoints);

    const int geographic_type = parseGeographicTypeCode(datum_value);
    if (geographic_type != KvUserDefined)
      GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, geographic_type);

    const int datum_code = parseDatumCode(datum_value);
    if (datum_code != KvUserDefined)
      GTIFKeySet(gtif, GeogGeodeticDatumGeoKey, TYPE_SHORT, 1, datum_code);
    else if (!datum_value.empty())
      std::cout << "warning: unknown datum: " << datum_value << std::endl;

    const int ellps_code = parseEllpsCode(ellps_value);
    if (ellps_code != KvUserDefined)
      GTIFKeySet(gtif, GeogEllipsoidGeoKey, TYPE_SHORT, 1, ellps_code);
    else if (!ellps_value.empty())
      std::cout << "warning: unknown ellipsoid: " << ellps_value << std::endl;

    GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, projected_cs_type);
    GTIFKeySet(gtif, ProjectionGeoKey, TYPE_SHORT, 1, KvUserDefined);

    const std::string unit_norm = toLower(units_value);
    if (!units_value.empty() && unit_norm != "m" && unit_norm != "meter" && unit_norm != "metre")
    {
      return fail("unknown unit type: " + units_value);
    }
    GTIFKeySet(gtif, ProjCenterLongGeoKey, TYPE_DOUBLE, 1, coord_long);
    GTIFKeySet(gtif, ProjCenterLatGeoKey, TYPE_DOUBLE, 1, coord_lat);

    // Store the keys into the TIFF Tags
    GTIFWriteKeys(gtif);
  }

  // get rid of the key parser
  GTIFFree(gtif);

  // save and close the TIFF file descriptor
  XTIFFClose(tif);

  return true;
}
#endif

namespace
{
/// Writes a value to a .csv file, as nan when it is NaN (whatever its sign)
std::string csvValue(double value)
{
  if (std::isnan(value))
  {
    return "nan";
  }
  std::stringstream ss;
  ss << value;
  return ss.str();
}

std::string densityLabel(const DensityGrid &grid)
{
  return grid.distribution() == ElementDistribution::VerticalCylinder ? "frontal_area_per_volume_Nd" :
                                                                        "leaf_area_per_volume";
}

/// Writes the density of each cell of size @c cell_size (tiling from @c min_bound) to a .csv file. This aggregates the
/// voxels (whose centres lie in each cell) into cells at a coarser resolution, such as for hydrodynamic models, with
/// aggregateDensityCells. The density is pooled over the cell's observed voxels (total corrected hits over total
/// G-weighted path length). With the beer_lambert estimator, this is identical to estimating with a single voxel the
/// size of the cell. The mean of the voxel densities and the patchiness ratio (mean / pooled) are also output, and the
/// fraction of the cell's voxels that were observed. Values that are undefined for a cell (e.g. without observed
/// voxels) are nan. @c dims are the (shifted) voxels to aggregate, and @c shift that of estimateDensityGrid
bool writeDensityCells(const std::string &file_name, const DensityGrid &grid, const Eigen::Vector3i &dims, int shift,
                       const Eigen::Vector3d &min_bound, const Eigen::Vector3d &cell_size)
{
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing" << std::endl;
    return false;
  }
  CellLayout layout;
  layout.min_bound = min_bound;
  layout.cell_size = cell_size;
  for (int i = 0; i < 3; i++)
  {
    layout.num_cells[i] =
      std::max(1, static_cast<int>(std::ceil(static_cast<double>(dims[i]) * grid.voxelWidth() / cell_size[i] - 1e-6)));
  }
  const std::vector<DensityCell> cells = aggregateDensityCells(grid, shift, Eigen::Vector3i(0, 0, 0), dims, layout);
  const std::string label = densityLabel(grid);
  ofs << "x_min,y_min,z_min,x_max,y_max,z_max,observed_fraction,hits,weighted_path,density_" << label
      << ",pooled_density_" << label << ",patchiness_ratio" << std::endl;
  for (int z = 0; z < layout.num_cells[2]; z++)
  {
    for (int y = 0; y < layout.num_cells[1]; y++)
    {
      for (int x = 0; x < layout.num_cells[0]; x++)
      {
        const DensityCell &cell = cells[layout.index(x, y, z)];
        const bool observed = cell.observed > 0;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const Eigen::Vector3d cell_min =
          min_bound + Eigen::Vector3d(x * cell_size[0], y * cell_size[1], z * cell_size[2]);
        const Eigen::Vector3d cell_max = cell_min + cell_size;
        ofs << cell_min[0] << "," << cell_min[1] << "," << cell_min[2] << "," << cell_max[0] << "," << cell_max[1]
            << "," << cell_max[2] << "," << csvValue(cell.observedFraction()) << ","
            << csvValue(observed ? cell.hits : nan) << "," << csvValue(observed ? cell.weighted_path : nan) << ","
            << csvValue(cell.density()) << "," << csvValue(cell.pooledDensity()) << ","
            << csvValue(cell.patchinessRatio()) << std::endl;
      }
    }
  }
  std::cout << "density cells written to " << file_name << std::endl;
  return true;
}

/// Writes the density of each horizontal layer of voxels to a .csv file: pooled across the layer's observed voxels
/// (total corrected hits over total G-weighted path length), and the mean of the voxel densities, as for
/// writeDensityCells. The mean of voxel densities is biased high when voxels contain few hits (around +10% at 2-3 hits
/// per voxel), so the mean hits per voxel is also output. The pooled estimate avoids this bias, but weights voxels by
/// path length, so it is biased low where the density varies within the layer (see the patchiness ratio). Only
/// voxels with sufficient rays are included, so unobserved space (e.g. below the ground) is excluded.
/// @c dims are the (shifted) voxels to aggregate, where voxel z spans heights min_bound[2] + voxel width * [z, z+1]
bool writeDensityProfile(const std::string &file_name, const DensityGrid &grid, const Eigen::Vector3i &dims, int shift,
                         const Eigen::Vector3d &min_bound)
{
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing" << std::endl;
    return false;
  }
  const double w = grid.voxelWidth();
  CellLayout layout;
  layout.min_bound = min_bound;
  layout.cell_size = Eigen::Vector3d(dims[0] * w, dims[1] * w, w);
  layout.num_cells = Eigen::Vector3i(1, 1, dims[2]);
  const std::vector<DensityCell> layers = aggregateDensityCells(grid, shift, Eigen::Vector3i(0, 0, 0), dims, layout);
  const std::string label = densityLabel(grid);
  ofs << "z_min,z_max,observed_voxels,total_voxels,mean_hits_per_voxel,density_" << label << ",pooled_density_" << label
      << ",patchiness_ratio" << std::endl;
  for (int z = 0; z < dims[2]; z++)
  {
    const DensityCell &layer = layers[static_cast<size_t>(z)];
    const double z_min = min_bound[2] + w * static_cast<double>(z);
    ofs << z_min << "," << z_min + w << "," << layer.observed << "," << layer.total << ","
        << csvValue(layer.meanHitsPerVoxel()) << "," << csvValue(layer.density()) << ","
        << csvValue(layer.pooledDensity()) << "," << csvValue(layer.patchinessRatio()) << std::endl;
  }
  std::cout << "density profile written to " << file_name << std::endl;
  return true;
}
}  // namespace

bool renderCloud(const std::string &cloud_file, const Cuboid &bounds, ViewDirection view_direction, RenderStyle style,
                 double pix_width, const std::string &image_file, const std::string &projection_file, bool mark_origin,
                 const std::string *const transform_file, const DensityOptions &density_options,
                 const DensityOutputs &density_outputs)
{
  // convert the view direction into useable parameters
  int axis = 0;
  if (view_direction == ViewDirection::Top)
    axis = 2;
  else if (view_direction == ViewDirection::Front || view_direction == ViewDirection::Back)
    axis = 1;
  double dir = 1;
  if (view_direction == ViewDirection::Left || view_direction == ViewDirection::Front)
    dir = -1;
  const bool flip_x = view_direction == ViewDirection::Left || view_direction == ViewDirection::Back;

  // pull out the main image axes (ax1,ax2 are the horiz,vertical axes)
  const Eigen::Vector3d extent = bounds.max_bound_ - bounds.min_bound_;
  // for each view axis (side,top,front = 0,1,2) we need to have an image x axis, and y axis.
  // e.g. x_axes[axis] is the 3D axis to use (x,y,z = 0,1,2) for the image horizontal direction
  const std::array<int, 3> x_axes = { 1, 0, 0 };
  const std::array<int, 3> y_axes = { 2, 2, 1 };
  const int ax1 = x_axes[axis];
  const int ax2 = y_axes[axis];
  const int width = 1 + static_cast<int>(extent[ax1] / pix_width);
  const int height = 1 + static_cast<int>(extent[ax2] / pix_width);
  const int depth = 1 + static_cast<int>(extent[axis] / pix_width);
  std::cout << "outputting " << width << "x" << height << " image" << std::endl;

  try  // there is a possibility of running out of memory here. So provide a helpful message rather than just asserting
  {
    // accumulated colour buffer
    std::vector<Eigen::Vector4d> pixels(width * height, Eigen::Vector4d(0, 0, 0, 0));
    // density calculation is a special case
    if (style == RenderStyle::Density || style == RenderStyle::Density_rgb)
    {
      const Eigen::Vector3i valid_grid_dims = (extent / pix_width).cast<int>() + Eigen::Vector3i(1, 1, 1);
      int shift = 0;
      std::unique_ptr<DensityGrid> grid_ptr =
        estimateDensityGrid(cloud_file, bounds, pix_width, valid_grid_dims, density_options, shift);
      if (!grid_ptr)
      {
        return false;
      }
      DensityGrid &grid = *grid_ptr;
      const bool cylinders = density_options.distribution == ElementDistribution::VerticalCylinder;
      // addNeighbourPriors shifts the voxel contents by -shift in each axis, so index 0 is at bounds.min_bound_. The
      // image has an extra voxel in each axis (beyond the bounds when the extent is a multiple of the voxel width),
      // so the profile and cells only use the voxels that the bounds overlap
      Eigen::Vector3i output_dims;
      for (int i = 0; i < 3; i++)
      {
        output_dims[i] =
          std::min(valid_grid_dims[i], std::max(1, static_cast<int>(std::ceil(extent[i] / pix_width - 1e-6))));
      }
      if (!density_outputs.profile_file.empty())
      {
        if (!writeDensityProfile(density_outputs.profile_file, grid, output_dims, shift, bounds.min_bound_))
        {
          return false;
        }
      }
      if (!density_outputs.cells_file.empty())
      {
        if (!writeDensityCells(density_outputs.cells_file, grid, output_dims, shift, bounds.min_bound_,
                               density_outputs.cell_size))
        {
          return false;
        }
      }
      double foliage_area = 0.0;
      for (int x = 0; x < width; x++)
      {
        for (int y = 0; y < height; y++)
        {
          double p =
            grid.peaks_[(x + shift) + grid.dimensions()[0] * (y + shift)];  // the location before it was shifted
          int top_z = p < -10000.0 ? -1 : (int)std::floor(p);

          double total_density = 0.0;
          for (int z = 0; z < depth; z++)
          {
            Eigen::Vector3i ind;
            ind[axis] = z;
            ind[ax1] = x;
            ind[ax2] = y;
            double d = grid.voxelDensity(grid.voxels()[grid.getIndex(ind)]);
#if defined FLAT_TOP_COMPENSATION
            if (z == top_z - shift)
            {
              d *= p - (double)top_z;
            }
#endif
            total_density += d * pix_width;
          }
          foliage_area += total_density * pix_width * pix_width;
          pixels[x + width * y] = Eigen::Vector4d(total_density, total_density, total_density, total_density);
        }
      }
      if (cylinders)
      {
        std::cout << "total frontal area (sum of N d h): " << foliage_area << " m^2" << std::endl;
      }
      else
      {
        std::cout << "total foliage area: " << foliage_area << " m^2" << std::endl;
      }
    }
    else  // otherwise we use a common algorithm, specialising on render style only per-ray
    {
      // this lambda expression lets us chunk load the ray cloud file, so we don't run out of RAM
      auto render = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                        std::vector<RGBA> &colours) {
        for (size_t i = 0; i < ends.size(); i++)
        {
          const RGBA &colour = colours[i];
          if (colour.alpha == 0)
            continue;
          const Eigen::Vector3d col = Eigen::Vector3d(colour.red, colour.green, colour.blue) / 255.0;
          const Eigen::Vector3d point = style == RenderStyle::Starts ? starts[i] : ends[i];
          const Eigen::Vector3d pos = (point - bounds.min_bound_) / pix_width;
          const Eigen::Vector3i p = (pos).cast<int>();
          const int x = p[ax1], y = p[ax2];
          // using 4 dimensions helps us to accumulate colours in a greater variety of ways
          Eigen::Vector4d &pix = pixels[x + width * y];
          switch (style)  // render the image according to the chosen style
          {
          case RenderStyle::Ends:
          case RenderStyle::Starts:
          case RenderStyle::Height:
            if (pos[axis] * dir > pix[3] * dir || pix[3] == 0.0)  // using 0.0 precisely as a flag here
            {
              pix = Eigen::Vector4d(col[0], col[1], col[2], pos[axis]);
            }
            break;
          case RenderStyle::Mean:
            pix += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
            break;
          case RenderStyle::Sum:
            pix += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
            break;
          case RenderStyle::Rays: {
            Eigen::Vector3d cloud_start = starts[i];
            Eigen::Vector3d cloud_end = ends[i];
            // clip to within the image (since we exclude unbounded rays from the image bounds)
            if (!bounds.clipRay(cloud_start, cloud_end))
            {
              continue;
            }
            Eigen::Vector3d start = (cloud_start - bounds.min_bound_) / pix_width;
            Eigen::Vector3d end = (cloud_end - bounds.min_bound_) / pix_width;
            const Eigen::Vector3d dir = cloud_end - cloud_start;

            // fast approximate 2D line rendering requires picking the long axis to iterate along
            const bool x_long = std::abs(dir[ax1]) > std::abs(dir[ax2]);
            const int axis_long = x_long ? ax1 : ax2;
            const int axis_short = x_long ? ax2 : ax1;
            const int width_long = x_long ? 1 : width;
            const int width_short = x_long ? width : 1;

            const double gradient = dir[axis_short] / dir[axis_long];
            if (dir[axis_long] < 0.0)
              std::swap(start, end);  // this lets us iterate from low up to high values
            const int start_long = static_cast<int>(start[axis_long]);
            const int end_long = static_cast<int>(end[axis_long]);
            // place a pixel at the height of each midpoint (of the pixel) in the long axis
            const double start_mid_point = 0.5 + static_cast<double>(start_long);
            double height = start[axis_short] + (start_mid_point - start[axis_long]) * gradient;
            for (int l = start_long; l <= end_long; l++, height += gradient)
            {
              const int s = static_cast<int>(height);
              pixels[width_long * l + width_short * s] += Eigen::Vector4d(col[0], col[1], col[2], 1.0);
            }
            break;
          }
          default:
            break;
          }
        }
      };
      if (!Cloud::read(cloud_file, render))
        return false;
    }

    double max_val = 1.0;
    double min_val = 0.0;
    const std::string image_ext = getFileNameExtension(image_file);
    const bool is_hdr = image_ext == "hdr" || image_ext == "tif";
    if (!is_hdr)  // limited range, so work out a sensible maximum value, I'm using mean + two standard deviations:
    {
      double sum = 0.0;
      double num = 0.0;
      for (auto &pixel : pixels)
      {
        sum += pixel[3];
        if (pixel[3] > 0.0)
          num++;
      }
      double mean = sum / num;
      double sum_sqr = 0.0;
      for (auto &pixel : pixels)
      {
        if (pixel[3] > 0.0)
          sum_sqr += sqr(pixel[3] - mean);
      }
      const double standard_deviation = std::sqrt(sum_sqr / num);
      max_val = mean + 2.0 * standard_deviation;
      min_val = mean - 2.0 * standard_deviation;
    }
    std::cout << "max val: " << max_val << std::endl;

    // The final pixel buffer
    std::vector<RGBA> pixel_colours;
    std::vector<float> float_pixel_colours;
    if (is_hdr)
      float_pixel_colours.resize(3 * width * height);
    else
      pixel_colours.resize(width * height);

    for (int x = 0; x < width; x++)
    {
      const int indx = flip_x ? width - 1 - x : x;  // possible horizontal flip, depending on view direction
      for (int y = 0; y < height; y++)
      {
        const Eigen::Vector4d colour = pixels[x + width * y];
        Eigen::Vector3d col3d(colour[0], colour[1], colour[2]);
        const uint8_t alpha = colour[3] == 0.0 ? 0 : 255;  // 'punch-through' alpha
        switch (style)  // convert to the colour data structure based on the chosen style
        {
        case RenderStyle::Mean:
        case RenderStyle::Rays:
          col3d /= colour[3];  // simple mean
          break;
        case RenderStyle::Height: {
          double shade =
            dir == 1.0 ? (colour[3] - min_val) / (max_val - min_val) : (colour[3] - max_val) / (min_val - max_val);
          col3d = Eigen::Vector3d(shade, shade, shade);
          break;
        }
        case RenderStyle::Sum:
        case RenderStyle::Density:
          if (!is_hdr)
            col3d /= max_val;  // rescale to within limited colour range
          break;
        case RenderStyle::Density_rgb: {
          if (is_hdr)
          {
            // brightness doubles and hue cycles every power of 10
            col3d = std::pow(colour[0], 0.3) * redGreenBlueSpectrum(std::log10(std::max(1e-6, colour[0])));
          }
          else
          {
            col3d = 0.7 * std::pow(colour[0], 0.3) * redGreenBlueSpectrum(std::log10(std::max(1e-6, colour[0])));

            // double shade = colour[0] / max_val;
            // col3d = redGreenBlueGradient(shade);
            // if (shade < 0.05)
            //   col3d *= 20.0 * shade;  // this blends the lowest densities down to black
          }
          break;
        }
        default:
          break;
        }
        const int ind = indx + width * y;
        if (is_hdr)
        {
          float_pixel_colours[3 * ind + 0] = (float)col3d[0];
          float_pixel_colours[3 * ind + 1] = (float)col3d[1];
          float_pixel_colours[3 * ind + 2] = (float)col3d[2];
        }
        else
        {
          RGBA col;
          col.red = uint8_t(std::max(0.0, std::min(255.0 * col3d[0], 255.0)));
          col.green = uint8_t(std::max(0.0, std::min(255.0 * col3d[1], 255.0)));
          col.blue = uint8_t(std::max(0.0, std::min(255.0 * col3d[2], 255.0)));
          col.alpha = alpha;
          pixel_colours[ind] = col;
        }
      }
    }
    if (mark_origin)  // an option to mark the lidar origin in the image
    {
      if (pixel_colours.empty())
      {
        std::cout << "warning: mark origin not implemented for hdr images" << std::endl;
      }
      else
      {
        const Eigen::Vector3d pos = -bounds.min_bound_ / pix_width;
        std::cout << "origin pixel location: " << pos[0] << ", " << pos[1] << std::endl;
        const Eigen::Vector3i p = pos.cast<int>();
        const int x = p[ax1], y = p[ax2];
#define DRAW_ARROW  // render the origin as an arrow, which therefore defines the x direction in the lidar frame
#if defined DRAW_ARROW
        for (int xx = x - 2; xx <= x + 10; xx++)
        {
          std::cout << "xx: " << xx << std::endl;
          for (int yy = y - 2; yy <= y + 2; yy++)
          {
            if (xx >= 6 && std::abs(yy - y) > ((x + 10) - xx) / 2)
              continue;
            if (xx >= 0 && xx < width && yy >= 0 && yy < height)
            {
              const int indx = flip_x ? width - 1 - xx : xx;
              RGBA &col = pixel_colours[indx + width * yy];
              col.red = col.green = 0;
              col.blue = col.alpha = 255;
            }
          }
        }
        std::cout << "done" << std::endl;
#endif
        if (x >= 0 && x < width && y >= 0 && y < height)
        {
          const int indx = flip_x ? width - 1 - x : x;  // possible horizontal flip, depending on view direction
          // using 4 dimensions helps us to accumulate colours in a greater variety of ways
          RGBA &col = pixel_colours[indx + width * y];
          col.red = col.blue = col.alpha = 255;
          col.green = 0;
          // we leave alpha alone as it might be needed to indicate the presence of points
        }
        else
        {
          std::cerr << "error: the origin cannot be marked on this image as it is not within the image bounds"
                    << std::endl;
        }
      }
    }
    // option to output the transformation of the image
    if (transform_file != nullptr)
    {
      // Compute transform
      const double scale = pix_width;
      const double translate_x = bounds.min_bound_.x();
      const double translate_y = bounds.max_bound_.y();
      const Eigen::Matrix3d transform =
        (Eigen::Translation2d(translate_x, translate_y) * Eigen::Scaling(scale, -scale)).matrix();

      // Write transform
      std::cout << "outputting transform: " << *transform_file << std::endl;
      std::ofstream ofs;
      ofs.open(*transform_file, std::ios::out);
      if (ofs.fail())
      {
        std::cerr << "Error: cannot open " << *transform_file << " for writing." << std::endl;
        return false;
      }
      ofs << "# Generated by rayrender." << std::endl;
      ofs << "# For a given pixel:" << std::endl;
      ofs << "#   P_pixel = [x_pixel, y_pixel, 1]" << std::endl;
      ofs << "# The position of the centre of the pixel in the ray cloud's frame can be" << std::endl;
      ofs << "# computed as:" << std::endl;
      ofs << "#   P_raycloud = [x_raycloud, y_raycloud, _]" << std::endl;
      ofs << "#   P_raycloud = T * P_pixel" << std::endl;
      ofs << "# Where T is the 3*3 transformation matrix defined in this file:" << std::endl;
      ofs << "#   T = [" << std::endl;
      ofs << "#     transform[0], transform[1], transform[2];" << std::endl;
      ofs << "#     transform[3], transform[4], transform[5];" << std::endl;
      ofs << "#     transform[6], transform[7], transform[8];" << std::endl;
      ofs << "#   ]" << std::endl;
      ofs << "# All z information is lost in rayrender." << std::endl;
      ofs << "transform: [" << std::endl;
      ofs << "  " << transform(0, 0) << ", " << transform(0, 1) << ", " << transform(0, 2) << "," << std::endl;
      ofs << "  " << transform(1, 0) << ", " << transform(1, 1) << ", " << transform(1, 2) << "," << std::endl;
      ofs << "  " << transform(2, 0) << ", " << transform(2, 1) << ", " << transform(2, 2) << "," << std::endl;
      ofs << "]" << std::endl;
      ofs.close();
    }
    std::cout << "outputting image: " << image_file << std::endl;

    // write the image depending on the file format
    const char *image_name = image_file.c_str();
    stbi_flip_vertically_on_write(1);
    if (image_ext == "png")
      stbi_write_png(image_name, width, height, 4, (void *)&pixel_colours[0], 4 * width);
    else if (image_ext == "bmp")
      stbi_write_bmp(image_name, width, height, 4, (void *)&pixel_colours[0]);
    else if (image_ext == "tga")
      stbi_write_tga(image_name, width, height, 4, (void *)&pixel_colours[0]);
    else if (image_ext == "jpg")
      stbi_write_jpg(image_name, width, height, 4, (void *)&pixel_colours[0], 100);  // 100 is maximal quality
    else if (image_ext == "hdr")
      stbi_write_hdr(image_name, width, height, 3, &float_pixel_colours[0]);
#if RAYLIB_WITH_TIFF
    else if (image_ext == "tif")
    {
      // obtain the origin offsets
      const Eigen::Vector3d origin(0, 0, 0);
      const Eigen::Vector3d pos = -(origin - bounds.min_bound_);
      const double x = pos[ax1], y = pos[ax2] + static_cast<double>(height - 1) * pix_width;

      std::cout << std::setprecision(3) << std::fixed << "x: " << x << ", y: " << y << ", absy: " << pos[ax2]
                << std::endl;
      // generate the geotiff file
      writeGeoTiffFloat(image_file, width, height, &float_pixel_colours[0], pix_width, false, projection_file, x, y);
    }
#endif
    else
    {
      std::cerr << "Error: image format " << image_ext << " not supported" << std::endl;
      return false;
    }
  }
  catch (std::bad_alloc const &)  // catch any memory allocation problems in generating large images
  {
    std::cout << "Not enough memory to process the " << width << "x" << height << " image." << std::endl;
    std::cout << "The --pixel_width option can be used to reduce the resolution." << std::endl;
  }
#if !RAYLIB_WITH_TIFF
  RAYLIB_UNUSED(projection_file);
#endif
  return true;
}

}  // namespace ray
