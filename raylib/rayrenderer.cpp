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
#include <cstdlib>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include "rayunused.h"

#define DENSITY_MIN_RAYS 10  // larger is more accurate but more blurred. 0 for no adaptive blending

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

double stripInterceptFraction(double width, double beam_diameter, double offset)
{
  if (beam_diameter <= 0.0)
  {
    return std::abs(offset) < 0.5 * width ? 1.0 : 0.0;
  }
  // the 1/e^2 diameter is 4 standard deviations. The power across the strip is the 1D marginal of the 2D Gaussian
  const double sigma = 0.25 * beam_diameter;
  const double scale = 1.0 / (std::sqrt(2.0) * sigma);
  return 0.5 * (std::erf((offset + 0.5 * width) * scale) - std::erf((offset - 0.5 * width) * scale));
}

double effectiveStemWidth(double diameter, double beam_diameter, double threshold)
{
  if (beam_diameter <= 0.0)
  {
    return threshold <= 1.0 ? diameter : 0.0;
  }
  if (stripInterceptFraction(diameter, beam_diameter, 0.0) < threshold)
  {
    return 0.0;  // never intercepts enough of the beam to be detected
  }
  // the intercepted fraction decreases monotonically with offset, so bisect for where it equals the threshold
  double low = 0.0;
  double high = 0.5 * diameter + 4.0 * beam_diameter;
  for (int i = 0; i < 60; i++)
  {
    const double mid = 0.5 * (low + high);
    if (stripInterceptFraction(diameter, beam_diameter, mid) >= threshold)
    {
      low = mid;
    }
    else
    {
      high = mid;
    }
  }
  return 2.0 * low;
}

StemWidthModel::StemWidthModel(double diameter, const BeamFootprint &footprint)
  : diameter_(diameter)
  , footprint_(footprint)
  , enabled_(diameter > 0.0 && footprint.enabled())
{
  if (!enabled_)
  {
    return;
  }
  const int num = static_cast<int>(table_max_range_ / table_step_) + 1;
  table_.resize(static_cast<size_t>(num));
  for (int i = 0; i < num; i++)
  {
    table_[static_cast<size_t>(i)] =
      effectiveStemWidth(diameter_, footprint_.diameter(table_step_ * static_cast<double>(i)), footprint_.threshold);
  }
}

double StemWidthModel::width(double range) const
{
  if (!enabled_)
  {
    return diameter_;
  }
  const double x = std::max(0.0, range) / table_step_;
  const size_t i = static_cast<size_t>(x);
  if (i + 1 >= table_.size())
  {
    return effectiveStemWidth(diameter_, footprint_.diameter(range), footprint_.threshold);
  }
  const double blend = x - static_cast<double>(i);
  return table_[i] * (1.0 - blend) + table_[i + 1] * blend;
}

double StemWidthModel::maxWidth(double max_range) const
{
  double max_width = diameter_;
  for (double range = 0.0; range <= max_range; range += table_step_)
  {
    max_width = std::max(max_width, width(range));
  }
  return max_width;
}

void DensityGrid::calculatePeaks(const std::string &file_name)
{
  auto calc_peaks = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                        std::vector<RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); ++i)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      if (!bounds_.clipRay(start, end, 1e-10))
      {
        continue;  // ray is outside of bounds
      }
      if (colours[i].alpha > 0)
      {
        const Eigen::Vector3d vox_end = (end - bounds_.min_bound_) / voxel_width_;
        const Eigen::Vector3i target =
          Eigen::Vector3d(std::floor(vox_end[0]), std::floor(vox_end[1]), std::floor(vox_end[2])).cast<int>();
        int peak_id = target[0] + target[1] * voxel_dims_[0];
        peaks_[peak_id] = std::max(peaks_[peak_id], vox_end[2]);
      }
    }
  };
  Cloud::read(file_name, calc_peaks);
}

/// Calculate the density within each voxel of the grid: by default the surface area per cubic metre, assuming an
/// unbiased distribution of surface angles. See @c ElementDistribution for the alternatives.
void DensityGrid::calculateDensities(const std::string &file_name)
{
  auto calculate = [&](std::vector<Eigen::Vector3d> &starts, std::vector<Eigen::Vector3d> &ends, std::vector<double> &,
                       std::vector<RGBA> &colours) {
    for (size_t i = 0; i < ends.size(); ++i)
    {
      Eigen::Vector3d start = starts[i];
      Eigen::Vector3d end = ends[i];
      if (!bounds_.clipRay(start, end, 1e-10))
      {
        continue;  // ray is outside of bounds
      }
      // a ray whose end point was clipped by the bounds has passed through the grid without hitting anything inside
      // it, so it is treated as a miss. Otherwise its hit would be wrongly counted in the boundary voxel.
      const double clip_tolerance = 1e-6;
      const bool end_clipped = (end - ends[i]).squaredNorm() > clip_tolerance * clip_tolerance;
      range_offset_ = (start - starts[i]).norm();
      bounded_ = colours[i].alpha > 0 && !end_clipped;
      const Eigen::Vector3d vox_start = (start - bounds_.min_bound_) / voxel_width_;
      const Eigen::Vector3d vox_end = (end - bounds_.min_bound_) / voxel_width_;
      source_ = vox_start;
      dir_ = (vox_end - vox_start).normalized();
      ray_g_ = gFunction(distribution_, dir_);
      if (ray_g_ < min_g_ || ray_g_ > max_g_ || ray_g_ <= 0.0)
      {
        num_ignored_rays_++;  // selecting rays on direction alone doesn't bias the estimate
        continue;
      }
      hit_voxel_ = Eigen::Vector3d(std::floor(vox_end[0]), std::floor(vox_end[1]), std::floor(vox_end[2])).cast<int>();
      hit_length_ = (vox_end - vox_start).norm();
      Eigen::Vector3d walk_end = vox_end;
      if (bounded_ && path_estimator_ == PathEstimator::Crossing && stem_diameter_ > 0.0)
      {
        // extend the ray through the stem's mean internal chord. The G-weighted length of this is pi/4 * width
        const double hit_range = range_offset_ + (end - start).norm();
        const double width = stem_width_.enabled() ? stem_width_.width(hit_range) : stem_diameter_;
        const double extension = std::min(1.0, 0.25 * kPi * width / ray_g_);
        Eigen::Vector3d ext_start = start;
        Eigen::Vector3d ext_end = end + (end - start).normalized() * extension;
        if (bounds_.clipRay(ext_start, ext_end, 1e-10))
        {
          walk_end = (ext_end - bounds_.min_bound_) / voxel_width_;
        }
      }
      walkGrid(vox_start, walk_end, *this);
    }
  };
  Cloud::read(file_name, calculate);
}

namespace
{
/// Deterministic uniform random number in (0, 1), independent of the standard library's distributions
inline double uniform01(std::mt19937 &gen)
{
  return (static_cast<double>(gen()) + 0.5) / 4294967296.0;
}
}  // namespace

namespace
{
/// One realisation of the free path simulation for freePathResponse, accumulating hits and free path length
void freePathRealisation(double solid_fraction, double width_ratio, double min_gap_ratio, unsigned int seed,
                         double &hits, double &path, double &density_out)
{
  // non-overlapping stems of unit diameter, placed by random sequential addition in a square
  const int num_stems = 1500;
  const int num_lines = 10000;
  const double diameter = 1.0;
  const double side = std::sqrt(num_stems * kPi * diameter * diameter / (4.0 * solid_fraction));
  density_out = num_stems / (side * side);
  std::mt19937 gen(seed);
  const double separation = diameter * (1.0 + min_gap_ratio);  // minimum distance between stem centres
  const double radius = 0.5 * width_ratio * diameter;          // of the widened stems
  const double cell = std::max(separation, radius);            // so stems within these distances are in adjacent cells
  const int cells = std::max(1, static_cast<int>(std::ceil(side / cell)));
  std::vector<std::vector<Eigen::Vector2d>> buckets(static_cast<size_t>(cells * cells));
  auto bucket = [&](double x) { return std::min(cells - 1, std::max(0, static_cast<int>(x / cell))); };
  // whether any stem centre is closer than dist to pos (dist no more than the cell width)
  auto near = [&](const Eigen::Vector2d &pos, double dist) {
    const int bx = bucket(pos[0]), by = bucket(pos[1]);
    for (int x = std::max(0, bx - 1); x <= std::min(cells - 1, bx + 1); x++)
    {
      for (int y = std::max(0, by - 1); y <= std::min(cells - 1, by + 1); y++)
      {
        for (const auto &c : buckets[static_cast<size_t>(x + cells * y)])
        {
          if ((pos - c).squaredNorm() < dist * dist)
          {
            return true;
          }
        }
      }
    }
    return false;
  };
  int placed = 0;
  for (long long attempt = 0; placed < num_stems && attempt < 100000000LL; attempt++)
  {
    const Eigen::Vector2d pos(side * uniform01(gen), side * uniform01(gen));
    if (!near(pos, separation))
    {
      buckets[static_cast<size_t>(bucket(pos[0]) + cells * bucket(pos[1]))].push_back(pos);
      placed++;
    }
  }
  // lines from random points that are free of the widened stems, until their first hit or the edge of the square
  for (int l = 0; l < num_lines; l++)
  {
    Eigen::Vector2d start;
    do
    {
      start = Eigen::Vector2d(side * (0.1 + 0.8 * uniform01(gen)), side * (0.1 + 0.8 * uniform01(gen)));
    } while (near(start, radius));
    const double angle = 2.0 * kPi * uniform01(gen);
    const Eigen::Vector2d dir(std::cos(angle), std::sin(angle));
    double t_max = std::numeric_limits<double>::max();
    for (int i = 0; i < 2; i++)
    {
      if (dir[i] > 1e-12)
      {
        t_max = std::min(t_max, (side - start[i]) / dir[i]);
      }
      else if (dir[i] < -1e-12)
      {
        t_max = std::min(t_max, -start[i] / dir[i]);
      }
    }
    // walk the cells along the line. Any stem whose widened disc reaches a cell has its centre in that cell or an
    // adjacent one, so once a hit lies within the current cell, no later cell can have an earlier hit
    double best = t_max;
    int cx = bucket(start[0]), cy = bucket(start[1]);
    const int step_x = dir[0] > 0.0 ? 1 : -1, step_y = dir[1] > 0.0 ? 1 : -1;
    double next_x = std::abs(dir[0]) > 1e-12 ? ((cx + (step_x > 0 ? 1 : 0)) * cell - start[0]) / dir[0] :
                                               std::numeric_limits<double>::max();
    double next_y = std::abs(dir[1]) > 1e-12 ? ((cy + (step_y > 0 ? 1 : 0)) * cell - start[1]) / dir[1] :
                                               std::numeric_limits<double>::max();
    const double delta_x = std::abs(dir[0]) > 1e-12 ? cell / std::abs(dir[0]) : std::numeric_limits<double>::max();
    const double delta_y = std::abs(dir[1]) > 1e-12 ? cell / std::abs(dir[1]) : std::numeric_limits<double>::max();
    while (true)
    {
      for (int x = std::max(0, cx - 1); x <= std::min(cells - 1, cx + 1); x++)
      {
        for (int y = std::max(0, cy - 1); y <= std::min(cells - 1, cy + 1); y++)
        {
          for (const auto &c : buckets[static_cast<size_t>(x + cells * y)])
          {
            const Eigen::Vector2d rel = c - start;
            const double b = rel.dot(dir);
            const double disc = b * b - (rel.squaredNorm() - radius * radius);
            if (disc >= 0.0)
            {
              const double t = b - std::sqrt(disc);
              if (t > 0.0 && t < best)
              {
                best = t;
              }
            }
          }
        }
      }
      const double t_exit = std::min(next_x, next_y);
      if (best <= t_exit || t_exit >= t_max)
      {
        break;
      }
      if (next_x < next_y)
      {
        cx += step_x;
        next_x += delta_x;
      }
      else
      {
        cy += step_y;
        next_y += delta_y;
      }
      if (cx < 0 || cy < 0 || cx >= cells || cy >= cells)
      {
        break;
      }
    }
    path += best;
    hits += best < t_max ? 1.0 : 0.0;
  }
}
}  // namespace

double freePathResponse(double solid_fraction, double width_ratio, double min_gap_ratio, unsigned int seed)
{
  if (solid_fraction <= 0.0)
  {
    return 1.0;
  }
  min_gap_ratio = std::max(0.0, min_gap_ratio);
  const double separation_scale = (1.0 + min_gap_ratio) * (1.0 + min_gap_ratio);
  solid_fraction = std::min(solid_fraction, 0.4 / separation_scale);
  width_ratio = std::max(width_ratio, 1.0);
  // pool several independent realisations, as one realisation's response varies by about 2-3% between seeds
  const int num_realisations = 64;
  double hits = 0.0, path = 0.0, density = 0.0;
  for (int i = 0; i < num_realisations; i++)
  {
    freePathRealisation(solid_fraction, width_ratio, min_gap_ratio, seed + 7919u * static_cast<unsigned int>(i), hits,
                        path, density);
  }
  return hits / (path * density * width_ratio);
}

double StemResponseTable::response(double solid_fraction, double width_ratio)
{
  // bilinear interpolation over a grid of solid fractions and width ratios, evaluated on demand
  static const std::vector<double> fractions = { 0.0, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4 };
  const double ratio_step = 0.25, max_ratio = 5.0;
  // the table is over the solid fraction of the exclusion discs (diameter d (1 + gap)), which have the packing limit
  const double separation_scale = (1.0 + min_gap_ratio_) * (1.0 + min_gap_ratio_);
  solid_fraction = std::max(0.0, std::min(solid_fraction * separation_scale, fractions.back()));
  const double r = (std::max(1.0, std::min(width_ratio, max_ratio)) - 1.0) / ratio_step;
  const int r0 = std::min(static_cast<int>(r), static_cast<int>((max_ratio - 1.0) / ratio_step) - 1);
  size_t f0 = 0;
  while (f0 + 2 < fractions.size() && fractions[f0 + 1] < solid_fraction)
  {
    f0++;
  }
  const double fb = (solid_fraction - fractions[f0]) / (fractions[f0 + 1] - fractions[f0]);
  const double rb = r - r0;
  auto value = [&](size_t f, int ri) {
    const auto key = std::make_pair(static_cast<int>(f), ri);
    auto it = table_.find(key);
    if (it == table_.end())
    {
      it = table_.emplace(key, freePathResponse(fractions[f] / separation_scale, 1.0 + ratio_step * ri, min_gap_ratio_))
             .first;
    }
    return it->second;
  };
  return (1.0 - rb) * ((1.0 - fb) * value(f0, r0) + fb * value(f0 + 1, r0)) +
         rb * ((1.0 - fb) * value(f0, r0 + 1) + fb * value(f0 + 1, r0 + 1));
}

void DensityGrid::enableOverlapCorrection()
{
  if (distribution_ != ElementDistribution::VerticalCylinder || stem_diameter_ <= 0.0)
  {
    return;
  }
  overlap_path_.assign(voxels_.size(), 0.0f);
}

void DensityGrid::applyOverlapCorrection()
{
  const bool pooled = path_estimator_ == PathEstimator::PooledBeerLambert;
  if (overlap_path_.empty() && !pooled)
  {
    return;
  }
  // neighbourhood of about +-0.1 m horizontally, to estimate the local density with less noise than a single voxel
  const int radius = std::max(0, static_cast<int>(std::round(pool_radius_ / voxel_width_)));
  std::vector<float> new_hits(voxels_.size());
  for (size_t i = 0; i < voxels_.size(); i++)
  {
    new_hits[i] = voxels_[i].numHits();
  }
  size_t num_corrected = 0, num_saturated = 0;
  double sum_factor = 0.0, sum_k_path = 0.0, sum_path = 0.0;
  for (int z = 0; z < voxel_dims_[2]; z++)
  {
    for (int y = 0; y < voxel_dims_[1]; y++)
    {
      for (int x = 0; x < voxel_dims_[0]; x++)
      {
        const int ind = getIndex(Eigen::Vector3i(x, y, z));
        if (pooled ? voxels_[ind].numRays() <= 0.0f : voxels_[ind].numHits() <= 0.0f)
        {
          continue;
        }
        double hits = 0.0, path = 0.0, k_path = 0.0;
        for (int i = std::max(0, x - radius); i <= std::min(voxel_dims_[0] - 1, x + radius); i++)
        {
          for (int j = std::max(0, y - radius); j <= std::min(voxel_dims_[1] - 1, y + radius); j++)
          {
            const int n = getIndex(Eigen::Vector3i(i, j, z));
            hits += voxels_[n].numHits();
            path += voxels_[n].pathLength();
            if (!overlap_path_.empty())
            {
              k_path += overlap_path_[static_cast<size_t>(n)];
            }
          }
        }
        if (path <= 0.0)
        {
          continue;
        }
        const double measured = hits / path;  // local measured density
        double factor = 1.0;
        if (k_path > 0.0 && measured > 0.0 && path_estimator_ != PathEstimator::Crossing)
        {
          // measured = a R(pi d a / 4, rho) for true density a = N d. Solve for a by bisection, as a R is increasing
          const double rho = k_path / path;  // path-weighted mean ratio of effective width to diameter
          const double max_density =
            max_solid_fraction_ * 4.0 / (kPi * stem_diameter_ * (1.0 + min_gap_ratio_) * (1.0 + min_gap_ratio_));
          auto response = [&](double a) { return a * response_.response(0.25 * kPi * stem_diameter_ * a, rho); };
          // as the response R >= 1, the true density is at most the measured density, so search below it. This
          // only evaluates the response up to the measured solid fraction (the response is slow to compute near
          // the packing limit)
          double a = max_density;
          const double high_bound = std::min(measured, max_density);
          if (measured <= response(high_bound))
          {
            double low = 0.0, high = high_bound;
            for (int it = 0; it < 40; it++)
            {
              const double mid = 0.5 * (low + high);
              (response(mid) < measured ? low : high) = mid;
            }
            a = 0.5 * (low + high);
          }
          else
          {
            num_saturated++;
          }
          factor = a / measured;
          num_corrected++;
          sum_factor += factor;
          sum_k_path += overlap_path_[static_cast<size_t>(ind)];
          sum_path += voxels_[ind].pathLength();
        }
        // for the pooled estimator, the voxel takes the neighbourhood's density, keeping its own path length
        new_hits[static_cast<size_t>(ind)] =
          static_cast<float>(factor * (pooled ? measured * voxels_[ind].pathLength() : voxels_[ind].numHits()));
      }
    }
  }
  for (size_t i = 0; i < voxels_.size(); i++)
  {
    voxels_[i].numHits() = new_hits[i];
  }
  if (num_corrected > 0 && sum_path > 0.0)
  {
    std::cout << "stem occlusion correction: mean factor " << sum_factor / static_cast<double>(num_corrected)
              << " over " << num_corrected << " voxels, with mean effective width / diameter of "
              << sum_k_path / sum_path << std::endl;
    if (num_saturated > 0)
    {
      const double scale = (1.0 + min_gap_ratio_) * (1.0 + min_gap_ratio_);
      std::cout << "Warning: " << num_saturated << " voxels are beyond the maximum solid fraction ("
                << max_solid_fraction_ / scale << ") for " << stem_diameter_ << " m stems with gaps of "
                << min_gap_ratio_
                << " diameters, of N*d = " << max_solid_fraction_ * 4.0 / (kPi * stem_diameter_ * scale)
                << " m^-1. Their densities are set to this limit. Check the stem diameter" << std::endl;
    }
  }
}

// When the cloud has a sharp change in density at the top (e.g. grass or wheat field) between air and the crop, then
// this function adjusts the density estimation based on this two-phase density, rather than assuming the top voxel is
// uniform density
void DensityGrid::flatTopCompensation()
{
  for (int x = 0; x < voxel_dims_[0]; x++)
  {
    for (int y = 0; y < voxel_dims_[1]; y++)
    {
      double p = peaks_[x + voxel_dims_[0] * y];
      if (p < -10000.0)
        continue;
      int z = (int)std::floor(p);
      if (z < 0 || z >= voxel_dims_[2])
        continue;
      int ind = getIndex(Eigen::Vector3i(x, y, z));

      float r = voxels_[ind].numRays();
      if (r > 2.0f)
        voxels_[ind].numHits() *=
          (r - 2.0f) / (r - 1.0f);  // unbiased estimator reduces density estimation slightly in a manner that is
                                    // equivalent to choosing a grass height slightly above the peak point height
    }
  }
}

// This is a form of windowed average over the Moore neighbourhood (3x3x3) window.
void DensityGrid::addNeighbourPriors()
{
#if DENSITY_MIN_RAYS > 0
  const int X = 1;
  const int Y = voxel_dims_[0];
  const int Z = voxel_dims_[0] * voxel_dims_[1];
  DensityGrid::Voxel neighbours;
  double num_hit_points = 0.0;
  double num_hit_points_unsatisfied = 0.0;

  // This simple 3x3x3 convolution needs to be a bit sneaky to avoid having to double the memory cost.
  // well, not that sneaky, we just shift the output -1,-1,-1 for each cell
  for (int x = 1; x < voxel_dims_[0] - 1; x++)
  {
    for (int y = 1; y < voxel_dims_[1] - 1; y++)
    {
      for (int z = 1; z < voxel_dims_[2] - 1; z++)
      {
        const int ind = getIndex(Eigen::Vector3i(x, y, z));
        if (voxels_[ind].numHits() > 0)
          num_hit_points++;
        float needed = DENSITY_MIN_RAYS - voxels_[ind].numRays();
        const DensityGrid::Voxel corner_vox = voxels_[ind - X - Y - Z];
        voxels_[ind - X - Y - Z] = voxels_[ind];  // move centre up to corner
        DensityGrid::Voxel &voxel = voxels_[ind - X - Y - Z];
        if (needed < 0.0)
          continue;
        neighbours = voxels_[ind - X];
        neighbours += voxels_[ind + X];
        neighbours += voxels_[ind - Y];
        neighbours += voxels_[ind + Y];
        neighbours += voxels_[ind - Z];
        neighbours += voxels_[ind + Z];
        if (neighbours.numRays() >= needed)
        {
          voxel += neighbours * (needed / neighbours.numRays());  // add minimal amount to reach DENSITY_MIN_RAYS
          continue;
        }
        voxel += neighbours;
        needed -= neighbours.numRays();

        neighbours = voxels_[ind - X - Y];
        neighbours += voxels_[ind - X + Y];
        neighbours += voxels_[ind + X - Y];
        neighbours += voxels_[ind + X + Y];

        neighbours += voxels_[ind - X - Z];
        neighbours += voxels_[ind - X + Z];
        neighbours += voxels_[ind + X - Z];
        neighbours += voxels_[ind + X + Z];

        neighbours += voxels_[ind - Y - Z];
        neighbours += voxels_[ind - Y + Z];
        neighbours += voxels_[ind + Y - Z];
        neighbours += voxels_[ind + Y + Z];
        if (neighbours.numRays() >= needed)
        {
          voxel += neighbours * (needed / neighbours.numRays());  // add minimal amount to reach DENSITY_MIN_RAYS
          continue;
        }
        voxel += neighbours;
        needed -= neighbours.numRays();

        neighbours = corner_vox;
        neighbours += voxels_[ind - X - Y + Z];
        neighbours += voxels_[ind - X + Y - Z];
        neighbours += voxels_[ind + X - Y - Z];
        neighbours += voxels_[ind - X + Y + Z];
        neighbours += voxels_[ind + X - Y + Z];
        neighbours += voxels_[ind + X + Y - Z];
        neighbours += voxels_[ind + X + Y + Z];
        if (neighbours.numRays() >= needed)
        {
          voxel += neighbours * (needed / neighbours.numRays());  // add minimal amount to reach DENSITY_MIN_RAYS
          continue;
        }
        voxel += neighbours;
        if (voxels_[ind].numHits() > 0)
          num_hit_points_unsatisfied++;
      }
    }
  }
  const double percentage = 100.0 * num_hit_points_unsatisfied / num_hit_points;
  std::cout << "Density calculation: " << percentage << "% of voxels had insufficient (<" << DENSITY_MIN_RAYS
            << ") rays within them" << std::endl;
  if (percentage > 50.0)
  {
    std::cout << "This is high. Consider using a larger pixel size, or a denser cloud, or reducing DENSITY_MIN_RAYS, "
                 "for consistent results"
              << std::endl;
  }
  else if (percentage < 1.0)
  {
    std::cout << "This is low enough that you could get more fidelity from using a smaller pixel size" << std::endl;
    std::cout << "or more accuracy by increasing DENSITY_MIN_RAYS" << std::endl;
  }
#endif
}

/// The density of a voxel at valid (shifted) index x,y,z, as an average over the voxel
double voxelAverageDensity(const DensityGrid &grid, const DensityGrid::Voxel &voxel, int x, int y, int z,
                           const Eigen::Vector3i &full_dims, int shift)
{
  double d = voxel.density();
#if defined FLAT_TOP_COMPENSATION
  // the top voxel's density is estimated for the part below the column's highest point, so scale it to
  // the voxel average, as is done for the rendered image
  const double p = grid.peaks_[static_cast<size_t>((x + shift) + full_dims[0] * (y + shift))];
  const int top_z = p < -10000.0 ? -1 : static_cast<int>(std::floor(p));
  if (z == top_z - shift)
  {
    d *= p - static_cast<double>(top_z);
  }
#else
  RAYLIB_UNUSED(grid);
  RAYLIB_UNUSED(x);
  RAYLIB_UNUSED(y);
  RAYLIB_UNUSED(z);
  RAYLIB_UNUSED(full_dims);
  RAYLIB_UNUSED(shift);
#endif
  return d;
}

/// Writes the density of each cell of size @c cell_size (tiling from @c min_bound) to a .csv file. This aggregates the
/// voxels (whose centres lie in each cell) into cells at a coarser resolution, such as for hydrodynamic models. Two
/// estimates are output: the mean of the voxel densities, and the pooled estimate (total hits over total G-weighted
/// path length in the cell). With the beer_lambert estimator, the pooled estimate is identical to estimating with a
/// single voxel the size of the cell. The fraction of the cell's voxels that were observed is also output.
static bool writeDensityCells(const std::string &file_name, const DensityGrid &grid, const Eigen::Vector3i &dims,
                              const Eigen::Vector3i &full_dims, int shift, const Eigen::Vector3d &min_bound,
                              double pix_width, const Eigen::Vector3d &cell_size)
{
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing" << std::endl;
    return false;
  }
  Eigen::Vector3i num_cells;
  for (int i = 0; i < 3; i++)
  {
    num_cells[i] =
      std::max(1, static_cast<int>(std::ceil(static_cast<double>(dims[i]) * pix_width / cell_size[i] - 1e-6)));
  }
  struct Cell
  {
    double sum = 0.0, hits = 0.0, path_length = 0.0;
    int observed = 0, total = 0;
  };
  std::vector<Cell> cells(static_cast<size_t>(num_cells[0] * num_cells[1] * num_cells[2]));
  for (int z = 0; z < dims[2]; z++)
  {
    for (int y = 0; y < dims[1]; y++)
    {
      for (int x = 0; x < dims[0]; x++)
      {
        const Eigen::Vector3d centre = pix_width * (Eigen::Vector3d(x, y, z) + Eigen::Vector3d(0.5, 0.5, 0.5));
        Eigen::Vector3i c;
        for (int i = 0; i < 3; i++)
        {
          c[i] = std::min(num_cells[i] - 1, static_cast<int>(centre[i] / cell_size[i]));
        }
        Cell &cell = cells[static_cast<size_t>(c[0] + num_cells[0] * (c[1] + num_cells[1] * c[2]))];
        cell.total++;
        const DensityGrid::Voxel &voxel = grid.voxels()[grid.getIndex(Eigen::Vector3i(x, y, z))];
        if (voxel.numRays() > DensityGrid::min_voxel_hits)
        {
          cell.sum += voxelAverageDensity(grid, voxel, x, y, z, full_dims, shift);
          cell.observed++;
          cell.hits += voxel.numHits();
          cell.path_length += voxel.pathLength();
        }
      }
    }
  }
  const bool cylinders = grid.distribution() == ElementDistribution::VerticalCylinder;
  const std::string label = cylinders ? "_frontal_area_per_volume_Nd" : "_leaf_area_per_volume";
  ofs << "x_min,y_min,z_min,x_max,y_max,z_max,observed_fraction,hits,mean_density" << label << ",pooled_density"
      << label << std::endl;
  for (int z = 0; z < num_cells[2]; z++)
  {
    for (int y = 0; y < num_cells[1]; y++)
    {
      for (int x = 0; x < num_cells[0]; x++)
      {
        const Cell &cell = cells[static_cast<size_t>(x + num_cells[0] * (y + num_cells[1] * z))];
        const Eigen::Vector3d cell_min =
          min_bound + Eigen::Vector3d(x * cell_size[0], y * cell_size[1], z * cell_size[2]);
        const Eigen::Vector3d cell_max = cell_min + cell_size;
        ofs << cell_min[0] << "," << cell_min[1] << "," << cell_min[2] << "," << cell_max[0] << "," << cell_max[1]
            << "," << cell_max[2] << ","
            << (cell.total > 0 ? static_cast<double>(cell.observed) / static_cast<double>(cell.total) : 0.0) << ","
            << cell.hits << "," << (cell.observed > 0 ? cell.sum / static_cast<double>(cell.observed) : 0.0) << ","
            << (cell.path_length > 0.0 ? cell.hits / cell.path_length : 0.0) << std::endl;
      }
    }
  }
  std::cout << "density cells written to " << file_name << std::endl;
  return true;
}

/// Writes the density of each horizontal layer of the grid to a .csv file, both as the mean of the voxel densities
/// and as the pooled estimate (total hits over total G-weighted path length) across the layer.
/// The mean of voxel densities is biased high when voxels contain few hits (around +10% at 2-3 hits per voxel), so
/// the mean hits per voxel is also output. The pooled estimate avoids this bias, but weights voxels by path length.
/// Only voxels with sufficient rays are included, so unobserved space (e.g. below the ground) is excluded.
/// @c dims are the valid dimensions, where voxel z spans heights base_height + pix_width * [z, z+1]
/// @c shift is the voxel shift applied by addNeighbourPriors, and @c full_dims the full grid dimensions
static bool writeDensityProfile(const std::string &file_name, const DensityGrid &grid, const Eigen::Vector3i &dims,
                                const Eigen::Vector3i &full_dims, int shift, double base_height, double pix_width)
{
  std::ofstream ofs(file_name.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << file_name << " for writing" << std::endl;
    return false;
  }
  const bool cylinders = grid.distribution() == ElementDistribution::VerticalCylinder;
  const std::string label = cylinders ? "_frontal_area_per_volume_Nd" : "_leaf_area_per_volume";
  ofs << "z_min,z_max,observed_voxels,total_voxels,mean_hits_per_voxel,mean_density" << label << ",pooled_density"
      << label << std::endl;
  for (int z = 0; z < dims[2]; z++)
  {
    double sum = 0.0;
    double hits = 0.0, path_length = 0.0;
    int count = 0;
    for (int x = 0; x < dims[0]; x++)
    {
      for (int y = 0; y < dims[1]; y++)
      {
        const DensityGrid::Voxel &voxel = grid.voxels()[grid.getIndex(Eigen::Vector3i(x, y, z))];
        if (voxel.numRays() > DensityGrid::min_voxel_hits)
        {
          sum += voxelAverageDensity(grid, voxel, x, y, z, full_dims, shift);
          count++;
          hits += voxel.numHits();
          path_length += voxel.pathLength();
        }
      }
    }
    const double z_min = base_height + pix_width * static_cast<double>(z);
    ofs << z_min << "," << z_min + pix_width << "," << count << "," << dims[0] * dims[1] << ","
        << (count > 0 ? hits / static_cast<double>(count) : 0.0) << ","
        << (count > 0 ? sum / static_cast<double>(count) : 0.0) << "," << (path_length > 0.0 ? hits / path_length : 0.0)
        << std::endl;
  }
  std::cout << "density profile written to " << file_name << std::endl;
  return true;
}

std::unique_ptr<DensityGrid> estimateDensityGrid(const std::string &cloud_file, const Cuboid &bounds,
                                                 double voxel_width, const Eigen::Vector3i &dims_in,
                                                 const DensityOptions &options, int &shift)
{
  Eigen::Vector3i dims = dims_in;
  Cuboid grid_bounds = bounds;
  shift = 0;
#if DENSITY_MIN_RAYS > 0
  shift = 1;
  dims += 2 * Eigen::Vector3i(shift, shift, shift);  // so that we have extra space to convolve
  grid_bounds.min_bound_ -= Eigen::Vector3d(voxel_width, voxel_width, voxel_width);
#endif
  const bool cylinders = options.distribution == ElementDistribution::VerticalCylinder;
  auto grid_ptr = std::make_unique<DensityGrid>(grid_bounds, voxel_width, dims, options.distribution,
                                                cylinders ? options.min_sin_zenith : 0.0);
  DensityGrid &grid = *grid_ptr;
  grid.setPathEstimator(options.path_estimator);
  grid.setPoolRadius(options.pool_radius);
  grid.setMinGapRatio(options.min_gap_ratio);
  if (cylinders)
  {
    grid.setMaxG(options.max_sin_zenith);
  }
  if (cylinders && options.stem_diameter > 0.0)
  {
    grid.setStemDiameter(options.stem_diameter);
  }
  if (cylinders && options.footprint.enabled())
  {
    if (options.stem_diameter <= 0.0)
    {
      std::cerr << "Error: the beam footprint correction requires the stem diameter" << std::endl;
      return nullptr;
    }
    grid.setFootprint(options.footprint, options.stem_diameter);
    std::cout << "correcting for beam footprint, effective stem width of " << options.stem_diameter << " m stems is "
              << effectiveStemWidth(options.stem_diameter, options.footprint.diameter(1.0), options.footprint.threshold)
              << " m at 1 m range, "
              << effectiveStemWidth(options.stem_diameter, options.footprint.diameter(5.0), options.footprint.threshold)
              << " m at 5 m range" << std::endl;
  }

#if defined FLAT_TOP_COMPENSATION
  if (options.flat_top_compensation)  // without the peaks, no compensation is applied
  {
    grid.calculatePeaks(cloud_file);
  }
#endif
  if (cylinders && options.overlap_correction && options.stem_diameter > 0.0)
  {
    grid.enableOverlapCorrection();
  }
  grid.calculateDensities(cloud_file);
  if (grid.numIgnoredRays() > 0)
  {
    std::cout << grid.numIgnoredRays() << " rays ignored outside of the sin zenith angle range "
              << options.min_sin_zenith << " to " << options.max_sin_zenith << std::endl;
  }
#if defined FLAT_TOP_COMPENSATION
  if (options.flat_top_compensation)
  {
    grid.flatTopCompensation();
  }
#endif
  grid.applyOverlapCorrection();
#if DENSITY_MIN_RAYS > 0
  grid.addNeighbourPriors();
#endif
  return grid_ptr;
}

bool renderCloud(const std::string &cloud_file, const Cuboid &bounds, ViewDirection view_direction, RenderStyle style,
                 double pix_width, const std::string &image_file, const std::string &projection_file, bool mark_origin,
                 const std::string *const transform_file, const DensityOptions &density_options)
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
      const Eigen::Vector3i dims = valid_grid_dims + 2 * Eigen::Vector3i(shift, shift, shift);
      const bool cylinders = density_options.distribution == ElementDistribution::VerticalCylinder;
      if (!density_options.profile_file.empty())
      {
        // addNeighbourPriors shifts the voxel contents by -shift in each axis, so these are the valid dimensions,
        // with index 0 at bounds.min_bound_
        const Eigen::Vector3i valid_dims = dims - 2 * Eigen::Vector3i(shift, shift, shift);
        if (!writeDensityProfile(density_options.profile_file, grid, valid_dims, dims, shift, bounds.min_bound_[2],
                                 pix_width))
        {
          return false;
        }
      }
      if (!density_options.cells_file.empty())
      {
        const Eigen::Vector3i valid_dims = dims - 2 * Eigen::Vector3i(shift, shift, shift);
        if (!writeDensityCells(density_options.cells_file, grid, valid_dims, dims, shift, bounds.min_bound_, pix_width,
                               density_options.cell_size))
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
            double d = grid.voxels()[grid.getIndex(ind)].density();
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
