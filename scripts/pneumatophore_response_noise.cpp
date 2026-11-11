// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Accuracy of the stem occlusion response, as in Section 5 of docs/pneumatophore_density_report.md: the
// seed-to-seed noise of ray::freePathResponse (64 realisations each, 6 seeds) at a few nodes, for both return
// conventions; the embedded, interpolated table (ray::StemResponseTable) against direct evaluation off the grid
// nodes; and direct values for minimum gaps of 0.5 and 1 diameter. Prints CSV to stdout.
//
// Build against raylib (in a raycloudtools build directory <build>, for the source directory <src>), e.g.:
//   g++ -O2 -std=c++17 -I<src> -I<src>/raylib -I<build>/raylib -I<build> -I/usr/include/eigen3
//     <src>/scripts/pneumatophore_response_noise.cpp -L<build>/raylib -lraylib -lpthread -o response_noise
//   ./response_noise > response.csv
#include "raylib/raystemmodel.h"
#include <cmath>
#include <cstdio>
#include <vector>
int main()
{
  struct P
  {
    double phi, rho;
  };
  const std::vector<P> nodes = { { 0.02, 2.25 }, { 0.1, 2.25 }, { 0.2, 2.25 }, { 0.1, 1.0 }, { 0.25, 3.0 } };
  printf("kind,phi,rho,convention,mean,sd,min,max,n_seeds\n");
  for (int conv = 0; conv < 2; conv++)
    for (const P &p : nodes)
    {
      std::vector<double> v;
      for (unsigned s = 0; s < 6; s++)
        v.push_back(ray::freePathResponse(p.phi, p.rho, 0.0, conv == 1, 1000u + 104729u * s));
      double m = 0, mn = 1e9, mx = -1e9;
      for (double x : v)
      {
        m += x;
        mn = std::fmin(mn, x);
        mx = std::fmax(mx, x);
      }
      m /= v.size();
      double sd = 0;
      for (double x : v) sd += (x - m) * (x - m);
      sd = std::sqrt(sd / (v.size() - 1));
      printf("seeds,%g,%g,%s,%.4f,%.4f,%.4f,%.4f,6\n", p.phi, p.rho, conv ? "outline" : "surface", m, sd, mn, mx);
      fflush(stdout);
    }
  const std::vector<P> off = { { 0.0236, 2.3 }, { 0.0785, 2.3 }, { 0.2, 2.3 }, { 0.15, 2.6 }, { 0.3, 2.1 } };
  for (int conv = 0; conv < 2; conv++)
  {
    ray::StemResponseTable table(0.0, conv == 1);
    for (const P &p : off)
      printf("table_vs_direct,%g,%g,%s,%.4f,%.4f,,,\n", p.phi, p.rho, conv ? "outline" : "surface",
             table.response(p.phi, p.rho), ray::freePathResponse(p.phi, p.rho, 0.0, conv == 1));
    fflush(stdout);
  }
  for (double gap : { 0.5, 1.0 })
  {
    for (const P &p : { P{ 0.0785, 2.3 }, P{ 0.0785, 1.0 } })
      printf("gap%g,%g,%g,surface,%.4f,,,,\n", gap, p.phi, p.rho, ray::freePathResponse(p.phi, p.rho, gap, false));
  }
  return 0;
}
