// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Build against raylib (in a raycloudtools build directory <build>, for the source directory <src>), e.g.:
//   g++ -O2 -std=c++17 -I<src> -I<src>/raylib -I<build>/raylib -I<build> -I/usr/include/eigen3
//     <src>/scripts/pneumatophore_self_shadow_mc2d.cpp -L<build>/raylib -lraylib -lpthread -o mc2d
// Run with key value pairs, e.g. ./mc2d N 300 d 0.01 deff 0.0234 D 0.22 w 0.05 k 1 M 9000000 seed 1
// (scripts/pneumatophore_self_shadow_study.py is the 3D study; this isolates the mechanism).
//
// 2D Monte Carlo of the pooled-window stem density estimator, for infinitely tall vertical stems (so the horizontal
// projection is exact). Periodic square domain (torus), no edges. Outline return convention: a line stops where it
// enters a stem's widened outline (width d_eff); lines start outside the widened outlines.
//
// Each line starts at a uniform random free point with a uniform random direction, travels a "pre-path" D without
// being recorded (it can be stopped there), then is recorded (path per voxel, hit in the voxel of the return) until
// its return or a maximum recorded length. D = 0 is the ensemble of the response table (paths start at random
// points); D >= sqrt(2) W makes every stem within a window upstream of a recorded point lie on the line's past, as for
// tall stems in 3D (the shadow condition is a property of the stem, not of the line).
//
// Per voxel window ((2k+1) x (2k+1) voxels of width w, W = (2k+1) w), with H hits and P path:
//   mean      : mean over windows of invert(H / P)                      (the library's `density`)
//   pooled    : invert(sum H / sum P)                                    (ratio of expectations)
//   selffree  : mean over windows of invert(sum_{i in v} h_i / (P + C_iv)) where C_iv is the counterfactual path
//               within the window of the lines that returned from stem i, continued through stem i (but not the
//               others) to their next return: removes each stem's own shadow from its own window only
//   credited  : mean over windows of invert(H / (P + sum_hits s)), s the geometric distance from the return to the
//               window edge (full crediting, as the crossing estimator does at voxel scale)
//   corr_a1   : mean over windows of invert((H/P) / (1 + c d_eff / W))
//   corr_a2   : mean over windows of a solving H/P = a R(a) (1 + S(phi(a)) c d_eff / W), S = (1-p)^3/(1+p) at the
//               solid fraction p of the exclusion discs
// Also the Fano factor of the window stem counts. Output: one line of key=value.
#include "raylib/raystemmodel.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
const double kPi = 3.14159265358979323846;
struct Rng
{
  std::mt19937_64 g;
  explicit Rng(unsigned long long s)
    : g(s)
  {}
  double u() { return (static_cast<double>(g() >> 11) + 0.5) / 9007199254740992.0; }
};

double structureFactor(double p)
{
  return std::pow(1.0 - p, 3.0) / (1.0 + p);
}
}  // namespace

int main(int argc, char **argv)
{
  // defaults
  double N = 300, d = 0.01, deff = 0.02, gap = 0.0, L = 3.0, w = 0.05, D = 0.25, maxrec = 1.0;
  int k = 1, poisson = 0;
  long long M = 2000000;
  unsigned long long seed = 1;
  for (int i = 1; i + 1 < argc; i += 2)
  {
    std::string a = argv[i];
    double v = std::atof(argv[i + 1]);
    if (a == "N")
      N = v;
    else if (a == "d")
      d = v;
    else if (a == "deff")
      deff = v;
    else if (a == "gap")
      gap = v;
    else if (a == "L")
      L = v;
    else if (a == "w")
      w = v;
    else if (a == "k")
      k = static_cast<int>(v);
    else if (a == "D")
      D = v;
    else if (a == "maxrec")
      maxrec = v;
    else if (a == "M")
      M = static_cast<long long>(v);
    else if (a == "seed")
      seed = static_cast<unsigned long long>(v);
    else if (a == "poisson")
      poisson = static_cast<int>(v);
  }
  Rng rng(seed * 7919ULL + 17ULL);
  const int nv = static_cast<int>(std::round(L / w));
  L = nv * w;
  const double r_eff = 0.5 * deff;
  const double sep = d * (1.0 + gap);
  // stems
  const double cw = std::max(std::max(r_eff, sep), 0.02);
  const int nc = std::max(3, static_cast<int>(std::floor(L / cw)));
  const double cwr = L / nc;
  std::vector<std::vector<int>> cells(static_cast<size_t>(nc * nc));
  std::vector<double> sx, sy;
  auto cell_of = [&](double x) {
    int c = static_cast<int>(std::floor(x / cwr));
    return ((c % nc) + nc) % nc;
  };
  auto wrapd = [&](double dx) { return dx - L * std::round(dx / L); };
  long long target;
  if (poisson)
  {
    // Poisson count by normal approximation of a large mean (exact enough)
    std::normal_distribution<double> nd(N * L * L, std::sqrt(N * L * L));
    target = static_cast<long long>(std::llround(std::max(0.0, nd(rng.g))));
  }
  else
  {
    target = static_cast<long long>(std::llround(N * L * L));
  }
  long long attempts = 0;
  while (static_cast<long long>(sx.size()) < target && attempts < 200000000LL)
  {
    attempts++;
    const double x = L * rng.u(), y = L * rng.u();
    bool ok = true;
    if (!poisson)
    {
      const int cx = cell_of(x), cy = cell_of(y);
      for (int dx = -1; dx <= 1 && ok; dx++)
        for (int dy = -1; dy <= 1 && ok; dy++)
        {
          const int bx = (cx + dx + nc) % nc, by = (cy + dy + nc) % nc;
          for (int j : cells[static_cast<size_t>(bx + nc * by)])
          {
            const double ex = wrapd(sx[j] - x), ey = wrapd(sy[j] - y);
            if (ex * ex + ey * ey < sep * sep)
            {
              ok = false;
              break;
            }
          }
        }
    }
    if (ok)
    {
      cells[static_cast<size_t>(cell_of(x) + nc * cell_of(y))].push_back(static_cast<int>(sx.size()));
      sx.push_back(x);
      sy.push_back(y);
    }
  }
  const size_t ns = sx.size();
  const double Nreal = ns / (L * L);
  const double phi = Nreal * kPi * d * d / 4.0;
  // whether a point is inside a widened outline
  auto inside = [&](double x, double y) {
    const int cx = cell_of(x), cy = cell_of(y);
    for (int dx = -1; dx <= 1; dx++)
      for (int dy = -1; dy <= 1; dy++)
      {
        const int bx = (cx + dx + nc) % nc, by = (cy + dy + nc) % nc;
        for (int j : cells[static_cast<size_t>(bx + nc * by)])
        {
          const double ex = wrapd(sx[j] - x), ey = wrapd(sy[j] - y);
          if (ex * ex + ey * ey < r_eff * r_eff)
            return true;
        }
      }
    return false;
  };
  // first return along a line from (x, y) in direction (ux, uy), from t > tmin to t < tmax, ignoring stem `skip`.
  // Returns t (or tmax) and the stem index (or -1)
  auto cast = [&](double x, double y, double ux, double uy, double tmin, double tmax, int skip, int &hit) {
    hit = -1;
    double best = tmax;
    // walk cells in unwrapped coordinates
    double px = x + ux * tmin, py = y + uy * tmin;
    int cx = static_cast<int>(std::floor(px / cwr)), cy = static_cast<int>(std::floor(py / cwr));
    const int stx = ux > 0 ? 1 : -1, sty = uy > 0 ? 1 : -1;
    double nx = std::abs(ux) > 1e-12 ? ((cx + (stx > 0)) * cwr - x) / ux : 1e300;
    double ny = std::abs(uy) > 1e-12 ? ((cy + (sty > 0)) * cwr - y) / uy : 1e300;
    const double dxs = std::abs(ux) > 1e-12 ? cwr / std::abs(ux) : 1e300;
    const double dys = std::abs(uy) > 1e-12 ? cwr / std::abs(uy) : 1e300;
    while (true)
    {
      for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++)
        {
          const int ux_c = cx + dx, uy_c = cy + dy;
          const int bx = ((ux_c % nc) + nc) % nc, by = ((uy_c % nc) + nc) % nc;
          const double offx = (ux_c - bx) / nc * L, offy = (uy_c - by) / nc * L;
          for (int j : cells[static_cast<size_t>(bx + nc * by)])
          {
            if (j == skip)
              continue;
            const double rx = sx[j] + offx - x, ry = sy[j] + offy - y;
            const double b = rx * ux + ry * uy;
            const double o2 = rx * rx + ry * ry - b * b;
            if (o2 < r_eff * r_eff)
            {
              const double t = b - std::sqrt(r_eff * r_eff - o2);
              if (t > tmin && t < best)
              {
                best = t;
                hit = j;
              }
            }
          }
        }
      const double te = std::min(nx, ny);
      if (best <= te || te >= tmax)
        break;
      if (nx < ny)
      {
        cx += stx;
        nx += dxs;
      }
      else
      {
        cy += sty;
        ny += dys;
      }
    }
    return best;
  };
  // per voxel path and hits
  const size_t nvox = static_cast<size_t>(nv) * nv;
  std::vector<double> path(nvox, 0.0), hits(nvox, 0.0);
  auto vox = [&](double x, double y) {
    int ix = static_cast<int>(std::floor(x / w)), iy = static_cast<int>(std::floor(y / w));
    ix = ((ix % nv) + nv) % nv;
    iy = ((iy % nv) + nv) % nv;
    return static_cast<size_t>(ix + nv * iy);
  };
  // add the path of segment t0..t1 of a line to the voxels (or to a callback)
  auto walk = [&](double x, double y, double ux, double uy, double t0, double t1, auto &&add) {
    double t = t0;
    double px = x + ux * t0, py = y + uy * t0;
    int ix = static_cast<int>(std::floor(px / w)), iy = static_cast<int>(std::floor(py / w));
    const int stx = ux > 0 ? 1 : -1, sty = uy > 0 ? 1 : -1;
    double nx = std::abs(ux) > 1e-12 ? ((ix + (stx > 0)) * w - x) / ux : 1e300;
    double ny = std::abs(uy) > 1e-12 ? ((iy + (sty > 0)) * w - y) / uy : 1e300;
    const double dxs = std::abs(ux) > 1e-12 ? w / std::abs(ux) : 1e300;
    const double dys = std::abs(uy) > 1e-12 ? w / std::abs(uy) : 1e300;
    while (t < t1)
    {
      const double te = std::min(std::min(nx, ny), t1);
      add(((ix % nv) + nv) % nv, ((iy % nv) + nv) % nv, te - t);
      t = te;
      if (nx < ny)
      {
        ix += stx;
        nx += dxs;
      }
      else
      {
        iy += sty;
        ny += dys;
      }
    }
  };
  // per stem, the counterfactual continuation path in each of its (2k+1)^2 windows, and its hits
  const int nwin = (2 * k + 1) * (2 * k + 1);
  std::vector<float> cf(ns * nwin, 0.0f);
  std::vector<double> stem_hits(ns, 0.0);
  std::vector<int> stem_vox_x(ns), stem_vox_y(ns);
  for (size_t j = 0; j < ns; j++)
  {
    stem_vox_x[j] = static_cast<int>(std::floor(sx[j] / w)) % nv;
    stem_vox_y[j] = static_cast<int>(std::floor(sy[j] / w)) % nv;
  }
  // per window (indexed by centre voxel), the sum over its hits of the geometric distance to the window edge
  std::vector<double> credit(nvox, 0.0);
  const double W = (2 * k + 1) * w;
  long long recorded = 0;
  for (long long r = 0; r < M; r++)
  {
    double x, y;
    do
    {
      x = L * rng.u();
      y = L * rng.u();
    } while (inside(x, y));
    const double ang = 2.0 * kPi * rng.u();
    const double ux = std::cos(ang), uy = std::sin(ang);
    int hit;
    const double tend = cast(x, y, ux, uy, 0.0, D + maxrec, -1, hit);
    if (tend <= D)
      continue;  // stopped before being recorded
    recorded++;
    walk(x, y, ux, uy, D, tend, [&](int ix, int iy, double len) { path[static_cast<size_t>(ix + nv * iy)] += len; });
    if (hit < 0)
      continue;
    const double hx = x + ux * tend, hy = y + uy * tend;
    const size_t hv = vox(hx, hy);
    hits[hv] += 1.0;
    stem_hits[static_cast<size_t>(hit)] += 1.0;
    // counterfactual continuation through stem `hit` to the next return (within reach of its windows)
    int hit2;
    const double reach = std::sqrt(2.0) * W + deff;
    const double t2 = cast(x, y, ux, uy, tend, tend + reach, hit, hit2);
    const int svx = stem_vox_x[static_cast<size_t>(hit)], svy = stem_vox_y[static_cast<size_t>(hit)];
    // accumulate per voxel of the continuation, then distribute to the stem's windows
    std::vector<std::pair<int, double>> seg;
    walk(x, y, ux, uy, tend, t2, [&](int ix, int iy, double len) { seg.emplace_back(ix + nv * iy, len); });
    for (int wy = -k; wy <= k; wy++)
      for (int wx = -k; wx <= k; wx++)
      {
        // window centred at voxel (svx + wx, svy + wy) contains the stem's voxel
        const int cxv = svx + wx, cyv = svy + wy;
        double sum = 0.0;
        for (const auto &s : seg)
        {
          const int ix = s.first % nv, iy = s.first / nv;
          int ddx = ix - cxv, ddy = iy - cyv;
          ddx -= nv * static_cast<int>(std::round(static_cast<double>(ddx) / nv));
          ddy -= nv * static_cast<int>(std::round(static_cast<double>(ddy) / nv));
          if (std::abs(ddx) <= k && std::abs(ddy) <= k)
            sum += s.second;
        }
        cf[static_cast<size_t>(hit) * nwin + static_cast<size_t>((wx + k) + (2 * k + 1) * (wy + k))] +=
          static_cast<float>(sum);
      }
    // geometric credit for the windows containing the hit voxel: distance from the return to the window edge
    const int hvx = static_cast<int>(hv % nv), hvy = static_cast<int>(hv / nv);
    for (int wy = -k; wy <= k; wy++)
      for (int wx = -k; wx <= k; wx++)
      {
        const int cxv = ((hvx + wx) % nv + nv) % nv, cyv = ((hvy + wy) % nv + nv) % nv;
        // window bounds in unwrapped coordinates around the hit point
        const double x0 = (hvx + wx - k) * w, x1 = (hvx + wx + k + 1) * w;
        const double y0 = (hvy + wy - k) * w, y1 = (hvy + wy + k + 1) * w;
        const double hxu = hx - L * std::floor(hx / L), hyu = hy - L * std::floor(hy / L);
        double tx = ux > 0 ? (x1 - hxu) / ux : (ux < 0 ? (x0 - hxu) / ux : 1e300);
        double ty = uy > 0 ? (y1 - hyu) / uy : (uy < 0 ? (y0 - hyu) / uy : 1e300);
        credit[static_cast<size_t>(cxv + nv * cyv)] += std::max(0.0, std::min(tx, ty));
      }
  }
  // windows
  const double rho = deff / d;
  ray::StemResponseTable table(gap, true);
  const double c = (4.0 * std::log(1.0 + std::sqrt(2.0)) - (4.0 / 3.0) * (std::sqrt(2.0) - 1.0)) / (2.0 * kPi);
  auto invert = [&](double measured_rate) {
    // measured_rate = hits / path, of rate N d_eff R. The library's measured density is hits / (rho path)
    bool sat;
    return ray::invertStemResponse(measured_rate / rho, rho, d, gap, 0.4, table, sat);
  };
  const double scale = (1.0 + gap) * (1.0 + gap);
  auto invert_a2 = [&](double measured_rate) {
    bool sat;
    const double max_density = 0.4 * 4.0 / (kPi * d * scale);
    auto response = [&](double a) {
      const double p = 0.25 * kPi * d * a;
      return table.response(p, rho) * (1.0 + structureFactor(p * scale) * c * deff / W);
    };
    return ray::invertDensityResponse(measured_rate / rho, max_density, response, sat);
  };
  double sum_mean = 0, sum_raw = 0, sum_self = 0, sum_cred = 0, sum_a1 = 0, sum_a2 = 0, sumH = 0, sumP = 0;
  double sum_n = 0, sum_n2 = 0;
  long long nwins = 0;
  // stems per voxel
  std::vector<int> nvoxstems(nvox, 0);
  for (size_t j = 0; j < ns; j++) nvoxstems[static_cast<size_t>(stem_vox_x[j] + nv * stem_vox_y[j])]++;
  // stems per voxel list
  std::vector<std::vector<int>> voxstems(nvox);
  for (size_t j = 0; j < ns; j++)
    voxstems[static_cast<size_t>(stem_vox_x[j] + nv * stem_vox_y[j])].push_back(static_cast<int>(j));
  for (int cy = 0; cy < nv; cy++)
    for (int cx = 0; cx < nv; cx++)
    {
      double H = 0, P = 0;
      int n = 0;
      double self_num_over = 0.0;                 // sum over stems of h_i / (P + C_iv), filled after P known
      std::vector<std::pair<int, int>> stems_in;  // stem, window slot
      for (int dy = -k; dy <= k; dy++)
        for (int dx = -k; dx <= k; dx++)
        {
          const size_t v = static_cast<size_t>(((cx + dx) % nv + nv) % nv + nv * (((cy + dy) % nv + nv) % nv));
          H += hits[v];
          P += path[v];
          n += nvoxstems[v];
          for (int j : voxstems[v])
          {
            // the stem's voxel is at offset (dx, dy) from this window's centre, so the window is centred at offset
            // (-dx, -dy) from the stem's voxel
            stems_in.emplace_back(j, (-dx + k) + (2 * k + 1) * (-dy + k));
          }
        }
      if (P <= 0)
        continue;
      nwins++;
      sumH += H;
      sumP += P;
      sum_n += n;
      sum_n2 += static_cast<double>(n) * n;
      const double m = H / P;
      sum_raw += m / rho;
      sum_mean += invert(m);
      for (const auto &s : stems_in)
        self_num_over +=
          stem_hits[static_cast<size_t>(s.first)] / (P + cf[static_cast<size_t>(s.first) * nwin + s.second]);
      sum_self += invert(self_num_over);
      sum_cred += invert(H / (P + credit[static_cast<size_t>(cx + nv * cy)]));
      sum_a1 += invert(m / (1.0 + c * deff / W));
      sum_a2 += invert_a2(m);
    }
  const double truth = Nreal * d;
  const double nmean = sum_n / nwins;
  const double fano = (sum_n2 / nwins - nmean * nmean) / nmean;
  const double globalR = (sumH / sumP) / (Nreal * deff);
  printf("N=%g d=%g deff=%g gap=%g poisson=%d w=%g k=%d W=%g D=%g M=%lld seed=%llu stems=%zu phi=%.4f recorded=%lld "
         "hits=%.0f fano=%.4f S_formula=%.4f R_global=%.4f R_table=%.4f mean=%.4f raw_mean=%.4f pooled=%.4f "
         "selffree=%.4f credited=%.4f corr_a1=%.4f corr_a2=%.4f pred_H1=%.4f pred_F=%.4f\n",
         N, d, deff, gap, poisson, w, k, W, D, M, seed, ns, phi, recorded, sumH, fano, structureFactor(phi * scale),
         globalR, table.response(phi, rho), sum_mean / nwins / truth, sum_raw / nwins / truth,
         invert(sumH / sumP) / truth, sum_self / nwins / truth, sum_cred / nwins / truth, sum_a1 / nwins / truth,
         sum_a2 / nwins / truth, 1.0 + c * deff / W, 1.0 + fano * c * deff / W);
  return 0;
}
