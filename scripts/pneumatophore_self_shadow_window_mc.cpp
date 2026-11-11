// Copyright (c) 2026
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// 2D Monte Carlo of the self-shadow factor of the pooled stem density estimator's neighbourhoods, for tall vertical
// stems (the horizontal projection is exact), in a periodic square, as used to fit ray::selfShadowFactor. Standalone:
//   g++ -O2 -std=c++17 scripts/pneumatophore_self_shadow_window_mc.cpp -o window_mc
//   ./window_mc N 300 d 0.01 deff 0.0225 w 0.05 conv 0 M 6e6 seed 1 D 0.08 L 3
// Stems of diameter d (non-overlapping, by random sequential addition, or Poisson with poisson 1) at N per m^2, seen
// with effective width deff. Lines start at random points outside the widened outlines with random directions, travel
// D unrecorded (D of at least 1.5 neighbourhood widths, so that recorded lines enter the neighbourhoods from outside),
// then are recorded to their first return or maxrec. The return is at the widened outline (conv 1), or at the
// physical surface or the closest approach of the line to the stem centre (conv 0, the library's default).
// Output, for neighbourhoods of (2k+1)^2 voxels of width w, k = 0..3:
//   blk  the mean over neighbourhoods of hits / path divided by the pooled sum of hits / sum of path: the
//        self-shadow factor 1 + beta that the library divides each neighbourhood by
//   crk  the same with the crossing chord (to the exit of the hit voxel) instead of the path to the hit
// The fit used N d_eff w <= 0.3 and d_eff / w from 0.1 to 1 (k = 0), with d / d_eff of 0.3, 0.44 and 1, two seeds
// each, and about 1e6 recorded hits per run.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

int main(int argc, char **argv)
{
  double N = 300, d = 0.01, deff = 0.023, w = 0.05, L = 3.0, D = 0.3, maxrec = 1.0, M = 2e6;
  int conv = 0, poisson = 0;
  unsigned long long seed = 1;
  for (int i = 1; i + 1 < argc; i += 2)
  {
    std::string a = argv[i];
    double v = std::atof(argv[i + 1]);
    if (a == "N") N = v;
    else if (a == "d") d = v;
    else if (a == "deff") deff = v;
    else if (a == "w") w = v;
    else if (a == "L") L = v;
    else if (a == "D") D = v;
    else if (a == "maxrec") maxrec = v;
    else if (a == "M") M = v;
    else if (a == "conv") conv = static_cast<int>(v);
    else if (a == "poisson") poisson = static_cast<int>(v);
    else if (a == "seed") seed = static_cast<unsigned long long>(v);
  }
  std::mt19937_64 g(seed * 1000003ULL + 7ULL);
  auto u = [&]() { return (static_cast<double>(g() >> 11) + 0.5) / 9007199254740992.0; };
  const int nv = static_cast<int>(std::round(L / w));
  L = nv * w;
  const double re = 0.5 * deff, rd = 0.5 * d;
  const double cw = std::max(std::max(re, d) * 1.0001, std::min(0.05, 1.5 / std::sqrt(std::max(N, 1.0))));
  const int nc = std::max(3, static_cast<int>(std::floor(L / cw)));
  const double cwr = L / nc;
  std::vector<std::vector<int>> cells(static_cast<size_t>(nc) * nc);
  std::vector<double> sx, sy;
  auto cell_of = [&](double x) {
    int c = static_cast<int>(std::floor(x / cwr));
    return ((c % nc) + nc) % nc;
  };
  auto wrapd = [&](double dx) { return dx - L * std::round(dx / L); };
  long long target = std::llround(N * L * L);
  if (poisson)
  {
    std::poisson_distribution<long long> pd(N * L * L);
    target = pd(g);
  }
  long long attempts = 0;
  while (static_cast<long long>(sx.size()) < target && attempts < 400000000LL)
  {
    attempts++;
    const double x = L * u(), y = L * u();
    bool ok = true;
    if (!poisson)
    {
      const int cx = cell_of(x), cy = cell_of(y);
      for (int a = -1; a <= 1 && ok; a++)
        for (int b = -1; b <= 1 && ok; b++)
          for (int j : cells[static_cast<size_t>(((cx + a + nc) % nc) + nc * ((cy + b + nc) % nc))])
          {
            const double ex = wrapd(sx[j] - x), ey = wrapd(sy[j] - y);
            if (ex * ex + ey * ey < d * d)
            {
              ok = false;
              break;
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
  const double Nreal = sx.size() / (L * L);
  auto inside = [&](double x, double y) {
    const int cx = cell_of(x), cy = cell_of(y);
    for (int a = -1; a <= 1; a++)
      for (int b = -1; b <= 1; b++)
        for (int j : cells[static_cast<size_t>(((cx + a + nc) % nc) + nc * ((cy + b + nc) % nc))])
        {
          const double ex = wrapd(sx[j] - x), ey = wrapd(sy[j] - y);
          if (ex * ex + ey * ey < re * re) return true;
        }
    return false;
  };
  // first return along the line from t > 0 to tmax
  auto cast = [&](double x, double y, double ux, double uy, double tmax) {
    double best = tmax;
    int cx = static_cast<int>(std::floor(x / cwr)), cy = static_cast<int>(std::floor(y / cwr));
    const int stx = ux > 0 ? 1 : -1, sty = uy > 0 ? 1 : -1;
    double nx = std::abs(ux) > 1e-12 ? ((cx + (stx > 0)) * cwr - x) / ux : 1e300;
    double ny = std::abs(uy) > 1e-12 ? ((cy + (sty > 0)) * cwr - y) / uy : 1e300;
    const double dxs = std::abs(ux) > 1e-12 ? cwr / std::abs(ux) : 1e300;
    const double dys = std::abs(uy) > 1e-12 ? cwr / std::abs(uy) : 1e300;
    while (true)
    {
      for (int a = -1; a <= 1; a++)
        for (int b = -1; b <= 1; b++)
        {
          const int ucx = cx + a, ucy = cy + b;
          const int bx = ((ucx % nc) + nc) % nc, by = ((ucy % nc) + nc) % nc;
          const double offx = static_cast<double>(ucx - bx) / nc * L, offy = static_cast<double>(ucy - by) / nc * L;
          for (int j : cells[static_cast<size_t>(bx + nc * by)])
          {
            const double rx = sx[j] + offx - x, ry = sy[j] + offy - y;
            const double bb = rx * ux + ry * uy;
            const double o2 = rx * rx + ry * ry - bb * bb;
            if (o2 < re * re)
            {
              double t;
              if (conv == 1)
                t = bb - std::sqrt(re * re - o2);
              else
                t = o2 < rd * rd ? bb - std::sqrt(rd * rd - o2) : bb;
              if (t > 1e-12 && t < best) best = t;
            }
          }
        }
      const double te = std::min(nx, ny);
      if (best <= te || te >= tmax) break;
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
  const size_t nvox = static_cast<size_t>(nv) * nv;
  std::vector<double> P(nvox, 0.0), C(nvox, 0.0), H(nvox, 0.0);
  // walk the segment t0..t1 adding path; returns the voxel of t1 and the distance from t1 to that voxel's exit
  auto walk = [&](double x, double y, double ux, double uy, double t0, double t1, std::vector<double> &acc,
                  double &rest) {
    double t = t0;
    const double px = x + ux * t0, py = y + uy * t0;
    int ix = static_cast<int>(std::floor(px / w)), iy = static_cast<int>(std::floor(py / w));
    const int stx = ux > 0 ? 1 : -1, sty = uy > 0 ? 1 : -1;
    double nx = std::abs(ux) > 1e-12 ? ((ix + (stx > 0)) * w - x) / ux : 1e300;
    double ny = std::abs(uy) > 1e-12 ? ((iy + (sty > 0)) * w - y) / uy : 1e300;
    const double dxs = std::abs(ux) > 1e-12 ? w / std::abs(ux) : 1e300;
    const double dys = std::abs(uy) > 1e-12 ? w / std::abs(uy) : 1e300;
    size_t last = 0;
    while (true)
    {
      const double te = std::min(std::min(nx, ny), t1);
      last = static_cast<size_t>(((ix % nv) + nv) % nv + nv * (((iy % nv) + nv) % nv));
      acc[last] += te - t;
      t = te;
      if (t >= t1)
      {
        rest = std::min(nx, ny) - t1;
        break;
      }
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
    return last;
  };
  const long long Mi = static_cast<long long>(M);
  long long recorded = 0;
  for (long long r = 0; r < Mi; r++)
  {
    double x, y;
    do
    {
      x = L * u();
      y = L * u();
    } while (inside(x, y));
    const double ang = 2.0 * M_PI * u();
    const double ux = std::cos(ang), uy = std::sin(ang);
    const double tend = cast(x, y, ux, uy, D + maxrec);
    if (tend <= D) continue;
    recorded++;
    double rest = 0.0;
    const size_t v = walk(x, y, ux, uy, D, tend, P, rest);
    double unused;
    walk(x, y, ux, uy, D, tend, C, unused);
    if (tend < D + maxrec)
    {
      H[v] += 1.0;
      C[v] += rest;
    }
  }
  double sumH = 0, sumP = 0, sumC = 0;
  for (size_t i = 0; i < nvox; i++)
  {
    sumH += H[i];
    sumP += P[i];
    sumC += C[i];
  }
  printf("N=%g Nreal=%g d=%g deff=%g w=%g conv=%d poisson=%d M=%g seed=%llu recorded=%lld hits=%.0f phi=%.5f "
         "R_global=%.4f",
         N, Nreal, d, deff, w, conv, poisson, M, seed, recorded, sumH, Nreal * M_PI * d * d / 4,
         (sumH / sumP) / (Nreal * deff));
  for (int k = 0; k <= 3; k++)
  {
    double sm = 0, sc = 0;
    long long nw = 0;
    for (int cy = 0; cy < nv; cy++)
      for (int cx = 0; cx < nv; cx++)
      {
        double h = 0, p = 0, c = 0;
        for (int b = -k; b <= k; b++)
          for (int a = -k; a <= k; a++)
          {
            const size_t i = static_cast<size_t>(((cx + a) % nv + nv) % nv + nv * (((cy + b) % nv + nv) % nv));
            h += H[i];
            p += P[i];
            c += C[i];
          }
        if (p <= 0 || c <= 0) continue;
        nw++;
        sm += h / p;
        sc += h / c;
      }
    printf(" W%d=%.4f bl%d=%.5f cr%d=%.5f", k, (2 * k + 1) * w, k, (sm / nw) / (sumH / sumP), k,
           (sc / nw) / (sumH / sumC));
  }
  printf("\n");
  return 0;
}
