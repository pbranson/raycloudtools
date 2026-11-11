#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Synthetic rod (dowel) test for calibrating the beam footprint (exit diameter w0, divergence g, threshold t) of
the stem density estimator, and the design analysis behind it.

design: for candidate rod tests (rod diameters x ranges), the Fisher information of the d_eff(r) measurements about
  (ln w0, ln g, ln t), the parameter standard errors and correlations, and the propagated error in the survey's mean
  rho = d_eff / d (1 cm stems at 1-3 m, and at 4 m), which the density error follows (d ln density ~ -d ln rho).
simulate: isolated vertical rods of one height (raycreate pneumatophores, about 0.3 rods per m^2, Poisson positions)
  scanned along one walk line. For every ray passing a rod at the rod's heights, the horizontal offset x of the ray
  axis from the rod axis at its closest approach, the range r there, and whether the ray returned from the rod. In
  each range bin, d_eff = 2 X (hit fraction of the rays with |x| < X). Then (w0, g, t) or a subset is fitted to
  d_eff(r) of one or several rod diameters by weighted least squares, and the fitted footprint's mean rho for the
  survey compared with the true one.
Requires raycreate on the PATH (or --bin) for simulate. Example:
  scripts/pneumatophore_rod_calibration.py design
  scripts/pneumatophore_rod_calibration.py simulate --out rodtest --bin build/raycloudtools/raycreate
"""
import argparse
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pneumatophore_validation_matrix as vm  # noqa: E402
import rayply_numpy  # noqa: E402

TRUE = (0.01, 0.003, 0.05)
SURVEY_RANGES = np.linspace(1.0, 3.0, 21)  # the ranges at which the survey's stems (0-0.2 m) are mostly seen


def deff(d, r, fp):
    return vm.effective_width(d, fp[0] + fp[1] * r, fp[2])


def mean_rho(fp, d=0.01, ranges=SURVEY_RANGES):
    return float(np.mean([deff(d, r, fp) / d for r in ranges]))


def jacobian(points, fp, names=("w0", "g", "t")):
    """d d_eff / d ln(parameter) at the points (d, r), by central differences"""
    idx = {"w0": 0, "g": 1, "t": 2}
    J = np.zeros((len(points), len(names)))
    for j, n in enumerate(names):
        for s, sign in ((1.002, 1), (1 / 1.002, -1)):
            p = list(fp)
            p[idx[n]] *= s
            J[:, j] += sign * np.array([deff(d, r, p) for d, r in points])
        J[:, j] /= 2 * math.log(1.002)
    return J


def survey_gradient(fp, names=("w0", "g", "t"), ranges=SURVEY_RANGES):
    idx = {"w0": 0, "g": 1, "t": 2}
    out = []
    for n in names:
        vals = []
        for s in (1.002, 1 / 1.002):
            p = list(fp)
            p[idx[n]] *= s
            vals.append(math.log(mean_rho(p, ranges=ranges)))
        out.append((vals[0] - vals[1]) / (2 * math.log(1.002)))
    return np.array(out)


def design(args):
    lines = ["Rod-test design: Fisher analysis of d_eff(r) measurements about (ln w0, ln g, ln t), true footprint",
             f"{TRUE[0] * 1000:g} mm, {TRUE[1] * 1000:g} mrad, t {TRUE[2]:g}; each (rod, range) point measured with",
             f"s.e. {args.sigma * 1000:g} mm in d_eff (about 400 rays across +-3 cm per point, see simulate).",
             "Columns: s.e. of ln w0, ln g, ln t (all three fitted); their correlations; the s.e. of ln(mean rho) for",
             "the survey (1 cm stems, 1-3 m) and at 4 m, all three fitted, and with w0, g fixed at the true values",
             "(t only fitted). 'inf' = not identifiable (the information matrix is singular).", ""]
    lines.append("design | se ln w0 | se ln g | se ln t | corr w0-g, w0-t, g-t | survey rho se (3 fitted) | 4 m | "
                 "survey rho se (t only)")
    designs = [
        ("1 cm at 2 m", [(0.01, 2.0)]),
        ("1 cm at 1, 2, 4 m", [(0.01, r) for r in (1, 2, 4)]),
        ("1 cm at 1-5 m (9 ranges)", [(0.01, r) for r in np.arange(1, 5.01, 0.5)]),
        ("0.5, 1, 2 cm at 2 m", [(d, 2.0) for d in (0.005, 0.01, 0.02)]),
        ("0.5, 1, 2 cm at 1, 2, 4 m", [(d, r) for d in (0.005, 0.01, 0.02) for r in (1, 2, 4)]),
        ("0.5, 1, 2, 4 cm at 1-5 m", [(d, r) for d in (0.005, 0.01, 0.02, 0.04) for r in np.arange(1, 5.01, 0.5)]),
        ("0.3, 0.5, 1 cm at 1, 3, 5 m", [(d, r) for d in (0.003, 0.005, 0.01) for r in (1, 3, 5)]),
        ("1 cm at 1-5 m + 0.5 cm at 1-5 m", [(d, r) for d in (0.005, 0.01) for r in np.arange(1, 5.01, 0.5)]),
        ("1 cm at 0.5-8 m (16 ranges)", [(0.01, r) for r in np.arange(0.5, 8.01, 0.5)]),
    ]
    for name, pts in designs:
        J = jacobian(pts, TRUE) / args.sigma
        F = J.T @ J
        g_s = survey_gradient(TRUE)
        g_4 = survey_gradient(TRUE, ranges=[4.0])
        try:
            if np.linalg.cond(F) > 1e12:
                raise np.linalg.LinAlgError
            C = np.linalg.inv(F)
            se = np.sqrt(np.diag(C))
            corr = C / np.outer(se, se)
            s_sur = math.sqrt(g_s @ C @ g_s)
            s_4 = math.sqrt(g_4 @ C @ g_4)
            txt = (f"{se[0]:.3f} | {se[1]:.3f} | {se[2]:.3f} | {corr[0, 1]:+.3f} {corr[0, 2]:+.3f} {corr[1, 2]:+.3f} | "
                   f"{s_sur:.4f} | {s_4:.4f}")
        except np.linalg.LinAlgError:
            txt = "inf | inf | inf | - | - | -"
        Jt = J[:, 2:3]
        Ct = np.linalg.inv(Jt.T @ Jt)
        s_t = math.sqrt(g_s[2] ** 2 * Ct[0, 0])
        lines.append(f"{name} | {txt} | {s_t:.4f}")
    # the degeneracy directions, at one range: d ln d_eff / d ln parameter of 1 cm stems, against range and diameter
    lines.append("")
    lines.append("Local sensitivities d ln d_eff / d ln parameter (rows: rod diameter, range); a direction is degenerate")
    lines.append("where two columns stay proportional over the design's points")
    lines.append("d (cm) | r (m) | d_eff (mm) | w0 | g | t | t-equivalent of a 10% w0 error (t ratio)")
    for d in (0.005, 0.01, 0.02, 0.04):
        for r in (1, 2, 4):
            J = jacobian([(d, r)], TRUE)[0] / deff(d, r, TRUE)
            lines.append(f"{d * 100:g} | {r} | {deff(d, r, TRUE) * 1000:.2f} | {J[0]:+.3f} | {J[1]:+.3f} | "
                         f"{J[2]:+.3f} | {math.exp(-J[0] * math.log(1.1) / J[2]):.3f}")
    text = "\n".join(lines) + "\n"
    print(text)
    if args.out:
        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, "rod_design.txt"), "w") as f:
            f.write(text)


def extract(work, d, height, half_width):
    starts, ends, _, colours = rayply_numpy.read_rays(os.path.join(work, "pneumatophores.ply"))
    stems = np.loadtxt(os.path.join(work, "pneumatophores_stems.csv"), delimiter=",", skiprows=1, ndmin=2)
    V = ends - starts
    Vh = V[:, :2]
    vh2 = np.einsum("ij,ij->i", Vh, Vh)
    length = np.linalg.norm(V, axis=1)
    bounded = colours[:, 3] > 0
    out = []
    for sx, sy, h, _ in stems:
        P = np.array([sx, sy])
        rel = P[None, :] - starts[:, :2]
        with np.errstate(invalid="ignore", divide="ignore"):
            t = np.einsum("ij,ij->i", rel, Vh) / vh2
        cross = np.abs(rel[:, 0] * Vh[:, 1] - rel[:, 1] * Vh[:, 0]) / np.sqrt(vh2)
        near = (cross < half_width) & (t > 0) & (vh2 > 1e-12)
        if not near.any():
            continue
        idx = np.nonzero(near)[0]
        z = starts[idx, 2] + t[idx] * V[idx, 2]
        ok = (z > 0.03) & (z < h - 0.03)
        idx, z = idx[ok], z[ok]
        end_off = np.linalg.norm(ends[idx, :2] - P[None, :], axis=1)
        hit = bounded[idx] & (end_off < 0.5 * d + 0.04) & (ends[idx, 2] < h + 0.01)
        # rays that end before reaching the rod (at another rod) carry no information about this one
        reach = hit | (t[idx] <= 1.0)
        idx, hit = idx[reach], hit[reach]
        r = t[idx] * length[idx]
        out.append(np.stack([r, cross[idx], hit.astype(float)], axis=1))
    return np.concatenate(out) if out else np.zeros((0, 3))


def bin_deff(rec, half_width, edges):
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (rec[:, 0] >= lo) & (rec[:, 0] < hi)
        n = int(m.sum())
        if n < 50:
            continue
        f = rec[m, 2].mean()
        rows.append((rec[m, 0].mean(), 2 * half_width * f, 2 * half_width * math.sqrt(max(f * (1 - f), 1e-4) / n),
                     n))
    return rows


def fit(data, start, names):
    """weighted least squares of (w0, g, t) (those named; the others fixed at start) to d_eff points
    [(d, r, deff, se)], Gauss-Newton with step halving in the log parameters"""
    idx = {"w0": 0, "g": 1, "t": 2}
    p = list(start)

    def resid(q):
        return np.array([(deff(d, r, q) - y) / s for d, r, y, s in data])

    res = resid(p)
    for _ in range(50):
        J = jacobian([(d, r) for d, r, _, _ in data], p, names) / np.array([[s] for _, _, _, s in data])
        F = J.T @ J
        step = np.linalg.solve(F + 1e-9 * np.eye(len(names)), -J.T @ res)
        lam = 1.0
        while lam > 1e-4:
            q = list(p)
            for j, n in enumerate(names):
                q[idx[n]] = p[idx[n]] * math.exp(lam * step[j])
            rq = resid(q)
            if rq @ rq < res @ res:
                break
            lam *= 0.5
        if lam <= 1e-4:
            break
        p, res = q, rq
        if np.max(np.abs(lam * step)) < 1e-6:
            break
    J = jacobian([(d, r) for d, r, _, _ in data], p, names) / np.array([[s] for _, _, _, s in data])
    try:
        C = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        C = np.full((len(names), len(names)), np.nan)
    return p, C, float(res @ res), len(data) - len(names)


def simulate(args):
    if args.bin:
        os.environ["PATH"] = os.pathsep.join(os.path.abspath(b) for b in args.bin.split(",")) + os.pathsep + \
            os.environ["PATH"]
    os.makedirs(args.out, exist_ok=True)
    fp = vm.fp_text(*TRUE)
    edges = np.arange(1.0, 5.51, 0.5)
    data = {}
    lines = [f"Synthetic rod test: rods of one height (1 m), true footprint {fp}, {args.rays} rays per scene along",
             f"one 8 m walk line, seeds {args.seeds}, range noise {args.range_noise} m. d_eff per 0.5 m range bin from",
             f"the hit fraction of rays with |offset| < {args.half_width * 100:g} cm (s.e. binomial), against the model.",
             ""]
    lines.append("rod d (cm) | range (m) | rays | d_eff measured (mm) | model (mm) | (meas - model) / se")
    for d in (0.005, 0.01, 0.02):
        recs = []
        for seed in args.seeds.split(","):
            work = tempfile.mkdtemp(prefix=f"rod_{d}_{seed}_", dir=args.out)
            noise = f"--range_noise {args.range_noise}" if args.range_noise > 0 else ""
            cmd = (f"cd {work} && raycreate pneumatophores {seed} --stems_per_m2 0.3 --overlapping --diameter {d} "
                   f"--height 1.0,0 --walk_spacing 8 --plot_width 8 --rays {args.rays} --footprint {fp} {noise}")
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)
            recs.append(extract(work, d, 1.0, args.half_width))
            if not args.keep:
                for name in os.listdir(work):
                    os.remove(os.path.join(work, name))
                os.rmdir(work)
        rec = np.concatenate(recs)
        rows = bin_deff(rec, args.half_width, edges)
        data[d] = [(d, r, y, s) for r, y, s, n in rows]
        for r, y, s, n in rows:
            m = deff(d, r, TRUE)
            lines.append(f"{d * 100:g} | {r:.2f} | {n} | {y * 1000:.2f} ± {s * 1000:.2f} | {m * 1000:.2f} | "
                         f"{(y - m) / s:+.2f}")
    lines.append("")
    lines.append("Fits (start: the specified footprint, or a wrong one where stated). rho err = fitted / true mean rho")
    lines.append("for the survey (1 cm stems at 1-3 m), and at 4 m; the density error is about minus these.")
    lines.append("case | w0 (mm) | g (mrad) | t | se ln w0, ln g, ln t | chi2/dof | rho err survey | rho err 4 m")
    true_rho, true_rho4 = mean_rho(TRUE), mean_rho(TRUE, ranges=[4.0])
    cases = [
        ("t only, 1 cm rod, w0 g as specified", data[0.01], TRUE, ("t",)),
        ("t only, 1 cm rod, w0 specified 25% high", data[0.01], (TRUE[0] * 1.25, TRUE[1], TRUE[2]), ("t",)),
        ("t only, 1 cm rod, g specified 50% high", data[0.01], (TRUE[0], TRUE[1] * 1.5, TRUE[2]), ("t",)),
        ("t only, 1 cm rod at 1.5-2.5 m only, g 50% high",
         [x for x in data[0.01] if 1.5 <= x[1] < 2.5], (TRUE[0], TRUE[1] * 1.5, TRUE[2]), ("t",)),
        ("w0, g, t, 1 cm rod", data[0.01], (TRUE[0] * 1.3, TRUE[1] * 0.7, TRUE[2] * 1.5), ("w0", "g", "t")),
        ("w0, g, t, 0.5 + 1 + 2 cm rods", data[0.005] + data[0.01] + data[0.02],
         (TRUE[0] * 1.3, TRUE[1] * 0.7, TRUE[2] * 1.5), ("w0", "g", "t")),
        ("g, t, 0.5 + 1 + 2 cm rods, w0 specified 25% high", data[0.005] + data[0.01] + data[0.02],
         (TRUE[0] * 1.25, TRUE[1], TRUE[2]), ("g", "t")),
    ]
    for name, pts, start, names in cases:
        if len(pts) <= len(names):
            continue
        p, C, chi2, dof = fit(pts, start, names)
        se = {n: math.sqrt(C[j, j]) if C[j, j] >= 0 else float("nan") for j, n in enumerate(names)}
        se_txt = ", ".join(f"{se[n]:.3f}" if n in se else "fixed" for n in ("w0", "g", "t"))
        lines.append(f"{name} | {p[0] * 1000:.2f} | {p[1] * 1000:.2f} | {p[2]:.4f} | {se_txt} | "
                     f"{chi2 / max(dof, 1):.2f} | {mean_rho(p) / true_rho:.4f} | "
                     f"{mean_rho(p, ranges=[4.0]) / true_rho4:.4f}")
    text = "\n".join(lines) + "\n"
    print(text)
    with open(os.path.join(args.out, "rod_simulation.txt"), "w") as f:
        f.write(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("design")
    p.add_argument("--sigma", type=float, default=0.0005, help="s.e. of each d_eff measurement (m)")
    p.add_argument("--out", default="")
    p = sub.add_parser("simulate")
    p.add_argument("--out", default="rodtest")
    p.add_argument("--bin", default="")
    p.add_argument("--seeds", default="201,202")
    p.add_argument("--rays", type=int, default=4000000)
    p.add_argument("--range_noise", type=float, default=0.0)
    p.add_argument("--half_width", type=float, default=0.04)
    p.add_argument("--keep", action="store_true")
    args = parser.parse_args()
    return design(args) if args.cmd == "design" else simulate(args)


if __name__ == "__main__":
    sys.exit(main())
