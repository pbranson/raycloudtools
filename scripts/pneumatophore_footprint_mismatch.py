#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Analysis of block L (footprint mismatch) of the validation matrix (scripts/pneumatophore_validation_matrix.py
--blocks L --out RUN): scenes generated with a true beam footprint (exit diameter w0, divergence g, threshold t) and
estimated assuming another footprint (and, for the default footprint, another stem diameter). Writes to RUN:
  footprint_tables.txt    d/t and d/t90 per layer, mean +- s.e. over seeds, and the ratio to the true-footprint
                          estimate (paired by seed), for every true footprint, density, estimator and assumed one
  footprint_sens.txt      d ln(density) / d ln(parameter) for w0, g, t (local: the 0.75/1.25 levels, and the
                          0.03/0.08 thresholds; global: a fit over all levels), and the parameter errors that keep
                          the density within +-2% and +-5%
  footprint_deff.txt      the d_eff mapping: density ratio against the ratio Q of the mean effective width
                          (weighted path of the assumed render / that of the true render, which is the path-weighted
                          mean of rho = d_eff / d over the cell's rays), one curve fitted over all assumed footprints
                          and diameters, with residuals per kind of error
  footprint_grid.txt      assumed diameter x assumed threshold, and the threshold that a counted-plot calibration
                          would pick for each diameter error, and how it transfers to the other density and range
  footprint_range.txt     the single-line 8 m plot: ratio to the true-footprint estimate by distance from the line
  footprint_rows.csv      all (true, density, estimator, assumed, layer) rows
Example: scripts/pneumatophore_footprint_mismatch.py RUN
"""
import argparse
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pneumatophore_validation_matrix as vm  # noqa: E402

ESTIMATORS = ("base", "ss", "sw")


def mean_se(v):
    v = [x for x in v if x is not None and not math.isnan(x)]
    if not v:
        return float("nan"), float("nan")
    m = sum(v) / len(v)
    if len(v) < 2:
        return m, float("nan")
    return m, math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1) / len(v))


def fmt(ms, digits=3):
    m, s = ms
    if math.isnan(m):
        return "  n/a  "
    return f"{m:.{digits}f}±{s:.{digits}f}" if not math.isnan(s) else f"{m:.{digits}f}"


def load(run):
    configs = [c for c in vm.build_configs() if c.block == "L"]
    by_ident = {c.ident(): c for c in configs}
    jobs = {}
    root = os.path.join(run, "jobs", "L")
    for name in sorted(os.listdir(root)):
        if not name.endswith(".json"):
            continue
        with open(os.path.join(root, name)) as f:
            r = json.load(f)
        if r.get("error") or r["ident"] not in by_ident:
            continue
        jobs.setdefault(r["ident"], []).append(r)
    return configs, jobs


def sums(cells, good_only=False):
    """sum truth, density, weighted path over the cells"""
    if good_only:
        cells = [c for c in cells if c[7] >= 0.9]
    t = sum(c[4] for c in cells)
    d = sum(vm.nan0(c[5]) for c in cells)
    wp = sum(c[9] for c in cells)
    return t, d, wp, len(cells)


def index_cells(job):
    out = {}
    for c in job["cells"]:
        out.setdefault((c[0], c[1]), []).append(c)
    return out


def seed_values(job_cells, variant, ref, z, good_only, select=None):
    """d/t of the variant, its ratio to the reference variant's estimate, and Q (weighted path ratio), over the
    cells (paired: the same cells for both)"""
    a = job_cells.get((variant, z), [])
    b = job_cells.get((ref, z), [])
    bk = {(c[2], c[3]): c for c in b}
    pairs = [(c, bk[(c[2], c[3])]) for c in a if (c[2], c[3]) in bk]
    if select:
        pairs = [p for p in pairs if select(p[1])]
    if good_only:
        pairs = [p for p in pairs if p[1][7] >= 0.9]  # cells observed >= 0.9 (by the true-footprint render)
    if not pairs:
        return None
    t = sum(p[0][4] for p in pairs)
    da = sum(vm.nan0(p[0][5]) for p in pairs)
    db = sum(vm.nan0(p[1][5]) for p in pairs)
    wa = sum(p[0][9] for p in pairs)
    wb = sum(p[1][9] for p in pairs)
    if t <= 0 or db <= 0 or wb <= 0:
        return None
    return dict(dt=da / t, ratio=da / db, q=wa / wb, n=len(pairs))


def true_name(cfg):
    return cfg.label.split()[1]


def density_of(cfg):
    return int(cfg.label.split(", ")[-1].split("/")[0])


def build_rows(configs, jobs):
    rows = []
    for cfg in configs:
        runs = jobs.get(cfg.ident(), [])
        if not runs or cfg.line_bins:
            continue
        indexed = [index_cells(r) for r in runs]
        for name, _ in cfg.variants:
            meta = cfg.meta[name]
            ref = f"true {meta['est']}"
            for z in cfg.layers:
                row = dict(true=true_name(cfg), stems=density_of(cfg), est=meta["est"], tag=meta["tag"],
                           kind=meta["kind"], level=meta["level"], exit=meta["exit"], divergence=meta["divergence"],
                           threshold=meta["threshold"], diameter=meta["diameter"], z=z)
                for good, suffix in ((False, ""), (True, "90")):
                    vals = [seed_values(ic, name, ref, z, good) for ic in indexed]
                    vals = [v for v in vals if v]
                    row["seeds"] = len(vals)
                    for k in ("dt", "ratio", "q"):
                        row[k + suffix], row[k + suffix + "_se"] = mean_se([v[k] for v in vals])
                    row["per_seed" + suffix] = vals
                rows.append(row)
    return rows


def find(rows, **kw):
    return [r for r in rows if all(r[k] == v for k, v in kw.items())]


def one(rows, **kw):
    f = find(rows, **kw)
    return f[0] if f else None


def tables(rows, out):
    lines = ["Block L: footprint mismatch. d/t (all cells) and d/t90 (cells observed >= 0.9), mean +- s.e. over seeds;",
             "r90 = the estimate / the estimate with the true footprint (same estimator, same cells observed >= 0.9),",
             "paired by seed; Q90 = mean rho assumed / mean rho true over those cells (weighted path ratio).", ""]
    for true in ("default", "narrow", "wide"):
        for n in (300, 1000):
            sub = find(rows, true=true, stems=n)
            if not sub:
                continue
            tags = []
            for r in sub:
                if r["tag"] not in tags and r["kind"] != "grid":
                    tags.append(r["tag"])
            fpt = one(sub, tag="true", est="base", z=0.0)
            lines.append(f"True {true} ({fpt['exit'] * 1000:g} mm, {fpt['divergence'] * 1000:g} mrad, "
                         f"t {fpt['threshold']:g}), {n}/m2, {fpt['seeds']} seeds")
            for z in (0.0, 0.1):
                head = f"  layer {z:.1f}-{z + 0.1:.1f} m | assumed | Q90 |"
                for e in ESTIMATORS:
                    head += f" {e} d/t | {e} d/t90 | {e} r90 |"
                lines.append(head)
                for tag in tags:
                    rb = one(sub, tag=tag, est="base", z=z)
                    if rb is None:
                        continue
                    if tag == "nofp":
                        desc = "no footprint"
                    else:
                        desc = (f"{rb['exit'] * 1000:.3g} mm, {rb['divergence'] * 1000:.3g} mrad, "
                                f"t {rb['threshold']:.3g}")
                    line = f"  {tag:8s} | {desc:28s} | {fmt((rb['q90'], rb['q90_se']))} |"
                    for e in ESTIMATORS:
                        r = one(sub, tag=tag, est=e, z=z)
                        line += (f" {fmt((r['dt'], r['dt_se']))} | {fmt((r['dt90'], r['dt90_se']))} | "
                                 f"{fmt((r['ratio90'], r['ratio90_se']))} |")
                    lines.append(line)
            lines.append("")
    with open(os.path.join(out, "footprint_tables.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return lines


def slope(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx if sxx > 0 else float("nan")


def crossing(points, target):
    """the parameter ratios (below and above 1) where the density ratio crosses target (points: (level, ratio),
    sorted), by linear interpolation in log-log; None where it is not crossed within the levels"""
    lo = hi = None
    pts = sorted(points)
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
        for t in target:
            if (y0 - t) * (y1 - t) <= 0 and y0 != y1:
                lx = math.log(x0) + (math.log(t) - math.log(y0)) * (math.log(x1) - math.log(x0)) / \
                    (math.log(y1) - math.log(y0))
                x = math.exp(lx)
                if x < 1 and (lo is None or x > lo):
                    lo = x
                if x > 1 and (hi is None or x < hi):
                    hi = x
    return lo, hi


def sensitivities(rows, out):
    lines = ["Sensitivity of the density to each footprint parameter: S = d ln(density) / d ln(parameter), from the",
             "paired ratio r90 (estimate / true-footprint estimate, cells observed >= 0.9; layer 0.1-0.2 m where 0-0.1",
             "has none). local: the x0.75 and x1.25 levels (threshold: 0.03 and 0.08); global: least squares over all",
             "levels and the true one. tol 2% / 5%: the range of assumed/true parameter ratios within which",
             "|r90 - 1| stays within 2% / 5% (interpolated; '<' / '>' = not reached within the levels).", ""]
    lines.append("true | stems | est | layer | S(w0) local global | S(g) local global | S(t) local global | "
                 "w0 2% | w0 5% | g 2% | g 5% | t 2% | t 5%")
    results = []
    for true in ("default", "narrow", "wide"):
        for n in (300, 1000):
            for e in ESTIMATORS:
                for z in (0.0, 0.1):
                    rt = one(rows, true=true, stems=n, est=e, tag="true", z=z)
                    if rt is None or math.isnan(rt["ratio90"]):
                        continue
                    t_true = rt["threshold"]
                    line = f"{true} | {n} | {e} | {z:.1f} |"
                    res = dict(true=true, stems=n, est=e, z=z)
                    tol = []
                    for kind, key in (("exit", "w0"), ("divergence", "g"), ("threshold", "t")):
                        pts = [(1.0, 1.0)]
                        for r in find(rows, true=true, stems=n, est=e, kind=kind, z=z):
                            if not math.isnan(r["ratio90"]):
                                pts.append((r["level"], r["ratio90"]))
                        pts.sort()
                        d = dict(pts)
                        if kind == "threshold":
                            a, b = 0.03 / t_true, 0.08 / t_true
                        else:
                            a, b = 0.75, 1.25
                        loc = (math.log(d[b]) - math.log(d[a])) / (math.log(b) - math.log(a)) \
                            if a in d and b in d else float("nan")
                        glob = slope([math.log(p[0]) for p in pts], [math.log(p[1]) for p in pts])
                        line += f" {loc:+.3f} {glob:+.3f} |"
                        res[key + "_local"], res[key + "_global"] = loc, glob
                        for tl in (0.02, 0.05):
                            lo, hi = crossing(pts, (1 - tl, 1 + tl))
                            mn, mx = min(p[0] for p in pts), max(p[0] for p in pts)
                            if kind == "threshold":
                                s_lo = f"{lo * t_true:.3f}" if lo else f"<{mn * t_true:.2f}"
                                s_hi = f"{hi * t_true:.3f}" if hi else f">{mx * t_true:.2f}"
                            else:
                                s_lo = f"{lo:.3f}" if lo else f"<{mn:.2f}"
                                s_hi = f"{hi:.3f}" if hi else f">{mx:.2f}"
                            tol.append(f"{s_lo}-{s_hi}")
                            res[f"{key}_{int(tl * 100)}"] = (lo, hi)
                    line += " " + " | ".join(tol)
                    lines.append(line)
                    results.append(res)
    lines.append("")
    lines.append("Threshold tolerances are absolute thresholds (true 0.05).")
    with open(os.path.join(out, "footprint_sens.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return lines, results


def fit_quadratic(xs, ys):
    """least squares y = b1 x + b2 x^2 (through the origin)"""
    s11 = sum(x * x for x in xs)
    s12 = sum(x ** 3 for x in xs)
    s22 = sum(x ** 4 for x in xs)
    t1 = sum(x * y for x, y in zip(xs, ys))
    t2 = sum(x * x * y for x, y in zip(xs, ys))
    det = s11 * s22 - s12 * s12
    return (t1 * s22 - t2 * s12) / det, (s11 * t2 - s12 * t1) / det


def deff_mapping(rows, out):
    lines = ["d_eff mapping. Per seed and variant, x = ln Q (Q = mean rho assumed / mean rho true over the same cells,",
             "from the weighted paths, i.e. the error in the path-weighted mean d_eff / d at the scanned ranges) and",
             "y = ln r (density ratio to the true-footprint estimate), cells observed >= 0.9 (all cells for layers",
             "without). One curve y = b1 x + b2 x^2 fitted per (true footprint, density, estimator, layer) over all",
             "single-parameter errors (exit, divergence, threshold); then the residuals of every kind, including the",
             "joint (range-degenerate) cases and the assumed-diameter grid, in % of density.", ""]
    lines.append("true | stems | est | layer | b1 | b2 | n | rms resid fit % | resid % by kind (mean, max abs): exit | "
                 "divergence | threshold | joint | grid")
    summary = []
    groups = {}
    for r in rows:
        if r["tag"] in ("true", "nofp"):
            continue
        key = (r["true"], r["stems"], r["est"], r["z"])
        vals = r["per_seed90"] if r["per_seed90"] else r["per_seed"]
        for v in vals:
            if v["q"] > 0 and v["ratio"] > 0:
                groups.setdefault(key, []).append((r["kind"], r["tag"], math.log(v["q"]), math.log(v["ratio"])))
    for key in sorted(groups, key=lambda k: (("default", "narrow", "wide").index(k[0]), k[1], k[2], k[3])):
        pts = groups[key]
        fit = [p for p in pts if p[0] in ("exit", "divergence", "threshold")]
        if len(fit) < 3:
            continue
        b1, b2 = fit_quadratic([p[2] for p in fit], [p[3] for p in fit])
        res = {}
        for kind, tag, x, y in pts:
            res.setdefault(kind, []).append(100.0 * (math.exp(y - (b1 * x + b2 * x * x)) - 1.0))
        rms = math.sqrt(sum(v * v for k in ("exit", "divergence", "threshold") for v in res.get(k, [])) /
                        sum(len(res.get(k, [])) for k in ("exit", "divergence", "threshold")))
        line = f"{key[0]} | {key[1]} | {key[2]} | {key[3]:.1f} | {b1:+.3f} | {b2:+.3f} | {len(fit)} | {rms:.2f} |"
        for kind in ("exit", "divergence", "threshold", "joint", "grid"):
            v = res.get(kind)
            line += f" {sum(v) / len(v):+.2f} {max(abs(x) for x in v):.2f} |" if v else " - |"
        lines.append(line)
        summary.append(dict(key=key, b1=b1, b2=b2, rms=rms, res=res))
    # the joint cases in detail
    lines.append("")
    lines.append("Joint cases (beam diameter or d_eff matched at 2 m): Q90 and r90 (mean over seeds), and r predicted by")
    lines.append("the curve from Q, for the ss estimator, layer 0-0.1 m (0.1-0.2 where 0-0.1 has no cells >= 0.9)")
    lines.append("true | stems | tag | assumed | Q90 | r90 | predicted | d_eff ratio at 1 / 2 / 4 m (model)")
    fits = {s["key"]: s for s in summary}
    for r in rows:
        if r["kind"] != "joint" or r["est"] != "ss":
            continue
        z = r["z"]
        if z != 0.0 and not math.isnan(one(rows, true=r["true"], stems=r["stems"], est="ss", tag=r["tag"],
                                               z=0.0)["ratio90"]):
            continue
        if z == 0.0 and math.isnan(r["ratio90"]):
            continue
        s = fits.get((r["true"], r["stems"], "ss", z))
        if s is None:
            continue
        x = math.log(r["q90"])
        pred = math.exp(s["b1"] * x + s["b2"] * x * x)
        rt = one(rows, true=r["true"], stems=r["stems"], est="ss", tag="true", z=z)
        model = [vm.effective_width(0.01, r["exit"] + r["divergence"] * rr, r["threshold"]) /
                 vm.effective_width(0.01, rt["exit"] + rt["divergence"] * rr, rt["threshold"]) for rr in (1, 2, 4)]
        lines.append(f"{r['true']} | {r['stems']} | {r['tag']} | {r['exit'] * 1000:.3g} mm {r['divergence'] * 1000:.3g}"
                     f" mrad t {r['threshold']:.4g} | {r['q90']:.4f} | {r['ratio90']:.4f} | {pred:.4f} | "
                     + " / ".join(f"{m:.3f}" for m in model))
    with open(os.path.join(out, "footprint_deff.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return lines, summary


def grid(rows, out):
    lines = ["Assumed diameter x assumed threshold (true default footprint 10 mm, 3 mrad, t 0.05; true d 1 cm).",
             "d/t90 per estimator, mean over seeds (layer 0-0.1 m at 300/m2, both layers listed).", ""]
    dts = (0.8, 0.9, 1.0, 1.1, 1.2)
    ths = vm.L_GRID_THRESHOLDS
    calib = {}
    for n in (300, 1000):
        for e in ESTIMATORS:
            for z in (0.0, 0.1):
                lines.append(f"{n}/m2, {e}, layer {z:.1f}: rows assumed d / true d, columns assumed threshold")
                lines.append("  d ratio | " + " | ".join(f"t {t:g}" for t in ths) + " | t* (d/t90 = 1)")
                for dr in dts:
                    vals = []
                    for t in ths:
                        if dr == 1.0:
                            tag = "true" if t == 0.05 else f"t{t:g}"
                        else:
                            tag = f"d{dr:g} t{t:g}"
                        r = one(rows, true="default", stems=n, est=e, tag=tag, z=z)
                        vals.append(r["dt90"] if r else float("nan"))
                    # the threshold a counted-plot calibration would choose: d/t90 = 1, interpolated in ln t
                    tstar = None
                    pairs = list(zip(ths, vals))
                    for (t0, v0), (t1, v1) in zip(pairs[:-1], pairs[1:]):
                        if not (math.isnan(v0) or math.isnan(v1)) and (v0 - 1) * (v1 - 1) <= 0 and v0 != v1:
                            lt = math.log(t0) + (1 - v0) * (math.log(t1) - math.log(t0)) / (v1 - v0)
                            tstar = math.exp(lt)
                            break
                    calib[(n, e, z, dr)] = (tstar, vals)
                    ts = f"{tstar:.4f}" if tstar else "-"
                    lines.append(f"  {dr:7g} | " + " | ".join(f"{v:.3f}" for v in vals) + f" | {ts}")
                lines.append("")
    # transfer: calibrate t at 300/m2 (layer 0-0.1, cells >= 0.9), apply at 1000/m2 and at other ranges
    lines.append("Transfer of a counted-plot threshold calibration (t* chosen so that d/t90 = 1 at 300/m2, layer 0-0.1):")
    lines.append("d/t90 at 1000/m2 with t* (interpolated in ln t), layer 0.1-0.2 (where 1000/m2 is observed >= 0.9),")
    lines.append("and the model ratio of the assumed rho to the true rho at 1, 2 and 4 m with (d', t*):")
    lines.append("est | d ratio | t* | 300/m2 0.1-0.2 | 1000/m2 0-0.1 | 1000/m2 0.1-0.2 | rho ratio 1 / 2 / 4 m | "
                 "rho ratio for 2 cm stems at 1 / 2 / 4 m")
    for e in ESTIMATORS:
        for dr in dts:
            tstar = calib.get((300, e, 0.0, dr), (None,))[0]
            if not tstar:
                continue
            vals = []
            for n, z in ((300, 0.1), (1000, 0.0), (1000, 0.1)):
                v = calib[(n, e, z, dr)][1]
                y = interp_log(ths, v, tstar)
                vals.append(y)
            rr = []
            for d_true in (0.01, 0.02):
                for rng in (1, 2, 4):
                    w = 0.01 + 0.003 * rng
                    rho_a = vm.effective_width(d_true * dr, w, tstar) / (d_true * dr)
                    rho_t = vm.effective_width(d_true, w, 0.05) / d_true
                    rr.append(rho_a / rho_t)
            lines.append(f"{e} | {dr:g} | {tstar:.4f} | " + " | ".join(f"{v:.3f}" for v in vals) + " | " +
                         " / ".join(f"{x:.3f}" for x in rr[:3]) + " | " + " / ".join(f"{x:.3f}" for x in rr[3:]))
    with open(os.path.join(out, "footprint_grid.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return lines


def interp_log(xs, ys, x):
    lx = math.log(x)
    for (x0, y0), (x1, y1) in zip(list(zip(xs, ys))[:-1], list(zip(xs, ys))[1:]):
        if x0 <= x <= x1:
            if math.isnan(y0) or math.isnan(y1):
                return float("nan")
            return y0 + (y1 - y0) * (lx - math.log(x0)) / (math.log(x1) - math.log(x0))
    return float("nan")


def range_table(configs, jobs, out):
    cfg = [c for c in configs if c.line_bins]
    lines = ["Range dependence: one walk line along the middle of an 8 m plot (true default footprint, 300/m2).",
             "Per distance of the cell centre from the line (m): r = estimate / true-footprint estimate (all cells of",
             "that distance, paired by seed; the far cells are poorly observed, so d/t itself is not the point), and",
             "Q = rho ratio from the weighted paths. Layer 0-0.1 m (0.1-0.2 m in brackets).", ""]
    if not cfg or cfg[0].ident() not in jobs:
        return lines
    cfg = cfg[0]
    runs = jobs[cfg.ident()]
    indexed = [index_cells(r) for r in runs]
    dists = sorted(set(round(c[8], 2) for c in runs[0]["cells"]))
    tags = []
    for name, _ in cfg.variants:
        if cfg.meta[name]["tag"] not in tags:
            tags.append(cfg.meta[name]["tag"])
    for e in ESTIMATORS:
        lines.append(f"estimator {e}: d/t of the true footprint, then r (Q) of each assumed footprint")
        lines.append("  assumed | " + " | ".join(f"{d:.1f} m" for d in dists))
        for tag in tags:
            name = f"{tag} {e}"
            line = f"  {tag:8s} |"
            for d in dists:
                parts = []
                for z in (0.0, 0.1):
                    vals = [seed_values(ic, name, f"true {e}", z, False, select=lambda c, d=d: round(c[8], 2) == d)
                            for ic in indexed]
                    vals = [v for v in vals if v]
                    if tag == "true":
                        parts.append(f"{mean_se([v['dt'] for v in vals])[0]:.3f}")
                    else:
                        parts.append(f"{mean_se([v['ratio'] for v in vals])[0]:.3f} "
                                     f"({mean_se([v['q'] for v in vals])[0]:.3f})")
                line += f" {parts[0]} [{parts[1]}] |"
            lines.append(line)
        lines.append("")
    obs = {}
    for r in runs:
        for c in r["cells"]:
            if c[0] == "true ss":
                obs.setdefault((round(c[8], 2), c[1]), []).append(c[7])
    lines.append("mean observed fraction by distance (layer 0-0.1 / 0.1-0.2): " +
                 ", ".join(f"{d:.1f} m {sum(obs[(d, 0.0)]) / len(obs[(d, 0.0)]):.2f}/"
                           f"{sum(obs[(d, 0.1)]) / len(obs[(d, 0.1)]):.2f}" for d in dists))
    with open(os.path.join(out, "footprint_range.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run", help="the --out directory of the validation matrix run of block L")
    args = parser.parse_args()
    configs, jobs = load(args.run)
    rows = build_rows(configs, jobs)
    with open(os.path.join(args.run, "footprint_rows.csv"), "w", newline="") as f:
        keys = [k for k in rows[0] if not k.startswith("per_seed")]
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    for lines in (tables(rows, args.run), sensitivities(rows, args.run)[0], deff_mapping(rows, args.run)[0],
                  grid(rows, args.run), range_table(configs, jobs, args.run)):
        print("\n".join(lines))
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
