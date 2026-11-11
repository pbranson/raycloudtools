#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Self-shadowing bias of the pooled stem density estimator, and its correction (rayrender --self_shadow_correction),
on small synthetic hand-held scans of pneumatophore fields. Requires raycreate and rayrender on the PATH (or --bin).

Each case is a scene (raycreate pneumatophores options) and a set of render variants (rayrender options), evaluated
over several seeds in 1 x 1 m cells per height layer. Per case, variant and layer it reports, as mean +- standard error
over the seeds:
  density/truth    the main output (mean of the voxel densities), summed over the plot's cells against the truth
  pooled/truth     the pooled cell density (hits / path over the cell's voxels)
  mean/pooled      their ratio (for uniform stems, the self-shadow bias without the occlusion correction's errors)
  observed         the mean observed fraction of the cells' voxels
  H1               the first-order prediction 1 + s d_eff / W^2 (s = 0.4732 W for the square window of width W, d_eff
                   the mean effective width reported by rayrender)
  obs>=0.9         density/truth over the cells whose observed fraction is at least 0.9 (and the mean number of them)
Writes <out>/rows.csv (one row per seed, with the observed-fraction bins), <out>/summary.txt and <out>/compact.txt
(case x variant tables). Example:
  scripts/pneumatophore_self_shadow_study.py --sets density,window --seeds 101,102,103 --out study
The small-window sets: K (voxel x pool radius, with and without the first-order correction, at 1/4, 1/2 and 1 of
the standard ray density), K2 (the same scenes for the small-window configuration SMALL_WINDOW and its ablations),
K3 (stem diameter and footprint) and K4 (assumed diameter +-20%), e.g.
  scripts/pneumatophore_self_shadow_study.py --sets K2,K3,K4 --seeds 101,102,103,104 --jobs 2 --out small_window
"""
import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

FOOTPRINT = "0.01,0.003,0.05"  # the default scanner model: 1 cm exit diameter, 3 mrad, 5% threshold
WIDE = "0.02,0.005,0.05"  # about 3 cm effective width for 1 cm stems at 1-2 m
EXIT_DISTANCE_SQUARE = (4.0 * math.log(1.0 + math.sqrt(2.0)) - (4.0 / 3.0) * (math.sqrt(2.0) - 1.0)) / (2.0 * math.pi)


def render_args(footprint, diameter, voxel=0.05, pool=0.05, extra=""):
    fp = f"--footprint {footprint} " if footprint else ""
    return (f"--estimator pooled --stem_diameter {diameter} {fp}--pixel_width {voxel} --pool_radius {pool} {extra}")


def variants(footprint=FOOTPRINT, diameter=0.01, voxel=0.05, pool=0.05, extra=""):
    """the uncorrected estimator and the self-shadow correction"""
    base = render_args(footprint, diameter, voxel, pool, extra)
    return [("base", base), ("ss", base + " --self_shadow_correction")]


# the small-window configuration (rayrender options added to --pool_radius 0)
SMALL_WINDOW = ("--self_shadow_correction --self_shadow_second_order --occlusion_radius 0.05 --ratio_bias_correction "
                "--no_neighbour_priors")


def small_window_variants(footprint=FOOTPRINT, diameter=0.01, short=False):
    """the default window (with and without the first-order correction), the first-order correction with single
    voxels, the small-window configuration (sw) and, unless short, its ablations"""
    default = render_args(footprint, diameter, 0.05, 0.05)
    single = render_args(footprint, diameter, 0.05, 0.0)
    v = [("default", default), ("default ss", default + " --self_shadow_correction"),
         ("p0 ss", single + " --self_shadow_correction"), ("sw", single + " " + SMALL_WINDOW)]
    if not short:
        v += [("sw first order", single + " " + SMALL_WINDOW.replace("--self_shadow_second_order ", "")),
              ("sw coupled", single + " " + SMALL_WINDOW.replace("--occlusion_radius 0.05 ", "")),
              ("sw no ratio", single + " " + SMALL_WINDOW.replace("--ratio_bias_correction ", "")),
              ("sw priors", single + " " + SMALL_WINDOW.replace(" --no_neighbour_priors", "")),
              ("sw v0.025p0.025", render_args(footprint, diameter, 0.025, 0.025) + " " + SMALL_WINDOW)]
    return v


def scene(stems, footprint=FOOTPRINT, extra=""):
    fp = f"--footprint {footprint} " if footprint else ""
    return f"--stems_per_m2 {stems} {fp}{extra}"


LOW = [0.0, 0.1]  # the lowest two 0.1 m layers
OBS_BINS = (0.9, 0.7, 0.5, 0.0)  # lower limits of the observed-fraction bins (cells in the first bin they reach)
MID = [0.1, 0.2, 0.3]  # layers below the tops of 0.6 m stems, away from the ground


def build_cases():
    cases = {}
    # density sweep, standard scene (stems 0.25 +- 0.08 m) and scan
    cases["density"] = [(f"random {n}/m2", scene(n), variants(), LOW, 0.01) for n in (30, 100, 300, 600, 1000, 1500)]
    # stems of one height (0.6 m), and stems taller than the scanner (no tops in view)
    cases["tall"] = [(f"height 0.6 m, {n}/m2", scene(n, extra="--height 0.6,0"), variants(), MID, 0.01)
                     for n in (100, 300, 1000)]
    cases["tall"] += [(f"height 2 m (no tops), {n}/m2", scene(n, extra="--height 2,0"), variants(), MID, 0.01)
                      for n in (100, 300)]
    # beam: none (d_eff = d) and wide (about 3 cm)
    cases["beam"] = [(f"no footprint, {n}/m2", scene(n, footprint=""), variants(footprint=""), LOW, 0.01)
                     for n in (300, 1000, 1500)]
    cases["beam"] += [(f"wide footprint, {n}/m2", scene(n, footprint=WIDE), variants(footprint=WIDE), LOW, 0.01)
                      for n in (100, 300, 1000)]
    # window width W = (2k+1) w, k = round(pool / w)
    cases["window"] = []
    for voxel, pool in ((0.05, 0.0), (0.025, 0.05), (0.025, 0.1), (0.05, 0.1), (0.05, 0.15), (0.1, 0.05),
                        (0.1, 0.15)):
        for n in (100, 300, 1000):
            cases["window"].append((f"voxel {voxel} pool {pool}, {n}/m2", scene(n),
                                    variants(voxel=voxel, pool=pool), LOW, 0.01))
    # stem diameter
    cases["diameter"] = []
    for d, n in ((0.005, 300), (0.005, 1500), (0.02, 30), (0.02, 300)):
        cases["diameter"].append((f"d {d} m, {n}/m2", scene(n, extra=f"--diameter {d}"), variants(diameter=d), LOW, d))
    # patches and cable roots
    cases["patches"] = [
        ("patches r0.2 30% 300/m2", scene(300, extra="--patchy 0.2,0.3,0"), variants(), LOW, 0.01),
        ("patches r0.3 30% 300/m2", scene(300, extra="--patchy 0.3,0.3,0"), variants(), LOW, 0.01),
        ("patches r0.5 30% 300/m2", scene(300, extra="--patchy 0.5,0.3,0"), variants(), LOW, 0.01),
        ("patches r0.3 30% 1000/m2", scene(1000, extra="--patchy 0.3,0.3,0"), variants(), LOW, 0.01),
        ("cable roots 150/m2", scene(150, extra="--clustered"), variants(), LOW, 0.01),
        ("cable roots 500/m2", scene(500, extra="--clustered"), variants(), LOW, 0.01),
        ("cable roots 150/m2 pool 0.1", scene(150, extra="--clustered"), variants(pool=0.1), LOW, 0.01),
    ]
    # smaller windows, which the self-shadow correction makes usable, in patches and along cable roots
    cases["small_window"] = []
    for label, extra, n in (("random", "", 300), ("patches r0.2 30%", "--patchy 0.2,0.3,0", 300),
                            ("patches r0.3 30%", "--patchy 0.3,0.3,0", 300), ("cable roots", "--clustered", 150)):
        for voxel, pool in ((0.05, 0.0), (0.025, 0.025)):
            cases["small_window"].append((f"{label} {n}/m2, voxel {voxel} pool {pool}", scene(n, extra=extra),
                                          variants(voxel=voxel, pool=pool), LOW, 0.01))
    # independence of the number of rays (run with --rays_per_m reduced, e.g. 14062 for a quarter)
    cases["rays"] = [(f"random {n}/m2, pool {pool}", scene(n), variants(pool=pool), LOW, 0.01)
                     for n in (30, 300) for pool in (0.0, 0.05)]
    # regular (min gap) and overlapping (Poisson) stems
    cases["spacing"] = [
        ("min gap 1, 1000/m2", scene(1000, extra="--min_gap 1"), variants(extra="--min_gap 1"), LOW, 0.01),
        ("overlapping 300/m2", scene(300, extra="--overlapping"), variants(), LOW, 0.01),
    ]
    # assumed diameter error (true 1 cm)
    cases["diameter_error"] = []
    for n in (300, 1000):
        v = []
        for d in (0.008, 0.012):
            v += [(f"{name} d={d}", args) for name, args in variants(diameter=d)]
        cases["diameter_error"].append((f"random {n}/m2, assumed d", scene(n), v, LOW, 0.01))
    # low observation: dense, tall, sparse scan lines, high scanner
    cases["lowobs"] = [
        ("1000/m2 height 0.6", scene(1000, extra="--height 0.6,0"), variants(), [0.0, 0.1, 0.2, 0.3, 0.4], 0.01),
        ("1500/m2 walk 2 m", scene(1500, extra="--walk_spacing 2"), variants(), LOW, 0.01),
        ("600/m2 walk 2 m height 0.4", scene(600, extra="--walk_spacing 2 --height 0.4,0.1"), variants(),
         [0.0, 0.1, 0.2], 0.01),
        ("300/m2 scanner 2 m walk 2 m", scene(300, extra="--walk_spacing 2 --scanner_height 2"), variants(), LOW,
         0.01),
    ]
    # measurement noise
    cases["noise"] = [
        ("300/m2 range noise 1 cm", scene(300, extra="--range_noise 0.01"), variants(), LOW, 0.01),
        ("300/m2 pose noise 2 cm", scene(300, extra="--pose_noise 0.02,1"), variants(), LOW, 0.01),
    ]
    # block K (small windows): voxel x pool radius, with and without the self-shadow correction, over stem densities,
    # patches and cable roots, at 1/4, 1/2 and the standard ray density
    k_windows = ((0.05, 0.0), (0.025, 0.025), (0.05, 0.05), (0.025, 0.0))
    k_variants = []
    for voxel, pool in k_windows:
        k_variants += [(f"v{voxel}p{pool} {name}", args) for name, args in variants(voxel=voxel, pool=pool)]
    k_scenes = [(f"random {n}/m2", scene(n)) for n in (30, 100, 300, 600, 1000, 1500)]
    k_scenes += [(f"patches r{r} 30% 300/m2", scene(300, extra=f"--patchy {r},0.3,0")) for r in (0.2, 0.3, 0.5)]
    k_scenes += [(f"cable roots {n}/m2", scene(n, extra="--clustered")) for n in (150, 500)]
    cases["K"] = [(f"{label}, rays x{scale}", args, k_variants, LOW, 0.01, scale)
                  for label, args in k_scenes for scale in (0.25, 0.5, 1.0)]
    # the same scenes for the small-window configuration (sw: single 0.05 m voxels with the second-order self-shadow
    # factor, the occlusion response inverted over 3 x 3 voxels, the finite-ray correction and no neighbour priors),
    # its ablations, and the default window with and without the first-order correction
    cases["K2"] = [(f"{label}, rays x{scale}", args, small_window_variants(), LOW, 0.01, scale)
                   for label, args in k_scenes for scale in (0.25, 0.5, 1.0)]
    # stem diameter and footprint, at the standard and a quarter of the ray density
    k3 = [(f"d 0.005 m, {n}/m2", scene(n, extra="--diameter 0.005"), 0.005, FOOTPRINT) for n in (300, 1500)]
    k3 += [(f"d 0.02 m, {n}/m2", scene(n, extra="--diameter 0.02"), 0.02, FOOTPRINT) for n in (30, 300)]
    k3 += [(f"no footprint, {n}/m2", scene(n, footprint=""), 0.01, "") for n in (300, 1000)]
    k3 += [(f"wide footprint, {n}/m2", scene(n, footprint=WIDE), 0.01, WIDE) for n in (100, 300, 1000)]
    cases["K3"] = [(f"{label}, rays x{scale}", args, small_window_variants(footprint=fp, diameter=d, short=True), LOW,
                    d, scale) for label, args, d, fp in k3 for scale in (0.25, 1.0)]
    # stems of one height (0.6 m, so no stem tops in the layers below 0.5 m) and overlapping (Poisson) stems, to judge
    # the second-order factor away from the stem tops
    k5 = [(f"height 0.6 m, {n}/m2", scene(n, extra="--height 0.6,0")) for n in (30, 300, 1000)]
    k5 += [("overlapping 300/m2", scene(300, extra="--overlapping"))]
    cases["K5"] = [(f"{label}, rays x{scale}", args, small_window_variants(short=True), [0.0, 0.1, 0.2, 0.3], 0.01,
                    scale) for label, args in k5 for scale in (0.25, 1.0)]
    # assumed diameter error of +-20% (true 1 cm), render only
    cases["K4"] = []
    for n in (300, 1000):
        v = []
        for d in (0.008, 0.012):
            v += [(f"{name} d={d}", a) for name, a in small_window_variants(diameter=d, short=True)]
        cases["K4"].append((f"random {n}/m2, assumed d", scene(n), v, LOW, 0.01))
    return cases


def run(cmd, log):
    with open(log, "a") as f:
        f.write("$ " + cmd + "\n")
        f.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT)


def read(file_name):
    with open(file_name) as f:
        return list(csv.DictReader(f))


def key(row):
    return tuple(round(float(row[k]), 4) for k in ("x_min", "y_min", "z_min"))


def run_case(job):
    set_name, label, scene_args, render_variants, layers, diameter, ray_scale, seed, args = job
    work = tempfile.mkdtemp(prefix="ss_", dir=args.tmp)
    rows = []
    try:
        log = os.path.join(work, "log.txt")
        half = 0.5 * args.plot_width
        walk = 1.0
        m = re.search(r"--walk_spacing ([0-9.]+)", scene_args)
        if m:
            walk = float(m.group(1))
        num_lines = max(1, int(args.plot_width / walk))
        rays = int(args.rays_per_m * ray_scale * args.plot_width * num_lines)
        top = max(layers) + 0.1
        run(f"cd {work} && raycreate pneumatophores {seed} {scene_args} --plot_width {args.plot_width} "
            f"--rays {rays} --truth_cells 1,1,0.1", log)
        truth = {key(r): float(r["true_frontal_area_per_volume_Nd"])
                 for r in read(os.path.join(work, "pneumatophores_truth_cells.csv"))}
        for name, render in render_variants:
            cells_file = os.path.join(work, "cells.csv")
            rlog = os.path.join(work, "render.txt")
            if os.path.exists(rlog):
                os.remove(rlog)
            start = time.time()
            run(f"cd {work} && rayrender pneumatophores.ply top density --density_model vertical_cylinder --no_flat_top "
                f"{render} --min_bound {-half},{-half},0 --max_bound {half},{half},{top} --cells cells.csv "
                f"--cell_size 1,1,0.1", rlog)
            render_time = time.time() - start
            text = open(rlog).read()
            m = re.search(r"mean effective width / diameter of ([0-9.e+-]+)", text)
            rho = float(m.group(1)) if m else 1.0
            voxel = float(re.search(r"--pixel_width ([0-9.]+)", render).group(1))
            pool = float(re.search(r"--pool_radius ([0-9.]+)", render).group(1))
            k = max(0, int(math.floor(pool / voxel + 0.5)))
            window = (2 * k + 1) * voxel
            d = float(re.search(r"--stem_diameter ([0-9.]+)", render).group(1))
            h1 = 1.0 + EXIT_DISTANCE_SQUARE * rho * d / window
            cells = read(cells_file)
            for z in layers:
                zr = round(z, 4)
                est = pooled = tr = obs = sq = 0.0
                count = 0
                # cells with observed fraction >= 0.9, and the sums per observed-fraction bin
                est9 = tr9 = 0.0
                bins = {b: [0.0, 0.0] for b in OBS_BINS}
                for c in cells:
                    kk = key(c)
                    if kk[2] != zr or kk not in truth:
                        continue
                    dens = float(c["density_frontal_area_per_volume_Nd"])
                    pool_d = float(c["pooled_density_frontal_area_per_volume_Nd"])
                    if math.isnan(dens):
                        dens = pool_d = 0.0
                    est += dens
                    pooled += pool_d
                    tr += truth[kk]
                    sq += (dens - truth[kk]) ** 2
                    frac = float(c["observed_fraction"])
                    frac = 0.0 if math.isnan(frac) else frac
                    obs += frac
                    count += 1
                    if frac >= 0.9:
                        est9 += dens
                        tr9 += truth[kk]
                    for b in OBS_BINS:
                        if frac >= b:
                            bins[b][0] += dens
                            bins[b][1] += truth[kk]
                            break
                if count == 0 or tr <= 0.0:
                    continue
                row = dict(set=set_name, case=label, variant=name, seed=seed, z_min=z, truth=tr / count,
                           density_ratio=est / tr, pooled_ratio=pooled / tr, observed=obs / count,
                           cell_rel_rmse=math.sqrt(sq / count) / (tr / count), rho=rho,
                           window=window, h1=h1, rays=rays, render_s=render_time,
                           density_ratio_obs90=est9 / tr9 if tr9 > 0 else float("nan"),
                           cells_obs90=sum(1 for c in cells if round(float(c["z_min"]), 4) == zr and
                                           not math.isnan(float(c["observed_fraction"])) and
                                           float(c["observed_fraction"]) >= 0.9))
                for b in OBS_BINS:
                    row[f"ratio_obs{b:g}"] = bins[b][0] / bins[b][1] if bins[b][1] > 0 else float("nan")
                rows.append(row)
    finally:
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)
    return rows


def mean_se(values):
    n = len(values)
    m = sum(values) / n
    if n < 2:
        return m, float("nan")
    return m, math.sqrt(sum((v - m) ** 2 for v in values) / (n - 1) / n)


def summarise(rows, out):
    groups = {}
    for r in rows:
        groups.setdefault((r["set"], r["case"], r["variant"], r["z_min"]), []).append(r)
    lines = ["set | case | variant | z | n | density/truth | obs>=0.9 density/truth (cells) | pooled/truth | "
             "mean/pooled | observed | cell rel. RMSE | H1 (W, rho) | render s"]
    for (set_name, case, variant, z), rs in groups.items():
        d = mean_se([r["density_ratio"] for r in rs])
        p = mean_se([r["pooled_ratio"] for r in rs])
        mp = mean_se([r["density_ratio"] / r["pooled_ratio"] for r in rs if r["pooled_ratio"] > 0])
        o = mean_se([r["observed"] for r in rs])
        e = mean_se([r["cell_rel_rmse"] for r in rs])
        v9 = [r["density_ratio_obs90"] for r in rs if not math.isnan(r["density_ratio_obs90"])]
        d9 = mean_se(v9) if v9 else (float("nan"), float("nan"))
        cells9 = sum(r["cells_obs90"] for r in rs) / len(rs)
        t = sum(r["render_s"] for r in rs) / len(rs)
        lines.append(f"{set_name} | {case} | {variant} | {z:.1f} | {len(rs)} | {d[0]:.3f} ± {d[1]:.3f} | "
                     f"{d9[0]:.3f} ± {d9[1]:.3f} ({cells9:.1f}) | "
                     f"{p[0]:.3f} ± {p[1]:.3f} | {mp[0]:.3f} ± {mp[1]:.3f} | {o[0]:.2f} | {e[0]:.3f} | "
                     f"{rs[0]['h1']:.3f} ({rs[0]['window']:.3f}, {rs[0]['rho']:.2f}) | {t:.2f}")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(out, "summary.txt"), "w") as f:
        f.write(text)
    print(text)
    compact(rows, out)


def compact(rows, out):
    """<out>/compact.txt: one row per case and layer, one column per variant, of density/truth (mean +- s.e. over the
    seeds) over all cells, and in brackets over the cells with observed fraction >= 0.9"""
    groups, names = {}, []
    for r in rows:
        groups.setdefault((r["set"], r["case"], r["z_min"]), {}).setdefault(r["variant"], []).append(r)
        if r["variant"] not in names:
            names.append(r["variant"])
    lines = ["set | case | z | observed | " + " | ".join(names)]
    for (set_name, case, z), vs in groups.items():
        cols = []
        for name in names:
            rs = vs.get(name)
            if not rs:
                cols.append("-")
                continue
            d = mean_se([r["density_ratio"] for r in rs])
            v9 = [r["density_ratio_obs90"] for r in rs if not math.isnan(r["density_ratio_obs90"])]
            cols.append(f"{d[0]:.3f} ± {d[1]:.3f} [{mean_se(v9)[0]:.3f}]" if v9 else f"{d[0]:.3f} ± {d[1]:.3f} [-]")
        observed = sum(r["observed"] for r in next(iter(vs.values()))) / len(next(iter(vs.values())))
        lines.append(f"{set_name} | {case} | {z:.1f} | {observed:.2f} | " + " | ".join(cols))
    with open(os.path.join(out, "compact.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    cases = build_cases()
    parser.add_argument("--sets", default="density", help="comma separated, of: " + ", ".join(cases) + ", or all")
    parser.add_argument("--cases", default="", help="only cases whose label contains this text")
    parser.add_argument("--seeds", default="101,102,103")
    parser.add_argument("--plot_width", type=float, default=4.0)
    parser.add_argument("--rays_per_m", type=float, default=56250.0,
                        help="rays per metre of walk (the default scene has 2 M rays over 36 m)")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--out", default="self_shadow_study")
    parser.add_argument("--bin", default="", help="directory added to the PATH, e.g. a build's tool directories")
    parser.add_argument("--keep", action="store_true", help="keep the scenes")
    args = parser.parse_args()
    if args.bin:
        os.environ["PATH"] = os.pathsep.join(args.bin.split(",")) + os.pathsep + os.environ["PATH"]
    os.makedirs(args.out, exist_ok=True)
    args.tmp = os.path.abspath(args.out)
    sets = list(cases) if args.sets == "all" else args.sets.split(",")
    jobs = []
    for s in sets:
        for case in cases[s]:
            label, scene_args, render_variants, layers, diameter = case[:5]
            ray_scale = case[5] if len(case) > 5 else 1.0
            if args.cases and args.cases not in label:
                continue
            for seed in (int(v) for v in args.seeds.split(",")):
                jobs.append((s, label, scene_args, render_variants, layers, diameter, ray_scale, seed, args))
    rows = []
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        for result in pool.map(run_case, jobs):
            rows += result
    rows.sort(key=lambda r: (sets.index(r["set"]), r["case"], r["variant"], r["z_min"], r["seed"]))
    with open(os.path.join(args.out, "rows.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    summarise(rows, args.out)


if __name__ == "__main__":
    sys.exit(main())
