#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Trunks alone, as in Section 7.8 of docs/pneumatophore_density_report.md: scenes with almost no pneumatophores
(0.01 per m2, the generator minimum) and trunks at 0.1 per m2 of mean diameter 0.2 m. Reports the estimate (the
`density` column) over the true trunk frontal area of the 6 x 6 m plot, layers 0-0.3 m, for three settings:
  (1) the pneumatophore settings (--stem_diameter 0.01, footprint)
  (2) --stem_diameter set to the trunks' mean diameter (0.2 m), footprint
  (3) no footprint correction and no occlusion correction
per seed (31..34 by default; a seed with no trunk in the plot gives nan) and as mean +- s.d. Requires raycreate and
rayrender on the PATH (or --bin), and numpy. Run in an empty working directory (or give --out). Example:
  scripts/pneumatophore_trunk_only_study.py --out trunks_alone
"""
import argparse
import csv
import os
import subprocess

import numpy as np

F = "0.01,0.003,0.05"
VARIANTS = {"stem settings (d 0.01, footprint)": f"--stem_diameter 0.01 --footprint {F}",
            "trunk diameter (d 0.2, footprint)": f"--stem_diameter 0.2 --footprint {F}",
            "no footprint, no occlusion correction": ""}


def run(cmd, log):
    with open(log, "a") as f:
        f.write(f"$ {cmd}\n")
        f.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT)


def key(r):
    return tuple(round(float(r[k]), 4) for k in ("x_min", "y_min", "z_min"))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seeds", default="31,32,33,34", help="comma separated seeds")
    p.add_argument("--out", default=".", help="output directory (created if needed)")
    p.add_argument("--bin", default="", help="comma separated directories added to the PATH (e.g. a build's tools)")
    a = p.parse_args()
    if a.bin:
        os.environ["PATH"] = os.pathsep.join(os.path.abspath(b) for b in a.bin.split(",")) + os.pathsep + \
            os.environ["PATH"]
    os.makedirs(a.out, exist_ok=True)
    os.chdir(a.out)
    seeds = [int(s) for s in a.seeds.split(",")]
    res = {v: [] for v in VARIANTS}
    top = os.getcwd()
    for seed in seeds:
        d = f"seed_{seed}"
        os.makedirs(d, exist_ok=True)
        os.chdir(d)
        run(f"raycreate pneumatophores {seed} --stems_per_m2 0.01 --trunks 0.1,0.2 --footprint {F} "
            "--truth_cells 1,1,0.1", "log.txt")
        with open("pneumatophores_truth_cells.csv") as f:
            truth = {key(r): float(r["true_trunk_frontal_area_per_volume"]) for r in csv.DictReader(f)}
        for i, (v, opts) in enumerate(VARIANTS.items()):
            run(f"rayrender pneumatophores.ply top density --density_model vertical_cylinder --estimator pooled "
                f"--no_flat_top --pixel_width 0.05 --min_bound -3,-3,0 --max_bound 3,3,0.3 --cells cells_{i}.csv "
                f"{opts}", "log.txt")
            with open(f"cells_{i}.csv") as f:
                cells = list(csv.DictReader(f))
            col = "density_frontal_area_per_volume_Nd"
            used = [c for c in cells if c[col] != "nan" and key(c) in truth]
            e = sum(float(c[col]) for c in used)
            t = sum(truth[key(c)] for c in used)
            res[v].append(e / t if t > 0 else float("nan"))
        os.remove("pneumatophores.ply")
        os.chdir(top)
    print(f"trunks alone (0.1 /m2, mean d 0.2 m), estimate (density column) / true trunk frontal area over the "
          f"6 x 6 m plot, 0-0.3 m; per seed {a.seeds} (seeds with no trunk in the plot give nan)")
    for v, r in res.items():
        vals = np.array(r)
        ok = ~np.isnan(vals)
        mean = vals[ok].mean() if ok.any() else float("nan")
        sd = vals[ok].std(ddof=1) if ok.sum() > 1 else 0.0
        print(f"  {v:40s} " + " ".join(f"{x:.3f}" for x in vals) + f"  mean {mean:.3f} ± {sd:.3f}")


if __name__ == "__main__":
    main()
