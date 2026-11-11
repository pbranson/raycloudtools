#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Sensitivity of the cell density estimates to the random seed of the synthetic scene (stem layout, walking path
and rays), with and without the stem occlusion correction. Requires raycreate and rayrender on the PATH, and numpy.
Run in an empty working directory.

For each scene and seed, reports estimate / truth summed over the 1 x 1 m cells of each 0.1 m layer, and the
relative RMSE of the cells; then the mean and standard deviation over seeds. Both cell estimates are evaluated: the
mean of the voxel densities (density_..., the main output) and the density pooled over the cell's voxels
(pooled_density_...), with their ratio (patchiness_ratio, the mean over the layer's cells).
"""
import argparse
import concurrent.futures
import csv
import os
import subprocess

import numpy as np

FOOTPRINT = "0.01,0.003,0.05"
# the cells file's estimates: pooled over the cell's voxels, and the mean of the voxel densities
COLUMNS = {"pooled": "pooled_density_", "mean": "density_"}
SCENES = {"random 100/m2": "--stems_per_m2 100", "random 300/m2": "--stems_per_m2 300",
          "random 1000/m2": "--stems_per_m2 1000", "patches r 0.3 m, 300/m2": "--stems_per_m2 300 --patchy 0.3,0.4,0"}


def run(cmd, cwd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=cwd)


def read_csv(name):
    with open(name) as f:
        return list(csv.DictReader(f))


def key(row):
    """a cell's identity, from its minimum corner rounded to 0.1 mm"""
    return tuple(round(float(row[k]), 4) for k in ("x_min", "y_min", "z_min"))


def one(scene, options, seed, layers, extra_args):
    d = f"{scene.replace(' ', '_').replace('/', '_').replace(',', '')}_{seed}"
    os.makedirs(d, exist_ok=True)
    run(f"raycreate pneumatophores {seed} {options} --footprint {FOOTPRINT} --truth_cells 1,1,0.1", d)
    truth = {key(r): float(r["true_frontal_area_per_volume_Nd"])
             for r in read_csv(f"{d}/pneumatophores_truth_cells.csv")}
    out = {}
    for name, extra in [("corrected", ""), ("uncorrected", "--no_overlap_correction")]:
        run(f"rayrender pneumatophores.ply top density --density_model vertical_cylinder --estimator pooled "
            f"--no_flat_top --stem_diameter 0.01 --footprint {FOOTPRINT} --pixel_width 0.05 --min_bound -3,-3,0 "
            f"--max_bound 3,3,{0.1 * layers:.1f} --cells cells_{name}.csv {extra} {extra_args}", d)
        cells = read_csv(f"{d}/cells_{name}.csv")
        for label, prefix in COLUMNS.items():
            col = [k for k in cells[0] if k.startswith(prefix)][0]
            result = []
            for z in range(layers):
                layer = [c for c in cells if abs(float(c["z_min"]) - 0.1 * z) < 1e-6 and key(c) in truth]
                pairs = [(float(c[col]), truth[key(c)]) for c in layer if c[col] != "nan"]
                if not pairs:
                    raise RuntimeError(f"{d}: no estimated cells in layer {z} match the truth cells")
                e = np.array([p[0] for p in pairs])
                t = np.array([p[1] for p in pairs])
                patchiness = np.nanmean([float(c["patchiness_ratio"]) for c in layer])
                result.append((e.sum() / t.sum(), np.sqrt(np.mean((e - t) ** 2)) / t.mean(), patchiness))
            out[(name, label)] = result
    return scene, seed, out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", type=int, default=8)
    parser.add_argument("--layers", type=int, default=3)
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    parser.add_argument("--scenes", default="", help="only scenes containing this text")
    parser.add_argument("--extra", default="", help="extra rayrender arguments, e.g. --pool_radius 0.1")
    args = parser.parse_args()
    jobs = [(scene, options, seed, args.layers, args.extra) for scene, options in SCENES.items()
            if args.scenes in scene for seed in range(101, 101 + args.seeds)]
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for scene, seed, out in pool.map(lambda j: one(*j), jobs):
            results.setdefault(scene, {})[seed] = out
    with open("seed_study.csv", "w") as f:
        f.write("scene,seed,correction,column,layer,ratio,relative_rmse,patchiness_ratio\n")
        for scene, by_seed in results.items():
            for seed, out in by_seed.items():
                for (name, label), layers in out.items():
                    for z, (ratio, rmse, patchiness) in enumerate(layers):
                        f.write(f"\"{scene}\",{seed},{name},{label},{z},{ratio:.4f},{rmse:.4f},{patchiness:.4f}\n")
    print(f"estimate / truth over 1 x 1 m cells, mean ± s.d. (min-max) over {args.seeds} seeds; cell relative RMSE")
    for scene, by_seed in results.items():
        print(f"== {scene}")
        for name in ["corrected", "uncorrected"]:
            for label in COLUMNS:
                parts = []
                for z in range(args.layers):
                    r = np.array([by_seed[s][(name, label)][z][0] for s in by_seed])
                    e = np.array([by_seed[s][(name, label)][z][1] for s in by_seed])
                    parts.append(f"{r.mean():.3f} ± {r.std(ddof=1):.3f} ({r.min():.3f}-{r.max():.3f}), "
                                 f"RMSE {e.mean():.3f}")
                print(f"  {name:12s} {label:6s} " + " | ".join(parts))
            p = [np.mean([by_seed[s][(name, "pooled")][z][2] for s in by_seed]) for z in range(args.layers)]
            print(f"  {name:12s} patchiness ratio (mean / pooled, mean over cells) " +
                  " | ".join(f"{v:.3f}" for v in p))


if __name__ == "__main__":
    main()
