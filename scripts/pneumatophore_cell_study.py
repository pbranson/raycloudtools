#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Sensitivity of hydrodynamic-cell (e.g. 1 x 1 x 0.1 m) frontal area density estimates to occlusion, ray angle
and estimation voxel size, on synthetic hand-held scans of pneumatophore fields. Requires raycreate and rayrender
on the PATH. Run in an empty working directory. Writes cell_study.csv and prints a summary.

For each scene, estimation method and band of ray angles, cells are compared to the truth per height layer:
  bias      = sum of estimates / sum of truth, over observed cells
  rel_rmse  = cell RMSE / mean true cell value
  coverage  = fraction of cells with at least min_hits hits (so a usable estimate)
"""
import argparse
import csv
import math
import subprocess


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def read(file_name):
    with open(file_name) as f:
        return list(csv.DictReader(f))


def key(row):
    return tuple(round(float(row[k]), 4) for k in ("x_min", "y_min", "z_min"))


SCENES = {
    "dense random 300/m2": "--stems_per_m2 300",
    "random 100/m2": "--stems_per_m2 100",
    "patches r=0.5m, 30% cover": "--stems_per_m2 300 --patchy 0.5,0.3,0",
    "patches r=0.2m, 30% cover": "--stems_per_m2 300 --patchy 0.2,0.3,0",
}
# the cells file's columns: the density pooled over the cell's voxels, and the mean of the voxel densities
POOLED = "density_"
MEAN = "mean_voxel_density_"
# name: (estimator, voxel width, column of the cells file)
METHODS = {
    "pooled 0.05m voxels": ("pooled", 0.05, POOLED),
    "pooled 0.05m voxels, mean": ("pooled", 0.05, MEAN),
    "pooled 0.1m voxels": ("pooled", 0.1, POOLED),
    "pooled 0.1m voxels, mean": ("pooled", 0.1, MEAN),
    "crossing 0.05m voxels, mean": ("crossing", 0.05, MEAN),
    "crossing 0.1m voxels, mean": ("crossing", 0.1, MEAN),
    "beer_lambert 0.05m voxels, mean": ("beer_lambert", 0.05, MEAN),
    "beer_lambert single 1x1x0.1m voxel": ("beer_lambert", 0.05, POOLED),
}
# bands of sin(zenith angle): 0.5 is 30 degrees from vertical, 0.9 is 64 degrees
BANDS = [(0.1, 1.0), (0.5, 1.0), (0.9, 1.0), (0.1, 0.5), (0.5, 0.9)]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--footprint", default="0.01,0.003,0.05", help="exit diameter,divergence,threshold")
    parser.add_argument("--diameter", type=float, default=0.01)
    parser.add_argument("--cell_size", default="1,1,0.1")
    parser.add_argument("--layers", type=int, default=3, help="number of cell layers from the ground to evaluate")
    parser.add_argument("--walk_spacing", type=float, default=1.0)
    parser.add_argument("--rays", type=int, default=2000000)
    parser.add_argument("--min_hits", type=float, default=10.0)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--scenes", default="", help="only scenes containing this text")
    parser.add_argument("--methods", default="", help="only methods containing this text")
    parser.add_argument("--bands", default="", help="sin zenith bands to test, e.g. 0.1-1,0.5-1")
    parser.add_argument("--crossing_widths", default="", help="instead test the crossing estimator at these widths")
    parser.add_argument("--pooled_widths", default="", help="instead test the pooled estimator at these widths")
    parser.add_argument("--extra", default="", help="extra rayrender arguments, e.g. --no_overlap_correction")
    parser.add_argument("--scene_extra", default="", help="extra raycreate arguments, e.g. --min_gap 0.5")
    parser.add_argument("--stems_per_m2", default="", help="instead test random scenes of these stem densities")
    args = parser.parse_args()
    methods = METHODS
    if args.crossing_widths:
        methods = {f"crossing {w}m voxels, mean": ("crossing", float(w), MEAN) for w in args.crossing_widths.split(",")}
    if args.pooled_widths:
        methods = {}
        for w in args.pooled_widths.split(","):
            methods[f"pooled {w}m voxels"] = ("pooled", float(w), POOLED)
            methods[f"pooled {w}m voxels, mean"] = ("pooled", float(w), MEAN)
    scenes = SCENES
    if args.stems_per_m2:
        scenes = {f"random {n}/m2": f"--stems_per_m2 {n}" for n in args.stems_per_m2.split(",")}
    bands = [tuple(float(v) for v in b.split("-")) for b in args.bands.split(",")] if args.bands else BANDS
    cell_z = float(args.cell_size.split(",")[2])

    out = open("cell_study.csv", "w", newline="")
    writer = csv.writer(out)
    writer.writerow(["scene", "method", "min_sin", "max_sin", "z_min", "true_mean", "bias", "rel_rmse", "coverage"])
    for scene, scene_options in scenes.items():
        if args.scenes not in scene:
            continue
        run(f"raycreate pneumatophores {args.seed} {scene_options} --diameter {args.diameter} --rays {args.rays} "
            f"--footprint {args.footprint} --walk_spacing {args.walk_spacing} --truth_cells {args.cell_size} "
            f"{args.scene_extra}")
        truth = {key(r): float(r["true_frontal_area_per_volume_Nd"]) for r in read("pneumatophores_truth_cells.csv")}
        print(f"\n== {scene}")
        print("method | sin zenith band | per layer from ground: bias, rel. RMSE, coverage")
        for method, (estimator, width, column) in methods.items():
            if args.methods not in method:
                continue
            for low, high in bands:
                run(f"rayrender pneumatophores.ply top density --density_model vertical_cylinder "
                    f"--estimator {estimator} --no_flat_top --min_sin_zenith {low} --max_sin_zenith {high} "
                    f"--stem_diameter {args.diameter} --footprint {args.footprint} --pixel_width {width} "
                    f"--min_bound -3,-3,0 --max_bound 3,3,{cell_z * args.layers} --cells cells.csv "
                    f"--cell_size {args.cell_size} {args.extra}")
                cells = read("cells.csv")
                col = [k for k in cells[0].keys() if k.startswith(column)][0]
                summary = []
                for layer in range(args.layers):
                    z = round(layer * cell_z, 4)
                    pairs, total = [], 0
                    for c in cells:
                        k = key(c)
                        if k[2] != z or k not in truth:
                            continue
                        total += 1
                        if float(c["hits"]) >= args.min_hits:
                            pairs.append((float(c[col]), truth[k]))
                    true_mean = sum(t for _, t in pairs) / len(pairs) if pairs else 0.0
                    bias = sum(e for e, _ in pairs) / sum(t for _, t in pairs) if pairs and true_mean > 0 else float("nan")
                    rmse = math.sqrt(sum((e - t) ** 2 for e, t in pairs) / len(pairs)) / true_mean if pairs and true_mean > 0 else float("nan")
                    coverage = len(pairs) / total if total else 0.0
                    writer.writerow([scene, method, low, high, z, f"{true_mean:.4f}", f"{bias:.3f}", f"{rmse:.3f}",
                                     f"{coverage:.2f}"])
                    summary.append(f"{bias:.2f}, {rmse:.2f}, {coverage:.2f}")
                print(f"{method} | {low}-{high} | " + " | ".join(summary), flush=True)
    out.close()


if __name__ == "__main__":
    main()
