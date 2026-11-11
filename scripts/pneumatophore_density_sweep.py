#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Sweep stem arrangements and voxel sizes on synthetic pneumatophore scans, comparing the estimated frontal area
per volume (N*d) with the truth. Requires raycreate and rayrender on the PATH. Run in an empty working directory.

Example:
  pneumatophore_density_sweep.py --widths 0.025,0.05,0.1 --footprint 0.01,0.003,0.05 --estimator crossing
"""
import argparse
import csv
import subprocess


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def read(file_name):
    with open(file_name) as f:
        return list(csv.DictReader(f))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--widths", default="0.025,0.05,0.1,0.2", help="voxel widths (m)")
    parser.add_argument("--footprint", default="", help="beam exit diameter,divergence,threshold. Empty for none")
    parser.add_argument("--estimator", default="crossing", choices=["crossing", "beer_lambert"])
    parser.add_argument("--diameter", type=float, default=0.01, help="stem diameter (m)")
    parser.add_argument("--stems_per_m2", type=float, default=100.0)
    parser.add_argument("--spreads", default="0.1,0.03,0.01", help="cable root cluster spreads (m) to test")
    parser.add_argument("--seed", type=int, default=3)
    args = parser.parse_args()

    create = f"--diameter {args.diameter} --stems_per_m2 {args.stems_per_m2}"
    render = (f"--density_model vertical_cylinder --estimator {args.estimator} --no_flat_top --min_sin_zenith 0.9 "
              "--min_bound -3,-3,0 --max_bound 3,3,0.8 --profile profile.csv")
    if args.footprint:
        create += f" --footprint {args.footprint}"
        render += f" --stem_diameter {args.diameter} --footprint {args.footprint}"
    cases = [("random", "")] + [(f"cable roots, spread {s} m", f"--clustered --cluster_spread {s}")
                                for s in args.spreads.split(",") if s]
    print("arrangement, voxel width, hits per voxel (0-0.2 m), mean-of-voxels FAI / true, pooled FAI / true")
    for name, options in cases:
        run(f"raycreate pneumatophores {args.seed} {create} {options}")
        truth = read("pneumatophores_truth.csv")
        layer = float(truth[0]["z_max"]) - float(truth[0]["z_min"])
        true_fai = sum(float(r["true_frontal_area_per_volume_Nd"]) * layer for r in truth)
        for width in [float(w) for w in args.widths.split(",")]:
            run(f"rayrender pneumatophores.ply top density --pixel_width {width} {render}")
            profile = read("profile.csv")
            pooled = [k for k in profile[0] if k.startswith("pooled_density_")][0]
            mean = [k for k in profile[0] if k.startswith("density_")][0]
            # unobserved layers are nan, and contribute nothing
            mean_fai = sum(float(r[mean]) * width for r in profile if r[mean] != "nan")
            pooled_fai = sum(float(r[pooled]) * width for r in profile if r[pooled] != "nan")
            low = [r for r in profile if float(r["z_min"]) < 0.1999 and r["mean_hits_per_voxel"] != "nan"]
            hits = sum(float(r["mean_hits_per_voxel"]) for r in low) / len(low)
            print(f"{name}, {width}, {hits:.1f}, {mean_fai / true_fai:.3f}, {pooled_fai / true_fai:.3f}", flush=True)


if __name__ == "__main__":
    main()
