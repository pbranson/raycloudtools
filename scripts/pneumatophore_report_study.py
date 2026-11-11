#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Multi-seed, multi-variant study of the pneumatophore cell density estimates, as used for the tables of
docs/pneumatophore_density_report.md (Sections 5 and 7). Requires raycreate and rayrender on the PATH (or --bin),
and numpy.

For each scene (raycreate pneumatophores options) and seed, the scene is generated once (6 x 6 m plot, the default
scan) and estimated with each variant (extra rayrender options), with --estimator pooled, 5 cm voxels, the scene's
footprint and a 1 cm stem diameter unless the variant sets them. It reports estimate / truth summed over the 36
1 x 1 m cells of each 0.1 m layer, for `density` (mean of the voxel densities, the main output) and
`pooled_density`, the cell relative RMSE and the mean patchiness_ratio; then the mean +- s.d. over the seeds. Writes
<out>/study.csv and one directory per scene and seed (with its log and cells files), and prints the summary.
Same method as scripts/pneumatophore_seed_study.py. Example (Table 2):

  scripts/pneumatophore_report_study.py --out table2 --layers 2 \\
    --scene "random 300/m2=--stems_per_m2 300" --scene "random 1000/m2=--stems_per_m2 1000" \\
    --variant "pooled corrected=" --variant "pooled uncorrected=--no_overlap_correction" \\
    --variant "crossing=--estimator crossing"
"""
import argparse
import concurrent.futures
import csv
import os
import subprocess

import numpy as np

FOOTPRINT = "0.01,0.003,0.05"
COLUMNS = {"mean": "density_", "pooled": "pooled_density_"}


def run(cmd, cwd, log):
    with open(os.path.join(cwd, log), "a") as f:
        f.write(f"$ {cmd}\n")
        f.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT, cwd=cwd)


def read_csv(name):
    with open(name) as f:
        return list(csv.DictReader(f))


def key(row):
    return tuple(round(float(row[k]), 4) for k in ("x_min", "y_min", "z_min"))


def slug(s):
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in s)


def one(scene, scene_opts, seed, variants, layers, keep):
    d = f"{slug(scene)}_{seed}"
    os.makedirs(d, exist_ok=True)
    if "--footprint" not in scene_opts:
        scene_opts += f" --footprint {FOOTPRINT}"
    run(f"raycreate pneumatophores {seed} {scene_opts} --truth_cells 1,1,0.1", d, "log.txt")
    truth = {key(r): float(r["true_frontal_area_per_volume_Nd"]) for r in read_csv(f"{d}/pneumatophores_truth_cells.csv")}
    out = {}
    for vname, vopts in variants:
        if "--footprint" not in vopts:
            vopts += f" --footprint {FOOTPRINT}"
        if "--stem_diameter" not in vopts:
            vopts += " --stem_diameter 0.01"
        est = "" if "--estimator" in vopts else "--estimator pooled "
        voxel = "" if "--pixel_width" in vopts else "--pixel_width 0.05 "
        cells_file = f"cells_{slug(vname)}.csv"
        run(f"rayrender pneumatophores.ply top density --density_model vertical_cylinder {est}--no_flat_top "
            f"{voxel}--min_bound -3,-3,0 --max_bound 3,3,{0.1 * layers:.1f} --cells {cells_file} {vopts}", d, "log.txt")
        cells = read_csv(f"{d}/{cells_file}")
        for label, prefix in COLUMNS.items():
            col = [k for k in cells[0] if k.startswith(prefix)][0]
            res = []
            for z in range(layers):
                layer = [c for c in cells if abs(float(c["z_min"]) - 0.1 * z) < 1e-6 and key(c) in truth]
                pairs = [(float(c[col]), truth[key(c)]) for c in layer if c[col] != "nan"]
                e = np.array([p[0] for p in pairs])
                t = np.array([p[1] for p in pairs])
                patch = np.nanmean([float(c["patchiness_ratio"]) for c in layer])
                res.append((e.sum() / t.sum(), np.sqrt(np.mean((e - t) ** 2)) / t.mean(), patch, len(pairs)))
            out[(vname, label)] = res
    if not keep:
        os.remove(f"{d}/pneumatophores.ply")
    return scene, seed, out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scene", action="append", required=True, help="'label=raycreate options' (repeatable)")
    p.add_argument("--variant", action="append", required=True, help="'label=rayrender options' (repeatable)")
    p.add_argument("--seeds", type=int, default=4)
    p.add_argument("--first_seed", type=int, default=101)
    p.add_argument("--layers", type=int, default=3, help="number of 0.1 m layers from the ground")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--out", default=".", help="output directory (created if needed)")
    p.add_argument("--bin", default="", help="comma separated directories added to the PATH (e.g. a build's tools)")
    p.add_argument("--keep_ply", action="store_true")
    a = p.parse_args()
    if a.bin:
        os.environ["PATH"] = os.pathsep.join(os.path.abspath(b) for b in a.bin.split(",")) + os.pathsep + \
            os.environ["PATH"]
    os.makedirs(a.out, exist_ok=True)
    os.chdir(a.out)
    scenes = [s.split("=", 1) for s in a.scene]
    variants = [v.split("=", 1) for v in a.variant]
    jobs = [(s, o, seed, variants, a.layers, a.keep_ply) for s, o in scenes
            for seed in range(a.first_seed, a.first_seed + a.seeds)]
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for scene, seed, out in pool.map(lambda j: one(*j), jobs):
            results.setdefault(scene, {})[seed] = out
    with open("study.csv", "w") as f:
        f.write("scene,seed,variant,column,layer,ratio,relative_rmse,patchiness_ratio,cells\n")
        for scene, by_seed in results.items():
            for seed, out in by_seed.items():
                for (v, label), res in out.items():
                    for z, (r, e, pr, n) in enumerate(res):
                        f.write(f"\"{scene}\",{seed},\"{v}\",{label},{z},{r:.4f},{e:.4f},{pr:.4f},{n}\n")
    print(f"estimate / truth over 1 x 1 m cells per 0.1 m layer, mean ± s.d. over {a.seeds} seeds "
          f"({a.first_seed}..{a.first_seed + a.seeds - 1}); [cell relative RMSE]; patchiness = mean density/pooled")
    for scene, by_seed in results.items():
        print(f"== {scene}")
        for v, _ in variants:
            for label in COLUMNS:
                parts = []
                for z in range(a.layers):
                    r = np.array([by_seed[s][(v, label)][z][0] for s in by_seed])
                    e = np.array([by_seed[s][(v, label)][z][1] for s in by_seed])
                    sd = r.std(ddof=1) if len(r) > 1 else 0.0
                    parts.append(f"{r.mean():.3f} ± {sd:.3f} [{e.mean():.3f}]")
                print(f"  {v:28s} {label:6s} " + " | ".join(parts))
            pr = [np.mean([by_seed[s][(v, 'mean')][z][2] for s in by_seed]) for z in range(a.layers)]
            print(f"  {v:28s} patch  " + " | ".join(f"{x:.3f}" for x in pr))


if __name__ == "__main__":
    main()
