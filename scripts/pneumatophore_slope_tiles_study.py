#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Sloped ground and tiling with raydensity, as in Section 7.7 of docs/pneumatophore_density_report.md. Requires
raycreate, raydensity, raysplit and raytranslate on the PATH (or --bin), and numpy. Run in an empty working directory
(or give --out).

For each seed (101, 102, ...): generate a 300 stems/m2 scene (flat ground), shear it vertically by a ground plane
z += sx x + sy y (keeping each ray's horizontal path, so its hits and G-weighted path), and estimate 1 x 1 x 0.1 m
cells (0-0.3 m) of the 6 x 6 m plot with raydensity:
  flat      : the original scene, no ground mesh
  ground    : the sheared scene, --ground plane mesh (heights above ground)
  absolute  : the sheared scene, no ground mesh (absolute heights)
  flattened : the sheared scene flattened with raytranslate subtract, no ground mesh
and compare each with the truth cells and with flat. Seed 101 also runs scripts/raydensity_tiles.sh on the flat
scene (3 x 3 tiles of 2 m) and compares it with the single (6 m tile) run, writing tiles_compare.txt. Writes
slope.csv and prints the summary. Example:
  scripts/pneumatophore_slope_tiles_study.py --seeds 4 --out slope
"""
import argparse
import csv
import json
import os
import struct
import subprocess
import sys

import numpy as np

FOOTPRINT = "0.01,0.003,0.05"


def run(cmd, log="log.txt"):
    with open(log, "a") as f:
        f.write(f"$ {cmd}\n")
        f.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT)


def zarr(store, name, chunk="0.0.0"):
    meta = json.load(open(f"{store}/{name}/.zarray"))
    a = np.fromfile(f"{store}/{name}/{chunk}", dtype="<f4").reshape(meta["chunks"])
    return a


def write_plane(name, sx, sy, lo=-40.0, hi=40.0, spacing=2.0):
    n = int(round((hi - lo) / spacing))
    verts = [(lo + i * spacing, lo + j * spacing) for j in range(n + 1) for i in range(n + 1)]
    faces = []
    for j in range(n):
        for i in range(n):
            v = i + (n + 1) * j
            faces += [(v, v + 1, v + n + 2), (v, v + n + 2, v + n + 1)]
    with open(name, "wb") as f:
        f.write(b"ply\nformat binary_little_endian 1.0\n")
        f.write(f"element vertex {len(verts)}\nproperty double x\nproperty double y\nproperty double z\n"
                "property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\n"
                f"element face {len(faces)}\nproperty list int int vertex_indices\nend_header\n".encode())
        for x, y in verts:
            f.write(struct.pack("<dddBBBB", x, y, sx * x + sy * y, 127, 127, 127, 255))
        for a, b, c in faces:
            f.write(struct.pack("<iiii", 3, a, b, c))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seeds", type=int, default=4)
    p.add_argument("--slope", default="0.1,-0.05", help="ground slope sx,sy")
    p.add_argument("--out", default=".", help="output directory (created if needed)")
    p.add_argument("--bin", default="", help="comma separated directories added to the PATH (e.g. a build's tools)")
    a = p.parse_args()
    a.scripts = os.path.dirname(os.path.abspath(__file__))
    if a.bin:
        os.environ["PATH"] = os.pathsep.join(os.path.abspath(b) for b in a.bin.split(",")) + os.pathsep + \
            os.environ["PATH"]
    os.makedirs(a.out, exist_ok=True)
    os.chdir(a.out)
    sys.path.insert(0, a.scripts)
    from rayply_numpy import read_ply, write_ply
    sx, sy = (float(v) for v in a.slope.split(","))
    common = (f"--origin -3,-3 --tile_width 6 --tiles 1,1 --tile 0,0 --heights 0,0.3 --stem_diameter 0.01 "
              f"--footprint {FOOTPRINT}")
    top = os.getcwd()
    rows = []
    for seed in range(101, 101 + a.seeds):
        d = f"seed_{seed}"
        os.makedirs(d, exist_ok=True)
        os.chdir(d)
        run(f"raycreate pneumatophores {seed} --stems_per_m2 300 --footprint {FOOTPRINT} --truth_cells 1,1,0.1")
        truth = np.zeros((3, 6, 6))
        with open("pneumatophores_truth_cells.csv") as f:
            for r in csv.DictReader(f):
                i, j, k = (int(round(float(r[c]) / s - o)) for c, s, o in
                           (("x_min", 1, -3), ("y_min", 1, -3), ("z_min", 0.1, 0)))
                if 0 <= i < 6 and 0 <= j < 6 and 0 <= k < 3:
                    truth[k, j, i] = float(r["true_frontal_area_per_volume_Nd"])
        header, data = read_ply("pneumatophores.ply")
        ray = ("nx", "ny", "nz") if "nx" in data.dtype.names else ("rayx", "rayy", "rayz")
        sheared = data.copy()
        sheared["z"] = data["z"] + sx * data["x"] + sy * data["y"]
        sheared[ray[2]] = data[ray[2]] + sx * data[ray[0]] + sy * data[ray[1]]
        write_ply("sloped.ply", header, sheared)
        write_plane("plane.ply", sx, sy)
        # raytranslate overwrites its input, so flatten a copy
        write_ply("flattened.ply", header, sheared)
        run("raytranslate flattened.ply subtract plane.ply")
        flattened = ["flattened.ply"]
        cases = {"flat": f"pneumatophores.ply flat.zarr {common}",
                 "ground": f"sloped.ply ground.zarr {common} --ground plane.ply",
                 "absolute": f"sloped.ply absolute.zarr {common}",
                 "flattened": f"{flattened[0]} flattened.zarr {common}"}
        est = {}
        for name, args in cases.items():
            run(f"raydensity {args}")
            est[name] = {c: zarr(f"{name}.zarr", c) for c in ("density", "pooled_density")}
        for name in cases:
            for col in ("density", "pooled_density"):
                for k in range(3):
                    e, t, fl = est[name][col][k], truth[k], est["flat"][col][k]
                    ok = ~np.isnan(e)
                    rows.append((seed, name, col, k, e[ok].sum() / t[ok].sum(), np.nanmax(np.abs(e / fl - 1)),
                                 e[ok].sum() / fl[ok].sum(), ok.sum()))
        if seed == 101:
            # (i) tiled vs single, on the flat scene
            run(f"bash {a.scripts}/raydensity_tiles.sh pneumatophores.ply tiled.zarr -3,-3 2 3,3 -1,1 4 -- "
                f"--heights 0,0.3 --stem_diameter 0.01 --footprint {FOOTPRINT}")
            with open("../tiles_compare.txt", "w") as f:
                for col in ("density", "pooled_density", "hits", "weighted_path", "observed_fraction",
                            "patchiness_ratio"):
                    single = zarr("flat.zarr", col)
                    worst = 0.0
                    for i in range(3):
                        for j in range(3):
                            t = zarr("tiled.zarr", col, f"0.{j}.{i}")
                            s = single[:, 2 * j:2 * j + 2, 2 * i:2 * i + 2]
                            both = ~np.isnan(s)
                            assert np.array_equal(both, ~np.isnan(t))
                            worst = max(worst, np.max(np.abs(t[both] / s[both] - 1)) if both.any() else 0)
                    f.write(f"{col}: max relative difference tiled (3 x 3 tiles of 2 m) vs single = {worst:.3g}\n")
        for f_ in ["pneumatophores.ply", "sloped.ply"] + flattened:
            os.remove(f_)
        os.chdir(top)
    with open("slope.csv", "w") as f:
        f.write("seed,case,column,layer,estimate_over_truth,max_cell_rel_diff_vs_flat,estimate_over_flat,cells\n")
        for r in rows:
            f.write(",".join(str(round(v, 5)) if isinstance(v, float) else str(v) for v in r) + "\n")
    print(f"slope {sx},{sy}: plot estimate / truth (and / flat), mean ± s.d. over {a.seeds} seeds; "
          "max |cell/flat - 1| over seeds")
    for name in ("flat", "ground", "absolute", "flattened"):
        for col in ("density", "pooled_density"):
            parts = []
            for k in range(3):
                sel = [r for r in rows if r[1] == name and r[2] == col and r[3] == k]
                v = np.array([r[4] for r in sel])
                fl = np.array([r[6] for r in sel])
                mx = max(r[5] for r in sel)
                parts.append(f"{v.mean():.3f} ± {v.std(ddof=1):.3f} (/flat {fl.mean():.3f}, max cell {mx:.3f})")
            print(f"  {name:9s} {col:15s} " + " | ".join(parts))


if __name__ == "__main__":
    main()
