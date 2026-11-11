#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Impact of sparse trunks on pneumatophore density estimates, and whether trunks are separable in the bottom
0.3 m. Requires raycreate and rayrender on the PATH, and numpy. Run in an empty working directory.

Three cases are estimated on the same synthetic scan (pooled estimator, 5 cm voxels, footprint corrected) and
compared to the pneumatophore truth in 1 x 1 x 0.1 m cells:
  A: all returns (trunks unseparated)
  B: trunk returns removed using the generator's labels (perfect separation)
  C: trunks detected from returns above the pneumatophores (a horizontal circle fit per trunk), and returns within
     each fitted circle (plus a margin) removed at all heights
Removed returns are set to zero alpha, so a ray stopped by a trunk counts as unobstructed up to the trunk.
"""
import argparse
import csv
import subprocess

import numpy as np

FOOTPRINT = "0.01,0.003,0.05"


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def read_csv(name):
    with open(name) as f:
        return list(csv.DictReader(f))


PLY_TYPES = {"float": "<f4", "float32": "<f4", "double": "<f8", "float64": "<f8", "uchar": "u1", "uint8": "u1",
             "int": "<i4", "uint": "<u4", "short": "<i2", "ushort": "<u2", "char": "i1"}


def ply_dtype(header):
    """the vertex record type from a binary little-endian .ply header (float or double ray clouds)"""
    fields = []
    for line in header.decode("ascii", "replace").splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0] == "property":
            fields.append((parts[2], PLY_TYPES[parts[1]]))
    return np.dtype(fields)


def read_ply(name):
    with open(name, "rb") as f:
        header = b""
        while not header.endswith(b"end_header\n"):
            header += f.readline()
        return header, np.fromfile(f, dtype=ply_dtype(header))


def write_ply(name, header, data):
    with open(name, "wb") as f:
        f.write(header)
        data.tofile(f)


def components(cells):
    """connected components (8-neighbour) of a set of integer grid cells"""
    remaining = set(cells)
    groups = []
    while remaining:
        seed = remaining.pop()
        group, stack = [seed], [seed]
        while stack:
            cx, cy = stack.pop()
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    n = (cx + dx, cy + dy)
                    if n in remaining:
                        remaining.remove(n)
                        group.append(n)
                        stack.append(n)
        groups.append(group)
    return groups


def fit_circle(x, y):
    """algebraic (Kasa) least squares circle fit"""
    a = np.stack([x, y, np.ones_like(x)], axis=1)
    b = x * x + y * y
    c, *_ = np.linalg.lstsq(a, b, rcond=None)
    cx, cy = c[0] / 2, c[1] / 2
    return cx, cy, np.sqrt(c[2] + cx * cx + cy * cy)


def detect_trunks(data, z_range, grid=0.05, min_points=200):
    """trunk circles (x, y, radius of returns) from bounded returns within the height range"""
    sel = (data["alpha"] > 0) & (data["z"] >= z_range[0]) & (data["z"] <= z_range[1])
    x, y = data["x"][sel].astype(float), data["y"][sel].astype(float)
    cells = {}
    for i, key in enumerate(zip(np.floor(x / grid).astype(int), np.floor(y / grid).astype(int))):
        cells.setdefault(key, []).append(i)
    trunks = []
    for group in components(list(cells.keys())):
        idx = np.concatenate([cells[c] for c in group])
        if len(idx) >= min_points:
            trunks.append(fit_circle(x[idx], y[idx]))
    return trunks


def estimate(cloud, out):
    run(f"rayrender {cloud} top density --density_model vertical_cylinder --estimator pooled --no_flat_top "
        f"--stem_diameter 0.01 --footprint {FOOTPRINT} --pixel_width 0.05 --min_bound -3,-3,0 --max_bound 3,3,0.3 "
        f"--cells {out}")
    rows = read_csv(out)
    col = [k for k in rows[0] if k.startswith("mean_density")][0]
    return {(round(float(r["x_min"])), round(float(r["y_min"])), round(float(r["z_min"]), 2)): float(r[col])
            for r in rows}


def summarise(name, est, truth, trunk_cells):
    print(f"{name}")
    for label, keep in [("all cells", lambda k: True), ("cells with a trunk", lambda k: k[:2] in trunk_cells),
                        ("cells without", lambda k: k[:2] not in trunk_cells)]:
        parts = []
        for z in (0.0, 0.1, 0.2):
            keys = [k for k in truth if k[2] == z and keep(k) and k in est]
            e = sum(est[k] for k in keys)
            t = sum(truth[k][0] for k in keys)
            tt = sum(truth[k][0] + truth[k][1] for k in keys)
            parts.append(f"{e / t:.3f} ({e / tt:.3f})" if t > 0 else "  -  ")
        print(f"  {label:20s} estimate / pneumatophore truth (/ incl. trunk) for layers 0-0.1, 0.1-0.2, 0.2-0.3 m: "
              + " | ".join(parts))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stems_per_m2", type=float, default=300)
    parser.add_argument("--trunks", default="0.1,0.2", help="trunks per m^2 and mean diameter (m)")
    parser.add_argument("--seed", type=int, default=31)
    parser.add_argument("--margin", type=float, default=0.015, help="beyond the fitted trunk circle (m)")
    args = parser.parse_args()
    run(f"raycreate pneumatophores {args.seed} --stems_per_m2 {args.stems_per_m2} --trunks {args.trunks} "
        f"--footprint {FOOTPRINT} --truth_cells 1,1,0.1")
    truth = {}
    for r in read_csv("pneumatophores_truth_cells.csv"):
        truth[(round(float(r["x_min"])), round(float(r["y_min"])), round(float(r["z_min"]), 2))] = (
            float(r["true_frontal_area_per_volume_Nd"]), float(r["true_trunk_frontal_area_per_volume"]))
    trunks = np.array([[float(r["x"]), float(r["y"]), float(r["diameter"])] for r in read_csv("pneumatophores_trunks.csv")])
    in_plot = (np.abs(trunks[:, 0]) < 3) & (np.abs(trunks[:, 1]) < 3)
    trunk_cells = {(int(np.floor(x)), int(np.floor(y))) for x, y, _ in trunks[in_plot]}
    print(f"{in_plot.sum()} trunks in the 6 x 6 m plot, diameters {np.round(trunks[in_plot, 2], 3)}; "
          f"trunk frontal area is {sum(t[1] for t in truth.values() if t[1] > 0) / max(1, sum(t[1] > 0 for t in truth.values())):.3f} m^-1 "
          f"in a cell holding a trunk")
    header, data = read_ply("pneumatophores.ply")
    is_trunk = (data["red"] == 150) & (data["green"] == 110)
    is_stem = (data["red"] == 40) & (data["green"] == 160) & (data["alpha"] > 0)

    est_a = estimate("pneumatophores.ply", "cells_a.csv")
    masked = data.copy()
    masked["alpha"][is_trunk] = 0
    write_ply("masked_b.ply", header, masked)
    est_b = estimate("masked_b.ply", "cells_b.csv")

    # C: detect trunks above the pneumatophores, then mask returns within the fitted circles at all heights
    top = float(np.max(data["z"][is_stem]))
    detected = detect_trunks(data, (top + 0.05, top + 0.65))
    near = (np.abs(trunks[:, 0]) < 4) & (np.abs(trunks[:, 1]) < 4)
    print(f"pneumatophore returns reach {top:.2f} m; trunks detected from returns at {top + 0.05:.2f}-{top + 0.65:.2f} m: "
          f"{len(detected)} (true trunks within 4 m of the plot centre: {near.sum()})")
    for tx, ty, d in trunks[near]:
        if not detected:
            break
        k = int(np.argmin([np.hypot(cx - tx, cy - ty) for cx, cy, _ in detected]))
        cx, cy, r = detected[k]
        print(f"  trunk at ({tx:6.2f},{ty:6.2f}) d {d:.3f}: centre error {100 * np.hypot(cx - tx, cy - ty):.1f} cm, "
              f"fitted diameter {2 * r:.3f} (+{100 * (2 * r - d):.1f} cm, the footprint widening)")
    x, y, z = data["x"].astype(float), data["y"].astype(float), data["z"].astype(float)
    remove = np.zeros(len(data), dtype=bool)
    for cx, cy, r in detected:
        remove |= np.hypot(x - cx, y - cy) < r + args.margin
    remove &= data["alpha"] > 0
    low = z < 0.3
    print(f"bottom 0.3 m: {is_trunk[low].sum()} trunk returns, {is_stem[low].sum()} pneumatophore returns. Separation: "
          f"{100 * (remove & is_trunk & low).sum() / max(1, (is_trunk & low).sum()):.1f}% of trunk returns removed, "
          f"{100 * (remove & is_stem & low).sum() / max(1, (is_stem & low).sum()):.2f}% of pneumatophore returns removed")
    masked = data.copy()
    masked["alpha"][remove] = 0
    write_ply("masked_c.ply", header, masked)
    est_c = estimate("masked_c.ply", "cells_c.csv")

    summarise("A: all returns", est_a, truth, trunk_cells)
    summarise("B: trunk returns removed (perfect separation)", est_b, truth, trunk_cells)
    summarise("C: trunks detected above the pneumatophores and removed", est_c, truth, trunk_cells)


if __name__ == "__main__":
    main()
