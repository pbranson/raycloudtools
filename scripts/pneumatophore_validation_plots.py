#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Validation plots for the pneumatophore density workflow, on synthetic hand-held scans. Requires raycreate and
rayrender on the PATH, and numpy and matplotlib. Run in an empty working directory; writes fig*.png.

Optionally pass a file of printed models (from a small program calling ray::effectiveStemWidth and
ray::freePathResponse, with columns kind,threshold_or_rho,x,value) with --models, for the model figure.
"""
import argparse
import csv
import os
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

# palette (validated reference palette): categorical slots, sequential blue ramp, blue-red diverging
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list("seq", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f0efec", "#ef8a8a", "#a32c2c"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": INK, "axes.titlecolor": INK, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8, "legend.frameon": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2, "font.family": "sans-serif",
})
FOOTPRINT = "0.01,0.003,0.05"
DIAMETER = 0.01
HALF = 3.0  # half the plot width (m)


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


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def read_csv(name):
    with open(name) as f:
        return list(csv.DictReader(f))


def read_rays(name):
    """Reads a ray cloud .ply (float or double precision), returning starts, ends, times and colours"""
    with open(name, "rb") as f:
        header = b""
        while not header.endswith(b"end_header\n"):
            header += f.readline()
        data = np.fromfile(f, dtype=ply_dtype(header))
    ray = ("nx", "ny", "nz") if "nx" in data.dtype.names else ("rayx", "rayy", "rayz")
    ends = np.stack([data["x"], data["y"], data["z"]], axis=1).astype(np.float64)
    starts = ends + np.stack([data[r] for r in ray], axis=1)  # the ray vector is start - end
    colours = np.stack([data["red"], data["green"], data["blue"], data["alpha"]], axis=1)
    return starts, ends, data["time"], colours

def scene(name, options, seed=11):
    """Generates a scene and its estimates (pooled estimator, 5 cm voxels, 1 x 1 x 0.1 m cells) in a subdirectory"""
    os.makedirs(name, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(name)
    if not os.path.exists("cells.csv"):
        run(f"raycreate pneumatophores {seed} {options} --footprint {FOOTPRINT} --truth_cells 1,1,0.1")
        run(f"rayrender pneumatophores.ply top density --density_model vertical_cylinder --estimator pooled "
            f"--no_flat_top --stem_diameter {DIAMETER} --footprint {FOOTPRINT} --pixel_width 0.05 "
            f"--min_bound -3,-3,0 --max_bound 3,3,0.6 --cells cells.csv --profile profile.csv")
    stems = np.array([[float(r["x"]), float(r["y"]), float(r["height"])] for r in read_csv("pneumatophores_stems.csv")])
    out = {"stems": stems, "rays": read_rays("pneumatophores.ply"), "cells": read_csv("cells.csv"),
           "truth_cells": read_csv("pneumatophores_truth_cells.csv"), "profile": read_csv("profile.csv"),
           "truth_profile": read_csv("pneumatophores_truth.csv")}
    os.chdir(cwd)
    return out


def classify(colours):
    """stem hits (bounded), ground returns (zero alpha, brown) and no return (zero alpha, grey)"""
    bounded = colours[:, 3] > 0
    brown = (colours[:, 0] == 120) & (colours[:, 1] == 80)
    return bounded, (~bounded) & brown, (~bounded) & (~brown)


def fig_scenes(scenes):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4), constrained_layout=True)
    for ax, (title, sc) in zip(axes, scenes.items()):
        s = sc["stems"]
        inside = (np.abs(s[:, 0]) < HALF + 1) & (np.abs(s[:, 1]) < HALF + 1)
        ax.scatter(s[inside, 0], s[inside, 1], s=0.4, color=AQUA, linewidths=0, label="stems", rasterized=True)
        starts = sc["rays"][0][::2000]
        ax.plot(starts[:, 0], starts[:, 1], color=ORANGE, lw=1.2, label="walking track")
        for v in np.arange(-HALF, HALF + 0.01, 1.0):
            ax.axvline(v, color=AXIS, lw=0.5, zorder=0)
            ax.axhline(v, color=AXIS, lw=0.5, zorder=0)
        ax.plot([-HALF, HALF, HALF, -HALF, -HALF], [-HALF, -HALF, HALF, HALF, -HALF], color=INK2, lw=1, ls="--",
                label="plot (1 m cells)")
        in_plot = (np.abs(s[:, 0]) < HALF) & (np.abs(s[:, 1]) < HALF)
        ax.set_title(f"{title}: {in_plot.sum() / (2 * HALF) ** 2:.0f} stems/m² in plot")
        ax.set_aspect("equal")
        ax.set_xlim(-HALF - 1, HALF + 1)
        ax.set_ylim(-HALF - 1, HALF + 1)
        ax.set_xlabel("x (m)")
        ax.grid(False)
    axes[0].set_ylabel("y (m)")
    axes[0].legend(loc="upper left", markerscale=12)
    fig.suptitle("Synthetic stem layouts and lidar walking track (stems are non-overlapping, 1 cm diameter)")
    fig.savefig("fig1_scenes.png", dpi=150)
    plt.close(fig)


def fig_side(sc):
    starts, ends, _, colours = sc["rays"]
    s = sc["stems"]
    y0 = 0.5 * (starts[0, 1] + starts[0, 1])  # the first walking line
    fig, axes = plt.subplots(2, 1, figsize=(14, 7), constrained_layout=True, height_ratios=[2, 1])
    ax = axes[0]
    # rays nearly in the vertical x-z plane of the first walking line, from the first 2 m of the line
    d = ends - starts
    horiz = np.linalg.norm(d[:, :2], axis=1) + 1e-12
    sel = (np.abs(starts[:, 1] - y0) < 0.02) & (np.abs(d[:, 1]) / horiz < 0.02) & (starts[:, 0] < -1.0)
    idx = np.nonzero(sel)[0][:400]
    hit, ground, none = classify(colours[idx])
    for mask, col, lab in [(none, MUTED, "no return (max range)"), (ground, ORANGE, "ground return (alpha 0)"),
                           (hit, BLUE, "stem return")]:
        first = True
        for i in idx[mask]:
            ax.plot([starts[i, 0], ends[i, 0]], [starts[i, 2], ends[i, 2]], color=col, lw=0.5, alpha=0.6,
                    label=lab if first else None)
            first = False
    ax.scatter(ends[idx[hit], 0], ends[idx[hit], 2], s=8, color=BLUE, zorder=3)
    slab = np.abs(s[:, 1] - y0) < 0.05
    for x, h in s[slab][:, [0, 2]]:
        ax.plot([x, x], [0, h], color=AQUA, lw=1.2)
    ax.plot([], [], color=AQUA, lw=1.2, label="stems within ±5 cm of the plane")
    ax.axhline(0, color=INK2, lw=1)
    ax.set_xlim(-3, 1)
    ax.set_ylim(-0.05, 1.35)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("height (m)")
    ax.set_title(f"Side view: rays close to the vertical plane y = {y0:.2f} m (a walking line)")
    ax.legend(loc="upper right", ncol=4)
    ax = axes[1]
    ax.hist(s[:, 2], bins=40, color=AQUA, edgecolor=SURFACE, linewidth=1)
    ax.set_xlabel("stem height (m)")
    ax.set_ylabel("stems")
    ax.set_title("Stem height distribution (mean 0.25 m, s.d. 0.08 m, truncated at 0.02 m)")
    fig.savefig("fig2_side_view.png", dpi=150)
    plt.close(fig)


def fig_sampling(sc):
    starts, ends, times, colours = sc["rays"]
    s = sc["stems"]
    hit, ground, none = classify(colours)
    d = ends - starts
    length = np.linalg.norm(d, axis=1)
    sin_zenith = np.linalg.norm(d[:, :2], axis=1) / np.maximum(length, 1e-12)
    zenith = np.degrees(np.arcsin(np.clip(sin_zenith, 0, 1)))
    zenith = np.where(d[:, 2] > 0, 180 - zenith, zenith)
    fig, axes = plt.subplots(2, 2, figsize=(13, 11), constrained_layout=True)
    ax = axes[0, 0]
    sub = slice(None, None, 500)
    sc_ = ax.scatter(starts[sub, 0], starts[sub, 1], c=times[sub], cmap=SEQ, s=4, linewidths=0)
    fig.colorbar(sc_, ax=ax, label="time (s)")
    ax.set_aspect("equal")
    ax.set_title("Scanner positions along the track, by time")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax = axes[0, 1]
    bins = np.linspace(0, 105, 43)
    ax.hist(zenith, bins=bins, color=MUTED, alpha=0.5, label="all rays", density=True)
    ax.hist(zenith[hit], bins=bins, color=BLUE, histtype="step", lw=2, label="stem returns", density=True)
    ax.hist(zenith[ground], bins=bins, color=ORANGE, histtype="step", lw=2, label="ground returns", density=True)
    ax.axvline(np.degrees(np.arcsin(0.1)), color=INK2, lw=1, ls="--")
    ax.text(np.degrees(np.arcsin(0.1)) + 1, ax.get_ylim()[1] * 0.9, "rays steeper than this\nare ignored",
            fontsize=8, color=INK2)
    ax.set_xlabel("ray zenith angle (degrees from straight down)")
    ax.set_ylabel("density")
    ax.set_title("Ray directions")
    ax.legend()
    ax = axes[1, 0]
    rb = np.linspace(0, 8, 41)
    ax.hist(length[hit], bins=rb, color=BLUE, edgecolor=SURFACE, linewidth=1)
    ax.set_xlabel("range to stem return (m)")
    ax.set_ylabel("returns")
    ax.set_title("Ranges of stem returns (the footprint grows with range)")
    ax = axes[1, 1]
    # close-up: hit points at 0-0.3 m height should lie on the effective (footprint-widened) stem outlines
    cx, cy, w = 0.0, starts[0, 1] + 0.6, 0.3
    near = (np.abs(ends[:, 0] - cx) < w) & (np.abs(ends[:, 1] - cy) < w) & hit & (ends[:, 2] < 0.3)
    ax.scatter(ends[near, 0], ends[near, 1], s=3, color=BLUE, linewidths=0, label="stem returns (z < 0.3 m)")
    stems_near = s[(np.abs(s[:, 0] - cx) < w + 0.05) & (np.abs(s[:, 1] - cy) < w + 0.05)]
    for x, y, _ in stems_near:
        ax.add_patch(plt.Circle((x, y), DIAMETER / 2, color=AQUA, zorder=3))
        for rng, ls in [(1.0, ":"), (3.0, "--")]:
            beam = 0.01 + 0.003 * rng
            eff = effective_width(DIAMETER, beam, 0.05)
            ax.add_patch(plt.Circle((x, y), eff / 2, fill=False, color=INK2, lw=0.6, ls=ls))
    ax.plot([], [], color=INK2, ls=":", lw=0.8, label="effective outline at 1 m range")
    ax.plot([], [], color=INK2, ls="--", lw=0.8, label="effective outline at 3 m range")
    ax.scatter([], [], color=AQUA, s=30, label="true stem (1 cm)")
    ax.set_xlim(cx - w, cx + w)
    ax.set_ylim(cy - w, cy + w)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title("Close-up: returns sit on the footprint-widened outlines")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7)
    fig.suptitle("Lidar sampling (dense random scene)")
    fig.savefig("fig3_sampling.png", dpi=150)
    plt.close(fig)


def effective_width(diameter, beam, threshold):
    """Python copy of ray::effectiveStemWidth (Gaussian beam, strip interception), for the close-up figure"""
    from math import erf, sqrt
    sigma = beam / 4.0

    def frac(off):
        k = 1.0 / (sqrt(2.0) * sigma)
        return 0.5 * (erf((off + diameter / 2) * k) - erf((off - diameter / 2) * k))

    if frac(0) < threshold:
        return 0.0
    lo, hi = 0.0, diameter / 2 + 4 * beam
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if frac(mid) >= threshold else (lo, mid)
    return 2 * lo


def fig_models(models_file):
    rows = read_csv(models_file)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    ax = axes[0]
    for t, col in zip(["0.02", "0.05", "0.15"], [VIOLET, BLUE, ORANGE]):
        r = [(float(x["x"]), float(x["value"]) * 100) for x in rows if x["kind"] == "width" and x["threshold_or_rho"] == t]
        ax.plot(*zip(*r), color=col, label=f"threshold {t}")
    ranges = np.linspace(0, 10, 50)
    ax.plot(ranges, (0.01 + 0.003 * ranges) * 100, color=MUTED, lw=1, ls="--", label="beam diameter (1/e²)")
    ax.axhline(1.0, color=AQUA, lw=1, ls=":", label="true stem diameter")
    ax.set_xlabel("range (m)")
    ax.set_ylabel("effective stem width (cm)")
    ax.set_title("Footprint model: 10 mm exit, 3 mrad divergence, 1 cm stem")
    ax.legend()
    ax = axes[1]
    for rho, col in zip(["1", "1.5", "2", "2.5", "3"], [AQUA, BLUE, VIOLET, ORANGE, "#e34948"]):
        r = [(float(x["x"]), float(x["value"])) for x in rows if x["kind"] == "response" and x["threshold_or_rho"] == rho]
        ax.plot(*zip(*r), color=col, marker="o", ms=4, label=f"d_eff / d = {rho}")
    f = np.linspace(0, 0.4, 50)
    ax.plot(f, 1 / (1 - f), color=MUTED, lw=1, ls="--", label="1 / (1 - solid fraction)")
    ax.set_xlabel("solid fraction  N π d² / 4")
    ax.set_ylabel("measured / true N·d")
    ax.set_title("Free-path response to non-overlapping stems (inverted by the correction)")
    ax.legend()
    fig.savefig("fig4_models.png", dpi=150)
    plt.close(fig)


def cell_arrays(sc, layer):
    """Maps of truth, mean-of-voxel estimate and observed fraction for one 0.1 m layer of 1 m cells"""
    truth = {}
    for r in sc["truth_cells"]:
        if abs(float(r["z_min"]) - 0.1 * layer) < 1e-6:
            truth[(round(float(r["x_min"])), round(float(r["y_min"])))] = float(r["true_frontal_area_per_volume_Nd"])
    est = np.full((6, 6), np.nan)
    tru = np.full((6, 6), np.nan)
    col = [k for k in sc["cells"][0].keys() if k.startswith("mean_density")][0]
    for r in sc["cells"]:
        if abs(float(r["z_min"]) - 0.1 * layer) < 1e-6:
            k = (round(float(r["x_min"])), round(float(r["y_min"])))
            if k in truth:
                i, j = k[1] + 3, k[0] + 3
                est[i, j] = float(r[col])
                tru[i, j] = truth[k]
    return tru, est


def fig_maps(scenes):
    names = list(scenes.keys())
    fig, axes = plt.subplots(len(names), 3, figsize=(13, 4.1 * len(names)), constrained_layout=True)
    for row, name in enumerate(names):
        tru, est = cell_arrays(scenes[name], 0)
        vmax = np.nanmax(np.concatenate([tru.ravel(), est.ravel()]))
        for ax, data, title in [(axes[row, 0], tru, "truth"), (axes[row, 1], est, "estimate")]:
            im = ax.imshow(data, origin="lower", extent=[-3, 3, -3, 3], cmap=SEQ, vmin=0, vmax=vmax)
            ax.set_title(f"{name}: {title}, N·d at 0-0.1 m (m⁻¹)")
            ax.grid(False)
            fig.colorbar(im, ax=ax, shrink=0.85)
        ratio = est / np.where(tru > 0.2, tru, np.nan)
        im = axes[row, 2].imshow(ratio, origin="lower", extent=[-3, 3, -3, 3], cmap=DIV,
                                 norm=TwoSlopeNorm(vcenter=1.0, vmin=0.7, vmax=1.3))
        axes[row, 2].set_title("estimate / truth (blank: truth < 0.2 m⁻¹)")
        axes[row, 2].grid(False)
        fig.colorbar(im, ax=axes[row, 2], shrink=0.85)
        for ax in axes[row]:
            ax.set_xlabel("x (m)")
            ax.set_ylabel("y (m)")
    fig.suptitle("1 × 1 m cells, lowest 0.1 m layer (pooled estimator, 5 cm voxels, footprint corrected)")
    fig.savefig("fig5_cell_maps.png", dpi=150)
    plt.close(fig)


def fig_profiles(scenes):
    fig, axes = plt.subplots(1, len(scenes) + 1, figsize=(5 * (len(scenes) + 1), 5), constrained_layout=True)
    for ax, (name, sc) in zip(axes, scenes.items()):
        tp = sc["truth_profile"]
        z = [0.5 * (float(r["z_min"]) + float(r["z_max"])) for r in tp]
        ax.plot([float(r["true_frontal_area_per_volume_Nd"]) for r in tp], z, color=INK2, lw=1.5, label="truth")
        pr = sc["profile"]
        col = [k for k in pr[0].keys() if k.startswith("mean_density")][0]
        zp = [0.5 * (float(r["z_min"]) + float(r["z_max"])) for r in pr]
        ax.plot([float(r[col]) for r in pr], zp, color=BLUE, marker="o", ms=5, label="estimate (5 cm layers)")
        ax.set_ylim(0, 0.6)
        ax.set_xlabel("N·d (m⁻¹)")
        ax.set_title(name)
        ax.legend()
    axes[0].set_ylabel("height (m)")
    ax = axes[-1]
    for (name, sc), col in zip(scenes.items(), [BLUE, ORANGE, AQUA]):
        xs, ys = [], []
        for layer in range(3):
            tru, est = cell_arrays(sc, layer)
            ok = ~np.isnan(tru) & ~np.isnan(est)
            xs += list(tru[ok])
            ys += list(est[ok])
        ax.scatter(xs, ys, s=12, color=col, label=name, alpha=0.8)
    lim = ax.get_xlim()[1]
    ax.plot([0, lim], [0, lim], color=MUTED, lw=1, ls="--", label="1:1")
    ax.set_xlabel("true cell N·d (m⁻¹)")
    ax.set_ylabel("estimated cell N·d (m⁻¹)")
    ax.set_title("All 1 × 1 × 0.1 m cells, 0-0.3 m")
    ax.legend()
    fig.suptitle("Plot-mean profiles and cell estimates")
    fig.savefig("fig6_profiles.png", dpi=150)
    plt.close(fig)


def fig_methods():
    densities = [30, 100, 300, 600, 1000, 1500]
    methods = {"pooled, occlusion corrected (recommended)": ("pooled", "", "mean", BLUE),
               "pooled, no occlusion correction": ("pooled", "--no_overlap_correction", "mean", VIOLET),
               "crossing (no correction)": ("crossing", "", "mean", ORANGE),
               "single 1×1×0.1 m beer_lambert voxel": ("beer_lambert", "--no_overlap_correction", "pooled", AQUA)}
    results = {m: [] for m in methods}
    os.makedirs("density_sweep", exist_ok=True)
    cwd = os.getcwd()
    os.chdir("density_sweep")
    for n in densities:
        run(f"raycreate pneumatophores 11 --stems_per_m2 {n} --footprint {FOOTPRINT} --truth_cells 1,1,0.1")
        truth = {(r["x_min"], r["y_min"]): float(r["true_frontal_area_per_volume_Nd"])
                 for r in read_csv("pneumatophores_truth_cells.csv") if float(r["z_min"]) < 0.05}
        for m, (est, extra, column, _) in methods.items():
            run(f"rayrender pneumatophores.ply top density --density_model vertical_cylinder --estimator {est} "
                f"--no_flat_top --stem_diameter {DIAMETER} --footprint {FOOTPRINT} --pixel_width 0.05 "
                f"--min_bound -3,-3,0 --max_bound 3,3,0.3 --cells cells.csv {extra}")
            cells = read_csv("cells.csv")
            col = [k for k in cells[0].keys() if k.startswith(column)][0]
            pairs = [(float(c[col]), truth[(c["x_min"], c["y_min"])]) for c in cells
                     if float(c["z_min"]) < 0.05 and (c["x_min"], c["y_min"]) in truth]
            results[m].append(sum(e for e, _ in pairs) / sum(t for _, t in pairs))
    os.chdir(cwd)
    fig, ax = plt.subplots(figsize=(9, 5.2), constrained_layout=True)
    ax.axhspan(0.95, 1.05, color=GRID, alpha=0.6, lw=0)
    ax.axhline(1.0, color=MUTED, lw=1)
    for m, (_, _, _, col) in methods.items():
        ax.plot(densities, results[m], color=col, marker="o", ms=6, label=m)
        ax.annotate(f"{results[m][-1]:.2f}", (densities[-1], results[m][-1]), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=8, color=INK2)
    ax.set_xscale("log")
    ax.set_xticks(densities)
    ax.set_xticklabels([str(n) for n in densities])
    ax.set_xlabel("stems per m² (1 cm diameter, random, non-overlapping)")
    ax.set_ylabel("estimate / truth, lowest 0.1 m")
    ax.set_title("Estimator bias against stem density, with a 10 mm, 3 mrad beam (shaded: ±5%)")
    ax.legend(loc="lower left")
    fig.savefig("fig7_methods.png", dpi=150)
    plt.close(fig)
    with open("density_sweep.csv", "w") as f:
        f.write("method," + ",".join(str(n) for n in densities) + "\n")
        for m in methods:
            f.write(f"\"{m}\"," + ",".join(f"{v:.4f}" for v in results[m]) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--models", default="", help="csv of printed effective width and response models")
    args = parser.parse_args()
    scenes = {"dense random": scene("dense_random", "--stems_per_m2 300"),
              "patches (r 0.5 m, 30% cover)": scene("patches", "--stems_per_m2 300 --patchy 0.5,0.3,0"),
              "cable roots": scene("cable_roots", "--stems_per_m2 150 --clustered --cluster_spread 0.03")}
    fig_scenes(scenes)
    fig_side(scenes["dense random"])
    fig_sampling(scenes["dense random"])
    if args.models:
        fig_models(args.models)
    fig_maps({k: scenes[k] for k in list(scenes)[:2]})
    fig_profiles(scenes)
    fig_methods()


if __name__ == "__main__":
    main()
