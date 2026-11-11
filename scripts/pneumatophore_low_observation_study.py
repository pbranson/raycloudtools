#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Poorly observed cells, as in Section 8.3 of docs/pneumatophore_density_report.md. Requires raycreate and
rayrender on the PATH (or --bin), and numpy.

Each scene (raycreate pneumatophores options, on a 4 x 4 m plot with 56 250 rays per metre of walk) is estimated
with 5 cm voxel cells (one cell per voxel) without and with --self_shadow_correction. Per 1 x 1 x 0.1 m cell with
observed voxels it separates:
  est/truth          the cell's density (mean of its observed voxels) against the cell truth
  sel                the mean truth of the observed voxels' 3 x 3 windows over the cell truth: the selection of
                     less dense (less occluded) places
  est/truth_obs      the estimate against the truth of the windows it samples
  filled/truth       with unobserved voxels given the mean of their observed 3 x 3 neighbours (if any)
for the variants base (uncorrected), ss (--self_shadow_correction) and cw (the self-shadow factor scaled by the
window's observed voxel fraction). It prints these per scene and layer, and binned by the cell's observed fraction
over the layers below the stem tops (all layers for the fixed 0.6 m and 2 m stems, else the lowest 0.2 m). Writes
<out>/cells.csv. With no --scene, runs the scenes of the report (seeds 101, and 102 for the two densest). Example:
  scripts/pneumatophore_low_observation_study.py --out low_observation
"""
import argparse
import collections
import csv
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np

FOOTPRINT = "0.01,0.003,0.05"
PLOT = 4.0
VOXEL = 0.05
ZMAX = 0.4
# label, seed, raycreate options
SCENES = [("h06_1000", 101, "--stems_per_m2 1000 --height 0.6,0"),
          ("h06_1000", 102, "--stems_per_m2 1000 --height 0.6,0"),
          ("def_1500", 101, "--stems_per_m2 1500"),
          ("def_1500", 102, "--stems_per_m2 1500"),
          ("notop_300", 101, "--stems_per_m2 300 --height 2,0"),
          ("walk2_1500", 101, "--stems_per_m2 1500 --walk_spacing 2"),
          ("def_300", 101, "--stems_per_m2 300"),
          ("walk2_600", 101, "--stems_per_m2 600 --walk_spacing 2 --height 0.4,0.1")]


def run(cmd, cwd):
    subprocess.run(cmd, shell=True, check=True, cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def box(a, n):
    """mean over the 3 x 3 horizontal window of each voxel (clipped at the plot edge), ignoring NaN"""
    p = np.pad(a, ((0, 0), (1, 1), (1, 1)), constant_values=np.nan)
    s = np.zeros_like(a)
    c = np.zeros_like(a)
    for dj in range(3):
        for di in range(3):
            v = p[:, dj:dj + n, di:di + n]
            s += np.nan_to_num(v)
            c += ~np.isnan(v)
    return np.where(c > 0, s / np.maximum(c, 1), np.nan)


def run_scene(job):
    label, seed, scene, rays_per_m, tmp = job
    work = tempfile.mkdtemp(prefix="lowobs_", dir=tmp)
    try:
        walk = 1.0
        if "--walk_spacing" in scene:
            walk = float(scene.split("--walk_spacing")[1].split()[0])
        rays = int(rays_per_m * PLOT * max(1, int(PLOT / walk)))
        if "--footprint" not in scene:
            scene += f" --footprint {FOOTPRINT}"
        run(f"raycreate pneumatophores {seed} {scene} --plot_width {PLOT} --rays {rays}", work)
        n = int(round(PLOT / VOXEL))
        nz = int(round(ZMAX / VOXEL))
        half = 0.5 * PLOT
        # voxel truth N d from the stems
        truth = np.zeros((nz, n, n))
        with open(os.path.join(work, "pneumatophores_stems.csv")) as f:
            for r in csv.DictReader(f):
                x, y, h = float(r["x"]), float(r["y"]), float(r["height"])
                if abs(x) >= half or abs(y) >= half:
                    continue
                i = int((x + half) / VOXEL)
                j = int((y + half) / VOXEL)
                for k in range(nz):
                    truth[k, j, i] += 0.01 * max(0.0, min(h, (k + 1) * VOXEL) - k * VOXEL) / (VOXEL ** 3)
        twin = box(truth, n)
        out = []
        base_dens = None
        for ss in ("", "--self_shadow_correction"):
            run(f"rayrender pneumatophores.ply top density --density_model vertical_cylinder --no_flat_top "
                f"--estimator pooled --stem_diameter 0.01 --footprint {FOOTPRINT} --pixel_width {VOXEL} "
                f"--min_bound {-half},{-half},0 --max_bound {half},{half},{ZMAX} --cells vox.csv "
                f"--cell_size {VOXEL},{VOXEL},{VOXEL} {ss}", work)
            dens = np.full((nz, n, n), np.nan)
            with open(os.path.join(work, "vox.csv")) as f:
                for r in csv.DictReader(f):
                    i = int(round((float(r["x_min"]) + half) / VOXEL))
                    j = int(round((float(r["y_min"]) + half) / VOXEL))
                    k = int(round(float(r["z_min"]) / VOXEL))
                    if k < nz and i < n and j < n:
                        dens[k, j, i] = float(r["density_frontal_area_per_volume_Nd"])
            obs = ~np.isnan(dens)
            filled = np.where(obs, dens, box(np.where(obs, dens, np.nan), n))
            variants = [("ss" if ss else "base", dens, filled)]
            if not ss:
                base_dens = dens.copy()
            else:
                # coverage-weighted self-shadow factor: scaled by the window's observed voxel fraction
                fobs = box(obs.astype(float), n)
                ratio = np.where(dens > 0, base_dens / np.where(dens > 0, dens, 1.0), 1.0)
                cw = base_dens / (1.0 + fobs * (ratio - 1.0))
                variants.append(("cw", cw, filled))
            cpv = int(round(1.0 / VOXEL))
            lpv = int(round(0.1 / VOXEL))
            for vname, dd, ff in variants:
                for layer in range(int(round(ZMAX / 0.1))):
                    for cj in range(int(PLOT)):
                        for ci in range(int(PLOT)):
                            sl = (slice(layer * lpv, (layer + 1) * lpv), slice(cj * cpv, (cj + 1) * cpv),
                                  slice(ci * cpv, (ci + 1) * cpv))
                            t = truth[sl].mean()
                            o = obs[sl]
                            if t <= 0 or o.sum() == 0:
                                continue
                            est = dd[sl][o].mean()
                            tobs = twin[sl][o].mean()
                            f = ff[sl]
                            fv = f[~np.isnan(f)].mean()
                            out.append(dict(label=label, seed=seed, var=vname, z=round(layer * 0.1, 1), ci=ci, cj=cj,
                                            obs=o.mean(), truth=t, est=est / t, sel=tobs / t, est_obs=est / tobs,
                                            filled=fv / t))
        return out
    finally:
        shutil.rmtree(work, ignore_errors=True)


def below_tops(r):
    # layers without stem tops: all for the fixed 0.6 and 2 m stems, else the lowest 0.2 m
    return r["label"].startswith(("h06", "notop")) or r["z"] < 0.15


def summarise(rows):
    g = collections.OrderedDict()
    for r in sorted(rows, key=lambda r: (r["label"], r["var"], r["z"])):
        g.setdefault((r["label"], r["var"], r["z"]), []).append(r)
    print("scene | variant | z | cells | observed | est/truth | sel (truth_obs/truth) | est/truth_obs | filled/truth")
    for k, rs in g.items():
        t = sum(r["truth"] for r in rs)

        def w(key):
            return sum(r[key] * r["truth"] for r in rs) / t
        print(f"{k[0]} | {k[1]} | {k[2]:.1f} | {len(rs)} | {sum(r['obs'] for r in rs) / len(rs):.2f} | "
              f"{w('est'):.3f} | {w('sel'):.3f} | {w('est') / w('sel'):.3f} | {w('filled'):.3f}")
    print("\nbinned by cell observed fraction (all scenes, layers below the stem tops): variant | bin | cells | "
          "est/truth | sel | est/truth_obs | filled")
    bins = [0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 1.01]
    for var in ("base", "ss", "cw"):
        for lo, hi in zip(bins[:-1], bins[1:]):
            rs = [r for r in rows if r["var"] == var and lo <= r["obs"] < hi and below_tops(r)]
            if not rs:
                continue
            t = sum(r["truth"] for r in rs)

            def w(key):
                return sum(r[key] * r["truth"] for r in rs) / t
            print(f"{var} | {lo:.2f}-{min(hi, 1):.2f} | {len(rs)} | {w('est'):.3f} | {w('sel'):.3f} | "
                  f"{w('est') / w('sel'):.3f} | {w('filled'):.3f}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scene", action="append", default=[],
                   help="'label,seed,raycreate options' (repeatable); default: the report's scenes")
    p.add_argument("--rays_per_m", type=float, default=56250.0, help="rays per metre of walk")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--out", default="low_observation")
    p.add_argument("--bin", default="", help="comma separated directories added to the PATH (e.g. a build's tools)")
    a = p.parse_args()
    if a.bin:
        os.environ["PATH"] = os.pathsep.join(os.path.abspath(b) for b in a.bin.split(",")) + os.pathsep + \
            os.environ["PATH"]
    os.makedirs(a.out, exist_ok=True)
    scenes = SCENES
    if a.scene:
        scenes = [(s.split(",", 2)[0], int(s.split(",", 2)[1]), s.split(",", 2)[2]) for s in a.scene]
    jobs = [(label, seed, scene, a.rays_per_m, os.path.abspath(a.out)) for label, seed, scene in scenes]
    rows = []
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        for result in pool.map(run_scene, jobs):
            rows += result
    with open(os.path.join(a.out, "cells.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    summarise(rows)


if __name__ == "__main__":
    main()
