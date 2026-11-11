#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Wider synthetic validation of the pooled stem density estimator (rayrender --estimator pooled), with and without
--self_shadow_correction, over scan geometries and error sources: blocks A-J of the stage-2 validation matrix.
Requires raycreate and rayrender on the PATH (or --bin).

Each configuration is a scene (raycreate pneumatophores options), a set of render variants (rayrender options) and
the height layers to evaluate, run over several seeds. Every scene variant runs the uncorrected estimator ("base") and
the self-shadow correction ("ss"). Fixed settings: 4 x 4 m plot, 8 m stem margin, 56 250 rays per metre of walk,
1 cm stems (hard-core RSA, heights 0.25 +- 0.08 m) at 300 /m^2, the default beam (10 mm exit, 3 mrad, 5% threshold)
with returns at the stem surface, 0.05 m voxels pooled over 0.05 m, 1 x 1 x 0.1 m evaluation cells.

Each (configuration, seed) is a job whose cells are saved to <out>/jobs/<block>/<config id>_<seed>.json. Existing
job files are skipped, so a run can be stopped and resumed, and further seeds added later. The summary (from all the
job files present) is written to <out>/summary.csv, <out>/bins.csv and <out>/summary.txt. Per configuration, variant
and layer, as mean +- standard error over the seeds:
  d/t, d/t90      density (mean of the voxel densities) / truth, summed over all the plot's cells, and over the
                  cells with observed fraction >= 0.9 (the pass metric)
  p/t, p/t90      the same for the pooled cell density (hits / path over the cell's voxels)
  d/p             density / pooled_density (the self-shadow signature for uniform stems; patchiness otherwise)
  rmse, in10      cell relative RMSE (all cells) and the fraction of cells within +-10% of the truth
  obs, obs_min    mean and minimum observed fraction of the cells
  ss              rayrender's mean self-shadow factor, and H1 = 1 + 0.4732 d_eff / W (W the window width)
  t_gen, t_est    run time (s) of the scene generation and of each estimate
bins.csv gives density / truth in bins of observed fraction (pooled over seeds and cells), and line_distance.csv by
distance of the cell centre from the nearest walk line (block A). Block X (run only when named) holds the diagnostic
runs of the failures. Block L (run only when named) is the footprint mismatch: scenes generated with a true beam
footprint (default 10 mm, 3 mrad, 0.05; narrow 5 mm, 0.5 mrad; wide 20 mm, 5 mrad) and estimated assuming another
(exit diameter and divergence x0.5-2, thresholds 0.02-0.15, joint errors that match d_eff at 2 m, an assumed
diameter x threshold grid, and no footprint correction), each with the default estimator, the self-shadow correction
and the small-window configuration (SW_EXTRA); analysed by scripts/pneumatophore_footprint_mismatch.py. Pass: |d/t90 - 1| <= 0.02 in 0-0.1 and 0.1-0.2 m for the pass-eligible configurations;
--rerun_failures adds seeds to the failing ones. Examples:
  scripts/pneumatophore_validation_matrix.py --blocks A,B --out matrix --bin build/raycloudtools/raycreate,...
  scripts/pneumatophore_validation_matrix.py --blocks all --summary_only --out matrix
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

FOOTPRINT = "0.01,0.003,0.05"
EXIT_DISTANCE_SQUARE = (4.0 * math.log(1.0 + math.sqrt(2.0)) - (4.0 / 3.0) * (math.sqrt(2.0) - 1.0)) / (2.0 * math.pi)
SEEDS = "101,102,103,104"
OBS_BINS = [0.0, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0001]
PASS_BLOCKS = set("ABCDEHIJ")
# the small-window configuration of stage 2 (remedy d), with --pool_radius 0
SW_EXTRA = ("--self_shadow_correction --self_shadow_second_order --occlusion_radius 0.05 --ratio_bias_correction "
            "--no_neighbour_priors")
PASS_TOLERANCE = 0.02
LOW = [0.0, 0.1]


def layers_to(top):
    return [round(0.1 * i, 1) for i in range(int(round(top / 0.1)))]


def strip_intercept_fraction(width, beam_diameter, offset):
    """as ray::stripInterceptFraction: the fraction of a Gaussian beam (1/e^2 diameter) intercepted by a strip"""
    if beam_diameter <= 0.0:
        return 1.0 if abs(offset) < 0.5 * width else 0.0
    scale = 1.0 / (math.sqrt(2.0) * 0.25 * beam_diameter)
    return 0.5 * (math.erf((offset + 0.5 * width) * scale) - math.erf((offset - 0.5 * width) * scale))


def effective_width(diameter, beam_diameter, threshold):
    """as ray::effectiveStemWidth: the width of the band of beam offsets whose intercepted fraction >= threshold"""
    if beam_diameter <= 0.0:
        return diameter if threshold <= 1.0 else 0.0
    if strip_intercept_fraction(diameter, beam_diameter, 0.0) < threshold:
        return 0.0
    low, high = 0.0, 0.5 * diameter + 4.0 * beam_diameter
    for _ in range(60):
        mid = 0.5 * (low + high)
        if strip_intercept_fraction(diameter, beam_diameter, mid) >= threshold:
            low = mid
        else:
            high = mid
    return 2.0 * low


def matching_threshold(diameter, beam_diameter, target_width):
    """the threshold for which a beam of this diameter gives the effective width @c target_width (None if none)"""
    lo, hi = 1e-6, strip_intercept_fraction(diameter, beam_diameter, 0.0)
    if effective_width(diameter, beam_diameter, lo) < target_width or hi <= lo:
        return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if effective_width(diameter, beam_diameter, mid) >= target_width:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fp_text(exit_diameter, divergence, threshold):
    return f"{round(exit_diameter, 6):g},{round(divergence, 7):g},{round(threshold, 6):g}"


def render(diameter=0.01, voxel=0.05, pool=0.05, footprint=FOOTPRINT, extra=""):
    fp = f"--footprint {footprint} " if footprint else ""
    return (f"--estimator pooled --stem_diameter {diameter} {fp}--pixel_width {voxel} --pool_radius {pool} {extra}"
            ).strip()


def pair(prefix="", **kw):
    """the default estimator and the self-shadow correction"""
    base = render(**kw)
    p = prefix + " " if prefix else ""
    return [(p + "base", base), (p + "ss", base + " --self_shadow_correction")]


def scene(stems=300, extra="", footprint=FOOTPRINT):
    fp = f"--footprint {footprint} " if footprint else ""
    return f"--stems_per_m2 {stems} {fp}{extra}".strip()


class Config:
    def __init__(self, block, label, scene_args, variants, layers=LOW, plot_width=4.0, eval_half=None,
                 line_bins=False, pass_eligible=None, note="", meta=None):
        self.block = block
        self.label = label
        self.scene_args = scene_args
        self.variants = variants
        self.layers = layers
        self.plot_width = plot_width
        self.eval_half = eval_half  # evaluate only the cells within this half-width of the plot centre
        self.line_bins = line_bins  # also summarise by distance of the cell centre from the nearest walk line
        self.pass_eligible = (block in PASS_BLOCKS) if pass_eligible is None else pass_eligible
        self.note = note
        self.meta = meta or {}  # per variant name: what the variant assumes (block L: footprint and diameter)

    def ident(self):
        text = json.dumps([self.scene_args, self.variants, self.layers, self.plot_width, self.eval_half])
        slug = re.sub(r"[^A-Za-z0-9]+", "_", self.label).strip("_")[:60]
        return f"{slug}_{hashlib.sha1(text.encode()).hexdigest()[:8]}"


def build_configs():
    c = []
    # ---- A. scan geometry
    for n in (300, 1000):
        for spacing in (0.5, 1, 2, 3):
            # a 3 m line spacing needs a 6 m plot for two lines; the central 4 x 4 m is evaluated
            plot = 6.0 if spacing == 3 else 4.0
            c.append(Config("A", f"walk {spacing} m, {n}/m2", scene(n, f"--walk_spacing {spacing}"), pair(),
                            plot_width=plot, eval_half=2.0 if plot > 4 else None, line_bins=True))
    for n in (300, 1000):
        for h in (0.8, 1.2, 2.0):
            c.append(Config("A", f"scanner {h} m, {n}/m2", scene(n, f"--scanner_height {h}"), pair()))
    for z0, z1 in ((0, 101.5), (30, 90), (60, 95)):
        c.append(Config("A", f"zenith {z0}-{z1}, 300/m2", scene(300, f"--zenith {z0},{z1}"), pair()))
    c.append(Config("A", "plot 8 m, 300/m2", scene(300), pair(), plot_width=8.0, line_bins=True))
    # range: a single walk line along the middle of an 8 m plot, cells binned by distance from the line
    c.append(Config("A", "single line, 8 m plot, 300/m2", scene(300, "--walk_spacing 8"), pair(), plot_width=8.0,
                    line_bins=True, pass_eligible=True))
    c.append(Config("A", "single line, 8 m plot, 100/m2", scene(100, "--walk_spacing 8"), pair(), plot_width=8.0,
                    line_bins=True, pass_eligible=True))
    # ---- B. range noise, both return conventions (generator and estimator matched)
    for n in (300, 1000):
        for conv in ("surface", "outline"):
            for sd in (0.0, 0.005, 0.01, 0.02, 0.03):
                g = "" if conv == "surface" else "--return_range outline"
                e = "" if conv == "surface" else "--return_range outline"
                noise = f"--range_noise {sd}" if sd > 0 else ""
                c.append(Config("B", f"range noise {sd}, {conv}, {n}/m2", scene(n, f"{g} {noise}".strip()),
                                pair(extra=e)))
    # ---- C. pose noise
    c.append(Config("C", "pose noise 0, 300/m2", scene(300), pair()))
    for sd in (0.01, 0.02, 0.05):
        for knot in (0.5, 2):
            c.append(Config("C", f"pose noise {sd} per {knot} m, 300/m2", scene(300, f"--pose_noise {sd},{knot}"),
                            pair()))
    # ---- D. diameter error (render only: one scene per seed, several assumed diameters)
    for d in (0.005, 0.01, 0.02):
        for n in (300, 1000):
            v = []
            for ratio in (0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3):
                v += pair(f"x{ratio}", diameter=round(d * ratio, 6))
            c.append(Config("D", f"true d {d}, {n}/m2", scene(n, f"--diameter {d}"), v,
                            note="pass only at ratio 1.0"))
    c.append(Config("D", "mixture 6/14 mm (50%), assumed 10 mm, 300/m2",
                    scene(300, "--diameter_mixture 0.006,0.014,0.5"), pair()))
    c.append(Config("D", "log-normal 10 +- 3 mm, assumed 10 mm, 300/m2", scene(300, "--diameter_sd 0.003"), pair()))
    # ---- E. stem heights
    for n in (300, 1000):
        for label, extra, top, ok in (("normal 0.25+-0.08", "", 0.5, True),
                                      ("normal 0.4+-0.15", "--height 0.4,0.15", 0.8, True),
                                      ("uniform 0.05-0.5", "--uniform_height 0.05,0.5", 0.5, False),
                                      ("fixed 0.6", "--height 0.6,0", 0.6, True),
                                      ("fixed 2 (no tops)", "--height 2,0", 1.0, True)):
            c.append(Config("E", f"heights {label}, {n}/m2", scene(n, extra), pair(), layers=layers_to(top),
                            pass_eligible=ok))
    # ---- F. low observation (reported, not pass/fail)
    for n in (600, 1000, 1500):
        for hl, he, top in (("fixed 0.6", "--height 0.6,0", 0.6), ("normal 0.25", "", 0.4)):
            for spacing in (1, 2):
                c.append(Config("F", f"{n}/m2, {hl}, walk {spacing} m", scene(n, f"{he} --walk_spacing {spacing}"),
                                pair(), layers=layers_to(top)))
    # ---- G. clustering (reported, not pass/fail)
    for n in (150, 500):
        for spread in (0.015, 0.03, 0.06):
            c.append(Config("G", f"cable roots {n}/m2, spread {spread}",
                            scene(n, f"--clustered --cluster_spread {spread}"),
                            pair("pool0.05") + pair("pool0.1", pool=0.1)))
    for radius in (0.2, 0.3, 0.5):
        for cover in (0.3, 0.6):
            for matrix in (0, 0.2):
                c.append(Config("G", f"patches r{radius} cover {cover} matrix {matrix}, 300/m2",
                                scene(300, f"--patchy {radius},{cover},{matrix}"), pair()))
    # ---- H. return convention: generator x estimator
    for n in (300, 1000):
        for hl, he, layers in (("normal", "", LOW), ("fixed 0.6", "--height 0.6,0", [0.0, 0.1, 0.2, 0.3])):
            for gen in ("surface", "outline"):
                g = "" if gen == "surface" else "--return_range outline"
                v = pair("est surface") + pair("est outline", extra="--return_range outline")
                c.append(Config("H", f"gen {gen}, heights {hl}, {n}/m2", scene(n, f"{g} {he}".strip()), v,
                                layers=layers, note="pass only for matched conventions"))
    # ---- I. spacing (estimator given the true gap; Poisson stems given gap 0)
    for n in (300, 1000):
        c.append(Config("I", f"Poisson (overlapping), {n}/m2", scene(n, "--overlapping"), pair()))
        for gap in (0, 0.5, 1):
            g = f"--min_gap {gap}" if gap > 0 else ""
            c.append(Config("I", f"min gap {gap}, {n}/m2", scene(n, g), pair(extra=g)))
    # ---- J. estimator window, over the A-C scenes at 300/m2
    windows = ((0.025, 0.05), (0.05, 0.05), (0.05, 0.1), (0.05, 0.15), (0.1, 0.1))
    scenes = set()
    for cfg in list(c):
        if cfg.block not in "ABC" or not re.search(r"(^|[^0-9])300/m2", cfg.label):
            continue
        # the default scene appears in several blocks; run it once
        canonical = (re.sub(r" --(zenith 0,101.5|scanner_height 1.2|walk_spacing 1)\b", "", cfg.scene_args),
                     cfg.plot_width)
        if canonical in scenes:
            continue
        scenes.add(canonical)
        e = "--return_range outline" if "outline" in cfg.label else ""
        v = []
        for voxel, pool in windows:
            v += pair(f"v{voxel} p{pool}", voxel=voxel, pool=pool, extra=e)
        c.append(Config("J", f"[{cfg.block}] {cfg.label}", cfg.scene_args, v, layers=cfg.layers,
                        plot_width=cfg.plot_width, eval_half=cfg.eval_half))
    c += footprint_mismatch_configs()
    # ---- X. diagnostics of the failures (not part of "all"): the errors that appear in 0.1-0.2 m, where stems of
    # the default heights end, repeated with stems of one height (0.6 m, no tops below 0.6 m)
    fixed = "--height 0.6,0"
    top5 = layers_to(0.5)
    for n in (300, 1000):
        c.append(Config("X", f"fixed 0.6, range noise 0.03, {n}/m2", scene(n, f"{fixed} --range_noise 0.03"), pair(),
                        layers=top5, pass_eligible=False))
        c.append(Config("X", f"fixed 0.6, outline, range noise 0.03, {n}/m2",
                        scene(n, f"{fixed} --return_range outline --range_noise 0.03"),
                        pair(extra="--return_range outline"), layers=top5, pass_eligible=False))
    for sd in (0.02, 0.05):
        c.append(Config("X", f"fixed 0.6, pose noise {sd} per 0.5 m, 300/m2",
                        scene(300, f"{fixed} --pose_noise {sd},0.5"), pair(), layers=top5, pass_eligible=False))
    # wider windows for wide stems (d_eff / W smaller, so less second-order self-shadowing)
    for d in (0.02,):
        v = pair("p0.1", diameter=d, pool=0.1) + pair("p0.15", diameter=d, pool=0.15)
        c.append(Config("X", f"true d {d}, 300/m2, wider windows", scene(300, f"--diameter {d}"), v,
                        pass_eligible=False))
        c.append(Config("X", f"true d {d}, fixed 0.6, 300/m2", scene(300, f"--diameter {d} {fixed}"),
                        pair(diameter=d), layers=top5, pass_eligible=False))
    c.append(Config("X", "mixture 6/14 mm, fixed 0.6, 300/m2",
                    scene(300, f"--diameter_mixture 0.006,0.014,0.5 {fixed}"), pair(), layers=top5,
                    pass_eligible=False))
    # the layer holding the tops of stems of one height: tops on a layer boundary (0.6 m) evaluated with the grid
    # extending above them, and tops within a layer (0.55 m)
    for h in (0.6, 0.55):
        c.append(Config("X", f"fixed {h}, grid to 0.7 m, 300/m2", scene(300, f"--height {h},0"), pair(),
                        layers=layers_to(0.7), pass_eligible=False))
    # no beam footprint: the stem tops' contribution without the footprint's widening of the top discs
    c.append(Config("X", "no footprint, heights normal 0.25, 300/m2", scene(300, footprint=""), pair(footprint=""),
                    layers=layers_to(0.5), pass_eligible=False))
    return c

# ---- L. footprint mismatch: scenes generated with a true beam footprint, estimated assuming another one
L_TRUE = (("default", 0.01, 0.003, 0.05), ("narrow", 0.005, 0.0005, 0.05), ("wide", 0.02, 0.005, 0.05))
L_RATIOS = (0.5, 0.75, 1.25, 1.5, 2.0)
L_THRESHOLDS = (0.02, 0.03, 0.08, 0.10, 0.15)
L_GRID_DIAMETERS = (0.8, 0.9, 1.1, 1.2)
L_GRID_THRESHOLDS = (0.02, 0.03, 0.05, 0.08, 0.10, 0.15)
L_MATCH_RANGE = 2.0  # the joint cases match d_eff (or the beam diameter) at this range, not at 1 or 4 m


def footprint_estimators(fp, diameter=0.01, extra=""):
    """the three estimators (default window, self-shadow correction, small window) with an assumed footprint"""
    base = render(diameter=diameter, footprint=fp, extra=extra)
    sw = render(diameter=diameter, footprint=fp, pool=0, extra=(extra + " " + SW_EXTRA).strip())
    return [("base", base), ("ss", base + " --self_shadow_correction"), ("sw", sw)]


def assumed_footprints(w0, g, t, d=0.01):
    """(tag, kind, level, exit diameter, divergence, threshold) of the assumed footprints for a true one"""
    out = [("true", "true", 1.0, w0, g, t)]
    out += [(f"e{r:g}", "exit", r, w0 * r, g, t) for r in L_RATIOS]
    out += [(f"g{r:g}", "divergence", r, w0, g * r, t) for r in L_RATIOS]
    out += [(f"t{x:g}", "threshold", x / t, w0, g, x) for x in L_THRESHOLDS]
    beam = w0 + g * L_MATCH_RANGE
    # the same beam diameter at 2 m: twice the divergence (narrower exit), or no divergence (wider exit)
    out.append(("j_steep", "joint", 1.0, w0 - g * L_MATCH_RANGE if w0 > g * L_MATCH_RANGE else 0.5 * w0,
                2.0 * g if w0 > g * L_MATCH_RANGE else (beam - 0.5 * w0) / L_MATCH_RANGE, t))
    out.append(("j_flat", "joint", 1.0, beam, 0.0, t))
    # the same d_eff at 2 m from another exit diameter and the threshold that compensates for it
    target = effective_width(d, beam, t)
    for r in (0.5, 1.5, 2.0):
        tm = matching_threshold(d, w0 * r + g * L_MATCH_RANGE, target)
        if tm is not None:
            out.append((f"jt_e{r:g}", "joint", r, w0 * r, g, tm))
    return out


def footprint_mismatch_configs():
    c = []
    for name, w0, g, t in L_TRUE:
        true_fp = fp_text(w0, g, t)
        for n in (300, 1000):
            variants, meta = [], {}
            for tag, kind, level, aw0, ag, at in assumed_footprints(w0, g, t):
                for est, args in footprint_estimators(fp_text(aw0, ag, at)):
                    variants.append((f"{tag} {est}", args))
                    meta[f"{tag} {est}"] = dict(tag=tag, kind=kind, level=level, exit=aw0, divergence=ag,
                                                threshold=at, diameter=0.01, est=est)
            for est, args in footprint_estimators(""):
                variants.append((f"nofp {est}", args))
                meta[f"nofp {est}"] = dict(tag="nofp", kind="none", level=0.0, exit=0.0, divergence=0.0,
                                           threshold=0.0, diameter=0.01, est=est)
            if name == "default":
                # assumed diameter x assumed threshold
                for dr in L_GRID_DIAMETERS:
                    for x in L_GRID_THRESHOLDS:
                        tag = f"d{dr:g} t{x:g}"
                        for est, args in footprint_estimators(fp_text(w0, g, x), diameter=round(0.01 * dr, 6)):
                            variants.append((f"{tag} {est}", args))
                            meta[f"{tag} {est}"] = dict(tag=tag, kind="grid", level=dr, exit=w0, divergence=g,
                                                        threshold=x, diameter=0.01 * dr, est=est)
            c.append(Config("L", f"true {name} {true_fp}, {n}/m2", scene(n, footprint=true_fp), variants,
                            pass_eligible=False, meta=meta,
                            note=f"true footprint {true_fp}; variants '<assumed> <estimator>'"))
    # range: one walk line along the middle of an 8 m plot, cells binned by distance from the line
    w0, g, t = L_TRUE[0][1:]
    variants, meta = [], {}
    for tag, kind, level, aw0, ag, at in assumed_footprints(w0, g, t):
        for est, args in footprint_estimators(fp_text(aw0, ag, at)):
            variants.append((f"{tag} {est}", args))
            meta[f"{tag} {est}"] = dict(tag=tag, kind=kind, level=level, exit=aw0, divergence=ag, threshold=at,
                                        diameter=0.01, est=est)
    c.append(Config("L", f"single line, 8 m plot, true default {fp_text(w0, g, t)}, 300/m2",
                    scene(300, "--walk_spacing 8", footprint=fp_text(w0, g, t)), variants, plot_width=8.0,
                    line_bins=True, pass_eligible=False, meta=meta))
    return c


def run(cmd, log):
    with open(log, "a") as f:
        f.write("$ " + cmd + "\n")
        f.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT)


def read_csv(file_name):
    with open(file_name) as f:
        return list(csv.DictReader(f))


def key(row):
    return tuple(round(float(row[k]), 4) for k in ("x_min", "y_min", "z_min"))


def walk_lines(scene_args, plot_width):
    m = re.search(r"--walk_spacing ([0-9.]+)", scene_args)
    spacing = float(m.group(1)) if m else 1.0
    num_lines = max(1, int(plot_width / spacing))
    return spacing, num_lines


def run_job(job):
    cfg, seed, args, out_file = job
    work = tempfile.mkdtemp(prefix="job_", dir=args.tmp)
    result = dict(block=cfg.block, label=cfg.label, ident=cfg.ident(), seed=seed, scene=cfg.scene_args,
                  plot_width=cfg.plot_width, variants={}, cells=[])
    try:
        log = os.path.join(work, "log.txt")
        half = 0.5 * cfg.plot_width
        spacing, num_lines = walk_lines(cfg.scene_args, cfg.plot_width)
        lines_y = [-half + spacing * (i + 0.5) for i in range(num_lines)]
        rays = int(args.rays_per_m * cfg.plot_width * num_lines)
        top = round(max(cfg.layers) + 0.1, 4)
        t0 = time.time()
        run(f"cd {work} && raycreate pneumatophores {seed} {cfg.scene_args} --plot_width {cfg.plot_width} "
            f"--rays {rays} --truth_cells 1,1,0.1", log)
        result["t_gen"] = time.time() - t0
        truth = {key(r): float(r["true_frontal_area_per_volume_Nd"])
                 for r in read_csv(os.path.join(work, "pneumatophores_truth_cells.csv"))}
        for name, render_args in cfg.variants:
            rlog = os.path.join(work, "render.txt")
            if os.path.exists(rlog):
                os.remove(rlog)
            t0 = time.time()
            run(f"cd {work} && rayrender pneumatophores.ply top density --density_model vertical_cylinder "
                f"--no_flat_top {render_args} --min_bound {-half},{-half},0 --max_bound {half},{half},{top} "
                f"--cells cells.csv --cell_size 1,1,0.1", rlog)
            t_est = time.time() - t0
            text = open(rlog).read()
            m = re.search(r"mean effective width / diameter of ([0-9.e+-]+)", text)
            rho = float(m.group(1)) if m else 1.0
            m = re.search(r"self-shadow correction(?: \(second order\))?: mean factor ([0-9.e+-]+)", text)
            ss = float(m.group(1)) if m else float("nan")
            voxel = float(re.search(r"--pixel_width ([0-9.]+)", render_args).group(1))
            pool = float(re.search(r"--pool_radius ([0-9.]+)", render_args).group(1))
            k = max(0, int(math.floor(pool / voxel + 0.5)))
            window = (2 * k + 1) * voxel
            d = float(re.search(r"--stem_diameter ([0-9.]+)", render_args).group(1))
            result["variants"][name] = dict(t_est=t_est, rho=rho, ss=ss, window=window,
                                            h1=1.0 + EXIT_DISTANCE_SQUARE * rho * d / window)
            for row in read_csv(os.path.join(work, "cells.csv")):
                kk = key(row)
                if kk not in truth or round(kk[2], 1) not in cfg.layers:
                    continue
                cx, cy = kk[0] + 0.5, kk[1] + 0.5
                if cfg.eval_half is not None and (abs(cx) > cfg.eval_half or abs(cy) > cfg.eval_half):
                    continue
                obs = float(row["observed_fraction"])
                result["cells"].append([name, round(kk[2], 1), cx, cy, truth[kk],
                                        float(row["density_frontal_area_per_volume_Nd"]),
                                        float(row["pooled_density_frontal_area_per_volume_Nd"]),
                                        0.0 if math.isnan(obs) else obs,
                                        min(abs(cy - y) for y in lines_y),
                                        float(row["weighted_path"])])
    except subprocess.CalledProcessError as e:
        result["error"] = str(e) + "\n" + open(os.path.join(work, "log.txt")).read()[-2000:]
    finally:
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    # a failed job is saved beside, not as, the job file, so that it runs again on resuming
    target = out_file + ".err" if result.get("error") else out_file
    with open(target + ".tmp", "w") as f:
        json.dump(result, f, allow_nan=True)
    os.replace(target + ".tmp", target)
    return out_file, result.get("error")


def mean_se(values):
    values = [v for v in values if not math.isnan(v)]
    n = len(values)
    if n == 0:
        return float("nan"), float("nan")
    m = sum(values) / n
    if n < 2:
        return m, float("nan")
    return m, math.sqrt(sum((v - m) ** 2 for v in values) / (n - 1) / n)


def ratio(a, b):
    return a / b if b > 0 else float("nan")


def nan0(v):
    return 0.0 if math.isnan(v) else v


def seed_metrics(cells):
    """metrics of one seed's cells (one variant and layer)"""
    t = sum(c[4] for c in cells)
    d = sum(nan0(c[5]) for c in cells)
    p = sum(nan0(c[6]) for c in cells)
    good = [c for c in cells if c[7] >= 0.9]
    t9 = sum(c[4] for c in good)
    d9 = sum(c[5] for c in good)
    p9 = sum(c[6] for c in good)
    n = len(cells)
    mt = t / n if n else float("nan")
    rmse = math.sqrt(sum((nan0(c[5]) - c[4]) ** 2 for c in cells) / n) / mt if n and mt > 0 else float("nan")
    in10 = sum(1 for c in cells if c[4] > 0 and abs(nan0(c[5]) / c[4] - 1.0) <= 0.1) / n if n else float("nan")
    rmse9 = float("nan")
    if good and t9 > 0:
        rmse9 = math.sqrt(sum((c[5] - c[4]) ** 2 for c in good) / len(good)) / (t9 / len(good))
    return dict(dt=ratio(d, t), pt=ratio(p, t), dp=ratio(d, p), dt90=ratio(d9, t9), pt90=ratio(p9, t9),
                rmse=rmse, rmse90=rmse9, in10=in10, obs=sum(c[7] for c in cells) / n if n else float("nan"),
                obs_min=min(c[7] for c in cells) if n else float("nan"), frac90=len(good) / n if n else 0.0,
                truth=mt)


def load_results(out, configs):
    by_ident = {cfg.ident(): cfg for cfg in configs}
    results = {}
    root = os.path.join(out, "jobs")
    if not os.path.isdir(root):
        return results
    for block in sorted(os.listdir(root)):
        for name in sorted(os.listdir(os.path.join(root, block))):
            if not name.endswith(".json"):
                continue
            with open(os.path.join(root, block, name)) as f:
                r = json.load(f)
            if r.get("error") or r["ident"] not in by_ident:
                continue
            results.setdefault(r["ident"], []).append(r)
    return results


def summarise(out, configs, seeds_filter=None):
    results = load_results(out, configs)
    rows, bins, line_rows = [], [], []
    for cfg in configs:
        runs = results.get(cfg.ident(), [])
        if seeds_filter:
            runs = [r for r in runs if r["seed"] in seeds_filter]
        if not runs:
            continue
        for name, _ in cfg.variants:
            for z in cfg.layers:
                per_seed = []
                pooled_cells = []
                for r in runs:
                    cells = [c for c in r["cells"] if c[0] == name and c[1] == z]
                    if cells:
                        per_seed.append(seed_metrics(cells))
                        pooled_cells += cells
                if not per_seed:
                    continue
                v = [r["variants"][name] for r in runs]
                row = dict(block=cfg.block, config=cfg.label, variant=name, z=z, seeds=len(per_seed),
                           pass_eligible=int(cfg.pass_eligible))
                for m in ("dt", "dt90", "pt", "pt90", "dp", "rmse", "rmse90", "in10", "obs", "frac90"):
                    row[m], row[m + "_se"] = mean_se([s[m] for s in per_seed])
                row["obs_min"] = min(s["obs_min"] for s in per_seed)
                row["truth"] = mean_se([s["truth"] for s in per_seed])[0]
                row["ss"] = mean_se([x["ss"] for x in v])[0]
                row["h1"] = mean_se([x["h1"] for x in v])[0]
                row["rho"] = mean_se([x["rho"] for x in v])[0]
                row["window"] = v[0]["window"]
                row["t_gen"] = mean_se([r["t_gen"] for r in runs])[0]
                row["t_est"] = mean_se([x["t_est"] for x in v])[0]
                rows.append(row)
                for lo, hi in zip(OBS_BINS[:-1], OBS_BINS[1:]):
                    inside = [c for c in pooled_cells if lo <= c[7] < hi]
                    t = sum(c[4] for c in inside)
                    if inside and t > 0:
                        bins.append(dict(block=cfg.block, config=cfg.label, variant=name, z=z, obs_lo=lo,
                                         obs_hi=min(hi, 1.0), cells=len(inside),
                                         dt=sum(nan0(c[5]) for c in inside) / t,
                                         pt=sum(nan0(c[6]) for c in inside) / t))
                if cfg.line_bins:
                    for dist in sorted(set(round(c[8], 2) for c in pooled_cells)):
                        per = []
                        for r in runs:
                            cells = [c for c in r["cells"] if c[0] == name and c[1] == z and round(c[8], 2) == dist]
                            if cells:
                                per.append(seed_metrics(cells))
                        dt = mean_se([s["dt"] for s in per])
                        dt90 = mean_se([s["dt90"] for s in per])
                        pt = mean_se([s["pt"] for s in per])
                        line_rows.append(dict(block=cfg.block, config=cfg.label, variant=name, z=z,
                                              line_distance=dist, seeds=len(per), dt=dt[0], dt_se=dt[1],
                                              dt90=dt90[0], dt90_se=dt90[1], pt=pt[0], pt_se=pt[1],
                                              obs=mean_se([s["obs"] for s in per])[0]))
    for name, data in (("summary.csv", rows), ("bins.csv", bins), ("line_distance.csv", line_rows)):
        if data:
            with open(os.path.join(out, name), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(data[0].keys()))
                w.writeheader()
                w.writerows(data)
    lines = ["block | config | variant | z | n | d/t | d/t90 | p/t | p/t90 | d/p | rmse | in10 | obs (min) | "
             "f90 | ss (H1) | t_gen t_est"]
    for r in rows:
        lines.append(f"{r['block']} | {r['config']} | {r['variant']} | {r['z']:.1f} | {r['seeds']} | "
                     f"{r['dt']:.3f} ± {r['dt_se']:.3f} | {r['dt90']:.3f} ± {r['dt90_se']:.3f} | "
                     f"{r['pt']:.3f} ± {r['pt_se']:.3f} | {r['pt90']:.3f} ± {r['pt90_se']:.3f} | "
                     f"{r['dp']:.3f} | {r['rmse']:.3f} | {r['in10']:.2f} | {r['obs']:.2f} ({r['obs_min']:.2f}) | "
                     f"{r['frac90']:.2f} | {r['ss']:.3f} ({r['h1']:.3f}) | {r['t_gen']:.1f} {r['t_est']:.1f}")
    with open(os.path.join(out, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return rows


def failures(rows, configs, variant="ss"):
    """pass-eligible configurations whose d/t90 in 0-0.2 m misses 1 +- PASS_TOLERANCE, for the variants ending with
    @c variant (all variants if empty). Layers without cells observed >= 0.9 are not tested"""
    eligible = {cfg.label: cfg for cfg in configs if cfg.pass_eligible}
    failed = set()
    for r in rows:
        if r["config"] not in eligible or r["z"] > 0.15 or not pass_variant(r) or not r["variant"].endswith(variant):
            continue
        if not math.isnan(r["dt90"]) and abs(r["dt90"] - 1.0) > PASS_TOLERANCE:
            failed.add(r["config"])
    return failed


def pass_variant(row):
    """whether a variant counts for pass/fail: blocks D (assumed = true diameter) and H (matched conventions) only
    count their matched variants"""
    v = row["variant"]
    if row["block"] == "D" and v.startswith("x") and not v.startswith("x1.0 "):
        return False
    if row["block"] == "H":
        gen_outline = "gen outline" in row["config"]
        return ("est outline" in v) == gen_outline
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    configs = build_configs()
    parser.add_argument("--blocks", default="all",
                        help="comma separated block letters (A-J, L for the footprint mismatch and X for the "
                        "diagnostics), or all (A-J)")
    parser.add_argument("--configs", default="", help="only configurations whose label matches this regex")
    parser.add_argument("--seeds", default=SEEDS)
    parser.add_argument("--rays_per_m", type=float, default=56250.0)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--out", default="validation_matrix")
    parser.add_argument("--bin", default="", help="directories added to the PATH, comma separated")
    parser.add_argument("--keep", action="store_true", help="keep the scenes")
    parser.add_argument("--summary_only", action="store_true")
    parser.add_argument("--list", action="store_true", help="list the configurations and exit")
    parser.add_argument("--rerun_failures", default="",
                        help="seeds to add for pass-eligible configurations that fail with the seeds present")
    parser.add_argument("--rerun_variant", default="ss",
                        help="with --rerun_failures: the variants (by name ending) whose failure counts")
    args = parser.parse_args()
    if args.bin:
        os.environ["PATH"] = os.pathsep.join(os.path.abspath(b) for b in args.bin.split(",")) + os.pathsep + \
            os.environ["PATH"]
    blocks = set("ABCDEFGHIJ") if args.blocks == "all" else set(args.blocks.split(","))  # L, X only on request
    selected = [c for c in configs if c.block in blocks and (not args.configs or re.search(args.configs, c.label))]
    if args.list:
        for c in selected:
            print(f"{c.block} | {c.label} | {c.scene_args} | {len(c.variants)} variants | layers {c.layers}")
        print(f"{len(selected)} configurations")
        return 0
    os.makedirs(args.out, exist_ok=True)
    args.tmp = os.path.abspath(args.out)
    seeds = [int(s) for s in args.seeds.split(",")]
    if args.rerun_failures:
        rows = summarise(args.out, configs)
        failed = failures(rows, selected, args.rerun_variant)
        selected = [c for c in selected if c.label in failed]
        seeds = [int(s) for s in args.rerun_failures.split(",")]
        print(f"re-running {len(selected)} failing configurations with seeds {seeds}")
    if not args.summary_only:
        jobs = []
        for cfg in selected:
            for seed in seeds:
                out_file = os.path.join(args.out, "jobs", cfg.block, f"{cfg.ident()}_{seed}.json")
                if not os.path.exists(out_file):
                    jobs.append((cfg, seed, args, out_file))
        print(f"{len(jobs)} jobs to run", flush=True)
        start = time.time()
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_job, j) for j in jobs]
            for i, fut in enumerate(as_completed(futures)):
                out_file, error = fut.result()
                status = "FAILED " + error.splitlines()[0] if error else "done"
                print(f"[{i + 1}/{len(jobs)} {time.time() - start:.0f} s] {status} {out_file}", flush=True)
    summarise(args.out, configs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
