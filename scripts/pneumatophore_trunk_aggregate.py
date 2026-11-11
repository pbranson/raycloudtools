#!/usr/bin/env python3
# Copyright (c) 2026
# Commonwealth Scientific and Industrial Research Organisation (CSIRO)
# ABN 41 687 119 230
"""Aggregates scripts/pneumatophore_trunk_study.py over seeds, for the table of Section 7.8 of
docs/pneumatophore_density_report.md. Expects <dir>/seed_<s>/summary.txt (the study's output for the density
column) and, if present, <dir>/seed_<s>/summary_pooled.txt (with --pooled), e.g. from the repository directory:
  for s in 31 32 33 34; do mkdir -p trunks/seed_$s && (cd trunks/seed_$s &&
    python3 ../../scripts/pneumatophore_trunk_study.py --seed $s > summary.txt &&
    python3 ../../scripts/pneumatophore_trunk_study.py --seed $s --pooled > summary_pooled.txt); done
  python3 scripts/pneumatophore_trunk_aggregate.py trunks
Prints, per case (A, B, C) and cell group, estimate / pneumatophore truth as mean +- s.d. over the seeds per 0.1 m
layer, with estimate / (pneumatophore + trunk truth) in brackets.
"""
import argparse
import os
import re

import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("dir", help="directory holding seed_<s>/summary.txt")
    p.add_argument("--seeds", default="31,32,33,34", help="comma separated seeds")
    a = p.parse_args()
    seeds = [int(s) for s in a.seeds.split(",")]
    for name, label in (("summary.txt", "density (mean, main)"), ("summary_pooled.txt", "pooled_density")):
        if not all(os.path.exists(os.path.join(a.dir, f"seed_{s}", name)) for s in seeds):
            continue
        acc = {}
        for s in seeds:
            case = None
            with open(os.path.join(a.dir, f"seed_{s}", name)) as f:
                lines = f.read().splitlines()
            for line in lines:
                if line[:2] in ("A:", "B:", "C:"):
                    case = line[0]
                m = re.match(r"\s+(all cells|cells with a trunk|cells without)\s+estimate.*?: (.*)", line)
                if m and case:
                    acc.setdefault((case, m.group(1)), []).append(
                        [tuple(map(float, re.findall(r"[\d.]+", part))) for part in m.group(2).split("|")])
        print(f"{label}: estimate / pneumatophore truth, mean ± s.d. over seeds {a.seeds} (in brackets: / "
              "pneumatophore + trunk truth), layers 0-0.1 | 0.1-0.2 | 0.2-0.3 m")
        for (case, group), v in acc.items():
            arr = np.array(v)
            sd = (lambda x: x.std(ddof=1)) if len(v) > 1 else (lambda x: 0.0)
            print(" ", case, f"{group:20s}", " | ".join(
                f"{arr[:, z, 0].mean():.3f}±{sd(arr[:, z, 0]):.3f} ({arr[:, z, 1].mean():.3f})"
                for z in range(arr.shape[1])))


if __name__ == "__main__":
    main()
