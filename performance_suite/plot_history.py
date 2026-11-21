#!/usr/bin/env python3
"""History of the run time of mcfost on one machine (history/<machine>.csv, record_history.py).

For each job, the median time of each binary relative to the reference binary of the same job
(ratio < 1: faster than the reference), against the date of the job; one line per workload,
binaries other than the reference only. Prints the table and writes history/<machine>.png and .pdf.

Usage: plot_history.py --machine <name> [--phase total] [--label <binary label>]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from datetime import datetime

from plot_scaling import COLORS


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--machine", required=True)
    parser.add_argument("--phase", default="total")
    parser.add_argument("--label", help="only this binary label (default: all but the reference)")
    args = parser.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    hist = os.path.join(here, "history", f"{args.machine}.csv")
    if not os.path.exists(hist):
        sys.exit(f"ERROR: {hist} not found (record_history.py)")
    series = defaultdict(list)  # workload -> [(date, ratio, label, sha)]
    with open(hist) as f:
        for r in csv.DictReader(f):
            if r["phase"] != args.phase or not r["ratio"] or r["label"] == r["ref_label"]:
                continue
            if args.label and r["label"] != args.label:
                continue
            date = datetime.fromisoformat(r["date"].replace("Z", "+00:00"))
            series[r["workload"]].append((date, float(r["ratio"]), r["label"], r["sha"][:7]))
    if not series:
        sys.exit(f"ERROR: no ratio to a reference binary for phase '{args.phase}' in {hist}")

    print(f"# {args.machine}, phase {args.phase}: time / time of the reference binary in the same job")
    print(f"{'date':>20} {'workload':>8} {'label':>20} {'sha':>8} {'ratio':>7}")
    for workload, pts in sorted(series.items()):
        for date, ratio, label, sha in sorted(pts):
            print(f"{date:%Y-%m-%d %H:%M} {workload:>8} {label:>20} {sha:>8} {ratio:7.3f}")

    import matplotlib
    matplotlib.use("agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    ax.axhline(1.0, color="0.4", lw=0.8)
    for workload, pts in sorted(series.items()):
        pts.sort()
        x, y = [p[0] for p in pts], [p[1] for p in pts]
        ax.plot(x, y, color=COLORS.get(workload, "k"), lw=2, marker="o", ms=4, label=workload)
        ax.annotate(workload, (x[-1], y[-1]), xytext=(5, 0), textcoords="offset points", va="center",
                    fontsize=8, color="0.15")
    ax.set_ylabel(f"time / reference ({args.phase})")
    ax.set_title(args.machine, fontsize=9)
    ax.grid(True, axis="y", color="0.9", lw=0.6)
    ax.legend(frameon=False, fontsize=8, ncol=len(series), loc="upper left")
    fig.autofmt_xdate()
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(here, "history", f"{args.machine}.{ext}"), dpi=150)
    print(f"plots: {os.path.relpath(os.path.join(here, 'history', args.machine + '.png'))}, .pdf")


if __name__ == "__main__":
    main()
