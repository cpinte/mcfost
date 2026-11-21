#!/usr/bin/env python3
"""Thread scaling of mcfost from the results of scaling.sbatch (one results file per thread count).

Prints, per workload, the median total time, the speed-up relative to the smallest thread count
(scaled by that count, so it reads as an equivalent number of cores) and the parallel efficiency,
and plots the speed-up against the number of threads (scaling.png and scaling.pdf in the directory).

Usage: plot_scaling.py <results/scaling_... directory> [--label <binary label>] [--phase total]
"""

import argparse
import glob
import os
import statistics
import sys

from record_history import read

# categorical slots in fixed order (validated palette, light surface): one per workload
COLORS = {"W1": "#2a78d6", "W2": "#eb6834", "W3": "#1baf7a", "W4": "#eda100", "W5": "#e87ba4"}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory")
    parser.add_argument("--label", help="binary label (default: the first one)")
    parser.add_argument("--phase", default="total", help="phase to plot (default: total)")
    args = parser.parse_args()

    data = {}  # workload -> {threads: median time}
    cpu, label = "", args.label
    for name in sorted(glob.glob(os.path.join(args.directory, "t*.tsv"))):
        meta, shas, times = read(name)
        cpu = meta.get("cpu", cpu)
        n = int(meta["omp_num_threads"])
        label = label or next(iter(shas))
        for (lab, workload, phase), ts in times.items():
            if lab == label and phase == args.phase:
                data.setdefault(workload, {})[n] = statistics.median(ts)
    if not data:
        sys.exit(f"ERROR: no '{args.phase}' times for {label} in {args.directory}")

    print(f"# {cpu}, binary {label}, phase {args.phase}")
    print(f"{'workload':>8} {'threads':>7} {'time[s]':>9} {'speed-up':>9} {'efficiency':>10}")
    for workload, d in sorted(data.items()):
        n0 = min(d)
        for n in sorted(d):
            s = n0 * d[n0] / d[n]
            print(f"{workload:>8} {n:7d} {d[n]:9.2f} {s:9.2f} {s / n:10.2f}")

    import matplotlib
    matplotlib.use("agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    nmax = max(max(d) for d in data.values())
    ax.plot([1, nmax], [1, nmax], color="0.6", lw=1, ls="--")
    ax.annotate("ideal", (nmax, nmax), xytext=(-6, 4), textcoords="offset points", ha="right", va="bottom",
                fontsize=8, color="0.35")
    for workload, d in sorted(data.items()):
        n0 = min(d)
        x = sorted(d)
        y = [n0 * d[n0] / d[n] for n in x]
        ax.plot(x, y, color=COLORS.get(workload, "k"), lw=2, marker="o", ms=4, label=workload)
        ax.annotate(f"{workload} ({y[-1] / x[-1]:.0%})", (x[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=8, color="0.15")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ticks = [2 ** k for k in range(0, 12) if 2 ** k <= nmax]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
        axis.set_major_formatter(matplotlib.ticker.FixedFormatter([str(t) for t in ticks]))
        axis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("OpenMP threads")
    ax.set_ylabel(f"speed-up ({args.phase})")
    ax.set_title(f"{cpu}\n{label}", fontsize=8)
    ax.grid(True, which="major", color="0.9", lw=0.6)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_xlim(right=nmax * 2.2)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.directory, f"scaling.{ext}"), dpi=150)
    print(f"plots: {os.path.join(args.directory, 'scaling.png')}, .pdf")


if __name__ == "__main__":
    main()
