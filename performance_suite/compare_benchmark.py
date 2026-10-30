#!/usr/bin/env python3
"""Compare mcfost benchmark results (run_benchmark.sh).

Median wall-clock time of each workload and phase for a baseline and a
candidate binary, from one or several results files. A phase is flagged when
the candidate is slower than the baseline by more than the threshold and by
more than the spread of the repeats. Phases shorter than --min-seconds are
reported but never flagged. Exit status 1 if any phase is flagged.

Usage: compare_benchmark.py results.tsv [more.tsv ...] [--baseline <label>] [--candidate <label>]
By default the baseline and candidate are the first two binaries of the files.
"""

import argparse
import statistics
import sys
from collections import defaultdict


# Phases that do not overlap (dust_prop and opacity are within init or
# step2_setup): "other" = total - their sum, the time no timer covers
TOP_LEVEL = ("setup", "init", "mc_temperature", "temperature", "temperature_nRE", "step2_setup",
             "mc_step2", "source_fct", "ray_tracing", "output", "mol", "atom")


def read_results(files):
    times = defaultdict(list)  # (label, workload, phase) -> [seconds]
    runs = defaultdict(dict)   # (file, label, workload, repeat) -> {phase: seconds}
    meta = defaultdict(set)    # key -> values, to check the runs are comparable
    labels = []                # binaries in order of appearance
    for name in files:
        with open(name) as f:
            for line in f:
                if line.startswith("#"):
                    key, _, value = line[1:].partition(":")
                    if key.strip() == "binary" and value.split()[0] not in labels:
                        labels.append(value.split()[0])
                    if key.strip() in ("host", "cpu", "omp_num_threads", "modules"):
                        meta[key.strip()].add(value.strip())
                    continue
                cols = line.split("\t")
                if cols[0] == "label" or len(cols) < 6:
                    continue
                times[(cols[0], cols[1], cols[3])].append(float(cols[4]))
                runs[(name, cols[0], cols[1], cols[2])][cols[3]] = float(cols[4])
    for (name, label, workload, repeat), phases in runs.items():
        if len(phases) > 1 and "total" in phases:  # binaries with -timing only
            times[(label, workload, "other")].append(
                phases["total"] - sum(t for p, t in phases.items() if p in TOP_LEVEL))
    return times, meta, labels


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+")
    parser.add_argument("--baseline", help="label of the reference binary (default: first binary)")
    parser.add_argument("--candidate", help="label of the binary under test (default: second binary)")
    parser.add_argument("--threshold", type=float, default=0.10, help="relative slowdown flagged (default 0.10)")
    parser.add_argument("--min-seconds", type=float, default=0.5, help="phases shorter than this are not flagged")
    args = parser.parse_args()

    times, meta, labels = read_results(args.files)
    if args.baseline is None or args.candidate is None:
        if len(labels) < 2:
            sys.exit("ERROR: fewer than two binaries in the results: give --baseline and --candidate")
        args.baseline = args.baseline or labels[0]
        args.candidate = args.candidate or labels[1]
    for key, values in meta.items():
        if len(values) > 1:
            print(f"WARNING: results from different {key}: {sorted(values)}")

    keys = sorted({(w, p) for (l, w, p) in times if l in (args.baseline, args.candidate)})
    flagged = []
    print(f"{'workload':8} {'phase':16} {args.baseline:>10} {args.candidate:>10} {'ratio':>7} {'spread':>7}")
    for w, p in keys:
        base = times.get((args.baseline, w, p))
        cand = times.get((args.candidate, w, p))
        if not base or not cand:
            print(f"{w:8} {p:16} {'-' if not base else statistics.median(base):>10} "
                  f"{'-' if not cand else statistics.median(cand):>10}   (missing)")
            continue
        mb, mc = statistics.median(base), statistics.median(cand)
        # relative spread of the repeats: noise level of this phase
        spread = max((max(t) - min(t)) / statistics.median(t) if statistics.median(t) > 0 else 0.0
                     for t in (base, cand))
        ratio = mc / mb if mb > 0 else float("inf")
        flag = ""
        if max(mb, mc) >= args.min_seconds and ratio - 1 > max(args.threshold, spread):
            flag = "  SLOWER"
            flagged.append((w, p))
        elif max(mb, mc) >= args.min_seconds and 1 - ratio > max(args.threshold, spread):
            flag = "  faster"
        print(f"{w:8} {p:16} {mb:10.3f} {mc:10.3f} {ratio:7.3f} {spread:7.3f}{flag}")

    if flagged:
        print(f"\n{len(flagged)} phase(s) slower than {args.baseline} by more than "
              f"{args.threshold:.0%} and the spread of the repeats")
        sys.exit(1)


if __name__ == "__main__":
    main()
