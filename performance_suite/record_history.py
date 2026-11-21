#!/usr/bin/env python3
"""Append benchmark results (run_benchmark.sh) to the history of a machine.

One line per binary, workload and phase in history/<machine>.csv: median, min and max over the
repeats, and the ratio to the reference binary of the same job. Timings of different jobs are only
compared through this ratio: a shared node varies from job to job, but the reference binary runs in
the same job, interleaved with the others, and sees the same node (benchmark.sbatch: REFERENCE_BIN,
label "ref"). Lines already in the history (same job, binary, workload, phase) are not added twice.

Usage: record_history.py <results.tsv> [more.tsv ...] --machine <name> [--reference <label>]
  --machine: node type, e.g. nt-milan or ozstar-skylake (one history file per node type)
  --reference: label of the reference binary (default: "ref" if present, otherwise none)
"""

import argparse
import csv
import os
import statistics
import sys
from collections import defaultdict

from compare_benchmark import TOP_LEVEL

FIELDS = ["date", "machine", "cpu", "threads", "proc_bind", "job", "label", "sha", "dirty", "workload", "phase",
          "median_s", "min_s", "max_s", "n", "ref_label", "ref_sha", "ref_median_s", "ratio", "results_file"]
PHASES = ("total", "other") + TOP_LEVEL + ("dust_prop", "opacity")


def read(name):
    """Metadata, binaries {label: sha, dirty} and times {(label, workload, phase): [s]} of a results file."""
    meta, shas, times = {}, {}, defaultdict(list)
    runs = defaultdict(dict)
    with open(name) as f:
        for line in f:
            if line.startswith("#"):
                key, _, value = line[1:].partition(":")
                key, value = key.strip(), value.strip()
                if key == "binary" and len(value.split()) >= 3:
                    label, sha = value.split()[0], value.split()[2]
                    shas.setdefault(label, (sha.replace("-dirty", ""), sha.endswith("-dirty")))
                elif key == "sha_binary" and len(value.split()) >= 2:
                    # SHA compiled into the binary: wins over the git directory of the binary
                    label, sha = value.split()[0], value.split()[1]
                    shas[label] = (sha, shas.get(label, ("", False))[1])
                else:
                    meta[key] = value
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] == "label" or len(cols) < 6:
                continue
            runs[(cols[0], cols[1], cols[2])][cols[3]] = float(cols[4])
    for (label, workload, _), phases in runs.items():
        for phase, t in phases.items():
            times[(label, workload, phase)].append(t)
        if len(phases) > 1 and "total" in phases:
            times[(label, workload, "other")].append(
                phases["total"] - sum(t for p, t in phases.items() if p in TOP_LEVEL))
    return meta, shas, times


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+")
    parser.add_argument("--machine", required=True, help="node type, e.g. nt-milan (history/<machine>.csv)")
    parser.add_argument("--reference", help='label of the reference binary (default: "ref" if present)')
    args = parser.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    hist = os.path.join(here, "history", f"{args.machine}.csv")
    os.makedirs(os.path.dirname(hist), exist_ok=True)
    seen = set()
    if os.path.exists(hist):
        with open(hist) as f:
            for row in csv.DictReader(f):
                seen.add((row["job"], row["date"], row["label"], row["workload"], row["phase"]))

    new_rows = []
    for name in args.files:
        meta, shas, times = read(name)
        ref = args.reference or ("ref" if "ref" in shas else None)
        if ref and ref not in shas:
            sys.exit(f"ERROR: reference binary '{ref}' not in {name} (binaries: {', '.join(shas)})")
        job = meta.get("slurm_job", "none").split()[0]
        for (label, workload, phase), ts in sorted(times.items()):
            if phase not in PHASES:
                continue
            key = (job, meta.get("date", ""), label, workload, phase)
            if key in seen:
                continue
            med = statistics.median(ts)
            rmed = statistics.median(times[(ref, workload, phase)]) if ref and (ref, workload, phase) in times else None
            new_rows.append({
                "date": meta.get("date", ""), "machine": args.machine, "cpu": meta.get("cpu", ""),
                "threads": meta.get("omp_num_threads", ""), "proc_bind": meta.get("omp_proc_bind", ""),
                "job": job, "label": label, "sha": shas.get(label, ("", False))[0],
                "dirty": int(shas.get(label, ("", False))[1]), "workload": workload, "phase": phase,
                "median_s": f"{med:.3f}", "min_s": f"{min(ts):.3f}", "max_s": f"{max(ts):.3f}", "n": len(ts),
                "ref_label": ref or "", "ref_sha": shas[ref][0] if ref else "",
                "ref_median_s": f"{rmed:.3f}" if rmed else "",
                "ratio": f"{med / rmed:.4f}" if rmed else "", "results_file": os.path.basename(name)})
            seen.add(key)

    write_header = not os.path.exists(hist)
    with open(hist, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if write_header:
            w.writeheader()
        w.writerows(new_rows)
    print(f"{len(new_rows)} lines added to {os.path.relpath(hist)}")


if __name__ == "__main__":
    main()
