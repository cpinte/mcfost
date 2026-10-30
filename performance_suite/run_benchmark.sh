#!/bin/bash
# Wall-clock benchmark of mcfost: runs a fixed set of workloads with one or more
# binaries and writes the time of each phase (mcfost -timing) to a results file,
# compared afterwards with compare_benchmark.py. Run on a full node with a fixed
# OMP_NUM_THREADS (see benchmark.sbatch). The repetitions of the binaries are
# interleaved so that a drift of the node affects them equally.
#
# Usage: run_benchmark.sh -b <label>=<mcfost binary> [-b <label>=<binary> ...]
#                         [-n <repeats>] [-w "<workloads>"] [-t <threads>] [-o <results file>]
# -t sets OMP_NUM_THREADS (e.g. -t 4 for a quick local test), otherwise taken
# from the environment.
#
# Workloads (models from test_suite/test_data, see README.md):
#   W1 disk temperature, W2 1 um image (after W1), W3 CO lines,
#   W4 Phantom dump on a Voronoi mesh (temperature), W5 PAH temperature and SED

set -e

script_dir=$(cd "$(dirname "$0")" && pwd)
test_data=$script_dir/../test_suite/test_data

binaries=()
n_repeats=3
workloads="W1 W2 W3 W4 W5"
results=""

while getopts "b:n:w:t:o:h" opt; do
    case $opt in
        b) binaries+=("$OPTARG") ;;
        n) n_repeats=$OPTARG ;;
        w) workloads=$OPTARG ;;
        t) export OMP_NUM_THREADS=$OPTARG ;;
        o) results=$OPTARG ;;
        *) sed -n '2,18p' "$0" ; exit 1 ;;
    esac
done

if [ ${#binaries[@]} -eq 0 ]; then
    echo "ERROR: at least one binary is needed: -b <label>=<mcfost binary>" ; exit 1
fi
if [ -z "$MCFOST_UTILS" ]; then
    echo "ERROR: MCFOST_UTILS is not set" ; exit 1
fi
if [ ! -d "$test_data" ]; then
    echo "ERROR: $test_data not found (test_suite/get_test_data.sh)" ; exit 1
fi
if [ -z "$OMP_NUM_THREADS" ]; then
    echo "ERROR: set OMP_NUM_THREADS (timings are only comparable at a fixed number of threads)" ; exit 1
fi

# the binaries under test must not replace themselves with the latest release
export MCFOST_AUTO_UPDATE=0
export OMP_STACKSIZE=${OMP_STACKSIZE:-512M}
ulimit -s unlimited 2>/dev/null || ulimit -s hard

host=$(hostname -s 2>/dev/null || hostname)
date=$(date -u +%Y-%m-%dT%H:%M:%SZ)
[ -z "$results" ] && results=$script_dir/results/$(date -u +%Y%m%d_%H%M%S)_${host}${SLURM_JOB_ID:+_$SLURM_JOB_ID}.tsv
mkdir -p "$(dirname "$results")"

work=$(mktemp -d "${TMPDIR:-/tmp}/mcfost_benchmark.XXXXXX")
trap 'rm -rf "$work"' EXIT

# Parameter files of the workloads, from the test data
sed -E -e 's/^[[:space:]]*[^[:space:]]+([[:space:]]+nbr_photons_eq_th)/  1.28e7\1/' \
       -e 's/^([[:space:]]*)T([[:space:]]+hydostatic equilibrium)/\1F\2/' \
       "$test_data/ref4.1_hydrostatic/ref4.1_hydrostatic.para" > "$work/disk.para"
# CO cube with 100 velocity channels
sed -E -e 's/^([[:space:]]*[^[:space:]]+[[:space:]]+)[0-9]+([[:space:]]+vmax \(km.s\^-1\), n_speed)/\1100\2/' \
       "$test_data/ref3.0/ref3.0.para" > "$work/mol.para"
cp "$test_data/discF_00500/discF_00500.para" "$work/phantom.para"
# 10 PAH sizes instead of 1, so that the stochastic heating has a significant cost
sed -E -e 's/^[[:space:]]*[^[:space:]]+([[:space:]]+nbr_photons_eq_th)/  1.28e6\1/' \
       -e 's/^[[:space:]]*5e-4[[:space:]]+5e-4[[:space:]]+3.5[[:space:]]+1([[:space:]]+amin)/  3.5e-4  1e-3 3.5 10\1/' \
       "$test_data/ref4.1_PAH/ref4.1_PAH.para" > "$work/pah.para"
for f in "$test_data"/ref4.1_hydrostatic/*.lambda "$test_data"/ref4.1_PAH/*.lambda; do
    [ -f "$f" ] && cp "$f" "$work/"
done

cpu_model() {
    if command -v lscpu > /dev/null; then
        lscpu | sed -n 's/^Model name:[[:space:]]*//p' | head -1
    else
        sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown
    fi
}

modules=""
if type module > /dev/null 2>&1; then
    modules=$(module -t list 2>&1 | grep -v ":$" | tr '\n' ' ')
fi

{
    echo "# mcfost wall-clock benchmark"
    echo "# date: $date"
    echo "# host: $host"
    echo "# cpu: $(cpu_model)"
    echo "# omp_num_threads: $OMP_NUM_THREADS"
    echo "# omp_proc_bind: ${OMP_PROC_BIND:-unset}"
    echo "# omp_places: ${OMP_PLACES:-unset}"
    echo "# slurm_job: ${SLURM_JOB_ID:-none} ${SLURM_JOB_PARTITION:-}"
    echo "# modules: $modules"
    echo "# repeats: $n_repeats"
    for b in "${binaries[@]}"; do
        label=${b%%=*} ; bin=${b#*=}
        dir=$(cd "$(dirname "$bin")" && pwd)
        sha=$(git -C "$dir" rev-parse HEAD 2>/dev/null || echo unknown)
        git -C "$dir" diff --quiet HEAD -- . 2>/dev/null || sha="$sha-dirty"
        echo "# binary: $label $dir/$(basename "$bin") $sha"
    done
    printf "label\tworkload\trepeat\tphase\tseconds\tcalls\n"
} > "$results"

now() { # wall-clock time [s]
    if [ -n "$EPOCHREALTIME" ]; then echo "$EPOCHREALTIME" ; else python3 -c 'import time; print(time.time())' ; fi
}

# Binaries older than the -timing option (October 2026): only the total wall time
# of the run, measured here
declare_timing() { "$1" -h 2>&1 | grep -q -- "-timing" ; }

run() { # label binary workload repeat directory para options
    local label=$1 bin=$2 w=$3 r=$4 dir=$5 para=$6 ; shift 6
    local timing="" t0 t1
    declare_timing "$bin" && timing=-timing
    # each run in its own directory: mcfost caches the dust properties and the
    # Voronoi tessellation in temporary files of the working directory, which
    # would otherwise be reused by the next binary
    mkdir -p "$dir"
    for f in "$work"/*.lambda; do [ -f "$f" ] && cp "$f" "$dir/" ; done
    t0=$(now)
    (cd "$dir" && "$bin" "$work/$para" -root_dir "$dir" $timing "$@" > "$dir.log" 2>&1) || \
        { echo "ERROR: $label $w failed, see $dir.log" ; tail -20 "$dir.log" ; exit 1 ; }
    t1=$(now)
    [ -z "$timing" ] && awk -v t0="$t0" -v t1="$t1" 'BEGIN {printf " TIMING total %.3f 1\n", t1-t0}' >> "$dir.log"
    # SHA compiled into the binary (MCFOST_GIT=1 builds): git may not be available on compute nodes
    if ! grep -q "^# sha_binary: $label " "$results"; then
        echo "# sha_binary: $label $(sed -n 's/^ *Git SHA = *//p' "$dir.log" | head -1)" >> "$results"
    fi
    grep "^ TIMING" "$dir.log" | awk -v l="$label" -v w="$w" -v r="$r" \
        '{printf "%s\t%s\t%s\t%s\t%s\t%s\n", l, w, r, $2, $3, $4}' >> "$results"
    # mcfost output of the first repeat, kept next to the results
    if [ "$r" -eq 1 ]; then
        mkdir -p "${results%.tsv}_logs"
        cp "$dir.log" "${results%.tsv}_logs/${label}_$w.log"
    fi
    echo "$label $w repeat $r: $(awk '$2=="total" {print $3}' "$dir.log") s"
}

for r in $(seq 1 "$n_repeats"); do
    for b in "${binaries[@]}"; do
        label=${b%%=*} ; bin=$(cd "$(dirname "${b#*=}")" && pwd)/$(basename "${b#*=}")
        d=$work/${label}_$r
        mkdir -p "$d"
        for w in $workloads; do
            case $w in
                W1) run "$label" "$bin" W1 "$r" "$d/disk" disk.para -no_SED ;;
                W2) [ -d "$d/disk/data_th" ] || run "$label" "$bin" W1 "$r" "$d/disk" disk.para -no_SED
                    run "$label" "$bin" W2 "$r" "$d/disk" disk.para -img 1.0 ;;
                W3) run "$label" "$bin" W3 "$r" "$d/mol" mol.para -mol -no_SED ;;
                W4) run "$label" "$bin" W4 "$r" "$d/phantom" phantom.para \
                        -phantom "$test_data/discF_00500/discF_00500" -not_random_Voronoi -no_SED ;;
                W5) run "$label" "$bin" W5 "$r" "$d/pah" pah.para ;;
                *) echo "ERROR: unknown workload $w" ; exit 1 ;;
            esac
        done
        rm -rf "$d"
    done
done

echo "Results written to $results"
