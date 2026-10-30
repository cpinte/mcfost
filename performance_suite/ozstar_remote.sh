#!/bin/bash -l
# Remote part of submit_ozstar.sh, run on the OzSTAR login node through ssh:
# clones/updates mcfost, builds each commit in its own worktree (build/<sha>),
# installs the test data (cached per archive), then submits benchmark.sbatch.
# The benchmark scripts are the ones copied by submit_ozstar.sh into performance_suite/
# (not those of the commits, which may predate them).
#
# Arguments: <remote base> <build setup> <build jobs> <job setup> <mcfost utils>
#            <sbatch options> <repeats> <workloads> <label>=<sha> [<label>=<sha> ...]
# The test data are those of the last commit (the candidate).

set -e

base=$HOME/$1 ; build_setup=$2 ; build_jobs=$3 ; job_setup=$4 ; utils=$5
sbatch_options=$6 ; repeats=$7 ; workloads=$8
shift 8

mkdir -p "$base/build" "$base/test_data_cache" "$base/results" "$base/logs"
cd "$base"
[ -d repo ] || git clone -q https://github.com/cpinte/mcfost.git repo
git -C repo fetch -q origin

eval "$build_setup"

args=()
for b in "$@"; do
    label=${b%%=*} ; sha=${b#*=}
    wt=$base/build/$sha
    if [ ! -d "$wt" ]; then
        git -C repo worktree add -q --detach "$wt" "$sha"
    fi
    if [ ! -x "$wt/src/mcfost" ]; then
        echo "Building $label ($sha)"
        (cd "$wt/src" && MCFOST_GIT=1 make static=yes -j"$build_jobs" > "$base/logs/build_$sha.log" 2>&1) || \
            { echo "ERROR: build of $sha failed, see $base/logs/build_$sha.log" ; tail -20 "$base/logs/build_$sha.log" ; exit 1 ; }
    fi
    args+=(-b "$label=$wt/src/mcfost")
done
echo "ldd of the candidate (libraries must be reachable from the compute nodes):"
ldd "$wt/src/mcfost" | grep -v -E "linux-vdso|ld-linux" || true

# test data of the candidate, downloaded once per archive, in ../test_suite/test_data of the scripts
name=$(sed -n 's/^TEST_DATA_FILE="\(.*\)\.tgz"/\1/p' "$wt/test_suite/get_test_data.sh")
if [ ! -d "$base/test_data_cache/$name/test_data" ]; then
    echo "Downloading the test data $name"
    mkdir -p "$base/test_data_cache/$name"
    (cd "$base/test_data_cache/$name" && bash "$wt/test_suite/get_test_data.sh" && rm -f ./*.tgz)
fi
mkdir -p "$base/test_suite"
rm -f "$base/test_suite/test_data"
ln -s "$base/test_data_cache/$name/test_data" "$base/test_suite/test_data"

# run-time environment of the job
{
    echo "$job_setup"
    echo "export MCFOST_UTILS=$utils"
} > "$base/job_setup.sh"

cd "$base/performance_suite"
results=$base/results/$(date -u +%Y%m%d_%H%M%S)_$(printf '%s_' "${@%%=*}" | sed 's/_$//').tsv
# shellcheck disable=SC2086
sbatch $sbatch_options --output="$base/logs/benchmark_%j.log" \
       --export=ALL,JOB_SETUP_FILE="$base/job_setup.sh" \
       benchmark.sbatch "${args[@]}" -n "$repeats" -w "$workloads" -o "$results"
echo "Results will be written to $results"
