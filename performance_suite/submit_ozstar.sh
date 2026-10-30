#!/bin/bash
# Submits the wall-clock benchmark on OzSTAR (run it yourself, needs ssh access):
# resolves the git references locally, then on the login node builds each commit
# in its own directory and submits benchmark.sbatch. The commits must be pushed
# to GitHub. Settings in ozstar.conf (copy ozstar.conf.example).
#
# Usage: submit_ozstar.sh [-n <repeats>] [-w "<workloads>"] <baseline ref> <candidate ref> [...]
#   e.g. submit_ozstar.sh main pah_nre
# Labels are <ref>-<short sha>; the first reference is the baseline. Binaries
# older than the -timing option only get the total wall time.
# Results: fetch_ozstar.sh, then compare_benchmark.py.

set -e

script_dir=$(cd "$(dirname "$0")" && pwd)
conf=$script_dir/ozstar.conf
[ -f "$conf" ] || { echo "ERROR: $conf not found (copy ozstar.conf.example)" ; exit 1 ; }
# shellcheck source=/dev/null
source "$conf"

repeats=3
workloads="W1 W2 W3 W4 W5"
while getopts "n:w:h" opt; do
    case $opt in
        n) repeats=$OPTARG ;;
        w) workloads=$OPTARG ;;
        *) sed -n '2,10p' "$0" ; exit 1 ;;
    esac
done
shift $((OPTIND-1))
[ $# -ge 1 ] || { sed -n '2,10p' "$0" ; exit 1 ; }

git -C "$script_dir" fetch -q origin
pairs=()
for ref in "$@"; do
    sha=$(git -C "$script_dir" rev-parse --verify "$ref^{commit}") || { echo "ERROR: unknown reference $ref" ; exit 1 ; }
    if [ -z "$(git -C "$script_dir" branch -r --contains "$sha" 2>/dev/null)" ]; then
        echo "ERROR: $ref ($sha) is not on GitHub: push it first" ; exit 1
    fi
    pairs+=("$(echo "$ref" | tr '/' '_')-${sha:0:7}=$sha")
done
echo "Submitting on $OZSTAR_HOST: ${pairs[*]}"
# benchmark scripts of this checkout (the commits may predate them)
ssh "$OZSTAR_HOST" "mkdir -p $REMOTE_BASE/performance_suite"
rsync -a "$script_dir/run_benchmark.sh" "$script_dir/benchmark.sbatch" "$script_dir/compare_benchmark.py" \
      "$OZSTAR_HOST:$REMOTE_BASE/performance_suite/"

remote_args=("$REMOTE_BASE" "$BUILD_SETUP" "$BUILD_JOBS" "$JOB_SETUP" "$REMOTE_MCFOST_UTILS"
             "$SBATCH_OPTIONS" "$repeats" "$workloads" "${pairs[@]}")
# arguments quoted for the remote shell
ssh "$OZSTAR_HOST" "bash -l -s -- $(printf '%q ' "${remote_args[@]}")" < "$script_dir/ozstar_remote.sh"
