#!/bin/bash
# Shows the benchmark jobs still queued on OzSTAR and copies the results files
# into results/ (run it yourself, needs ssh access). Settings in ozstar.conf.
# Then: ./compare_benchmark.py results/<file>.tsv

set -e

script_dir=$(cd "$(dirname "$0")" && pwd)
conf=$script_dir/ozstar.conf
[ -f "$conf" ] || { echo "ERROR: $conf not found (copy ozstar.conf.example)" ; exit 1 ; }
# shellcheck source=/dev/null
source "$conf"

echo "Benchmark jobs on $OZSTAR_HOST:"
ssh "$OZSTAR_HOST" "squeue -u \$USER -n mcfost-benchmark -o '%.10i %.8T %.10M %R'" || true

rsync -av "$OZSTAR_HOST:$REMOTE_BASE/results/" "$script_dir/results/"
rsync -a "$OZSTAR_HOST:$REMOTE_BASE/logs/" "$script_dir/results/logs/"
echo "Job and build logs in $script_dir/results/logs (not under version control)"
