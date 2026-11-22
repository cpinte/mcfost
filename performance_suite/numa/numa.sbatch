#!/bin/bash -l
# NUMA test of mcfost on a full node (no code change): the same binary
#   default    : OMP_PROC_BIND=close, OMP_PLACES=cores (as benchmark.sbatch)
#   interleave : numactl --interleave=all (pages spread over the NUMA nodes)
#   spread     : OMP_PROC_BIND=spread
# with the repeats interleaved (run_benchmark.sh). Workloads: WORKLOADS (default "W1L W2 W4L": W1 and W4 with
# 10 x the packets, so that the Monte Carlo part lasts long enough on a full node). If interleave is clearly faster, the shared
# arrays are bandwidth-limited and the per-thread accumulators should be first-touched in parallel
# (performance review P-O10, P-S11; tracker R4).
#
# From performance_suite/:
#   sbatch --cpus-per-task=<cores of a node> --export=ALL,MCFOST_BIN=$HOME/mcfost/src/mcfost numa/numa.sbatch
# Results: results/<date>_<host>_<job>.tsv, compared with
#   python3 compare_benchmark.py results/<file>.tsv --baseline default --candidate interleave
#   python3 compare_benchmark.py results/<file>.tsv --baseline default --candidate spread
#
#SBATCH --job-name=mcfost-numa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=mcfost-numa_%j.log
#SBATCH --mem-per-cpu=2g
#SBATCH --time=4:00:00
#SBATCH --tmp=20G

ulimit -s unlimited 2>/dev/null || ulimit -s hard
export MCFOST_UTILS=${MCFOST_UTILS:-$HOME/mcfost/utils}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-$(nproc)}
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=1024m
cd "${SLURM_SUBMIT_DIR:-.}"

[ -x "$MCFOST_BIN" ] || { echo "ERROR: MCFOST_BIN ($MCFOST_BIN) is not an executable" ; exit 1 ; }
command -v numactl > /dev/null || { echo "ERROR: numactl not found on $(hostname)" ; exit 1 ; }
# node topology, for the interpretation
lscpu | grep -E "Model name|^CPU\(s\)|Socket|NUMA"
numactl --hardware
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

./run_benchmark.sh -n 3 -w "${WORKLOADS:-W1L W2 W4L}" \
    -b default="$MCFOST_BIN" \
    -b interleave="$PWD/numa/mcfost_interleave.sh" \
    -b spread="$PWD/numa/mcfost_spread.sh"
