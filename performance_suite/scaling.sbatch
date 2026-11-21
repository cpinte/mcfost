#!/bin/bash -l
# Thread scaling of one mcfost binary on a full node: run_benchmark.sh at each number of threads,
# one results file per thread count in results/scaling_<date>_<host>_<job>/t<n>.tsv, plotted by
# plot_scaling.py. Threads pinned to cores, packed (close): the first sockets fill first.
# From performance_suite/:
#   sbatch --cpus-per-task=<cores of a node> --export=ALL,MCFOST_BIN=<path>/mcfost scaling.sbatch
# Optional (in --export): THREADS ("1 2 4 8 16 32 64", counts above the allocation are skipped),
# WORKLOADS ("W1 W4"), REPEATS (2), LABEL (mcfost).
#
#SBATCH --job-name=mcfost-scaling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=mcfost-scaling_%j.log
#SBATCH --mem-per-cpu=2g
#SBATCH --time=3:00:00

ulimit -s unlimited 2>/dev/null || ulimit -s hard
export MCFOST_UTILS=${MCFOST_UTILS:-$HOME/mcfost/utils}
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=1024m
cd "${SLURM_SUBMIT_DIR:-.}"

[ -x "$MCFOST_BIN" ] || { echo "ERROR: MCFOST_BIN ($MCFOST_BIN) is not an executable" ; exit 1 ; }
lscpu | grep -E "Model name|^CPU\(s\)|Socket|NUMA"
command -v numactl > /dev/null && numactl --hardware

max=${SLURM_CPUS_PER_TASK:-$(nproc)}
out=results/scaling_$(date -u +%Y%m%d_%H%M%S)_$(hostname -s)_${SLURM_JOB_ID:-local}
mkdir -p "$out"
for n in ${THREADS:-1 2 4 8 16 32 64}; do
    [ "$n" -le "$max" ] || { echo "skipping $n threads (allocation: $max)" ; continue ; }
    ./run_benchmark.sh -t "$n" -n "${REPEATS:-2}" -w "${WORKLOADS:-W1 W4}" \
        -b "${LABEL:-mcfost}=$MCFOST_BIN" -o "$out/t$n.tsv"
done
echo "Scaling results in $out (python3 plot_scaling.py $out)"
