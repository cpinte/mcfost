#!/bin/bash -l
# Slurm job for the mcfost wall-clock benchmark (run_benchmark.sh).
# Submitted by submit_ozstar.sh, or by hand from this directory:
#   sbatch --account=<account> benchmark.sbatch -b main=<path>/mcfost -b new=<path>/mcfost [-n 3] [-w "W1 W5"]
# The arguments are passed to run_benchmark.sh. Always use the same node type:
# timings are only compared between runs on the same CPU model, compiler and
# number of threads. Account, partition, node type and, where the site allows
# it, --exclusive are given on the sbatch command line (SBATCH_OPTIONS in
# ozstar.conf). OzSTAR does not allow --exclusive: request all the cores of a
# node with --cpus-per-task to avoid sharing it with other jobs.
#
#SBATCH --job-name=mcfost-benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=mcfost-benchmark_%j.log
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4g
#SBATCH --time=1:00:00

# run-time environment of the binaries (e.g. LD_LIBRARY_PATH), written by submit_ozstar.sh
[ -n "$JOB_SETUP_FILE" ] && source "$JOB_SETUP_FILE"

ulimit -s unlimited 2>/dev/null || ulimit -s hard

export MCFOST_UTILS=${MCFOST_UTILS:-$HOME/mcfost/utils}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-$(nproc)}
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_STACKSIZE=1024m

cd "${SLURM_SUBMIT_DIR:-.}"
./run_benchmark.sh "$@"
