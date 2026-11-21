#!/bin/bash
# mcfost with its memory pages interleaved over all NUMA nodes (numa.sbatch)
exec numactl --interleave=all "$MCFOST_BIN" "$@"
