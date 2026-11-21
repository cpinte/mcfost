#!/bin/bash
# mcfost with the threads spread over the node instead of packed (numa.sbatch)
OMP_PROC_BIND=spread exec "$MCFOST_BIN" "$@"
