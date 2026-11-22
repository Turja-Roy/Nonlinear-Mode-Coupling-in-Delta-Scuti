#!/usr/bin/env bash
#SBATCH -J network_runs
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8              # cluster admin, 2026-11: at most 8 cores a node
#SBATCH --mem=16000       # over half a node (31868 MB): no two 8-core tasks of ours share a node
#SBATCH -t 12:00:00
#SBATCH -o logs/network_runs_%A_%a.log
# One task of a sweep's run array, submitted by shell-scripts/network_sweep.sh (UTA HPC): every
# SLURM_ARRAY_TASK_COUNT-th run folder of the list, from this task's index, all at once, each plotted.
#   network_runs.sbatch <list of run folders>
# WALL [min] (from the sweep, default 660): a run stops itself at 90% of it and keeps what it has
# (stop = wall); `timeout` at WALL is only the backstop, logged as TIMEOUT with no summary. Keep WALL below
# --time, or Slurm kills the task first and its plots never run.

set -uo pipefail
source shell-scripts/config.sh
setup_python
JOBS=${SLURM_CPUS_PER_TASK:-8}
WALL=${WALL:-660}
source shell-scripts/network_lib.sh
one() {   # run folder
    timeout "${WALL}m" ./build/network_run "$1/network.data" --out "$1" $RUN_ARGS \
        --max-seconds $((WALL * 54)) > "$1/log.txt" 2>&1
    case $? in
        0) python scripts/network_plot.py --data "$1" --title "$TAG ${1#"$OUT/$TAG/"}" >> "$1/log.txt" 2>&1 \
               || echo "PLOT FAILED $1" ;;
        124) echo "TIMEOUT $1 (> ${WALL} min)" ;;
        *) echo "FAILED $1 (see $1/log.txt)" ;;
    esac
}
export -f one
export RUN_ARGS WALL OUT TAG
echo "=== $TAG, task $SLURM_ARRAY_TASK_ID of $SLURM_ARRAY_TASK_COUNT on $(hostname)"
awk -v k="$SLURM_ARRAY_TASK_ID" -v n="$SLURM_ARRAY_TASK_COUNT" '(NR - 1) % n == k' "$1" \
    | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
