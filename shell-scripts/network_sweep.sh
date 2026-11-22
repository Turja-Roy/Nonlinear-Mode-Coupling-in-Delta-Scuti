#!/usr/bin/env bash
#SBATCH -J network_sweep
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8              # cluster admin, 2026-11: at most 8 cores a node
#SBATCH --mem=16000       # over half a node (31868 MB): no two 8-core tasks of ours share a node
#SBATCH -t 2:00:00
#SBATCH -a 0-4
#SBATCH -o logs/network_sweep_%A_%a.log
# Network ladders on the UTA HPC cluster, one array task per model in config.sh's MODELS:
#
#   SETS=2p_sum sbatch --array=1 shell-scripts/network_sweep.sh
#   BUILD_ONLY=1 SETS=1p sbatch --array=1 shell-scripts/network_sweep.sh     # build and place, no runs
#
# Sets, settings and outputs: shell-scripts/network_lib.sh. This job compiles network_build/network_run
# (array tasks take turns on a lock), builds the sets, and lists the networks still to run; then it
# submits shell-scripts/network_runs.sbatch as an array of 8-core tasks, one node each, up to 8 networks
# a task (WALL [min] a run, default 660), so every network runs at once on as many nodes as it takes,
# each task starting as soon as a node is free; then a job after them that plots any run without a
# plot, the convergence plots and networks.md (shell-scripts/network_post.sh).

set -euo pipefail
source shell-scripts/config.sh
setup_python
mkdir -p build
flock build/.compile.lock ./compile.sh -j "${SLURM_CPUS_PER_TASK:-4}" network_build network_run
JOBS=${SLURM_CPUS_PER_TASK:-4}
export WALL=${WALL:-660}
source shell-scripts/network_lib.sh
echo "=== $MODEL (array task ${SLURM_ARRAY_TASK_ID:-none})"
[[ -z ${SETS// } ]] && exit 0
build_sets
[[ ${BUILD_ONLY:-0} == 1 ]] && exit 0

L=$NET_DATA/$TAG/todo.txt
todo > "$L"
n=$(wc -l < "$L")
if (( n == 0 )); then shell-scripts/network_post.sh "$PWD" "$OUT/$TAG"; exit 0; fi
k=$(( (n + 7) / 8 ))
export MODEL
a=$(sbatch --parsable --array=0-$((k - 1)) shell-scripts/network_runs.sbatch "$L")
sbatch --dependency=afterany:${a%%;*} -J network_post -p normal -c 1 --mem=4000 -t 2:00:00 \
    -o logs/network_post_%j.log --wrap "shell-scripts/network_post.sh $PWD $OUT/$TAG"
echo "$n networks: run array ${a%%;*} on $k nodes, then network_post"
