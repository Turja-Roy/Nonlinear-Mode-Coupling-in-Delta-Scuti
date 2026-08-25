#!/bin/bash
#SBATCH -J l-conv-figs
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -o logs/l-conv-figs_%j.out
#SBATCH -e logs/l-conv-figs_%j.err
#
# Every figure for every model, after the l-conv array has produced
# out/lconv/observables_<tag>.csv. Picks up fig_lconv and fig_l_ingredients,
# and regenerates fig_gamma and fig_detuning_pg for the two merged changes.
# Figures whose CSV is missing report SKIPPED rather than failing the job.
#
# Submitted with afterany, not afterok: one model failing its sweep must not
# block the figures for the four that succeeded.

set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORKDIR"
mkdir -p logs out

source shell-scripts/config.sh
setup_python
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

python -u scripts/make_figures.py --all-tags
echo "=== figures done  $(date -Is) ==="
