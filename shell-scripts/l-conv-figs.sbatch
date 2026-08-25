#!/bin/bash
#SBATCH -J l-conv-figs
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 08:00:00
#SBATCH -o logs/l-conv-figs_%j.out
#SBATCH -e logs/l-conv-figs_%j.err
#
# Figures for the l-convergence study, after the l-conv array has produced
# out/lconv/observables_<tag>.csv.
#
#   1. every figure for every model, from the narrow l <= 3 tables as usual
#   2. a ladder of figure sets, one per l_max, from the wide sweep truncated
#      at that daughter degree -- the same familiar figures side by side, so a
#      difference between rungs is visible without a summary statistic
#
# Submitted with afterany, not afterok: one model failing its sweep must not
# block the figures for the four that succeeded. Figures whose CSV is missing
# report SKIPPED rather than failing the job.

set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORKDIR"
mkdir -p logs out

source shell-scripts/config.sh
setup_python
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

# 1. The standard set, narrow tables.
python -u scripts/make_figures.py --all-tags

# 2. The ladder. Only the figures that read a coupling table -- the rest
#    (propagation, gamma, xi_kappa, kappa_cum) do not depend on l_max at all.
LADDER="mu eth eth_pg channels channels_pg_types detuning detuning_pg
        pg_census lowfreq spectrum lconv l_ingredients amplitude_sweep
        daughter_energy"

for N in 3 6 10 15 25; do
    echo "=== l_max = $N  $(date -Is) ==="
    python -u scripts/make_figures.py --all-tags \
        --source lconv --l-max-cut "$N" \
        --figs "out/plots_lmax/lmax$(printf '%02d' "$N")" \
        --only $LADDER \
      || echo "FAILED: l_max=$N"
done

echo "=== figures done  $(date -Is) ==="
