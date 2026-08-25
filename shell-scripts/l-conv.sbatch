#!/bin/bash
#SBATCH -J l-conv
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 24:00:00
#SBATCH -a 0-4
#SBATCH -o logs/l-conv_%A_%a.out
#SBATCH -e logs/l-conv_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=txr7972@mavs.uta.edu

set -euo pipefail

WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORKDIR"
mkdir -p logs out

source shell-scripts/config.sh
setup_python

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

TAGS=(dsct_M2.0 dsct_M2.2_T7888 dsct_M2.0_T7696 dsct_M1.85_T7555 dsct_M1.7_T7750)
TAG="${TAGS[$SLURM_ARRAY_TASK_ID]}"
JOBS="${SLURM_CPUS_PER_TASK:-4}"

echo "=== $TAG on $JOBS core(s)  $(date -Is) ==="
python -u scripts/three_mode_search.py \
    --model "models/$TAG" \
    --detail-dir detail_wide --inlist gyre_ad_wide.in \
    --nad summary_nad_wide.h5 summary_nad_wide_hi.h5 \
    --l-max 25 --parents-only --gamma tot --m000 \
    --out out/lconv --tag "$TAG" --force \
    -j "$JOBS"
    echo "=== $TAG done  $(date -Is) ==="
