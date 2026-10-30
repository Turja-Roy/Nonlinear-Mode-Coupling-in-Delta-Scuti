#!/usr/bin/env bash
#SBATCH --job-name=network_checks
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=logs/network_checks_%j.log
# The checks of Plans/network-growth.md sec. 8 that need runs of their own (the
# rest read network_sweep.sh's output). Independent of the sweep: submit both.
#
#   sbatch shell-scripts/network_checks.sh
#
# out/checks/mw25/             MW25 Figs 1-6 under Adams, beside the old rk8pd out/mw25_cpp.
# out/checks/permode/<net>/    the survey's 9 networks (data/atol_bench = out/atol_bench/nets)
#                              at the survey's settings, per-mode error weight; vs
#                              out/atol_bench/runs/adams/<net>/r1e-8_a1e-300 (per component).
# out/checks/ladder/<run>/     dsct_M2.0 ladder_1p d4 (the sweep's file), f_p 0.9: base
#                              (rtol 1e-8, f_d 1e-3), rtol 1e-10, f_d 1e-2, f_d 1e-4.

set -euo pipefail
source shell-scripts/config.sh
setup_python
mkdir -p build
flock build/.compile.lock ./compile.sh -j "${SLURM_CPUS_PER_TASK:-4}" network_build network_run mw25
export OMP_NUM_THREADS=1
OUT=out/checks
rm -rf "$OUT" data/checks
mkdir -p "$OUT"

(./build/mw25 --out "$OUT/mw25" && python scripts/mw25_plot.py --data "$OUT/mw25") \
    > "$OUT/mw25.log" 2>&1 || echo "mw25 FAILED (see $OUT/mw25.log)" &

./build/network_build --model models/dsct_M2.0 --kappa-cache out/four_mode_dsct_M2.0.kappa_cache.tsv \
    -j "${SLURM_CPUS_PER_TASK:-4}" --parent "$(model_field dsct_M2.0 4)" -N 4 --closure-cut 0.005 \
    --out data/checks
L=data/checks/cut0.15_d4.data
BASE="--seed 1 --e-max 1e-2 --t-end 1e13 --max-seconds 36000"

run() {   # run name, network file, network_run args...
    local d=$OUT/$1; shift
    mkdir -p "$d"
    ./build/network_run "$@" --out "$d" > "$d/log.txt" 2>&1 || echo "FAILED $d"
    python scripts/network_plot.py --data "$d" --title "${d#"$OUT"/}" >> "$d/log.txt" 2>&1 || echo "PLOT FAILED $d"
}
export -f run
export OUT
{
    for f in data/atol_bench/*.data; do
        echo "permode/$(basename "$f" .data) $f --q0 1e-6 --q0-daughter 1e-9 --eth-parent 0 --eth-daughter 0 $BASE --rtol 1e-8"
    done
    echo "ladder/base $L $BASE --rtol 1e-8"
    echo "ladder/rtol1e-10 $L $BASE --rtol 1e-10"
    echo "ladder/fd1e-2 $L $BASE --rtol 1e-8 --eth-daughter 1e-2"
    echo "ladder/fd1e-4 $L $BASE --rtol 1e-8 --eth-daughter 1e-4"
} | xargs -P "${SLURM_CPUS_PER_TASK:-4}" -L 1 bash -c 'run "$@"' _
wait

python - <<'EOF'
import pathlib, pandas as pd
rows = [pd.read_csv(p).assign(run=str(p.parent)) for p in sorted(pathlib.Path("out/checks").glob("*/*/summary.csv"))]
print(pd.concat(rows)[["run", "stop", "seconds", "steps", "failed_steps", "E_parents", "swing"]].to_string(index=False))
EOF
