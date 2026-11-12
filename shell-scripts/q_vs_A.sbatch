#!/usr/bin/env bash
#SBATCH -J model_checks
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16000
#SBATCH -t 12:00:00
#SBATCH -o logs/model_checks_%j.log
# One model: Fig. 4 rates, direct-triplet rankings, dq_dt vs dA_dt, energy-budget redraws.
# --clean: build/ here holds binaries copied from the HEP cluster.

set -euo pipefail
source shell-scripts/config.sh
setup_python
mkdir -p build logs
flock build/.compile.lock ./compile.sh --clean -j 4 --test
export OMP_NUM_THREADS=1

TAG=dsct_M2.0_T7696
MODEL=models/$TAG
TITLE='$M = 2.0\,M_\odot$, $\log g = 3.93$, $T_{\rm eff} = 7695$ K'

echo "=== damping"
./build/damping --model "$MODEL"
python scripts/damping_plot.py --data "out/damping_$TAG.csv" --title "$TITLE"

echo "=== direct triplets"
./build/network_build --model "$MODEL" --kappa-cache "out/four_mode_$TAG.kappa_cache.tsv" -j 4 --list-direct \
    > "out/direct_candidates_$TAG.txt"
head -3 "out/direct_candidates_$TAG.txt"

echo "=== slow_amplitude (background)"
SA=out/slow_amplitude
mkdir -p "$SA"
./build/slow_amplitude --out "$SA" > "$SA/log_toys.txt" 2>&1 &
for n in 1p2d 1p4d; do
    ./build/slow_amplitude --net "data/atol_bench/${n}_M2.0T_eth.data" --out "$SA" > "$SA/log_$n.txt" 2>&1 &
done

echo "=== network_plot over existing runs"
shopt -s nullglob
for f in out/networks/$TAG/*/*/energies.csv out/networks_mixed/$TAG/*/*/energies.csv; do
    d=$(dirname "$f")
    python scripts/network_plot.py --data "$d" --title "$TAG $(basename "$(dirname "$d")") $(basename "$d")" \
        > /dev/null || echo "PLOT FAILED $d"
done

wait
cat "$SA"/log_*.txt
python scripts/slow_amplitude_plot.py --data "$SA"
