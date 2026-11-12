#!/usr/bin/env bash
#SBATCH -J q_vs_A
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --mem=30000
#SBATCH -t 12:00:00
#SBATCH -o logs/q_vs_A_%j.log
# dq_dt against dA_dt on networks of real modes (Plans/merged-ladders.md sec. 5), one background
# run per network. Each q run stops at 5 h (slow_amplitude's default), so a network takes <= ~10 h.

set -euo pipefail
source shell-scripts/config.sh
setup_python
mkdir -p build logs
flock build/.compile.lock ./compile.sh -j 4 slow_amplitude
export OMP_NUM_THREADS=1

qa() {   # run folder under out/networks, t_end [s]
    mkdir -p "out/q_vs_A/$1"
    ./build/slow_amplitude --net "out/networks/$1/network.data" --out "out/q_vs_A/$1" --t-end "$2" \
        > "out/q_vs_A/$1/log.txt" 2>&1
    cat "out/q_vs_A/$1/log.txt"
}
qa dsct_M2.0_T7696/ladder_2p/l2_n-2_l1_n+7/d4 3e10 &
qa dsct_M2.0_T7696/ladder_2p/l2_n-3_l2_n-2/d3 3e11 &
qa dsct_M2.2_T7888/ladder_1p/l1_n+3/d4 1e11 &
qa dsct_M2.2_T7888/ladder_1p/l1_n+3/d2 1e10 &
qa dsct_M2.0/ladder_1p/l1_n-1/d6 3e11 &
wait
python scripts/slow_amplitude_plot.py --data out/q_vs_A
