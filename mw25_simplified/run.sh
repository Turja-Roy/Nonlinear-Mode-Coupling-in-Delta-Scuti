#!/usr/bin/env bash

set -euo pipefail
cd "$(dirname "$0")"

OUT=${1:-out/}

./build.sh
./mw25_simplified "$OUT"
python3 plot.py --data "$OUT"
