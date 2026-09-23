#!/usr/bin/env bash

# Only needs GSL
set -euo pipefail
cd "$(dirname "$0")"

GSL_C=$(pkg-config --cflags gsl 2>/dev/null || gsl-config --cflags)
GSL_L=$(pkg-config --libs gsl 2>/dev/null || gsl-config --libs)

g++ -std=c++17 -O2 -Wall -Wextra $GSL_C mw25_simplified.cpp -o mw25_simplified $GSL_L
