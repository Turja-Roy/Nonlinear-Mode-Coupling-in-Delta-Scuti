#!/usr/bin/env bash
# Build the C++ coupling pipeline.
#
#   ./compile.sh                 build everything into build/
#   ./compile.sh -j 16           more parallel jobs (default 4)
#   ./compile.sh --test          build, then run the unit tests
#   ./compile.sh --clean         remove build/ first
#   ./compile.sh three_mode_search        build one target only
#
# CXXFLAGS overrides the optimisation flags (default -O2). Do not add
# -march=native for a cluster build unless you compile on the same node type
# you run on -- a login node is often a different microarchitecture and the
# binary dies with an illegal instruction on the compute nodes.
#
# Header-only dependencies are vendored under third_party/ and are used in
# preference to anything installed, so the build does not depend on the
# cluster's Eigen or Boost version. GSL and HDF5 are compiled libraries and
# cannot be vendored as headers; the script looks for them and says what to
# load if they are missing.

set -euo pipefail
cd "$(dirname "$0")"

JOBS=4
TARGETS=()
RUN_TESTS=0
CLEAN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -j)        JOBS="$2"; shift 2 ;;
        -j*)       JOBS="${1#-j}"; shift ;;
        --test)    RUN_TESTS=1; shift ;;
        --clean)   CLEAN=1; shift ;;
        -h|--help) sed -n '2,16p' "$0" | sed 's/^# \?//'; exit 0 ;;
        *)         TARGETS+=("$1"); shift ;;
    esac
done

CXX=${CXX:-g++}
BUILD=build
[[ $CLEAN -eq 1 ]] && rm -rf "$BUILD"
mkdir -p "$BUILD/obj"

# ---------------------------------------------------------------- dependencies
say() { printf '%s\n' "$*"; }
die() { printf 'error: %s\n' "$*" >&2; exit 1; }

need_pkg() {   # name, pkg-config name, hint
    if pkg-config --exists "$2" 2>/dev/null; then
        printf "  %-9s %s\n" "$1:" "$(pkg-config --modversion "$2")"
        return 0
    fi
    return 1
}

say "dependencies:"
say "  Eigen:    third_party/eigen3 (vendored)"
say "  Boost:    third_party/boost (vendored, odeint subset)"
say "  CLI11:    third_party/CLI (vendored)"
say "  HighFive: third_party/highfive (vendored)"

GSL_C=""; GSL_L="-lgsl -lgslcblas"
if need_pkg "GSL" gsl; then
    GSL_C=$(pkg-config --cflags gsl); GSL_L=$(pkg-config --libs gsl)
elif command -v gsl-config >/dev/null 2>&1; then
    GSL_C=$(gsl-config --cflags); GSL_L=$(gsl-config --libs)
    say "  GSL:      $(gsl-config --version) (via gsl-config)"
else
    die "GSL not found. On a cluster try:  module load gsl
       Only gsl_sf_coupling_3j is used, in src/angular.cpp."
fi

HDF5_C=""; HDF5_L="-lhdf5"; HDF5_FOUND=0
for p in hdf5 hdf5-serial; do
    if pkg-config --exists "$p" 2>/dev/null; then
        # --cflags is legitimately empty when the headers are already on the
        # default include path, so track the hit separately.
        HDF5_C=$(pkg-config --cflags "$p"); HDF5_L=$(pkg-config --libs "$p")
        HDF5_FOUND=1
        say "  HDF5:     $(pkg-config --modversion "$p")"
        break
    fi
done
if [[ $HDF5_FOUND -eq 0 ]]; then
    if [[ -n ${HDF5_DIR:-} && -d $HDF5_DIR/include ]]; then
        HDF5_C="-I$HDF5_DIR/include"; HDF5_L="-L$HDF5_DIR/lib -lhdf5"
        say "  HDF5:     \$HDF5_DIR=$HDF5_DIR"
    elif [[ -f /usr/include/hdf5.h ]]; then
        say "  HDF5:     /usr/include (no pkg-config)"
    else
        die "HDF5 not found. On a cluster try:  module load hdf5
       or set HDF5_DIR to its install prefix."
    fi
fi

OMP=""
if echo 'int main(){}' | $CXX -fopenmp -x c++ - -o /dev/null 2>/dev/null; then
    OMP="-fopenmp"
    say "  OpenMP:   enabled"
else
    say "  OpenMP:   not available -- --jobs will run serially"
fi

# ---------------------------------------------------------------------- build
INC="-Isrc -Ithird_party -Ithird_party/eigen3 $GSL_C $HDF5_C"
CXXFLAGS=${CXXFLAGS:--O2}   # -march=native is opt-in: see the note in --help
FLAGS="-std=c++17 $CXXFLAGS -Wall -Wextra $OMP $INC"
LIBS="$GSL_L $HDF5_L $OMP -lm"

LIB_SRC=(numeric angular model kappa triplets amplitude stability csv)
APPS=(three_mode_search channels four_mode_search mixed_network make_inlists mw25)
TESTS=(numeric kappa amplitude)

compile_one() {   # src, obj, extra
    local s=$1 o=$2
    if [[ -f $o && $o -nt $s ]]; then return 0; fi
    $CXX $FLAGS ${3:-} -c "$s" -o "$o" || { printf 'FAILED: %s\n' "$s" >&2; exit 1; }
    printf '  cc %s\n' "$s"
}

say ""
say "compiling (-j $JOBS, $CXXFLAGS):"
pids=()
for m in "${LIB_SRC[@]}"; do
    compile_one "src/$m.cpp" "$BUILD/obj/$m.o" &
    pids+=($!)
    while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n; done
done
for p in "${pids[@]}"; do wait "$p" || die "compile failed"; done

# Rebuild the archive only when an object is newer than it. `ar` otherwise
# rewrites it every run, which makes it newer than every binary and forces a
# full re-link even when nothing changed.
LIB="$BUILD/libcoupling.a"
stale=0
[[ -f $LIB ]] || stale=1
for m in "${LIB_SRC[@]}"; do [[ $BUILD/obj/$m.o -nt $LIB ]] && stale=1; done
if [[ $stale -eq 1 ]]; then
    ar rcs "$LIB" "$BUILD"/obj/*.o
    say "  ar libcoupling.a"
fi

# Apps and tests compile to their own objects, then link. Doing both in one
# step meant that touching any library source recompiled every app from
# scratch -- 150 s where the link alone is a couple.
link_one() {   # name, source, extra-flags
    local out="$BUILD/$1" obj="$BUILD/obj/$1.o"
    if [[ ! -f $obj || $2 -nt $obj ]]; then
        $CXX $FLAGS ${3:-} -c "$2" -o "$obj" \
            || { printf 'FAILED to compile %s\n' "$2" >&2; exit 1; }
        printf '  cc %s\n' "$2"
    fi
    if [[ -f $out && $out -nt $obj && $out -nt "$LIB" ]]; then return 0; fi
    $CXX $FLAGS "$obj" "$LIB" $LIBS -o "$out" \
        || { printf 'FAILED to link %s\n' "$1" >&2; exit 1; }
    printf '  ld %s\n' "$1"
}

want() {   # is this target requested?
    [[ ${#TARGETS[@]} -eq 0 ]] && return 0
    local w; for w in "${TARGETS[@]}"; do [[ $w == "$1" ]] && return 0; done
    return 1
}

for a in "${APPS[@]}";  do want "$a" && link_one "$a" "apps/$a.cpp"; done
# -UNDEBUG: the tests are assert()s, and -O2 builds would otherwise define
# NDEBUG and compile every check out, so they would pass without testing.
for x in "${TESTS[@]}"; do
    want "test_$x" && link_one "test_$x" "tests/test_$x.cpp" "-UNDEBUG"
done

say ""
say "-> $BUILD/"

if [[ $RUN_TESTS -eq 1 ]]; then
    say ""
    for x in "${TESTS[@]}"; do "$BUILD/test_$x"; done
fi
