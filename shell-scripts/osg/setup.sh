#!/usr/bin/env bash
# Once, on the OSPool access point (ap40): what network_build/network_run need there and the access point
# lacks. SUNDIALS (its source in $S/sundials-7.9.0) and GSL, built static, so network_run carries them to
# the execute points; HDF5 comes from the MESA SDK. ~/miniforge -> the data area's miniforge, for the
# Python of shell-scripts/config.sh.
#   bash shell-scripts/osg/setup.sh
set -euo pipefail
S=/ospool/ap40/data/$USER/soft
[[ -e ~/miniforge ]] || ln -s "/ospool/ap40/data/$USER/miniforge" ~/miniforge

cmake -S "$S/sundials-7.9.0" -B "$S/sundials-7.9.0/build" -DCMAKE_INSTALL_PREFIX="$S/sundials" \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON \
    -DSUNDIALS_ENABLE_OPENMP=ON -DENABLE_OPENMP=ON \
    -DSUNDIALS_ENABLE_C_EXAMPLES=OFF -DSUNDIALS_ENABLE_CXX_EXAMPLES=OFF -DSUNDIALS_ENABLE_EXAMPLES_INSTALL=OFF \
    -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_CXX=OFF -DEXAMPLES_INSTALL=OFF
cmake --build "$S/sundials-7.9.0/build" -j 4
cmake --install "$S/sundials-7.9.0/build"

cd "$S/sources"
[[ -f gsl-2.8.tar.gz ]] || wget -nv https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz
rm -rf gsl-2.8 && tar xf gsl-2.8.tar.gz
(cd gsl-2.8 && ./configure -q --prefix="$S/gsl" --disable-shared && make -j 4 > /dev/null && make install > /dev/null)
rm -rf gsl-2.8
ls "$S/sundials/lib"*/libsundials_nvecopenmp.a "$S/gsl/lib/libgsl.a"
