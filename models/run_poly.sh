#!/bin/sh
# Polytrope harness (M0/M1): an independent structure for the kappa machinery,
# with no MESA input at all. poly_to_fgong scales to 1 M_sun / 1 R_sun in cgs,
# so coupling/ reads these exactly as it reads the delta Sct model.
#
# --drop-outer: build_poly integrates to theta_s = 0, where rho and P vanish
# and V diverges. FGONG cannot carry that point.
set -e
cd "$(dirname "$0")"

for n in 3 0; do
  d="poly_n${n}/gyre"
  echo "=== n = ${n} ==="
  mkdir -p "${d}/detail"
  (
    cd "${d}"
    "$GYRE_DIR/bin/build_poly" "poly_n${n}.h5" -n "${n}" -G 1.6666666666666667 -d 0.005
    "$GYRE_DIR/bin/poly_to_fgong" "poly_n${n}.h5" "poly_n${n}.fgong" --drop-outer
    "$GYRE_DIR/bin/gyre" gyre_ad.in
  )
done

echo "done: $(ls poly_n3/gyre/detail | wc -l) + $(ls poly_n0/gyre/detail | wc -l) detail files"
