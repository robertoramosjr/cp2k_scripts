#!/bin/bash
# Build the CUTOFF scan (step 0 of azo_au_bridge/README.md) for one junction xyz.
# Usage: make_grid_scan.sh XYZ OUTDIR REL_CUTOFF CUTOFF1 [CUTOFF2 ...]
set -e
xyz=$(readlink -f "${1:?XYZ}"); out=${2:?OUTDIR}; rel=${3:?REL_CUTOFF}; shift 3
for c in "$@"; do
    python3 "$(dirname "$(readlink -f "$0")")/make_junction_input.py" --xyz "$xyz" \
        --project "cut${c}" --outdir "$out/cut${c}" --cutoff "$c" --rel-cutoff "$rel" \
        --run-type ENERGY --scf-probe
done
