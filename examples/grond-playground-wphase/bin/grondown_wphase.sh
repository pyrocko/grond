#!/bin/bash

if [ -z "$1" ]; then
    echo "grondown_w_phase.sh <catalog-name>"
    exit 1
fi
bindir="$(dirname "$(realpath "$0")")"

"$bindir/grondown" \
    --radius-min=3000 \
    --nstations-wanted=40 \
    --padding-factor=10 \
    "$@" 11000 0.001 1.0 "$1"
