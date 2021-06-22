#!/bin/bash

if [ -z "$1" ]; then
    echo "grondown_w_phase.sh <catalog-name>"
    exit 1
fi

bindir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

"$bindir/grondown" \
    --radius-min=3000 \
    --nstations-wanted=40 \
    --padding-factor=10 \
    --network=G,GE,II,IU \
    "$@" 11000 0.001 1.0 "$1"
