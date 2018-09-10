#!/bin/bash

if [ -z "$1" ]; then
    echo "grondown_regional.sh <catalog-name>"
    exit 1
fi

bindir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

"$bindir/grondown" \
    --nstations-wanted=40 \
    --padding-factor=10 \
    "$@" 1000 .01 2. "$1"
