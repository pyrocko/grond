#!/bin/bash

set -e
DATA_PATH="data/events/2009laquila/insar/"
mkdir -p $DATA_PATH

wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/README -P $DATA_PATH
wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/asc_insar.npz -P $DATA_PATH
wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/asc_insar.yml -P $DATA_PATH
wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/dsc_insar.npz -P $DATA_PATH
wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/dsc_insar.yml -P $DATA_PATH
