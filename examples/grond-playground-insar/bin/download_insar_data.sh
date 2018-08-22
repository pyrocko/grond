#!/bin/bash

set -e
EVENT_PATH="data/events/2009laquila"
INSAR_PATH="$EVENT_PATH/insar/"
mkdir -p $EVENT_PATH
mkdir -p $INSAR_PATH

wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/README -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/asc_insar.npz -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/asc_insar.yml -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/dsc_insar.npz -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/dsc_insar.yml -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/dsc_insar.yml -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2009-LAquila-InSAR/event.txt -P $EVENT_PATH
