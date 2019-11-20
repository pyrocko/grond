#!/bin/bash
set -e
EVENT_PATH="data/events/2019-ridgecrest"
INSAR_PATH="$EVENT_PATH/insar/"
GNSS_PATH="$EVENT_PATH/gnss/"

mkdir -p $EVENT_PATH
mkdir -p $INSAR_PATH
mkdir -p $GNSS_PATH

# wget http://data.pyrocko.org/examples/2019-Ridgecrest-insar_gnss/README -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2019-Ridgecrest-insar_gnss/insar/ascending.npz -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2019-Ridgecrest-insar_gnss/insar/ascending.yml -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2019-Ridgecrest-insar_gnss/insar/descending.npz -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2019-Ridgecrest-insar_gnss/insar/descending.yml -P $INSAR_PATH
wget http://data.pyrocko.org/examples/2019-Ridgecrest-insar_gnss/gnss/2019-ridgecrest-gnss.yml -P $GNSS_PATH
wget http://data.pyrocko.org/examples/2019-Ridgecrest-insar_gnss/event.txt -P $EVENT_PATH

