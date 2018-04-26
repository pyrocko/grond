# Grond

## Introduction

is a battering ram to non-linearly explore earthquake source parameters based on geophysical data sets. Grond uses as foundation the pyrocko toolbox. 

Enabled source types so far are centroid moment tensors, double-couple point-sources and finite
rectangular sources. Enabled data sets to be used for the inverse modelling are seismic waveforms, InSAR coseismic near-field displacement maps and GNSS coseismic displacement vectors. These data sets can used also combined. 
For seismic waveforms tne standard data formats can be used. 
For using coseismic InSAR displacement maps, a [Kite] InSAR data container needs to be prepared. [Kite] is an interactive InSAR-data post-processor, which imports a lot of standard InSAR processor output formats. 
For using GNSS campaign coseismic displacement vectors, a simple data file needs to be prepared, containing the GNSS station positions, displacement values and measurement uncertainties.

For the forward modelling precalculated and stored Green's functions are used. These stores can be obtained from an online repository or are user-defined and user-calculated using pyrockos fomosto module.

## Installation

First, install [Pyrocko](http://pyrocko.org/docs/current/install/),
then install Grond:

```bash
git clone https://gitext.gfz-potsdam.de/heimann/grond.git
cd grond
sudo python setup.py install [Pyrocko](http://pyrocko.org/docs/current/install/),
```

For support of InSAR data modelling install also [Kite](https://github.com/pyrocko/kite)

## Updating an existing installation

```bash
cd grond  # change to the directory to where you cloned grond initially
git pull origin master
sudo python setup.py install
```

# Documentation

Find the documentation at https://pyrocko.github.io/grond/.


