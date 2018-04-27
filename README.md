# Grond

## Introduction

Grond explores characteristics of earthquake and deformation sources. It is a
framework to search model spaces in non-linear inversion problems. A variety of
source models and geophysical observations can be incorporated. Grond is not
limited to a specific optimisation algorithm.  The default optimizer in Grond
is a global Monte Carlo sampler, capable of retrieving multimodal and
irregularly shaped solution spaces. It is a multi-objective function optimizer
and can efficiently apply the Bayesian bootstrap to estimate earthquake source
model uncertainties. As foundation the [Pyrocko](https://pyrocko.org) toolbox
is used.

Enabled source types so far are centroid moment tensors, double-couple
point-sources and finite rectangular sources. Enabled data sets to be used for
the inverse modelling are seismic waveforms, InSAR coseismic near-field
displacement maps and GNSS coseismic displacement vectors. These data sets can
used also combined. For seismic waveforms standard data formats can be used.
For using coseismic InSAR displacement maps, a [Kite] InSAR data container
needs to be prepared. [Kite](https://github.com/pyrocko/kite) is an interactive
InSAR-data post-processor, which imports a lot of standard InSAR processor
output formats. For using GNSS campaign coseismic displacement vectors, a
simple data file needs to be prepared, containing the GNSS station positions,
displacement values and measurement uncertainties.

For the forward modelling precalculated and stored Green's functions are used.
These stores can be obtained from an online repository or are user-defined and
user-calculated using Pyrocko's fomosto module.


## Installation

First, install [Pyrocko](http://pyrocko.org/docs/current/install/),
then install Grond:

```bash
git clone https://gitext.gfz-potsdam.de/heimann/grond.git
cd grond
sudo python setup.py install
```

For support of InSAR data modelling install also
[Kite](https://github.com/pyrocko/kite)


## Updating an existing installation

```bash
cd grond  # change to the directory to where you cloned grond initially
git pull origin master
sudo python setup.py install
```

## License

GNU General Public License, Version 3, 29 June 2007

Copyright © 2018 Helmholtz Centre Potsdam GFZ German Research Centre for
Geosciences, Potsdam, Germany and University of Kiel, Kiel, Germany.

Grond is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. Grond is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.


# Documentation

Find the documentation at https://pyrocko.github.io/grond/.


# Community and support

Please use the [Pyrocko Hive](https://hive.pyrocko.org) for discussion,
feedback and bug reports.


