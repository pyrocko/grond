
Grond
=====

A probabilistic earthquake source inversion framework. Designed and crafted in Mordor.

Introduction
------------

Grond is a software tool and framework for robust characterization of
earthquake sources from diverse observations. It employs efficient forward
modelling as well as elaborate inversion strategies and is able to deliver
meaningful model uncertainties.

The software implements a Bayesian bootstrap-based probabilistic joint
inversion scheme for heterogeneous data. It explores the full model space and
helps to map potential model parameter trade-offs. The program is highly
flexible in how it can be adapted to specific dislocation problems, the design
of objective functions, and the diversity of measurements.

Pre-computed Green's function databases handled by the Pyrocko software
library are the backone of the forward modelling. They serve synthetic
near-field geodetic data (InSAR and GNSS) and synthetic seismic waveforms at
all distances and receiver depths for arbitrary earthquake source models.

Contents
--------

.. toctree::
   :maxdepth: 2

   install/index
   quickstart/index
   overview/index
   examples/index
   method/index
   config/index
   report/index
   cli/index
   library/index


Recommended citation
--------------------

TODO DOI

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
