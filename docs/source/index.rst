
Grond
=====

A probabilistic earthquake source inversion framework. Designed and crafted in Mordor.

Abstract
--------

Grond is an open source software tool for robust characterization of earthquake sources. Moment tensors and finite fault rupture models can be estimated from seismic, InSAR, and GNSS observations. It serves as a modular framework for the analysis of diverse magmatic, tectonic, and other geophysical processes.

In Grond the optimisation explores the full model space and maps model parameter trade-offs. Meaningful model uncertainties are delivered through a Bayesian bootstrap-based probabilistic joint inversion scheme.
Various data-fitting and weighting options provide easy customization of the objective function.

Rapid forward modelling is enabled by using pre-computed Green's function databases, handled through the Pyrocko software library. They serve synthetic near-field geodetic displacements (InSAR and GNSS) and synthetic seismic waveforms for arbitrary earthquake source models and geometries.


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


Citation
--------

The recommended citation for Grond is:

    Heimann, Sebastian; Isken, Marius; Kühn, Daniela; Sudhaus, Henriette;
    Steinberg, Andreas; Vasyura-Bathke, Hannes; Daout, Simon; Cesca, Simone; Dahm, Torsten (2018): Grond - A probabilistic earthquake source inversion framework. V. 1.0. GFZ Data Services.
    https://doi.org/10.5880/GFZ.2.1.2018.003

Download BibTeX citatation file: :download:`CITATION.bib <../../CITATION.bib>`


License
-------

Grond is licensed under the :download:`GPLv3 <../../LICENSE>`.


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
