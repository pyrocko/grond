
Grond
=====

A probabilistic earthquake source inversion framework. Designed and crafted in Mordor.

Abstract
--------

Grond is an open source software tool for robust characterization of earthquake sources. Moment tensors and finite fault rupture models can be estimated from a combination of seismic waveforms, waveform attributes and geodetic observations like InSAR and GNSS. It helps you to investigate diverse magmatic, tectonic, and other geophysical processes at all scales.

It delivers meaningful model uncertainties through a Bayesian bootstrap-based probabilistic joint inversion scheme. The optimisation explores the full model space and maps model parameter trade-offs with a flexible design of objective functions.

Rapid forward modelling is enabled by using pre-computed Green's function databases, handled through the Pyrocko software library. They serve synthetic near-field surface displacements and synthetic seismic waveforms for arbitrary earthquake source models and geometries.


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


Community and support
---------------------

Please use the `Grond channel
<https://hive.pyrocko.org/pyrocko-support/channels/grond/>`_ in the `Pyrocko Hive
<https://hive.pyrocko.org/>`_ for discussion and feedback. To file a bug
report, please use the `Issues page in our code repository <https://git.pyrocko.org/pyrocko/grond/issues>`_.


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


Source code
-----------

The source code of Grond is hosted on `git.pyrocko.org <https://git.pyrocko.org/pyrocko/grond/>`_.


Contact
-------

* Sebastian Heimann - sebastian.heimann@gfz-potsdam.de
* Marius Isken - marius.isken@gfz-potsdam.de


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
