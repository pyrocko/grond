`Grond <https://pyrocko.org/grond/>`_ is an open source software tool for
robust characterization of earthquake sources. Moment tensors and finite fault
rupture models can be estimated from a combination of seismic waveforms,
waveform attributes and geodetic observations like InSAR and GNSS. It helps you
to investigate diverse magmatic, tectonic, and other geophysical processes at
all scales.

It delivers meaningful model uncertainties through a `Bayesian bootstrap-based
probabilistic joint inversion scheme
<https://pyrocko.org/grond/docs/current/method/>`_. The optimisation explores
the full model space and maps model parameter trade-offs with a flexible design
of objective functions.

Rapid forward modelling is enabled by using pre-computed `Green’s function
databases <https://greens-mill.pyrocko.org/>`_, handled through the `Pyrocko
<https://pyrocko.org/docs>`_ software library. They serve synthetic near-field
surface displacements and synthetic seismic waveforms for arbitrary earthquake
source models and geometries.

Installation with pip
---------------------

*See also:* `Grond Manual: Installation
<https://pyrocko.org/grond/docs/current/install>`_

Grond and all its dependencies can be installed by running 

.. code-block:: bash

   pip install grond  # read below...

**but**, we recommend to make a conscious decision about how its main
dependency `Pyrocko <https://pyrocko.org/docs>`_ and especially Pyrocko's own
dependencies are installed. The `Pyrocko Installation Manual
<https://pyrocko.org/docs/current/install/>`_ describes different installation
schemes.

As a general advice, we recommend to exclusively use either, (1) the system's
native package manager, (2) Anaconda, or (3) pip only. In (1) and (2), only
resort to use pip for those few packages which are not available as native
packages. Otherwise, competing package managers will ruin your day!

To prevent pip from automatically resolving dependencies run

.. code-block:: bash

   pip install --no-deps grond

This assumes that `Pyrocko <https://pyrocko.org/docs>`_ and `Kite
<https://pyrocko.org/kite/>`_ have been installed beforehand.

Documentation
--------------

Documentation and examples can be found in the `Grond Manual
<https://pyrocko.org/grond/>`_.

Community
---------

Meet people from all over the world doing awesome research with Grond in our
community chat: use the *Grond* channel in the friendly `Pyrocko Hive
<https://hive.pyrocko.org>`_. This is the best place to talk about new features,
special techniques or to get help on setting up your first inversion with
Grond.

Development
-----------

Grond is open source.

Join us at our `Git repository <https://git.pyrocko.org/pyrocko/grond/>`_ and
read the `Contribution guide
<https://git.pyrocko.org/pyrocko/grond/src/branch/master/CONTRIBUTING.md>`_.

-- The Grond Developers
