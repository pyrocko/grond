Problems
========

The problem is the optimisation task - the specific :term:`source` model we invert.
The search ranges and inversion conditions can be configured flexibly. So far defined problems in Grond are the optimisations of different source types, which are derived from `Pyrocko Sources`_.

The optimisation problems which can be solved with Grond are:

* :ref:`CMTProblem <cmt>`
    A problem that solves for a centroid moment tensor point source (derived from Pyrocko's ``CMTSource``). This problem fits the very general earthquake source analysis based on far-field seismic waveforms.

* :ref:`DoubleDCProblem <double_dc>`
    A problem that solves for two double-couple point sources (derived from ``DoubleDCSource``). This problem can be used to solve for somewhat complex, *segmented earthquake sources* to better fit far-field seismic data.
    
* :ref:`RectangularProblem <rectangular>`
    A problem that solves for a rectangular finite source (derived from Pyrocko's ``RectangularSource``). This problem fits well to large earthquakes and/or problems for which near-field surface displacement data (InSAR, GNSS, etc.) are available.

* :ref:`VolumePointProblem <volume_point>`
    A problem that solves for a spherical volume point (infinite) to model magmatic or volcanic processes. Only static targets (GNSS or InSAR) are supported.

To define and configure a problem the part called ``problem_config`` in the configuration is set up.


General configuration
---------------------

.. code-block :: yaml
  :caption: Generic example ``ProblemConfig``

  problem_config: !grond.ProblemConfig

    # Name used to identify the output
    name_template: '${event_name}'

    # How to combine the target misfits. For L1 norm: 1, L2 norm: 2, etc.
    norm_exponent: 2

    # Define the ranges of the solution space.
    ranges:
      ..
      # Problem specific, see below ...


The following problem parameters are shared by all problems and are part of all
problem configurations:

.. glossary::

  ``name_template``
    can be any string and provides a stem for the result folders `runs` and `report` to identify different optimisations. Meaningful is to use short event and problem identifications in this string.

  ``norm_exponent``
    defines the norm of combining several `normalization_family` in the global misfit. This integer value is 1 or larger. Please find here more information on the global `misfit calculation in Grond`_.

  ``ranges``
    defines the bounds of individual and specific source model parameters. See the details for the source ranges of different problems in the sections below.

An example for the configuration of a rectangular fault problem is given here:


.. _cmt:

``CMTProblem`` configuration
----------------------------

.. literalinclude :: ../../../../src/data/snippets/problem_cmt.gronf
  :language: yaml
  :caption: Example ``CMTProblemConfig``

The Grond CMTProblem represents one of most popular problems in seismology.  
Sought-after are moment tensor source models that well fit the observed seismic waveforms. The 
waveforms are usually far-field observations such that this point-source approximation
is well suited. The ``CMTProblem`` can be configured to a pure double-couple problem. 
The source time function is fixed to half-sinusoid (see the documentation of 
`Pyrocko Sources`_ for details on the source time function).

Configuration parameters that common for all problems are listed above and following are parameters 
particularly for the ``CMTProblem``.

**Non-general problem parameters**:


.. glossary::

  ``distance_min``
    is a minimum target distance to the source used to exclude targets closer than this. Tailored to the problem, too close targets will not be considered in the misfit evaluation. Finite-rupture effects on near targets may be excluded efficiently with a meaning setting
    for this parameter. 
    
  ``distance_max``
    is a maximum target distance to the source used to exclude targets farther than this. Tailored to the problem, too far away targets will not be considered in the misfit evaluation. Like this certain phase interferences may be efficiently excluded.

  ``mt_type``
    configures the type of moment tensor. The source model can be set to be a full moment tensor (``mt_type: full``) or can be constrained to a deviatoric moment tensor (``mt_type: deviatoric``) or even to a pure double couple source (``mt_type: dc``).

**Source parameters**:

(Please check for more details the description of the `Pyrocko Sources`_.)


.. rubric :: ``ranges``:

.. glossary::

    ``east_shift``
      is a relative position in east direction to the reference location given in :file:`event.txt`. It is given in meters.

    ``north_shift``
      is a relative position in north direction to the reference location given in 'event.txt'. It is given in meters.

    ``depth``
      is the depth of the point source in meters.

    ``time``
      is the relative time to the origin time given in 'event.txt' in seconds.

    ``magnitude``
      is the earthquake moment magnitude.
    
    ``rmnn`` & ``rmee`` & ``rmdd`` & ``rmne`` & ``rmnd`` & ``rmed``
      are the moment tensor components.

    ``duration``
      is the duration of the source time function in seconds.

.. _double_dc:

``DoubleDCProblem`` configuration
---------------------------------

.. literalinclude :: ../../../../src/data/snippets/problem_double_dc.gronf
  :language: yaml
  :caption: Example ``DoubleDCProblemConfig``


This problem has two double-couple point sources (derived from ``DoubleDCSource``). They are dependent in location and relative timing to avoid overlapping in either space or time. The mechanisms, the durations and the moments of the two sources are independent. Using this model more complex earthquakes with two prominent rupture phases or with a change of mechanism along the rupture plane can be studied. Or simply the potential of a major source
complexity of an earthquake can be tested.

Configuration parameters that common for all problems are listed above and following are parameters 
particularly for the ``DoubleDCProblem``.

**Non-general problem parameters**:

``distance_min``
  is a minimum target distance to the source used to exclude targets closer than this. Tailored to the problem, too close targets will not be considered in the misfit evaluation. Finite-rupture effects on near targets may be excluded efficiently with a meaning setting for this parameter.
    
    
**Source parameters**:

(Please check for more details the description of the `Pyrocko Sources`_.)

.. rubric :: ``ranges``:

.. glossary ::
    ``east_shift``
      is a relative position in east direction to the reference location given in 'event.txt'. It is given in meters.

    ``north_shift``
      is a relative position in north direction to the reference location given in 'event.txt'. It is given in meters.

    ``depth``
      is the depth of the starting point source in meters.
    
    ``time``
      is the relative time to the origin time given in 'event.txt' in seconds.

    ``magnitude``
      is the total earthquake moment magnitude.

    ``strike1`` & ``dip1`` & ``rake1``
      constrain the mechanism of the first double-couple source.
    
    ``strike2`` & ``dip2`` & ``rake2``
      constrain the mechanism of the second double-couple source.
    
    ``delta_time``
      is the time difference between the two sources in seconds. Needs to be larger than zero to separate the sources in time and to make source 2 the later source.
        
    ``delta_depth``
      is the depth difference of the two sources in meters.
    
    ``azimuth``
      the azimuth of source 2 with respect to source 1 (clockwise from north) in degrees.
    
    ``distance``
      is the distance between the two sources in meters. Needs to be larger than zero to separate the sources in space.
    
    ``mix``
      is a value between ``0`` and ``1`` that defines the relative moment contributions of the sources to the total moment. In the extreme, with ``mix=0`` all the moment is in the first source and none in the second or else ``mix=1`` put all moment in the second source which leaves none for the first source. ``mix=0.25`` defines three quarters of the total moment on the first source and one quarter on the second, while obviously ``mix=0.5`` gives two sources of the same strength.
        
    ``duration1`` & ``duration2``
      are the durations of the first and second source's source time functions, respectively, in seconds.


.. _rectangular:

``RectangularProblem`` configuration
------------------------------------

.. literalinclude :: ../../../../src/data/snippets/problem_rectangular.gronf
  :language: yaml
  :caption: Example ``RectangularProblemConfig``

The rectangular source is a simple finite source model with a rectangular shape and uniform moment or slip across the rupture plane. It resembles the source model defined by `Haskell (1964)`_, but has a nucleation point from which spreads a circular rupture. The position of the nucleation point on the rupture plane can be part of the problem. Uniform and bilateral ruptures are therefore possible. With the ``RectangularProblem`` also directivity effects in the observations of large earthquake may be predicted.

The static rectangular source is very similar to the analytical rectangular dislocation source as described by `Okada (1985)`_, which is embedded in an isotropic elastic half-space. The ``RectangularProblem`` is therefore well suited to predict near-field static surface displacements observed at GNSS stations or with InSAR. For a joint optimisation of seismic waveforms and near-field static surface displacements a ``RectangularProblem`` is the appropriate choice.

Configuration parameters that common for all problems are listed above and following are parameters particularly for the ``RectangularProblem``.

**Non-general problem parameters**:

``decimation_factor``
  is only valid for finite sources. It defines a reduced number of sub-sources that build the finite source rectangle. A reduced number speeds up the forward modelling but may lead to artefacts in the source near-field. Default is no decimation (``decimation_factor: 1``)

**Source parameters**:

(Please check for more details the description of the `Pyrocko Sources`_.)

For the source parameter configuration, please note that the last three parameters ``nucleation_x``, ``nucleation_y`` and ``time`` are needed to define the Rectangular Source for the forward modelling of seismic waveforms. If they are missing waveform targets are ignored in the optimisation. If only static targets are defined, the source parameters for the nucleation point and origin time, if given, are ignored.


``ranges``:

.. glossary ::
    ``east_shift``
      is a relative position in east direction to the reference location given in :file:`event.txt`. It is given in meters.

    ``north_shift``
      is a relative position in north direction to the reference location given in 'event.txt'. It is given in meters.

    ``depth``
      is the depth of upper fault edge (not centroid!) in meters.

    ``length``
      is the along-strike length of the fault in meters.

    ``width``
      is the along-dip width of the fault in meters.

    ``strike``
      is the strike angle of fault against north in degrees.

    ``dip``
      is the dip angle of fault against horizontal in degrees.

    ``rake``
      is the rake angle of slip in degrees.

    ``time``
      is the relative time to the origin time given in 'event.txt' in seconds.

    ``nucleation_x``
      relative horizontal position of the rupture nucleation point on the fault to the centre location. This parameter may range from -1 to 1. With 0 being in the centre, -1 being at the left-side fault edge, 1 at the right-side fault edge, and 0.5 is half-way between centroid and right-side fault edge.

    ``nucleation_y``
      relative along-dip position of the rupture nucleation point on the fault to the centre location. This parameter may range from -1 to 1. With 0 being in the centre, -1 being at the top fault edge, 1 at the bottom fault edge, and 0.5 is half-way between centroid and bottom fault edge.


.. _volume_point:

``VolumePointProblem`` configuration
------------------------------------

.. literalinclude :: ../../../../src/data/snippets/problem_volume_point.gronf
  :language: yaml
  :caption: Example ``VolumePointProblemConfig``


The volume point problem is a simple dislocation source with a spherical shape, as a model for magmatic or volcanic processes.

This problem can be used for static targets only (InSAR or GNSS), and is similar to the Mogi model, with the difference that the ``VolumePointProblem`` is using Pyrocko's Green's function databases and can act in layered media.

**Source parameters**:

``ranges:``

  ``volume_change``
    Allowed volume change of the point source, given in :math:`\mathrm{m}^3`. Positive for inflational volumes, negative values for deflation.

.. _misfit calculation in Grond: ../method/index.html#Misfit calculation
.. _Pyrocko Sources: _https://pyrocko.org/docs/current/library/reference/gf.html#module-pyrocko.gf.seismosizer
.. _Haskell (1964): https://pubs.geoscienceworld.org/ssa/bssa/article/54/6A/1811/116295/total-energy-and-energy-spectral-density-of 
.. _Okada (1985): https://pubs.geoscienceworld.org/ssa/bssa/article/75/4/1135/118782/surface-deformation-due-to-shear-and-tensile
