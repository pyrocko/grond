Problems and problem configuration
==================================

(unfinished)

A problem is a task to optimise a specific source in Grond given certain conditions that can be configured
rather flexibly. So far defined problems in Grond are the optimisations of different 
source types, which are derived from `Pyrocko Sources`_. 

These problems are in short (more details in the corresponding sections below):

* ``CMTProblem``: 
    A problem that solves for a centroid moment tensor point source (derived from ``CMTSource``). This
    problem fits the very general earthquake source analysis based on far-field seismic waveforms.

* ``DoubleDCProblem``: 
    A problem that solves for two double-couple point sources (derived from ``DoubleDCSource``). This 
    problem can be used to solve for somewhat complex, segmented earthquake sources to better fit 
    far-field seismic data. 
    
* ``RectangularProblem``:
    A problem that solves for a rectangular finite source (derived from ``RectangularSource``). This
    problem fits well to large earthquakes and/or problems for which near-field data (InSAR, GNSS, etc.)
    are available.

To define and configure a problem the part called ``problem_config`` in the configuration is set up.


General configuration
---------------------

The general problem configuration contains the following parameters:

``name_template``:
    can be any string and provides a stem for the result folders 
    `runs` and `reports` to identify different optimisations. 
    Meaningful is to use short event and problem identifications 
    in this string.

``norm_exponent``:
    defines the norm of combining several `normalization_family` in the
    global misfit. This integer value is 1 or larger. Please find here
    more information on the global `misfit calculation in Grond`_.   

``ranges``
    defines the bounds of individual and specific source model parameters.
    See the details for the source ranges of different problems in the sections below. 
 
An example for the configuration of a rectangular fault problem is given here:

.. code-block :: sh

  problem_config: !grond.RectangularProblemConfig
    name_template: '${event_name}'
    norm_exponent: 0
    decimation_factor: 4
    ranges:
      north_shift: '-2000 .. 20000'
      east_shift: '-2000 .. 20000'
      depth: '5000 .. 30000'
      length: '12000 .. 18000'
      width: '4000 .. 14000'
      slip: '0.2 .. 2.'
      strike: '80 .. 330'
      dip: '0 .. 60'
      rake: '60 .. 90'

CMTProblem source configuration
-------------------------------

description: general earthquake source problem.  

- full moment tensor

Here are non-general parameters and their description for the ``RectangularProblem`` in particular
with a concise example at the end.

**Non-general problem parameters**:

``distance_min``: is a minimum target distance to the source used to exclude targets closer than this. 
    Tailored to the problem, too close targets will not be considered in the misfit evaluation. 
    Finite-rupture effects on near targets may be excluded efficiently with a meaning setting 
    for this parameter. 
    
``distance_min``: is a maximum target distance to the source used to exclude targets farther than this. 
    Tailored to the problem, too far away targets will not be considered in the misfit evaluation. 
    Like this certain phase interferences may be efficiently excluded. 

``mt_type``: configures the type of moment tensor. The source can be a full moment tensor (`` mtfull`, or a 'deviatoric'

**Source parameters**:

``ranges``:
    ``east_shift`` is a relative position in east direction to the reference location 
        given in 'event.txt'. It is given in meters.

    ``north_shift`` is a relative position in north direction to the reference location
        given in 'event.txt'. It is given in meters.

    ``depth``: is the depth of the point source in meters.

    ``time``: is the relative time to the origin time given in 'event.txt' in seconds.

    ``magnitude``: is the earthquake moment magnitude.



.. code-block :: sh

  problem_config: !grond.CMTProblemConfig
    name_template: '${event_name}_regional_mt'
    norm_exponent: 1
    distance_min: 0
    mt_type: 'deviatoric'
    ranges:
      time: '-10 .. 10 | add'
      north_shift: '-40e3 .. 40e3'
      east_shift: '-40e3 .. 40e3'
      depth: '4e3 .. 50e3'
      magnitude: '4.0 .. 7.0'
      rmnn: '-1.41421 .. 1.41421'
      rmee: '-1.41421 .. 1.41421'
      rmdd: '-1.41421 .. 1.41421'
      rmne: '-1 .. 1'
      rmnd: '-1 .. 1'
      rmed: '-1 .. 1'
      duration: '0. .. 0.'




DoubleDCProblem configuration
-----------------------------

TODO: alles, problem description

Here are non-general parameters and their description for the ``DoubleDCProblem`` 
in particular with a concise example at the end.

``DoubleDC`` stands for two double-couple sources. This problem 


**Non-general problem parameters**:

``distance_min``: 0

``mt_type``: 'deviatoric'


TODO: example

RectangularProblem configuration
--------------------------------

TODO: problem description

Here are non-general parameters and their description for 
the ``RectangularProblem`` in particular with a concise 
example at the end.


**Non-general problem parameters**:

``decimation_factor``
    is only valid for finite sources. It defines a reduced number
    of sub-sources that build the finite source rectangle. A reduced
    number speeds up the forward modelling but may lead to artefacts 
    in the source near-field. Default is no decimation 
    (``decimation_factor: 1``)

**Source parameters**:

For the source parameter configuration, please note that the 
last three parameters
``nucleation_x``, ``nucleation_y`` and ``time`` are needed
to define the Rectangular Source for the forward modelling
of seismic waveforms. If they are missing waveform targets are
ignored in the optimisation. If only static targets are
defined, the source parameters for the nucleation point and
origin time, if given,
are ignored.


``ranges``:
    ``east_shift`` is a relative position in east direction to the reference location 
        given in 'event.txt'. It is given in meters.

    ``north_shift`` is a relative position in north direction to the reference location
        given in 'event.txt'. It is given in meters.

    ``depth``: is the depth of upper fault edge (not centroid!) in meters.

    ``length``: is the along-strike length of the fault in meters.

    ``width``: is the along-dip width of the fault in meters.

    ``strike``: is the strike angle of fault against north in degrees.

    ``dip``: is the dip angle of fault against horizontal in degrees.

    ``rake``: is the rake angle of slip in degrees.

    ``time``: is the relative time to the origin time given in 'event.txt' in seconds.

    ``nucleation_x``: relative horizontal position of the rupture nucleation point 
        on the fault to the centre location. This parameter may range from -1 to 1.
        With 0 being in the centre, -1 being at the left-side fault edge, 1 at the 
        right-side fault edge, and 0.5 is half-way between centroid and right-side 
        fault edge.

    ``nucleation_y``: relative along-dip position of the rupture nucleation point
        on the fault to the centre location. This parameter may range from -1 to 1. 
        With 0 being in the centre, -1 being at the top fault edge, 1 at the bottom 
        fault edge, and 0.5 is half-way between centroid and bottom fault edge.

**Example configuration**:

.. code-block :: sh

  problem_config: !grond.RectangularProblemConfig
    name_template: '${event_name}_static'
    norm_exponent: 1
    decimation_factor: 4
    ranges:
      north_shift: '-2000 .. 20000'
      east_shift: '-2000 .. 20000'
      depth: '5000 .. 30000'
      length: '12000 .. 18000'
      width: '4000 .. 14000'
      slip: '0.2 .. 2.'
      strike: '80 .. 330'
      dip: '0 .. 60'
      rake: '60 .. 90'
      time: '-15. .. 10. | add'
      nucleation_x: '-1. .. 1.'
      nucleation_y: '-1. .. 1.'


.. _misfit calculation in Grond: ../method/index.html#Misfit calculation
.. _Pyrocko Sources: _https://pyrocko.org/docs/current/library/reference/gf.html#module-pyrocko.gf.seismosizer
