Targets
=======

TODO: describe the different target types available in Grond, and when to use
what. Special things like alignment with picks, cross correlation based fitting
log amp spec fitting, overall amplitude fitting, envelope fitting, etc can be
described here (with examples). Technical details should go to the library
reference.

**Targets** are generic data representations, derived or postprocessed 
from observables or synthesised data. A target can be, a filtered waveform, a spectrum or InSAR displacement. Each target has properties which can be tuned. These can be frequency filters, selected observed components and essentially a Green's functions store which will modell the synthetics for a particular target.

The waveform target groups can be combined to solve the inverse problem, leading to a joint optimisation from different data sources and observeables.


.. note ::

    Weighting between targets is described in the :doc:`../method/index` section.


Waveform targets
----------------

**Tapering**:
    Time-windows around phase arrivals of interest are cut out and tapered. Only these windows will be compared:
    :math:`{\bf d}_{raw, synth}` and the restituted observed waveforms. Only these parts are used in the misfit calculation. The taper window duration is configured with the waveform `targets`_ as well as the phase.

    The tapering is source-model dependent, since the tapering time is given 
    with respect to the theoretic phase arrival time. This arrival time depends on the source location, which is often part of the optimisation itself and therefore may change continuously with each iteration.
    Therefore, restitution, tapering and filtering are done for each misfit calculation anew. Grond uses the pyrocko `CosTaper`_ taper, which is a taper-function with Cosine-shaped fade-outs. The `fade-out` time can be configured or it is calculated as the inverse of the minimum frequency of the chosen bandpass filter.


**Filtering**: 
    Filtering to the desired frequency band is part of the 
    restitution. Minimum and maximum (``fmin, fmax``) frequencies can be configured.

**Lp Normalisation**:
    The `Lp normalisation <https://en.wikipedia.org/wiki/Lp_space>`_ for calculating the waveform misfit.

**Domain**:
    Can be selection from

    * ``time_domain``
        Misfit calculated in time domain, here it is useful to configure the ``tautoshift_max`` and ``autoshift_penalty_max`` to allow for small time shifts of the synthetic data.

    * ``frequency_domain``
        Waveform misfit is calculated in the frequency domain.

    * ``log_frequency_domain``
        Waveform misfit is calculated in the logarithmic frequency domain.

    * ``envelope``
        Waveform envelops are compared.

    * ``absolute``
        The absolute amplitudes are used to calculate the misfit

    * ``cc_max_norm``
        Misfit is calculated from cross-correlation of the traces.



Example :class:`~grond.targets.waveform.WaveformTargetGroup` configuration block:

.. code-block :: yaml

  - !grond.WaveformTargetGroup
      enabled: true
      normalisation_family: time_domain
      path: all
      weight: 1.0
      distance_min: 10000.0
      distance_max: 1000000.0
      channels: [Z, R, T]
      misfit_config: !grond.WaveformMisfitConfig
        fmin: 0.01
        fmax: 0.1
        ffactor: 1.5
        tmin: vel_surface:5.5
        tmax: vel_surface:3.0
        domain: time_domain
        norm_exponent: 2
        tautoshift_max: 0.0
        autoshift_penalty_max: 0.0
      interpolation: multilinear
      store_id: crust2_ib



Satellite targets
-----------------

Observations of spatial surface displacements as derived from unwrapped InSAR data. These data must be hold in a special container format and prepared using the `kite <https://pyrocko.org/#kite>`_ software package.

Prior to optimisation we have to parametrise a quadtree of the surface displacements (spatial sub-sampling) and pre-calculate the data's covariance matrix with kite's ``spool`` tool:

.. code-block :: bash

    spool events/<event_name>/data/insar/scene_ascending.yml

Please see `kite's documentation <https://pyrocko.org/docs/kite/current/>`_ for insights into the pre-processing methods.

**Scene ID**:
    The InSAR scenes are identified by their kite ``scene_id``. Scenes can be explicitly selected, or the wildcard ``*all`` can be used.

**Optimise Orbital Ramps**:
    Optimisation for a 2D offset plane in each InSAR scene. This will compensate tradeoffs between the earthquake signal and uncorrected trends in the unwrapped surface displacements.
    The slopes of ``ramp_north`` and ``ramp_east`` are given in :math:`\frac{m}{m}`, the offset in :math:`m` - these parameters have to be tuned with touch.


Example :class:`~grond.targets.satellite.SatelliteTargetGroup` configuration block:

.. code-block :: yaml

    - !grond.SatelliteTargetGroup
      enabled: true
      normalisation_family: insar_target
      path: all
      weight: 1.0
      kite_scenes: ['*all']
      misfit_config: !grond.SatelliteMisfitConfig
        optimise_orbital_ramp: true
        ranges:
          offset: -0.5 .. 0.5
          ramp_east: -1e-4 .. 1e-4
          ramp_north: -1e-4 .. 1e-4
      interpolation: multilinear
      store_id: crust2_ib_static


GNSS campaign targets
---------------------

True 3D surface displacement as measured by GNSS stations can be included in the inversion process by defining a :class:`~grond.targets.gnss_campaign.GNSSCampaignTargetGroup`. The station's displacement data has to be stored according to :mod:`~pyrocko.model.gnss_campaign`. Please refer to pyrocko's documentation of the GNSS model (`See example <https://pyrocko.org/docs/current/library/examples/gnss_data.html>`_)

**GNSS Campaigns Name**:
    The campaigns are identified by their ``campaign_name``. Campaigns can be explicitly selected, or the wildcard ``*all`` can be used.

Example :class:`~grond.targets.gnss_campaign.GNSSCampaignTargetGroup` configuration block:

.. code-block :: yaml

    - !grond.GNSSCampaignTargetGroup
      enabled: true
      normalisation_family: gnss_target
      path: all
      weight: 1.0
      gnss_campaigns: ['*all']
      misfit_config: !grond.GNSSCampaignMisfitConfig {}
      interpolation: multilinear
      store_id: crust2_ib_static

