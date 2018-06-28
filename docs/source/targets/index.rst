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

**Lp Normalisation**
    The `Lp normalisation <https://en.wikipedia.org/wiki/Lp_space>`_ factor for calculating the waveform misfit.

**Domain**
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



Example :class:`~grond.targets.waveform.WaveformTargetGroup` configuration

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

Observations of spatial surface displacements as derived from unwrapped InSAR data. These data must be prepared using the `kite <https://pyrocko.org>`_ software package.

Prior to optimisation we have to define a subsamples Quadtree and Covariance matrix.

**Scene ID**
    The InSAR scenes are identified by their kite ``scene_id``, they can be explicitly selected or a wildcard ``*all*`` can be used.

**Optimise Orbital Ramps**
    Activates the optimisation of a 2D offset plane in the InSAR. This will compensate tradeoffs between the earthquake signal and trends in the unwrapped surface displacements.
    The ranges are given in :math:`\frac{m}{m}`, these parameters have to be tuned with touch.

Example :class:`~grond.targets.satellite.SatelliteTargetGroup` configuration

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


Example :class:`~grond.targets.gnss_campaign.GNSSCampaignTargetGroup` configuration

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

