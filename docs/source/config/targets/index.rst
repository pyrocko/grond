Targets
=======

**Targets** is what Grond tries to match: The misfit between an *observed* target and forward model is minimized. Targets are derived from observables or synthesised data. A target can be, a filtered waveform, a spectrum, InSAR or GNSS displacement. Each target has properties which can be tuned. These can be frequency filters, selected observed components and essentially a Green's functions store which is responsible for the synthetics at a particular target.

Different ``TargetGroups`` can be combined to solve the inverse problem, leading to a joint optimisation from different data sources and observables.

.. note ::

    Weighting between different targets is described in the :doc:`/method/index` section.


General ``Target`` configuration
--------------------------------

Parameters valid for all types of ``MisfitTargets`` are:

.. glossary ::

  ``normalisation_family``
    Normalisation family (see the Grond documentation for how it works). Use distinct normalisation families when mixing misfit contributors with different magnitude scaling, like e.g. cross-correlation based misfit and time-domain :math:`L^p` norm.

  ``weight``
    How to weight contributions from this group in the global misfit.

  ``path``
    Just a name used to identify targets from this group. Use dot-separated path notation to group related contributors.

  ``interpolation``
      Green's function store interpolation, choose from:

      * ``multilinear``
          Performs a linear interpolation between discrete Green's function for improved resolution of synthetic data. *This option is computationally more expensive.*

      * ``nearest_neighbor``
          Uses the Green's function calculation for the forward model.

      Choices other than 'nearest_neighbor' may require dense GF stores to avoid aliasing artefacts in the forward modelling.

  ``store_id``
      Name of the GF Store to use.


Waveform targets
----------------

See the :doc:`dataset configuration <../dataset/index>` for loading waveforms and response information.

.. code-block :: yaml
  :caption: Example ``WaveformTarget`` configuration

  - !grond.WaveformTargetGroup
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


.. glossary ::

  **Tapering**
      ``tmin`` and ``tmax`` define time-windows around phase arrivals of interest, those are cut out and tapered.

      :math:`{\bf d}_{raw, synth}` and the restituted observed waveforms. Only these parts are used in the misfit calculation. The taper window duration is configured for each seismic station individually by phase arrivals.

      The tapering is source-model dependent, since the tapering time is given with respect to the theoretic phase arrival time. This arrival time depends on the source location, which is often part of the optimisation itself and therefore may change continuously with each iteration. Therefore, restitution, tapering and filtering are done for each misfit calculation anew. Grond uses the Pyrocko `CosTaper`_ taper. The ``fade_out`` time can be configured or it is calculated as the inverse of the minimum frequency of the chosen bandpass filter.

  **Frequency filtering**
      ``fmin`` and ``fmax`` in Hz define the desired bandpass filter.

  ``norm_exponent``
      The `Lp normalisation <https://en.wikipedia.org/wiki/Lp_space>`_ for calculating the waveform misfit.

  ``domain``
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

  ``tautoshift_max``
      defines the maximum allowed time uin seconds the observed and synthetic trace may be shifted during the inversion.

  ``autoshift_penalty_max``
      is the misfit penalty for autoshifting seismic traces.

Example :class:`~grond.targets.waveform.WaveformTargetGroup` configuration section:


Satellite targets
-----------------

.. code-block :: yaml
    :caption: Example ``SatelliteTarget`` configuration

    - !grond.SatelliteTargetGroup
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

Observations of spatial surface displacements as derived from unwrapped InSAR data. These data must be hold in a special container format and prepared using the `kite <https://pyrocko.org/#kite>`_ software package.

Prior to optimisation we have to parametrise a quadtree of the surface displacements (spatial sub-sampling) and pre-calculate the data's covariance matrix with kite's ``spool`` tool:

.. code-block :: bash

    spool events/<event_name>/data/insar/scene_ascending.yml

Please see `kite's documentation <https://pyrocko.org/docs/kite/current/>`_ for insights into the pre-processing methods.

.. glossary::

  ``kite_scenes``
    The InSAR scenes are identified by their kite ``scene_id``. Scenes can be explicitly selected, or the wildcard ``*all`` can be used.

  ``optimise_orbital_ramp``:
    Optimisation for a 2D offset plane in each InSAR scene. This will compensate tradeoffs between the earthquake signal and uncorrected trends in the unwrapped surface displacements.
    The slopes of ``ramp_north`` and ``ramp_east`` are given in :math:`\frac{m}{m}`, the offset in :math:`m` - these parameters have to be tuned with touch.


Example :class:`~grond.targets.satellite.SatelliteTargetGroup` configuration section:


GNSS campaign targets
---------------------

.. code-block :: yaml
    :caption: Example ``GNSSTarget`` configuration

    - !grond.GNSSCampaignTargetGroup
      normalisation_family: gnss_target
      path: all
      weight: 1.0
      gnss_campaigns: ['*all']
      misfit_config: !grond.GNSSCampaignMisfitConfig {}
      interpolation: multilinear
      store_id: crust2_ib_static

True 3D surface displacement as measured by GNSS stations can be included in the inversion process by defining a :class:`~grond.targets.gnss_campaign.GNSSCampaignTargetGroup`. The station's displacement data has to be stored according to :mod:`~pyrocko.model.gnss_campaign`. Please refer to Pyrocko's documentation of the GNSS model (`See example <https://pyrocko.org/docs/current/library/examples/gnss_data.html>`_)

.. glossary ::

  ``gnss_campaigns``
    The campaigns are identified by their ``campaign_name``. Campaigns can be explicitly selected, or the wildcard ``*all`` can be used.

Example :class:`~grond.targets.gnss_campaign.GNSSCampaignTargetGroup` configuration section:


.. _CosTaper: https://pyrocko.org/docs/current/library/reference/trace.html#module-pyrocko.trace
