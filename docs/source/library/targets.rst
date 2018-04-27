The :mod:`targets` Module
=========================

The :class:`Target` Base
------------------------

This is the abstract base implementation all implemented `targets` inherit from.
See :class:`~grond.targets.WaveformMisfitTarget`, :class:`~grond.targets.SatelliteMisfitTarget`, :class:`~grond.targets.GNSSCampaignMisfitTarget` for implementations of :class:`~grond.targets.MisfitTarget`


.. automodule :: grond.targets.base
    :members: TargetGroup, MisfitTarget, MisfitConfig,
        MisfitResult


The Waveform Target
-------------------

.. automodule :: grond.targets.waveform
    :members: WaveformTargetGroup, WaveformMisfitTarget, WaveformMisfitConfig,
        WaveformMisfitResult


The Satellite Target
--------------------

.. automodule :: grond.targets.satellite
    :members: SatelliteTargetGroup, SatelliteMisfitTarget, SatelliteMisfitConfig,
        SatelliteMisfitResult


The GNSS Campaign Target
-------------------------

.. automodule :: grond.targets.gnss_campaign
    :members: GNSSCampaignTargetGroup, GNSSCampaignMisfitTarget,
        GNSSCampaignMisfitConfig, GNSSCampaignMisfitResult
