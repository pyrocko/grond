Analysers
=========

Analysers do their work before the actual optimisations, they analyse the observed data to weigh the :term:`targets <target>` in the following optimisation.

Grond come with two analysers, which work on waveform targets:

``TargetBalancingAnalyser``:
    This analyser balances the waveform target contributions to the misfit, based on estimated signal amplitudes.

    Synthetic waveforms from random source models are drawn from the model space defined in the problem configuration. Waveforms at stations far away from the source will have small signal strengths, while stations close to the sources will show larger signal amplitudes on average. To balance the importance of waveforms the inverse of the mean signal strength is used as a ``balancing_weight``. This compensates dampening effects from geometrical spreading and averages the radiation pattern towards the misfit evaluation.
    
    
``NoiseAnalyser``:
    This analyser evaluates the pre-event noise for the waveform targets to define empirical target weights. 
    
    This analyser assumes that the pre-event noise is random background noise and stationary in the time. However, a check for significant earthquakes in the critical time frame can be enabled.
    If this analyser is enabled, high pre-event noise variance optionally leads to either a small :term:`target` weight for the corresponding target in the misfit calculation 
    (``mode='weighting'``) or to zero weight (``mode='weeding'``). The latter case is realized if the noise of a trace exceeds a noise level that is above the median noise by a configurable amount (``cutoff``).
    
    
``TargetBalancingAnalyser`` configuration
-----------------------------------------

Balancing of target weights depends on the configured problem search space. This means the configuration of the problem also configures the ``target_balancing``. The only extra parameter needed is

``niterations``

    an integer number defining how many random waveform predictions are modelled. A meaningful number will depend on the problem. Large model spaces (loose model parameter bounds) may require a larger number of predictions for stable average values. Tightly constrained problems can draw from a small model space and may not vary much - here a smaller number is suffice.

    The computational effort increases linearly with the number of ``niterations``.

.. code-block :: yaml
 
  analyser_configs:
    - !grond.TargetBalancingAnalyserConfig
      niterations: 1000
      

``NoiseAnalyser`` configuration
-------------------------------

This analyser is independent from the problem. Only the reference time given in :file:`event.txt` is used to estimate the phase arrivals. These phase arrivals define where `pre-event` is.

.. glossary ::

    ``nwindows``
        is an integer number defining the number of sub-windows of the trace. If larger than ``1``, the noise variance in each sub-window and the total average noise variance is calculated.

    ``pre_event_noise_duration``
        defines the trace length of pre-event noise in seconds.

    ``check_events``
        is a boolean value. If ``True`` the IRIS global earthquake catalogue is searched for phase arrivals of other events, which may interfere with the pre-event noise.

    ``phase_def``
        is a string that defines the reference phase for the pre-event time window. See `Pyrocko's definition of phases <https://pyrocko.org/docs/current/apps/cake/manual.html>`_.
        
    ``statistics``
        Set weight to inverse of noise variance (var) or standard deviation (std).
        
    ``mode``
        Generate weights based on inverse of noise measure (weighting), or discrete on/off 
        style in combination with cutoff value (weeding).
        
    ``cutoff``
        Set weight to zero, when noise level exceeds median by the given cutoff factor.
      
.. code-block :: yaml
      
  analyser_configs:
    - !grond.NoiseAnalyserConfig
      nwindows: 1
      pre_event_noise_duration: 500.
      check_events: False
      phase_def: P
      statistics: 'var'
      mode: 'weeding'
      cutoff: 2.
