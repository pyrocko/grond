Analysers
=========

Analysers do their work before the actual optimisations. They analyse the data 
with respect to the problem-dependent expected signal strength for each target 
and with respect to the pre-event noise power. 
The analyser results are used to weight the data in the following optimisation.

TODO: review, link to the reference, link to papers?

There are two analysers implemented so far, which work on waveform targets:

``TargetBalancingAnalyser``:
    This analyser balances, based on estimated signal amplitudes, the waveform 
    target contributions to the misfit. 
    This is done with synthetic waveforms from random source models drawn from
    the model space defined in the problem configuration. Waveforms from 
    stations that are far away from the source will on average have small 
    signal strengths, while targets of stations close to the sources will have 
    larger signal amplitudes. To balance the importance of waveforms the 
    inverse of the mean signal strength is used as a ``balancing_weight``. 
    Like this the effect of simple geometrical spreading or through the 
    radiation pattern is lessened in the misfit evaluation. 
    
    
``NoiseAnalyser``:
    This analyser evaluates the pre-event noise for the waveform targets to 
    form empirical target weights. High pre-event noise variance leads to a 
    small target weight for the corresponding target in the misfit calculation.
    It is assumed that the pre-event noise is random background noise and 
    stationary in the time from the first pre-event sample to the last analysed 
    phase sample. A check for significant earthquakes in the critical time 
    frame can be enabled.
    
    
``TargetBalancingAnalyser`` configuration
-----------------------------------------

The balancing of target is problem-dependent. This means the configuration of 
the problem also configures the ``target_balancing``. The only extra parameter 
needed is 

``niterations``: an integer number that defines how many random waveform 
    predictions the target-balancing weights are based. 
    A meaningful number will depepend on the problem. Large model spaces 
    (loose model parameter bounds) may require a larger number of predicitons 
    for more stable average values. For more tightly constrained problems the 
    signal strength at the targets for the not-so-different source models, 
    which are drawn from a small model space, may not vary much. A smaller 
    number may suffice. Of course, the computional effort increases linearly 
    with the number of ``niterations``.

.. code-block :: sh
 
  analyser_configs:
    - !grond.TargetBalancingAnalyserConfig
      niterations: 1000
      

``NoiseAnalyser`` configuration
-------------------------------

This analyser is not strongly model-independent. Only the reference time given 
in 'event.txt' is used to estimate the phase arrivals. These phase arrivals 
define what `pre-event` is.

``nwindows``: is an integer number. It gives the number of sub-windows in which
    the noise trace of duration ``pre_event_noise_duration`` is split. If 
    larger than ``1``, the noise variance is in each sub-window and the average
    noise variance is returned.

``pre_event_noise_duration``: defines how long (in seconds) the pre-event noise
    trace is. 

``check_events``: is a boolean value. If ``True`` the iris global earthquake 
    catalog is looked up to find if phase arrivals of other events may disturb 
    the statistics of the noise.

``phase_def``: is a string that defines the reference phase for the pre-event 
    time window.
      
.. code-block :: sh
      
  analyser_configs:
    - !grond.NoiseAnalyserConfig
      nwindows: 1
      pre_event_noise_duration: 500.
      check_events: False
      phase_def: P
