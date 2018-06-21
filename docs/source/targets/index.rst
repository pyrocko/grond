Targets
=======

TODO: describe the different target types available in Grond, and when to use
what. Special things like alignment with picks, cross correlation based fitting
log amp spec fitting, overall amplitude fitting, envelope fitting, etc can be
described here (with examples). Technical details should go to the library
reference.

**Targets** are data representations, derived or postprocessed 
from observables or synthesised data. For example, a filtered waveform can be
a target. It is defined through the filter, the component, a Green's functions
store for the synthetics and else. A `target` can also be a waveform spectrum.


Waveform targets
----------------

**Tapering**:
    Tapering is done on the waveform targets and means that parts of the 
    recorded waveforms are cut out. Here the taper applies to specific seismic
    phases in both the raw synthetic waveforms calculated 
    :math:`{\bf d}_{raw, synth}` and 
    the recorded, restituted waveforms. Only these parts are then used in the 
    misfit calculation. 
    The taper window duration is configured with the waveform `targets`_ as 
    well as the time. 

    The tapering is source-model dependent, since the tapering time is given 
    with respect to the theoretic phase arrival
    time. This arrival time depends on the source location, which is often part of 
    the optimisation itself and therefore may change continuously. Therefore, 
    restitution, tapering and filtering are done for each misfit calculation anew.
    Grond uses the pyrocko `CosTaper`_ taper, which is a taper-function with 
    Cosine-shaped fade-outs. The `fade-out` time can be configured or it is 
    calculated as the inverse of the minimum frequency of the chosen bandpass 
    filter.


**Filtering**: 
    Filtering to the desired frequency band is part of the 
    restitution. Minimum and maximum frequencies are configured.

    TODO Explain the ffactor: filter factor fmin/ffactor  fmax ffactor
    
Satellite targets
-----------------




GNSS campaign targets
---------------------
