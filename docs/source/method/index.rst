Method
======

TODO: (general) this section should be as self-contained as possible, describe 
the method in general - give references to other sections how things are
implemented in Grond.

The very core of the optimisation is the data-point-wise calculation of the 
difference between observed and predicted data: 
:math:`{\bf d}_{obs} - {\bf d}_{synth}`.
Here described is the method how Grond handles this difference to define
misfit and objective functions and how the optimisation is set up to find the 
optimum models and to estimate model uncertainties.

`Observed data` here means full waveforms that are tapered to the defined 
phases, restituted and filtered. `Synthetic waveforms` are the forward-
modelled waveforms that are tapered and filtered in the same way as the 
observed waveforms. 

From the difference :math:`{\bf d}_{obs} - {\bf d}_{synth}` the 
`misfit` is defined based
on a certain L-norm and by the use of certain data weights. Then there are 
the objective functions defined in Grond that determine the ...





Forward modelling with precalculated Green's functions
------------------------------------------------------

The forward modelling of synthetic data for earthquake source models requires
the calculation of the Green's functions between all source points and 
receiver positions. In the general source problem, the positions of the sources
change during the optimisation because the misfit is calculated for many
different source receiver configurations. The calculation of the Green's
functions for each specific source-receiver pair would be a significant part 
of the computational effort in the optimisation and make is rather slow.
Therefore, in Grond pre-calculated Green's functions used that have been 
created with the `Pyrocko fomosto module`_.

...



Objective function design
-------------------------

multiple objective functions bootstrapping

....


Tapering
........

Tapering is done on the waveform targets and means that parts of the waveforms
containing the to-be-modelled, specific seismic phases are cut out. Only these
parts are then compared to forward modelled phases. 
The taper window duration is configured ( target config link) as well as the 
time. 

TODO: set link (line above)

However, the tapering time is given with respect to the theoretic phase arrival
time. This arrival time depends on the source location, which is often part of 
the optimisation itself and therefore may change continuously. Therefore, 
restitution, tapering and filtering are done for each misfit calculation anew.
Grond uses the pyrocko `CosTaper`_ taper, which is a taper-function with 
Cosine-shaped fade-outs. The `fade-out` time can be configured or it is 
calculated as the inverse of the minimum frequency of the chosen bandbass 
filter.


Filtering
.........

filtering to the desired fraquency band is part of the restitution.

TODO: Explain the ffactor
filter factor fmin/ffactor  fmax *ffactor

Misfit calculation
..................

The core of the optimisation is the data-point-wise  calculation of the 
difference between observed and predicted data: 
:math:`|{\bf d}_{obs} - {\bf d}_{synth}|`.

The misfit is based on the configurable :math:`L_x`-norm with 
:math:`x \quad \epsilon \quad [1, 2, 3, ...]`. For the often used norms 
:math:`L_1` and :math:`L_1`

.. math::
  :nowrap:

  \begin{align*}
    \mathrm{misfit}& &= \lVert {\bf{d}}_{obs} - {{\bf d}}_{synth} \rVert_x  &= \
        [\sum{({ d}_{i, obs} - {d}_{i, synth})^x}]^{\frac{1}{x}}\\
    \mathrm{misfit}&_{norm} &= \lVert {\bf{d}}_{obs}  \rVert_x  &= \
        [\sum{{ d}_{i, obs}^x}]^{\frac{1}{x}}
  \end{align*}
 
In case data weights are applied to the data points, the misfit calculation 
changes in the general case of uncorrelated data errors to:

.. math::

  \mathrm{misfit} = (\sum{ ({\bf W}^{\frac{1}{x}}({\bf{d}}_{obs} - \
  {{\bf d}}_{synth}))^{x}})^{\frac{1}{x}}

  \mathrm{misfit_{norm}} = (\sum{ ({\bf W}^{\frac{1}{x}}{\bf{d}}_{obs} )^{x}})^{\frac{1}{x}}

TODO: weights as factors (balancing manual bootstrap), 
Normalization fa(obs - synth) / obs), 
balancing_weight 1/misfit misfit is the synthetic waveform
normalization families


Weighting
.........



Grond implements several different kinds of weights. There are data weight 
factors or weight matrices based on empirical data error estimates - sometimes 
including  data error correlation (``station_noise_weight``, 
``variance_covariance_matrix``). There are data weight factors to balance the 
signal amplitudes of waveforms (``target_balancing_weight``). And there are 
arbitrary, user-defined data weights (``manual_weight``).
The different types of weights often apply to certain target groups.

**data weights** 


**Weights for waveform targets**:

``balancing weights``: from ``TargetBalancingAnalyser``


``station_noise_weight``: from ``Noise_Analyser``


**satellite targets**

**GNSS targets**





The bootstrap method
--------------------

``bootstrap_weights``

Optimisation 
------------

The BABO optimiser
..................

.. _Pyrocko fomosto module: https://pyrocko.org/docs/current/apps/fomosto/index.html
.. _CosTaper: https://pyrocko.org/docs/current/library/reference/trace.html#module-pyrocko.trace