Method
======

TODO: (general) this section should be as self-contained as possible, describe 
the method in general - give references to other sections how things are
implemented in Grond.

The very core of the optimisation is the data-point-wise calculation of the 
difference between observed and predicted data: 
:math:`{\bf d}_{obs} - {\bf d}_{synth}`.
Here described is the method on:

1. what exactly are the observed data :math:`{\bf d}_{obs}` and the synthetic 
   data :math:`{\bf d}_{synth}`
2. how handles Grond the difference :math:`{\bf d}_{obs} - {\bf d}_{synth}` 
   with respect to defining the objective functions through misfits and data
   weighting,
3. how the optimisation is set up to search the model space to find the 
   optimum models and 
4. which methods are used to estimate model uncertainties.

`Observed data` here means full waveforms that are tapered to the defined 
phases, restituted and filtered. `Synthetic waveforms` are the forward-
modelled waveforms that are tapered and filtered in the same way as the 
observed waveforms. 

From the difference :math:`{\bf d}_{obs} - {\bf d}_{synth}` the 
`misfit` is defined based
on a certain L-norm and by the use of certain data weights. Then there are 
the objective functions defined in Grond that determine the ...





Forward modelling with pre-calculated Green's functions
-------------------------------------------------------

The forward modelling of raw synthetic data  :math:`{\bf d}_{raw, synth}` for 
earthquake source models requires the calculation of the Green's functions
(GF) between all source points and 
receiver positions based on a medium model. In the general source problem, 
the positions of the sources change during the optimisation because the 
misfit is calculated for many different source receiver configurations. 
The calculation of the GFs for each specific source-receiver 
pair would be a significant part of the computational effort in the 
optimisation and make is rather slow.
Therefore, in Grond pre-calculated GFs, stored in a database called 'GF store`,
are used that have been created with the `Pyrocko fomosto module`_. 

Generally, we distinguish different types of GF stores (for detail see the 
`Pyrocko fomosto module`_ documentation). For the options possible in Grond
at the moment
we need only to distinguish between GFs that have been calculated based on 
different GF methods to either allow for the fast forward
calculation of dynamic seismic waveforms or static near-field displacements.


GF stores can be searched and downloaded on our `GF store database`_, for the 
some general global seismic waveform analyses and/or for InSAR 
and GNSS data analyses based, e.g. on the global 1d `PREM model`_,.
For more specific analyses, based on an individual choice of the medium, the
GF store can be created - usually very easy - with the
`Pyrocko fomosto module`_ for different GF methods.


**GFs for seismic waveforms**

For regional data analyses with optional near-field terms the ``QSEIS`` method 
by for layered media by `Wang et al.`_ (1999) is appropriate. For global data 
the ``QSSP`` method also by `Wang et al.`_ (2017) is more suited. 
 
 

**GFs for static near-field displacements** (measured by using GNSS or InSAR)

For the calculation of purely static coseismic displacements the use of the 
``PSGRN/PSCMP`` method by `Wang et al.`_ (2006) is suggested for fast 
forward modelling.


Objective function design
-------------------------

The `objective function` gives a scalar value based on which a source model is
evaluated to be better (smaller values) or worse (larger value) than other
source models. It is often called `misfit function`. The source model that 
results in the smallest values of the `objective function` is the global 
minimum of the misfit function optimum model



Tapering
........

Tapering is done on the waveform targets and means that parts of the waveforms
containing the to-be-modelled, specific seismic phases are cut out of raw 
synthetic waveforms calculated :math:`{\bf d}_{raw, synth}`. Only these
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

filtering to the desired frequency band is part of the restitution.

TODO Explain the ffactor
filter factor fmin/ffactor  fmax ffactor

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
.. _GF store database: http://kinherd.org/gfs.html

.. _PREM model: http://ds.iris.edu/spud/earthmodel/9991844
.. _Wang et al.: https://www.gfz-potsdam.de/en/section/physics-of-earthquakes-and-volcanoes/data-products-services/downloads-software/
