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



Misfit calculation
..................


The core of the optimisation is the data-point-wise  calculation of the 
difference between observed and predicted data: 
:math:`|{\bf d}_{obs} - {\bf d}_{synth}|`. Not entire traces and and not the
full spectrum of a trace are modelled. Before, observed and synthetic data 
are tapered and filtered.

**Tapering**:
    Tapering is done on the waveform targets and means that parts of the recorded 
    waveforms are cut out. Here the taper applies to specific seismic phases in
    both the raw synthetic waveforms calculated :math:`{\bf d}_{raw, synth}` and 
    the recorded, restituted waveforms. Only these parts are then used in the 
    misfit calculation. 
    The taper window duration is configured ( target config link) as well as the 
    time. 

    TODO: set link (line above)

    The tapering is source-model dependent, since the tapering time is given with 
    respect to the theoretic phase arrival
    time. This arrival time depends on the source location, which is often part of 
    the optimisation itself and therefore may change continuously. Therefore, 
    restitution, tapering and filtering are done for each misfit calculation anew.
    Grond uses the pyrocko `CosTaper`_ taper, which is a taper-function with 
    Cosine-shaped fade-outs. The `fade-out` time can be configured or it is 
    calculated as the inverse of the minimum frequency of the chosen bandbass 
    filter.


**Filtering**: 
    Filtering to the desired frequency band is part of the 
    restitution. Minimum and maximum frequencies are configured.

    TODO Explain the ffactor: filter factor fmin/ffactor  fmax ffactor


The misfit is based on the configurable :math:`L_x`-norm with 
:math:`x \,\, \epsilon \,\, [1, 2, 3, ...]`:

.. math::

    \lVert e \rVert_x = \lVert {\bf{d}}_{obs} - {{\bf d}}_{synth} \rVert_x  = \
        (\sum{|{ d}_{i, obs} - {d}_{i, synth}|^x})^{\frac{1}{x}}.
        
Also the norm of the data is associated with each misfit. This measure will be 
used to normalize the misfit values:
        
.. math::
        
    \lVert e_{\mathrm{0}} \rVert_x = \lVert {\bf{d}}_{obs}  \rVert_x  = \
        (\sum{|{d}_{i, obs}|^x})^{\frac{1}{x}}.

The normalized misfit

.. math::
 
    \lVert e_{\mathrm{norm}} \rVert_x = \
    \frac{\lVert e \rVert_x}{ \lVert e_{\mathrm{0}} \rVert_x}.

is a useful measure to evaluate the data fit at a glance. Only for model
predicitions that manage to explain parts of the observed data holds
:math:`\lVert e_{\mathrm{norm}} \rVert_x <1`. Furthermore, the data norm 
:math:`\lVert e_{\mathrm{0}} \rVert_x` is used in the normalization of data
groups.



Weighting
.........


Grond implements several different kinds of weights:

* target balancing (for waveforms ony)
* noise-based data weights
* user-defined, manual weights
* normalisation within data groups

These weights are applied as factors to the misfits, optionally as a products
of a combination of several of these weights.
Generally, the misfit and data norm calculations with a data weights vector 
:math:`\bf{w}` change to:

.. math::
  :nowrap:
  
  \begin{align*}
    \lVert e \rVert_x &= (\sum{ ({w_i}|{{d}}_{i,obs} - \
  {{ d}}_{i,synth}|)^{x}})^{\frac{1}{x}}\\
    \lVert e_{\mathrm{0}} \rVert_x  &= (\sum{ ({w_i} \ 
       |{{d}}_{i,obs} |)^{x}})^{\frac{1}{x}}
  \end{align*}
  
**Target balancing weights**

With these weights waveforms are `balanced` with respect to the expected signal
amplitude. 
Signal amplitudes in a trace :math:`|s_j|` depend on the source-receiver 
distance, on the phase type and the taper used. The problem tackled 
with this weight is that
large signal amplitude have higher contributions to the misfit than smaller
signal amplitudes, without carrying more information. From synthetic 
waveforms of `N` forward models that have been randomly drawn from the defined 
model space the mean signal amplitude of the traces is derived. The weight 
for each trace is simply the inverse of these mean signal amplitudes:

.. math::
    
    w_j = 1 / \sum_{i=1}^{N}|s_{ji}|.


Like this small 
signal are enhanced in the
objective function and large signals supressed. This is described as 
`adaptive station weighting` in the PhD `thesis by Heimann`_ (2011) (page 23).
In Grond they are called ``balancing weights`` and are received from the
``TargetBalancingAnalyser`` before the optimization.
  
**Data weights based on data error statistics**

There are direct data weight vectors :math:`\bf{w}` or weight matrices
:math:`\bf{W}` based on empirical data error variance estimates. Partly, e.g. 
for InSAR and GNSS data, these can include data error 
correlations expressed in the data error variance-covariance matrix 
:math:`\bf{\Sigma}`: 

.. math::
  {\bf w} = \frac{1}{{\bf \sigma}}, \quad  \bf{W} = \sqrt{{\bf \Sigma}^{-1}}.

  
For a ``WaveformTarget``  the data error statistics stem from real recordings 
of noise before the first phase arrival as described e.g. in 
`Duputel et al.`_ (2012). From the noise traces the inverse of their
standard deviation is used. In Grond they are called `station_noise_weights`` 
and are received from the ``Noise_Analyser`` before the optimization.

For a ``SatelliteTarget`` the data error statistics are loaded with the data 
sets. The estimation of the noise statistics has to be done before Grond
by using `kite`_.
In `kite`_ the noise estimation can be done in areas of the displacement map
that are not affected by coseismic deformation by using spatial sampling
methods and semi-variogram and covariogram formation, described e.g. in
`Sudhaus and Jonsson`_ (2009).

For a ``GNSSCampaignTarget`` the data error statistics are also loaded with
the data set. They have to be estimated before and given in the GNSS data 
`YAML`-file describing the data set. For details visit the corresponding 
chapter in the `Pyrocko tutorial`_. 

**manual data weight**

User-defined manual data weights enable an arbitrary weighting of data sets. 
No rules apply other from the user's rationale. In Grond they are called 
``manual_weight`` and are given in the configuration file.


TODO link to the target sheet

**normalization of data groups** 


TODO: weights as factors (balancing manual bootstrap), 
Normalization fa(obs - synth) / obs), 
balancing_weight 1/misfit misfit is the synthetic waveform
normalization families


The bootstrap method
--------------------

``bootstrap_weights``

Bayesian bootstrap

Residual bootstrap is a  computationally more light version of the 
`randomize-then-optimize` procedure. With empirical estimates of the data 
error statistics we add synthetic random noise to all residuals to evaluate the
misfit anew. 
keeping 


Optimisation 
------------

The BABO optimiser
..................

.. _Pyrocko fomosto module: https://pyrocko.org/docs/current/apps/fomosto/index.html
.. _CosTaper: https://pyrocko.org/docs/current/library/reference/trace.html#module-pyrocko.trace
.. _GF store database: http://kinherd.org/gfs.html
.. _kite: https://pyrocko.org/docs/kite/current/

.. _PREM model: http://ds.iris.edu/spud/earthmodel/9991844
.. _Wang et al.: https://www.gfz-potsdam.de/en/section/physics-of-earthquakes-and-volcanoes/data-products-services/downloads-software/
.. _Duputel et al.: https://academic.oup.com/gji/article/190/2/1243/645429
.. _Sudhaus and Jonsson: https://academic.oup.com/gji/article/176/2/389/2024820
.. _YAML: http://yaml.org/
.. _Pyrocko tutorial: https://pyrocko.org/docs/current/library/examples/gnss_data.html
.. _thesis by Heimann: http://ediss.sub.uni-hamburg.de/volltexte/2011/5357/pdf/Dissertation.pdf
