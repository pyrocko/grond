Method
======

TODO: REVIEW: This section should be as self-contained as possible, describe 
the method in general - give references to other sections how things are
implemented in Grond.

The very core of any optimisation is the evaluation of a misfit value between
observed :math:`{\bf d}_{obs}` and predicted data :math:`{\bf d}_{synth}`. This
is most often based on the difference  :math:`{\bf d}_{obs} - {\bf d}_{synth}`,
but can also be any other comparison, like a correlation measure for example.


`Observed data` here means post-processed data and not the `raw` measurements.
E.g. full waveforms are usually tapered to the defined 
phases, restituted and filtered. `Synthetic waveforms` are the forward-
modelled waveforms that are tapered and filtered in the same way as the 
observed waveforms. Find details on the post-processing in the `targets`_ 
section. The `targets` are derived from data defined in the `dataset`_.

This sheet describes the method of Grond on:

1. how Grond implements the differences between :math:`{\bf d}_{obs}` and
   :math:`{\bf d}_{synth}` 
   with respect to the definition of objective functions and data
   weighting,
2. how the optimisation is set up to search the model space to find the 
   optimum models and 
3. which methods are used to estimate model uncertainties.


Forward modelling with pre-calculated Green's functions
-------------------------------------------------------

The forward modelling of raw synthetic data :math:`{\bf d}_{raw, synth}` for
earthquake source models requires the calculation of the Green's function (GF)
between all source points and receiver positions involved, based on a medium
model. In the general source problem, the positions of the sources change
during the optimisation because the misfit is calculated for many different
source receiver configurations. The calculation of the GFs for each specific
source-receiver pair is computationally demanding and would be a significant
contribution to the total computational cost of an optimisation. Therefore, in
Grond pre-calculated GFs, stored in a database called 'GF store`, are used that
have been created with the `Pyrocko fomosto module`_. 

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
minimum of the misfit function optimum model.

The objective function defines what a `model fit` is and how `good` or
`poor` models are scaled with respect to others. Furthermore, the
objective function has rules how different data sets are handled, which 
L-norm is applied and how data 
errors are considered in optimisations. 


.. figure:: ../images/illu_combi_weights.svg
    :name: Fig. 1
    :height: 500px
    :align: center
    :alt: alternate text
    
    **Figure 1**: GROND objective function design in an overview illustration. 
    Details on how each function and weight vectors are formed follow below.

    
Misfit calculation
..................


The usual core of an optimisation is the data-point-wise calculation of the 
difference between observed and predicted data: 
:math:`|{\bf d}_{obs} - {\bf d}_{synth}|`. 

In Grond :math:`{\bf d}_{obs}` and :math:`{\bf d}_{synth}` can be

* seismic waveforms traces in time domain
* seismic waveforms in spectral domain
* seismic waveforms in logarithmic spectral domain
* static surface displacements measured by using InSAR or from pixel offsets
* static surface displacements measured by using GNSS sensors

TODO: add spectral phase ratio and more?

The misfit in Grond can further be based on the maximum waveform correlation. 

Not entire traces and and not the
full spectrum of a trace are compared for the misfit evaluation. 
Before, observed and synthetic data are tapered and filtered (see above).

The misfit is based on the configurable :math:`L_x`-norm with 
:math:`x \,\, \epsilon \,\, [1, 2, 3, ...]`:

.. math::
  :label: eq:ms

    \lVert e \rVert_x = \lVert {\bf{d}}_{obs} - {{\bf d}}_{synth} \rVert_x  = \
        (\sum{|{ d}_{i, obs} - {d}_{i, synth}|^x})^{\frac{1}{x}}.
        
Also the norm of the data is associated with each misfit. This measure will be 
used to normalise the misfit values:
        
.. math::
  :label: ns
        
    \lVert e_{\mathrm{0}} \rVert_x = \lVert {\bf{d}}_{obs}  \rVert_x  = \
        (\sum{|{d}_{i, obs}|^x})^{\frac{1}{x}}.

The normalised misfit

.. math::
  :label: ms_ns
 
    \lVert e_{\mathrm{norm}} \rVert_x = \
    \frac{\lVert e \rVert_x}{ \lVert e_{\mathrm{0}} \rVert_x}.

is a useful measure to evaluate the data fit at a glance. Only for model
predictions that manage to explain parts of the observed data holds
:math:`\lVert e_{\mathrm{norm}} \rVert_x <1`. Furthermore, the data norm 
:math:`\lVert e_{\mathrm{0}} \rVert_x` is used in the normalisation of data
groups.

For waveform data correlation the misfit function is based on the maximum
correlation :math:`\mathrm{max}(C)` of :math:`{\bf d}_{obs}` and 
:math:`{\bf d}_{synth}` defined as:

.. math::
  :nowrap:
  :label: cor
  
  \begin{align*}
    e_{\mathrm{cc}} &= \frac{1}{2} - \frac{1}{2}\, \mathrm{max}(C), \, \
    \mathrm{with} \,\,\,
    e_{\mathrm{0, cc}} = \frac{1}{2} \,\, , \mathrm{such\,\, that}  \\
    e_{\mathrm{norm}} &= 1 - \mathrm{max}(C).
  \end{align*}  


Weighting
.........

Grond implements several different kinds of weights:

* :math:`w_{\mathrm{tba},i}` - target balancing (for waveforms only)
* :math:`w_{\mathrm{noi},i}` - noise-based data weights
* :math:`w_{\mathrm{man},i}` - user-defined, manual weights of data groups
* normalisation within data groups (leads to balancing of data groups)

These weights are applied as factors to the misfits, optionally as a product
of weight combinations. E.g. for a waveform all data weights combined means:

.. math::
  :label: wcomb
  
   w_{\mathrm{comb},i} = w_{\mathrm{tba},i} \cdot w_{\mathrm{noi},i} \
   \cdot w_{\mathrm{man},i}.

The misfit and data norm calculations with data weights 
:math:`w_{\mathrm{comb},i}` change to:

.. math::
  :nowrap:
  :label: wms_wns

  \begin{align*}
    \lVert e \rVert_x &= (\sum{ ({w_{\mathrm{comb},i}} \cdot |{{d}}_{i,obs} - \
  {{ d}}_{i,synth}|)^{x}})^{\frac{1}{x}}\\
    \lVert e_{\mathrm{0}} \rVert_x  &= (\sum{ ({w_{\mathrm{comb},i}} \cdot \ 
       |{{d}}_{i,obs} |)^{x}})^{\frac{1}{x}}
  \end{align*}
  
**Target balancing weights**:
    With these weights waveforms are `balanced` with respect to the expected 
    signal amplitude. 
    Signal amplitudes in a trace :math:`|{\bf{d}}_{synth}|` depend on the 
    source-receiver distance, on the phase type and the taper used. The problem 
    tackled with this weight is that
    large signal amplitude have higher contributions to the misfit than smaller
    signal amplitudes, without carrying more information. From synthetic 
    waveforms of `N` forward models that have been randomly drawn from the 
    defined model space the mean signal amplitude of the traces is derived. 
    The weight for each trace is simply the inverse of these mean signal 
    amplitudes:

    .. math::
      :label: wtba
        
      {\bf w}_{\mathrm{tba}} = 1/ \lVert {\bf{d}}_{synth}  \rVert_x  = \
            (\sum^{N}{|{d}_{i, synth}|^x})^{\frac{1}{x}}.

    Like this small 
    signal are enhanced in the
    objective function and large signals suppressed. This is described as 
    `adaptive station weighting` in the PhD `thesis by Heimann`_ (2011) (page 23).
    In Grond they are called ``balancing weights`` and are received from the
    ``TargetBalancingAnalyser`` before the optimisation.

    .. figure:: ../images/illu_target_balancing.svg
        :name: Fig. 2
        :width: 300px
        :align: left
        :alt: alternate text
        :figclass: align-center
        
        **Figure 2**: Qualitative sketch how target balancing weight increases with 
        source distance to balance amplitude decrease caused by geometrical 
        spreading. 

**Data weights based on data error statistics**:
    There are direct data weight vectors :math:`\bf{w}` or weight matrices
    :math:`\bf{W}` based on empirical data error variance estimates. Partly,
    e.g. for InSAR and GNSS data, these can include data error 
    correlations expressed in the data error variance-covariance matrix 
    :math:`\bf{\Sigma}`: 
    
    .. math::
      :label: wnoi

      {\bf w} = \frac{1}{{\bf \sigma}}, \quad  \bf{W} = \sqrt{{\bf \Sigma}^{-1}}.

    For a ``WaveformTarget``  the data error statistics stem from real recordings 
    of noise before the first phase arrival as described e.g. in 
    `Duputel et al.`_ (2012). From the noise traces the inverse of their
    standard deviation is used. In Grond they are called `station_noise_weights`` 
    and are received from the ``Noise_Analyser`` before the optimisation.

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

**manual data weighting**:
    User-defined manual data weights enable an arbitrary weighting of data sets 
    in contrast to balancing of single observations through target balancing and 
    noise-based data weights. 
    No rules apply other than from the user's rationale. In Grond they are called 
    ``manual_weight`` and are given in the configuration file of the `targets`_.

**Normalisation of data and data groups**:
    The normalisation in Grond is applied to data groups that are member of the
    so called ``normalisation_family``. A `normalisation family` in Grond can 
    be composed in many ways. However, it is often meaningful to put data of 
    the same kind and with similar weighting schemes into the same 
    `normalisation family` (see also Fig. 1). 
    This could be P and S waves, or two InSAR data sets. As an explanation some 
    examples are given here:

**Example 1:** Fitting waveforms of P and S waves to solve 
for a source model 

    Let's say we use the waveform fit in time domain and in spectral domain 
    combined. We then have weighted misfits as 
    in Equation :eq:`wms_wns` for P waves with
    :math:`{\bf d}_{obs,\mathrm{Pt}}` 
    and :math:`{\bf d}_{synth,\mathrm{Pt}}` in time domain and 
    :math:`{\bf d}_{obs,\mathrm{Ps}}` and :math:`{\bf d}_{synth,\mathrm{Ps}}` 
    in spectral domain. We have also the corresponding weighted misfit norms 
    (see Equation :eq:`wms_wns`) and the same for S waveforms in time and 
    spectral domain. 
    Let's also say we are using the :math:`L_{\mathrm{2}}\,`-norm. 
    
    The waveforms of P and S waves in time domain are of a similar and kind 
    and can, maybe even should, be normalised together. The same may be 
    meaningful for the normalisation of the P and S waves in spectral domain.  
    
    In Grond we say the time-
    domain data and the spectral-domain data each 
    belong to a different ``normalisation_family``.

    The **global misfit** for two normalisations families will read:


.. math::
  :label: norm_ex1
  
    \lVert e_{\mathrm{norm,\,global}} \rVert_{2} = \sqrt{
       \frac{ ( \lVert e_{\mathrm{Pt}} \rVert_2)^2 + \
    (\lVert e_{\mathrm{St}} \rVert_2)^2 }{\
        (\lVert e_{\mathrm{0,Pt}} \rVert_2)^2 + \ 
    (\lVert e_{\mathrm{0,St}} \rVert_2)^2 } \
    +  \frac{ ( \lVert e_{\mathrm{Ps}} \rVert_2)^2 + \
    (\lVert e_{\mathrm{Ss}} \rVert_2)^2 }{\
     (\lVert e_{\mathrm{0,Ps}} \rVert_2)^2 + \ 
    (\lVert e_{\mathrm{0,Ss}} \rVert_2)^2 } \
    }

    
**Example 2:** Fitting waveforms of P waves and static surface displacements
to solve for a source model 
    
    Let's say we use P waveforms in the time domain 
    :math:`{\bf d}_{obs,\mathrm{Pt}}`. We combine the waveform
    misfit defined in Equation :eq:`wms_wns` with the misfit of the 
    maximum waveform defined in Equation :eq:`cor`
    correlation. Furthermore we use InSAR-measured
    static surface displacements  :math:`{\bf d}_{obs,\mathrm{insar}}` and 
    GNSS-measured static surface displacements 
    :math:`{\bf d}_{obs,\mathrm{gnss}}`.
    The static surface displacement misfit is defined as in 
    Equation :eq:`wms_wns`. 
    
    The waveform misfits and the correlations, even if the same weights are
    applied, are measures of a different nature. Also the dynamic waveforms
    and the static near-field displacements have different relationships to
    the source parameters. Different normalisation is meaningful. The static
    surface displacement data themselves should be comparable, even though
    InSAR and GNSS positing are very different measuring techniques. 
    
    The **global misfit** in this example is then:
    
.. math::
  :label: norm_ex2
  
    \lVert e_{\mathrm{norm,\,global}} \rVert_{2} = \sqrt{
       \frac{ ( \lVert e_{\mathrm{Pt}} \rVert_2)^2}{\
        (\lVert e_{\mathrm{0,Pt}} \rVert_2)^2 } \
    +  \frac{ ( \lVert e_{\mathrm{Ptcor}} \rVert_2)^2 }{\
     (\lVert e_{\mathrm{0,Ptcor}} \rVert_2)^2  } \
      +  \frac{ ( \lVert e_{\mathrm{insar}} \rVert_2)^2 + \
    (\lVert e_{\mathrm{gnss}} \rVert_2)^2 }{\
     (\lVert e_{\mathrm{0,insar}} \rVert_2)^2 + \ 
    (\lVert e_{\mathrm{0,gnss}} \rVert_2)^2 } \
    }   


The bootstrap method
--------------------

`Bootstrapping` in Grond (see also `Bootstrapping in wikipedia`_)  enables to 
suppress some types of bias in the 
optimization results. Observations that are affected by signals other than 
from the analysed source process often show a high misfits. Also observations
for which the Green's functions based on a medium model, which is at this 
particular site not a good approximation of the underground, can result in 
high misfit values. Already a few high misfit values may pull the optimisation 
to a biased optimum. With bootstrapping we can further estimate model 
parameter uncertainties in an efficient way, which include the propagation of
the data error, but also modelling errors are assessed to some extent.  

In Grond the bootstrapping is applied in a 
number of parallel `bootstrapping chains` where individual bootstrap weights
or bootstrap noise is applied to the model misfits. Basically, individual 
optimization are carried out in each bootstrap chain. Find more below for the 
`BABO Optimiser`.

In Grond **two** different bootstrapping types are implemented. There is 
bootstrapping realised through misfit weights, called `Classic` and `Bayesian
bootstrapping`, and there is bootstrapping realised adding noise to the 
residuals, which is the so-called  `Residual bootstrapping` 
(Fig. 1).

Classic and Bayesian bootstrap
..............................

These bootstrap types are based on weighting. We 
divert from the physics-related and noise-related target weights and create
additional random weight factors for each target. Virtually equal weights 
of 1 for each target are redistributed to new random weights, which add up
to equal the number of targets. In this way the 
final misfit values are comparable even without normalisation.
   
**Classic weights**:
    For `classic` bootstrap weights we draw :math:`N_{\mathrm{targets}}` 
    random integer numbers 
    :math:`{\bf r} \, \epsilon \, [0 \,\, N_{\mathrm{targets}}]`
    from a uniform distribution (Fig. 2, left). 
    We then sort these in :math:`N_{\mathrm{targets}}` bins (Fig. 2, right).
    The frequency in each bin forms the bootstrap target weights.


.. figure:: ../images/classic_bootstrap_weights.svg
    :name: Fig. 3
    :width: 1600px
    :align: center
    :alt: alternate text
    :figclass: align-center
    
    **Figure 3**: Formation of `classical` bootstrap weights. Uniformly random
    samples (left) and the corresponding histogram (right) with the frequencies
    being used as bootstrap weights.  

**Bayesian weights**
    For `Bayesian` bootstrap weights we draw :math:`N_{\mathrm{targets}}+1` 
    random real numbers :math:`{\bf r} \, \epsilon \, [0 \,\, N_{\mathrm{targets}}]`
    from a uniform distribution (Fig. 4, left). 
    We then sort the obtained random values in an ascending order (Fig. 4, 
    middle) 
    and calculate the bootstrap weights as the differences 
    :math:`w_{\mathrm{bootstr},\,i}=r_{i+1}-r_i`.

.. figure:: ../images/bayesian_bootstrap_weights.svg
    :name: Fig. 4
    :width: 1600px
    :align: center
    :alt: alternate text
    :figclass: align-center

    **Figure 4**: Formation of `Bayesian` bootstrap weights. Uniformly random
    samples (left) are sorted (middle) and the differences of neighbouring 
    points (right) are being used as bootstrap weights.  
    
Residual bootstrap
..................
    
Residual bootstrap actually is a computationally more efficient version of the 
`Randomize-then-optimize`_ procedure. The name of the latter method describes
the procedure - with empirical estimates of the data 
error statistics individual realisations of synthetic correlated random noise 
are added to the data for many slightly differing optimisations (Fig. 5). 
Source 
parameter distributions retrieved with the `Randomize-then-optimize`_ method 
based on the data error variance-covariance matrix have been shown to match the 
model parameter distributions obtained from `Marcov Chain Monte Carlo` sampling
of the model spaces by `Jonsson et al.`_ (2014).
In our `residual bootstrap` we add such individual realisations of synthetic 
correlated random noise (Fig. 5C) to the misfits to evaluate individual 
`global misfits`
(Fig. 1). Like this we save the calculation of many forward models compared to 
`Randomize-then-optimize`_, while obtaining the same result.

To generate random noise we use functions of the `kite`_ module. From the 
noise estimation region defined in the `kite`_ scenes (Fig. 5A), the 
noise power spectrum
is used directly with a randomised phase spectrum to create new random noise
with common characteristics in the spatial domain (Fig. 5B). The noise is 
then subsampled
exactly like the data to be used on the model residuals (Fig. 5C).

.. figure:: ../images/illu_residual_bootstrap_realisation.svg
    :name: Fig. 5
    :width: 1400px
    :align: center
    :alt: alternate text
    :figclass: align-center

    **Figure 5**: Residual bootstrap realisation in grond. From data noise (A)
    we synthesise random correlated data noise (B), which is then subsampled
    like the data (C) to be added to the residuals.  


Optimisation 
------------

Grond is open for many different optimisation schemes. So far implemented is 
the so-called `Bayesian Bootstrap Optimisation` (BABO). The `Optimiser` defines
the particular objective function or objective functions and options for them. 
The optimiser also defines the model space sampling schemes. Multiple objective
functions are realized in parallel running optimisation chains. So far these
are the bootstrap chains (see below).

The BABO optimiser
..................

BABO stands for `Bayesian Bootstrap Optimisation` that is done if the 
optimiser is configured to the full extent. As the name says, BABO allows for 
a source optimisation while providing the full information in the results for 
a fully Bayesian analysis. BABO is based on `Direct Search`, meaning model
parameters are drawn in a randomised way from the defined model space 
and synthetic data are then calculated to be compared with the observed data. 
This needs no assumptions on the topology of
the misfit space and is appropriate also for highly non-linear problems.

BABO can turn into a simple Monte-Carlo random direct search if some options 
are switched off. It can also resemble a simulated annealing optimisation 
approach using a certain problem configuration. Last but not least BABO
enables fully probabilistic bootstrapping of the optimisation results. This is 
realised in parallel with optimisation chains to which bootstrapping weights
are applied.

Note:
*Weights* are explained above. The specific
weighting is configured with the `targets`_ used and also with the `problem`_.
The *model space* in which the optimisation takes place is 
defined with the `problem`_.
Here described is the sampling and in the context of the multiple objective 
functions given by the bootstrapping.


Sampling scheme and sampling phases
...................................

Like in any `direct search` optimisation models are drawn from the model space.
From all visited and evaluated models we form and keep a so-called `highscore` 
list. The sampling is set up to progressively converge to the low-misfit 
regions efficiently.
However, for multi-modal model parameters distributions an 
efficient sampling can loose sight of multiple minima with significantly
low misfits. In Grond we can use measures to nurse these multiples.   

**highscore list**: 
    This list contains a defined number of the current best (lowest misfit)
    models. It is continuously updated. The `highscore` list length 
    :math:`L_{hs}` (i.e. number of member models) is `problem`_ dependend:
    :math:`L_{hs} = f_{\mathrm{len}} \cdot (N_{\mathrm{par}} -1)`, 
    with
    :math:`N_{\mathrm{par}}` being the number of model paramters.
    :math:`f_{\mathrm{len}}` is configurable
    (``chain_length_factor``, default is 8).

There are three sampling phases defined, based on which models are drawn from
the model space:

* ``UniformSamplerPhase`` - models are drawn randomly
* ``InjectionSamplerPhase`` - allows to inject specific models 
* ``DirectedSamplerPhase`` - existing low-misfit models `direct` the sampling

.. figure:: ../images/illu_sampling_phases.svg
    :name: Fig. 6
    :height: 300px
    :align: center
    :alt: alternate text
    :figclass: align-center

    **Figure 7**: Sketch of model parameter sampling 
    
    
**UniformSamplerPhase**:
    This is a starting sampler phase of the optimisation. A configurable number
    of models are drawn 
    randomly from the entire model space based on a uniform distribution.

**InjectionSamplerPhase**:
    This is a starting sampler phase of the 
    optimisation in case it should not start blind. It allows to inject 
    specific models at the start of the optimisation. These models could 
    stem from a previous optimisation.

**DirectedSamplerPhase**: 
    This sampler phase follows any starting phase. Using the positions and/or
    the distribution of the
    current `highscore` models the `directed` sampler draws a configurable 
    number of new models. 
    Like this convergence to low-misfit regions is enabled. There are quite 
    some noteworthy details to this sampler phase.
    
    **sampling distributions**: For drawing new models normal distributions
    are used. The standard deviations for the sampler are derived from the 
    `highscore` model parameter standard deviations by using a configurable 
    value (`scatter scale`, see below). Optionally, the covariance of model 
    parameter distributions is
    taken into account by configuring a ``multivariate_normal`` sampler
    distribtion instead of a ``normal`` sampler distribution. 
    The center points for the sampling distribution is configurable to be 
    the ``mean`` of the `highscore`` model parameter distributions, 
    to a ``random`` model of the `highscore` models or an 
    ``excentricity_compensated`` draw (see below). 
    
    **scatter scale**: This scale defines the search radius around the current
    `highscore` models. With a scatter scale of 2 the search for new models
    has a distribution with twice the standard deviation as estimated for the 
    current `highscore` models. It is possible to define a beginning scatter
    scale and an ending scatter scale. When defining a larger value for the 
    beginning scatter scale and a smaller value for the ending scatter scale,
    during the progressing optimisation, the search gets more and more 
    confined. In other words, the sampling evolves from being more explorative 
    to being more exploitive.

    **excentricity compensation**: This method applies to the center value of 
    the sampler distribution. Taking this option, the center point of the 
    sampler distribution is with an increased likelihood a `highscore` member 
    model off-center to the `highscore` model mean value compared to a random
    choice. The probability of drawing a model from the 
    `highscore` list is derived from distances the `highscore` models have
    to other `highscore` models in the model parameter space. 
    Excentricity is therefore compensated, because models with few neighbours 
    at larger distances have an increased likelihood to be drawn. 
    
    What's the use? Convergence is slowed down, yes, but to the benefit of 
    low-misfit region represented by only a few models drawn up 
    to the current point. 
    
    Let's say there are two separated groups of 
    low-misfit models in our `highscore` list, with one group forming the 75%
    majority. 
    In the directed sampler phase the choices of a mean center point
    for the distribution as well as a random starting point for the sampler 
    distribution would favour new samples in the region of the 
    `highscore` model majority. Models in the low-misfit region may be dying
    out in the `highscore` list due to favorism and related sparse sampling.
    `excentricity compensations` can help is these cases and keep models with 
    not significantly higher misfits in the game and in sight.
    
    TODO: correct? too many explanations? Sebastian,
    here is the perfect place for one of your movies.
 

Bootstrap chains
................

A `bootstrap chain` is set up with individual target bootstrap weights and/or 
target bootstrap residuals (Fig. 7A). Therefore each bootstrap chain has 
an individual objective function. With one 
forward model :math:`N_{\mathrm{bootstrap}}` 
different `global 
misfits` are calculated (Fig. 7B). Like this for each bootstrap chain we can 
run an individual optimisation, even though all bootstrap chains share the same 
forward models. 

The highscore list member models in each bootstrap chain (Fig. 7B) will differ 
to some
extent and therefore different bootstrap chains may converge to different 
places within the model space (Fig. 7C, Fig. 8). These differences mark the 
uncertainty of the models with respect to data errors.

.. figure:: ../images/illu_bootstrap_weights.svg
    :name: Fig. 7
    :height: 400px
    :align: center
    :alt: alternate text
    :figclass: align-center
    
    **Figure 7**:  Bootstrap chain graph. (A) Illustration of bootstrap 
    weights, (B) bootstrap chain highscore lists and  (C) their influence 
    on the convergence in the model parameter space due to the 
    individual objective function of each bootstrap chain.

The convergence of model parameters for the models within each bootstrap chain 
is dependent on the settings of the optimisation, e.g. the setup of parameter
bounds, `scatter scale` settings of the `directive sampling phase` and else.
With very `exploitive` settings convergence can be forced. However, if the 
convergence within each bootstrap chain starts to form individual solar systems
in the model space, further optimisation will not provide significantly better
models. In Fig. 8 the area of the `highscore` models of the three bootstrap
chains has only little overlap compared to an earlier stage visualised in 
Fig. 7C.



.. figure:: ../images/illu_babo_chains.svg
    :name: Fig. 8
    :height: 300px
    :align: left
    :alt: alternate text
    :figclass: align-left
    
    **Figure 8**: Drawing new candidate models based on the existing solution 
    space. (...)

    
    
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
.. _Bootstrapping in wikipedia: https://en.wikipedia.org/wiki/Bootstrapping_(statistics)
.. _Randomize-then-optimize: https://epubs.siam.org/doi/abs/10.1137/140964023
.. _Jonsson et al.: http://adsabs.harvard.edu/abs/2014AGUFM.S51C..05J

.. _dataset: ../dataset/index.html
.. _targets: ../targets/index.html
.. _problem: problems/index.html
