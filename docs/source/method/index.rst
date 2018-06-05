Method
======

TODO: (general) this section should be as self-contained as possible, describe 
the method in general - give references to other sections how things are
implemented in Grond.

The very core of the optimisation is the evaluation of a misfit between 
observed and predicted data. This is most often the difference  
:math:`{\bf d}_{obs} - {\bf d}_{synth}`, but can also be another comparision,
like a correlation measure for examples.

This sheet describes the method on:

1. what exactly are the observed data :math:`{\bf d}_{obs}` and the synthetic 
   data :math:`{\bf d}_{synth}`
2. how handles Grond the differences between :math:`{\bf d}_{obs}` and
   :math:`{\bf d}_{synth}` 
   with respect to defining the objective functions through misfits and data
   weighting,
3. how the optimisation is set up to search the model space to find the 
   optimum models and 
4. which methods are used to estimate model uncertainties.

`Observed data` here means full waveforms that are tapered to the defined 
phases, restituted and filtered. `Synthetic waveforms` are the forward-
modelled waveforms that are tapered and filtered in the same way as the 
observed waveforms. 


**Tapering**:
    Tapering is done on the waveform targets and means that parts of the 
    recorded waveforms are cut out. Here the taper applies to specific seismic
    phases in both the raw synthetic waveforms calculated 
    :math:`{\bf d}_{raw, synth}` and 
    the recorded, restituted waveforms. Only these parts are then used in the 
    misfit calculation. 
    The taper window duration is configured ( target config link) as well as 
    the time. 

    TODO: set link (line above)

    The tapering is source-model dependent, since the tapering time is given 
    with respect to the theoretic phase arrival
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


From the difference :math:`{\bf d}_{obs} - {\bf d}_{synth}` or any other
quantitative comparison (e.g. correlation) the `misfit` is defined based
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
minimum of the misfit function optimum model.

The objective function defines what a `model fit` is and how `good` or
`poor` models are scaled with respect to others. Furthermore, the
objective function has rules how different data sets are handled and how data 
errors are considered in optimisations. 


.. figure:: ../images/illu_combi_weights.svg
    :name: Fig. 1
    :height: 500px
    :align: center
    :alt: alternate text
    
    **Figure 1**: Overview illustration of the GROND objective function. 
    Details on how each function and weight vectors are formed follow below.

TODO: there is nothing about the L-norms in the grafic, Is that misleading or 
formation the sake of generalization ok or is there a nice way for that?  
    
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
predicitions that manage to explain parts of the observed data holds
:math:`\lVert e_{\mathrm{norm}} \rVert_x <1`. Furthermore, the data norm 
:math:`\lVert e_{\mathrm{0}} \rVert_x` is used in the normalisation of data
groups.

For waveform data correlation the misift function is based on the maximum
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
  
**Target balancing weights**

.. figure:: ../images/illu_target_balancing.svg
    :name: Fig. 2
    :width: 300px
    :align: left
    :alt: alternate text
    :figclass: align-center
    
    **Figure 2**: Qualitative sketch how target balancing weight increasse with 
    source distance to balance amplitude decrease caused by geometrical 
    spreading. 

With these weights waveforms are `balanced` with respect to the expected signal
amplitude. 
Signal amplitudes in a trace :math:`|{\bf{d}}_{synth}|` depend on the source-receiver 
distance, on the phase type and the taper used. The problem tackled 
with this weight is that
large signal amplitude have higher contributions to the misfit than smaller
signal amplitudes, without carrying more information. From synthetic 
waveforms of `N` forward models that have been randomly drawn from the defined 
model space the mean signal amplitude of the traces is derived. The weight 
for each trace is simply the inverse of these mean signal amplitudes:

.. math::
  :label: wtba
     
    {\bf w}_{\mathrm{tba}} = 1/ \lVert {\bf{d}}_{synth}  \rVert_x  = \
        (\sum^{N}{|{d}_{i, synth}|^x})^{\frac{1}{x}}.


Like this small 
signal are enhanced in the
objective function and large signals supressed. This is described as 
`adaptive station weighting` in the PhD `thesis by Heimann`_ (2011) (page 23).
In Grond they are called ``balancing weights`` and are received from the
``TargetBalancingAnalyser`` before the optimisation.
  
**Data weights based on data error statistics**

There are direct data weight vectors :math:`\bf{w}` or weight matrices
:math:`\bf{W}` based on empirical data error variance estimates. Partly, e.g. 
for InSAR and GNSS data, these can include data error 
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

**manual data weighting**

User-defined manual data weights enable an arbitrary weighting of data sets in contrast to balancing of single observations through target balancing and 
noise-based data weights. 
No rules apply other than from the user's rationale. In Grond they are called 
``manual_weight`` and are given in the configuration file.


TODO link to the target sheet

**Normalisation of data and data groups**

The normalisation in Grond is applied to data groups that are member of the so
called ``normalisation_family``. A `normalisation family` in Grond can be 
composed in many ways. However, it is often meaningful to put data of the same 
kind and with similar weights into the same `normalisation family`. This 
could be P and S waves, or two InSAR data sets. As an explanation some 
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

    
*Example 2:* Fitting waveforms of P waves and static surface displacements
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
surpress some types of bias in the 
optimization results. Observations that are affected by signals other than 
from the analysed source process often show a high misfits. Also observations
for which the Green's functions based on a medium model, which is at this 
particular site not a good approximation of the underground, can result in 
high misfit values. Already a few high misfit values may pull the optimisation 
to a biased optimium. With bootstrapping we can further estimate model 
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
   
**Classic weights**

For `classic` bootstrap weights we draw :math:`N_{\mathrm{targets}}` 
random integer numbers :math:`{\bf r} \, \epsilon \, [0 \quad N_{\mathrm{targets}}]`
from a uniform distribution (Fig. 2, left). 
We then sort these in :math:`N_{\mathrm{targets}}` bins (Fig. 2, right).
The frequency in each bin forms the bootstrap target weights.


.. figure:: ../images/classic_bootstrap_weights.svg
    :name: Fig. 3
    :width: 1400px
    :align: center
    :alt: alternate text
    :figclass: align-center
    
    **Figure 3**: Formation of `classical` bootstrap weights. Uniformely random
    samples (left) and the corresponding histogram (right) with the frequencies
    being used as bootstrap weights.  

**Bayesian weights**

For `Bayesian` bootstrap weights we draw :math:`N_{\mathrm{targets}}+1` 
random real numbers :math:`{\bf r} \, \epsilon \, [0 N_{\mathrm{targets}}]`
from a uniform distribution (Fig. 4, left). 
We then sort the obtained random values in an ascending order (Fig. 4, middle) 
and calculate the bootstrap weights as the differences 
:math:`w_{\mathrm{bootstr},\,i}=r_{i+1}-r_i`.

.. figure:: ../images/bayesian_bootstrap_weights.svg
    :name: Fig. 4
    :width: 1400px
    :align: center
    :alt: alternate text
    :figclass: align-center

    **Figure 4**: Formation of `Bayesian` bootstrap weights. Uniformely random
    samples (left) are sorted (middle) and the differences of neighbouring 
    points (right) are being used as bootstrap weights.  
    
Residual bootstrap
..................
    
Residual bootstrap actually is a computationally more efficient version of the 
`Randomize-then-optimize`_ procedure. The name of the latter method describes
the procedure - with empirical estimates of the data 
error statistics individual realisations of synthetic correlated random noise 
are added to the data for many slightly differing optimisations. Source 
parameter distributions retrieved with the `Randomize-then-optimize`_ method 
based on the data error variance-covariance matrix have been shown to match the 
model parameter distributions obtained from Marcov Chain Monte Carlo sampling
of the model spaces by `Jonsson et al.`_ (2014).
In `residual bootstrap` we add such individual realisations of synthetic 
correlated random noise to the misfits to evaluate individual `global misfits`
(Fig. 1). Like this we save the calculation of many forward models compared to 
`Randomize-then-optimize`_, while obtaining the same result.

To generate random noise we use functions of the `kite`_ module. From the 
noise estimation region defined in the `kite`_ scenes, the noise power spectrum
is used directly with a randomised phase spectrum to create new random noise
with common characteristics in the spatial domain. The noise is then subsampled
exactly like the data to be used on the model residuals.

TODO: Fig.5 residual bootstrap

Optimisation 
------------

high score 


The BABO optimiser
..................

 
.. figure:: ../images/illu_bootstrap_weights.svg
    :name: Fig. 6
    :height: 300px
    :align: center
    :alt: alternate text
    :figclass: align-center
    
    **Figure 6**:  Illustration of bootrap weights and their influence on the 
    convergence in the model parameter space due to the differing misfit 
    function for each bootstrap chain.


One history of smapling the model space. N misfit spaces

All bootstrap chains share the misfits of the 
same forward models. The misfit is re-evaluated without the sampling of new
models and new forward modelling. This makes the bootstrapping setup in Grond computationally very effcient.
 
.. figure:: ../images/illu_babo_chains.svg
    :name: Fig. 7
    :height: 300px
    :align: center
    :alt: alternate text
    :figclass: align-center
    
    **Figure 7**: Drawing new candiate models based on the existing solution 
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
