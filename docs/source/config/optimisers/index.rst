Optimisers
==========

The optimiser is the heart of the inversion, it samples the model space striving to lower the :term:`target's <target>` misfit.

Bayesian Bootstrap Optimiser (BABO)
-----------------------------------

BABO allows for a source optimisation while providing the full information in the results for a fully Bayesian analysis. BABO is based on Direct Search where model parameters are drawn in a randomised way from the defined model space. Those models are then calculated and compared with the :term:`target's <target>` observed data. This needs no assumptions on the topology of the misfit space and is appropriate also for highly non-linear problems.

BABO can be configured for a simple Monte-Carlo random direct search. It can also resemble a simulated annealing optimisation approach. Ultimately BABO enables fully probabilistic bootstrapping of the optimisation results. 

See the :doc:`/method/index` documentation for detailed information about the :ref:`optimisation <optimisation>`.


Example Configuration
~~~~~~~~~~~~~~~~~~~~~

Basic configuration section for the BABO optimiser:

.. code-block :: yaml

  optimiser_config: !grond.HighScoreOptimiserConfig

  # Number of bootstrap realisations to be tracked simultaneously in the
  # optimisation
  nbootstrap: 100

  # stages of the sampler. Start with uniform sampling of the model space
  # (!grond.UniformSamplerPhase), then narrow down to the interesting regions
  # (!grond.DirectedSamplerPhase).
  sampler_phases:

  - !grond.UniformSamplerPhase

      # Number of iterations
      niterations: 1000

  - !grond.DirectedSamplerPhase

      # Number of iterations
      niterations: 20000

      # Multiplicator for width of sampler distribution at start of this phase
      scatter_scale_begin: 2.0

      # Multiplicator for width of sampler distribution at end of this phase
      scatter_scale_end: 0.5


General configuration
~~~~~~~~~~~~~~~~~~~~~

General parameters are:

``nbootstrap``
  Number of bootstrap realisations to be tracked simultaneously during the optimisation.

``sampler_phase``
  List of sampling stages: Start with uniform sampling of the model model space and narrow down through directed sampling.


``UniformSamplerPhase`` configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
At the beginning of the optimisation this sampler phase explores the solution space uniformly. A configurable number of models are drawn randomly from the entire model space based on a uniform distribution.

``niterations``
    Number of iterations for this phase.

``DirectedSamplerPhase`` configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This sampler is used for the second phase and follows any of starting samplers above: Using existing models of the current `highscore` models the `directed` sampler draws a configurable number of new models. Like this convergence to low-misfit regions is enabled. There are quite some noteworthy configurable details to this sampler phase:

.. glossary ::

  ``niterations``
    Number of iterations for this phase.

  ``sampling_distributions``
    New models are drawn from normal distribution. The standard deviations are derived from the `highscore` models parameter's standard deviation and scaled by ``scatter_scale`` (see below). Optionally, the covariance of model parameters is taken into account by configuring when ``multivariate_normal`` is enabled (default is ``normal`` distribution). The distribution is centred around

      1. ``mean`` of the `highscore` model parameter distributions
      2. a ``random`` model from the `highscore` list or
      3. an ``excentricity_compensated`` draw (see below).

  ``scatter_scale``
    This scales search radius around the current `highscore` models. With a scatter scale of 2 the search for new models has a distribution with twice the standard deviation as estimated from the current `highscore` list. It is possible to define a beginning scatter scale and an ending scatter scale. This leads to a confining directed search. In other words, the sampling evolves from being more explorative to being more exploitive in the end.

  ``starting_point``
    This method tunes to the center value of the sampler distribution: This option, will increase the likelihood to draw a `highscore` member model off-center to the mean value. The probability of drawing a model from the `highscore` list is derived from distances the `highscore` models have to other `highscore` models in the model parameter space. Eccentricity is therefore compensated, because models with few neighbours at larger distances have an increased likelihood to be drawn.

What's the use? Convergence is slowed down, yes, but to the benefit of low-misfit region represented by only a few models drawn up to the current point.

Let's assume there are two separated groups of low-misfit models in our `highscore` list, with one group forming the 75% majority. In the directed sampler phase the choices of a mean center point for the distribution as well as a random starting point for the sampler distribution would favour new samples in the region of the `highscore` model majority. Models in the low-misfit region may be dying out in the `highscore` list due to favour and related sparse sampling. `eccentricity compensations` can help is these cases and keep models with not significantly higher misfits in the game and in sight.

``InjectionSamplerPhase`` configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This starting phase allows to inject pre-defined models at the start of the optimisation. These models could originate from a previous optimisation.

``xs_inject``
    Array with the reference model.

TODO: correct? too many explanations? Sebastian, here is the perfect place for one of your movies.
