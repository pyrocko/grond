Configuration
=============

Grond is configured by YAML text files. Read on :doc:`structure` to get started with the configuration.

.. toctree::
    :caption: Contents
    :name: Contents

    structure
    dataset/index
    targets/index
    problems/index
    analysers/index
    optimisers/index

.. rubric:: Glossary

.. glossary ::

    ``Dataset``
        The core of data input, selection and data configuration.

    ``Targets``
        are the observable (seismic or surface displacements) whos misfit is minimised.


    ``Problem``
        is specific source model (CMT, rectangular fault), which is optimised.


    ``Analysers``
        are responsible for data analysis before the optimisation.

    ``Optimiser``
        deliver the optimisation strategy.
