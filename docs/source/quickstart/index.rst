.. role:: bash(code)
   :language: bash

Quickstart
==========

This is a type-along introductory example to get your feet wet in operating
Grond.

We will estimate a centroid moment tensor (CMT) solution given a set of
simulated regional full waveform seismic records. Such synthetic tests are
useful e.g. for resolution analysis.

A synthetic test
----------------

With the following few commands, we will let Grond

* forward model a fully synthetic scenario of regional seismic waveforms.
* setup a processing environment (project folder).
* perform a probabilistic optimisation to retrieve a CMT solution for the scenario event.
* visualise the results and various details about the optimisation.

**Online help**

.. code-block :: sh

    grond --help           # get list of available subcommands
    grond scenario --help  # all subcommands come with built-in help

**Forward-model a random scenario and create project folder**

.. code-block :: sh

    grond scenario --targets=waveforms my_first_project

*Note: precomputed Green's functions (GF) needed to run this example will be
downloaded from the internet.*

**Checking data setup and configuration**

.. code-block :: sh

    cd my_first_project
    grond check config/scenario.gronf

**Starting optimisation**

.. code-block :: sh

    grond go config/scenario.gronf

.. note ::

    Interrupted or prolonged optimisation can be continued with:

    .. code-block :: sh

        grond go config/scenario.gronf

**Continuing an aborted/modified optimisation**

.. code-block :: sh

    grond continue config/scenario.gronf

**Plotting the results in a report**

.. code-block :: sh

    grond report runs/cmt_scenario_ev001.grun

**Open the report in your web browser**

.. code-block :: sh

    grond report -so

`See the generated report online. <https://pyrocko.org/grond/reports/quickstart>`_

Next steps
----------

Read the :doc:`overview section <../overview/index>` to see how to set up a
project folder with your own data or take a tour through our 
:doc:`example projects <../examples/index>`.
