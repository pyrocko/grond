Report
======

TODO: describe how to create and configure Grond's reports. Describe what
output plots are available.

Grond's reports are presented in interactive HTML webpages where you can browse and compare different events and inversion runs.

Generate reports
----------------

When an inversion is finished you can create and open the report with:

.. code-block:: sh
    
    grond report -so <rundir>

You can look up the meaning of ``-so`` with ``grond report --help``


Plot types
----------

This section will briefly describe the plots which are generated from the different modules. If a module was not activated during the run, it will not generate any plots.

To see which plots are available, check out ``grond plot list <rundir>``.

Optimizer
.........

Waveform Target
...............

Satellite Target
................

GNSS Campaign Target
....................
