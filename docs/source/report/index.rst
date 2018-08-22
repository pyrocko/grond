Report
======

Grond's reports are presented in interactive HTML webpages where you can browse and compare different events and inversion runs.

.. figure :: ../images/report_webpage.png
    :name: Grond Webpage
    :width: 80%
    :alt: Grond report webpage
    :figclass: align-center

    **Figure 1**: Example of a Grond report, here the bootstrapped misfit evolution of a static inversion is shown. Use the left navigation panel to navigate the plots.

Online example report
---------------------

Explore the `online example reports <https://pyrocko.org/grond/reports>`_ to see what information the inversion reveals.

Generate reports
----------------

When an inversion is finished and with a set-up project structure you can create and open the report with:

.. code-block:: sh
    
    grond report -so <rundir>

The flags ``-s`` will spin up a webserver and ``-o`` will open the browser; more information is given from ``grond report --help``.


Plot types
----------

To see which plots are available for a particular configuration, check out the subcommand

.. code :: bash

    grond plot list <rundir>
