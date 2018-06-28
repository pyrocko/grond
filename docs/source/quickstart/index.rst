.. role:: bash(code)
   :language: bash

Quickstart
==========

Grond ships with two options to quickstart a new project structure (see
:ref:`project-layout`), including Grond's YAML configuration files. Here
we present a quickstart solution employing forward-modelled data.

.. tip::
    
    All subcommands come with a help, e.g. ``grond go --help``!

Forward-model a scenario:

.. code-block:: sh

    grond scenario --targets=waveforms my_first_project
    In which cache directory shall the GF store be downloaded to? 
    Default 1, (C)ancel: 1

Check the configuration:

.. code-block:: sh

    cd my_first_project
    grond check config.yml

Start the optimisation:

.. code-block:: sh

    grond go config.yml

Plot the results in a report:

.. code-block:: sh

    grond report <rundir>


Custom Grond projects
---------------------

After initialising a new project you can add your own data and customise the configration file.

1. Run ``grond init <project-folder>`` to initialise an empty project
2. Add your own data to folders :file:`events/<name>/`
3. Add your :file:`events.txt`
4. Copy :file:`config.yml` to :file:`my-event.yml`
5. Customize values in :class:`~grond.dataset.DatasetConfig`
6. Run a ``grond check`` to check your configuration and input data
7. Run ``grond go my-event.yml`` to start the optimisation
8. Run ``grond report <rundir>`` to obtain a report on the results.

See the :doc:`../examples/index` for detailed information.
