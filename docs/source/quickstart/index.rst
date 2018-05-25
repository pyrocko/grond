.. role:: bash(code)
   :language: bash

Quickstart
==========

Grond ships with a two options to quickstart a new project structure (see :ref:`project-layout`), including Grond's YAML configuration files.

.. tip::
    
    All subcommands come with a help, e.g. ``grond go --help``!

Option 1: Forward-model a scenario
----------------------------------

The subcommand ``grond scenario`` will forward model observations for a modelled earthquake and create a ready-to-go Grond project. Different observations and source problems can be added by flags - see ``grond scenario --help`` for possible combinations and options.

The scenario can contain the following synthetic observations:

* Seismic waveforms
* InSAR surface displacements
* GNSS surface displacements

.. code-block:: sh
    
    grond scenario <project-folder>

A map of the random scenario is plotted in :file:`scenario_map.pdf`.

Option 2: Initialise an empty project
-------------------------------------

An empty project structure can be created with the subcommand ``grond init``. Different configurations can be added by flags, see ``grond init --help`` for more information.

.. code-block:: sh
    
    grond init <project-folder>

Start the optimisation
----------------------

To start optimising the bootstrapped Grond scenario use the subcommand `go`

.. code-block:: sh

    cd <project-folder>
    grond go config.yml

Optionally we can check the configuration beforehand with

.. code-block:: sh

    grond check config.yml

Custom Grond projects
---------------------

After initialising a new project you can add your own data and customise the configration file.

1. Add your own data to folders :file:`events/<name>/`
2. Add your :file:`events.txt`
3. Copy :file:`config.yml` to :file:`my-event.yml`
4. Customize values in :class:`~grond.dataset.DatasetConfig`
5. Run a ``grond check`` to check your configuration and input data
6. Run ``grond go my-event.yml`` to start the optimisation

See the :doc:`../examples/index` for detailed information.
