Configuration file structure and format (YAML)
----------------------------------------------

Grond is configured with plain text files in `YAML`_ format. The YAML format
has been chosen because it can represent arbitrarily nested data structures
built from mappings, lists, and scalar values. It also provides an excellent
balance between human and machine readability. When working with such files, it
is good to know that the **indentation is part of the syntax** and that
comments can be introduced with the ``#`` symbol. The configuration file can be
made to work with multiple events. A basic configuration template can be
generated using the :doc:`grond init </commands/init>` subcommand (see section
:ref:`project-init`).

The following commented listing presents the overall structure of a Grond
configuration file. It has a top level container (mapping), introduced with the
line ``--- !grond.Config`` with several child elements: ``path_prefix``,
``rundir_template``, ``dataset_config``, etc. Some of these entries may again
contain their own child elements (indented blocks of lines) or lists (lines
introduced with dashes). The type markers, like e.g. ``!grond.DatasetConfig``,
select the Grond object type of the following mapping and their documentation
can likely be found in the :doc:`/library/index`.

.. code-block :: sh

    %YAML 1.1
    --- !grond.Config

    # All file paths referenced below are treated relative to the location of this
    # configuration file, here we may give a common prefix. E.g. setting it to '..'
    # if the configuration file is in the sub-directory '${project_root}/config'
    # allows us to give the paths below relative to '${project_root}'.

    path_prefix: '..'

    # Path, where to store output (run directories). The placeholder
    # '${problem_name}' will be expanded to a name configured below in
    # problem_config.name_template and will typically include a config identifier
    # and the event name.

    rundir_template: 'runs/${problem_name}.run'

    # Configuration section for dataset (input data)

    dataset_config: !grond.DatasetConfig
      ...

    # Configuration section for the forward modelling engine (configures where
    # to look for GF stores)

    engine_config: !grond.EngineConfig
      ...

    # Configuration section selecting data to be included in the data optimisation. 
    # Amongst other parameters, the objective function for the optimisation is 
    # defined for each target group. The targets can be composed of one or more 
    # contributions, each represented by a !grond.TargetConfig section.

    target_groups:

    - !grond.WaveformTargetGroup       # Setup for seismic waveforms
      ...

    - !grond.SatelliteTargetGroup      # Setup for InSAR
      ...

    - !grond.GNSSCampaignTargetGroup   # Setup for coseismic GNSS displacements
      ...

    # Definition of the problem to be solved - source model, parameter space, and
    # global misfit configuration settings. Only one problem can be defined per 
    # configuration file.

    #problem_config: !grond.RectangularProblemConfig  # setup for an extended source
    #problem_config: !grond.DoubleDCProblemConfig     # setup for combination of two double-couples
    problem_config: !grond.CMTProblemConfig           # setup for a general moment tensor

      # Name used to identify the output the placeholder '${event_name}' will
      # be replaced with the current event name

      name_template: 'cmt_${event_name}'
      ...
      
    # Configuration of pre-optimisation analysis phase; e.g. balancing weights are
    # determined during this phase. Analysers can be combined.

    analyser_configs: 


    - !grond.TargetBalancingAnalyserConfig  # balancing weights
      ...

    - !grond.NoiseAnalyserConfig            # pre-event noise based weights
      ...

    # Configuration of the optimisation procedure.

    optimiser_config: !grond.HighScoreOptimiserConfig
      ...


.. _YAML: https://en.wikipedia.org/wiki/YAML
