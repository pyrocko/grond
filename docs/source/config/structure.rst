Configuration file structure and format (YAML)
==============================================

Grond is configured with plain text files in `YAML`_ format. The YAML format has been chosen because it can represent arbitrarily nested data structures built from mappings, lists, and scalar values. It also provides an excellent balance between human and machine readability. When working with such files, it is good to know that the **indentation is part of the syntax** and that comments can be introduced with the ``#`` symbol.

The configuration file can be made to work with multiple events. A basic configuration of some example types can be generated using the subcommand (see section :ref:`project-init`)


.. code-block:: sh

    grond init <example>

The following commented configuration presents the overall structure of a Grond configuration file. It has a top level container (mapping), introduced with the line ``--- !grond.Config`` with several child elements: ``path_prefix``, ``rundir_template``, ``dataset_config``, etc. Some of these entries may again contain their own child elements (indented blocks of lines) or lists (lines introduced with dashes). The type markers, like e.g. ``!grond.DatasetConfig``, select the Grond object type of the following mapping and their documentation can likely be found in the :doc:`/library/index`.

.. literalinclude :: ./config_example.yaml
    :caption: Example Grond YAML config

.. _YAML: https://en.wikipedia.org/wiki/YAML
