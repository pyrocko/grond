Command line interface
======================

Grond is a normally invoked as a command line program by calling ``grond`` from
the Unix/Linux shell. Its command line interface (CLI) uses standard Unix/Linux
conventions for its options and arguments. To get a brief summary on any Grond
subcommand, add the ``--help`` option. Options are recognized by their leading
double-dashes (``--help``). Some options have single character aliases names
accessible with single dash notation (``-h``). In the documentation,
placeholders for required arguments are denoted with angle brackets
(``<placeholder>``), placeholders for optional arguments with square brackets,
e.g. (``[placeholder]``).

.. command-output:: grond --help

The following pages list the self-documentation strings of all Grond
subcommands.

.. toctree::
    :maxdepth: 1

    scenario
    init
    events
    check
    go
    continue
    forward
    harvest
    plot
    movie
    export
    report
    diff
    qc-polarization
    upgrade-config
    version
