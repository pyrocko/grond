Command line interface
======================

Grond is a normally invoked as a command line program by calling ``grond`` from
the Unix/Linux shell. It's command line interface (CLI) uses the standard
Unix/Linux conventions for its options and arguments.  To get a brief summary
on any of its subcommands, add the ``--help`` option. Options are recognized by
their leading double-dashes, e.g. ``--help``. Some options in Grond have short aliases
with single dashes.

.. command-output:: grond --help

Grond subcommands
-----------------

.. toctree::
    :caption: Subcommands:

    scenario
    init
    events
    check
    go
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
