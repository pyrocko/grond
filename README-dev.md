# Grond development environment and contribution guide

## Language

Grond is written in the Python programming language (versions ==2.7 and >=3.4).

## Deployment

Grond uses Setuptools for its installation script. See ``setup.py`` and
``setup.cfg`` in the project root directory.

## Versioning and releases

Git is used for version control. Use development branches for new features.
Master branch should always point to a stable version.

*Commit message conventions:*

* start with lower case
* colon-prepend affected component 
* try to use imperative form
* examples:
  - ``docs: add section about weighting``
  - ``waveform targets: correct typo in component names``
  - ``waveform targets: fix issues with misaligned traces``

## Testing

* TODO

## Code style

Grond source code must follow the PEP8 coding standards. It must pass the
code style check provided by the ``flake8`` tool.

Additionally,

* use i/n convention for indices and counts
  - e.g. ``for istation in range(nstations):``
* do not abbreviate words unless this would result in ridiculously long names
* use British english, e.g.
  - 'modelling' rather than 'modeling'
  - 'analyser' rather than 'analyzer'
  - 'optimiser' rather than 'optimizer'
* log messages: TODO
* docstrings: TODO

## Documentation

Grond's documentation is built using the ``Sphinx`` tool. See the ``docs``
in the project root directory. Build with ``make html`` in ``docs``.

*Text style rules:*

* titles: only capitalize first word
* use British english
