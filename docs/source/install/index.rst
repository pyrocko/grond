Installation
============

Grond is compatible with Python 2.7 as well as Python 3.4 and above.

.. _manual-install:

Manual installation (from source)
---------------------------------

First install Grond's prerequisites. It depends on

* `Pyrocko`_ - required
* `Kite`_ - optional, needed for InSAR data modelling

Next, get Grond's source code using :program:`git clone`:

.. code-block :: sh

    git clone https://git.pyrocko.org/pyrocko/grond.git
    cd grond

Finally, decide if you want to install Grond system-wide or just for a single
user.

System wide installation::

    sudo python3 setup.py install


Installation in your home directory::

    python3 setup.py install --user


Updating a manual installation
------------------------------

For updating an existing manual installation of Grond, update the source code
using :program:`git pull`:

.. code-block :: sh

    cd grond  # change to the directory to where you cloned Grond initially
    git pull origin master

Then just reinstall as described in section :ref:`manual-install` above.


Installation under Anaconda
---------------------------

Pyrocko's pre-built Anaconda packages are available online, Grond can be
installed with Anaconda's pip into the Anaconda environment.


.. code-block :: sh

    conda install -c pyrocko pyrocko
    pip install grond


Installation with pip
---------------------

Grond and all its dependencies can be installed by running

.. code-block:: bash

   pip install grond  # read below...

**but**, we recommend to make a conscious decision about how its main
dependency `Pyrocko`_ and especially Pyrocko's own dependencies are installed.
The `Pyrocko Installation Manual <https://pyrocko.org/docs/current/install/>`_
describes different installation schemes.

As a general advice, we recommend to exclusively use either, (1) the system's
native package manager, (2) Anaconda, or (3) pip only. In (1) and (2), only
resort to use pip for those few packages which are not available as native
packages. Otherwise, competing package managers will ruin your day!

To prevent pip from automatically resolving dependencies run

.. code-block:: bash

   pip install --no-deps grond

This assumes that `Pyrocko`_ and `Kite`_ have been installed beforehand.

.. _kite: https://pyrocko.org/kite/
.. _pyrocko: https://pyrocko.org/docs/current/install/
