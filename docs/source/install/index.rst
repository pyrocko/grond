Installation
============

Grond is compatible with Python 2.7 as well as Python 3.4 and above.

Installation from source (manual install)
-----------------------------------------

First install Grond's prerequisites. It depends on

* `Pyrocko`_ - required
* `Kite`_ - optional, needed for InSAR data modelling

Next, get Grond's source code using `git clone`:

.. code-block :: sh

    git clone https://github.com/pyrocko/grond.git
    cd grond

Finally, decide if you want to install system-wide or just for a single user.

System wide installation::

    sudo python3 setup.py install


Installation in your home directory::

    python3 setup.py install --user


Updating a manual installation
------------------------------

For updating an existing manual installation of Grond, update the source code
using 'git pull':

.. code-block :: sh

    cd grond  # change to the directory to where you cloned Grond initially
    git pull origin master

Then just over-install as described in the section above.


Installation under Anaconda
---------------------------

Pyrocko's pre-built Anaconda packages are available online, Grond can be
installed with Anaconda's pip into the Anaconda environment.


.. code-block :: sh
    conda install -c pyrocko pyrocko
    pip install git@https://github.com/pyrocko/grond.git


Installation with pip
---------------------

If you want to install Grond's prerequisite Pyrocko with pip, carefully read
`Installation with pip`_ in the Pyrocko manual first. Grond's prerequisites
must be installed beforehand, we do not allow pip to automatically resolve
dependencies.

.. code-block :: sh

    pip3 install git@https://github.com/pyrocko/grond.git


.. _kite: https://pyrocko.org/docs/kite/current/
.. _pyrocko: https://pyrocko.org/docs/current/install/
.. _Installation with pip: https://pyrocko.org/docs/current/install/packages/pip.html
