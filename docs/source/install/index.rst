Installation
============

To install Grond, first install its prerequisites

* `Pyrocko`_ - required
* `Kite`_ - optional, needed for InSAR data modelling

Next, install the Grond library and command line program with

.. code-block :: sh

    git clone https://github.com/pyrocko/grond.git
    cd grond
    sudo python3 setup.py install


Updating
--------

For updating an existing installation of Grond:

.. code-block :: sh

    cd grond  # change to the directory to where you cloned grond initially
    git pull origin master
    sudo python3 setup.py install



.. _kite: https://pyrocko.org/docs/kite/current/
.. _pyrocko: https://pyrocko.org/docs/current/install/
