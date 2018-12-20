Installation
============

Grond is compatible with Python 2 and 3.

The installtion depends on:

* `Pyrocko`_ - required, waveform and GNSS modelling
* `Kite`_ - optional, needed for InSAR data modelling

Other dependencies are:

* `NumPy <https://www.numpy.org/>`_
* `SciPy <https://scipy.org/>`_
* `Matplotlib <https://matplotlib.org/>`_

Next, install the Grond library and command line program with

.. code-block :: sh

    git clone https://github.com/pyrocko/grond.git
    cd grond
    sudo python3 setup.py install


With PIP
--------

Installation through pip:

.. code-block :: sh

    pip install pyrocko
    pip install git@https://github.com/pyrocko/grond.git


With Anaconda
-------------

Pyrocko's pre-build Anaconda packages are available online,
grond is installed through pip.


.. code-block :: sh
    conda install -c pyrocko pyrocko
    pip install git@https://github.com/pyrocko/grond.git


Updating
--------

For updating an existing installation of Grond:

.. code-block :: sh

    cd grond  # change to the directory to where you cloned grond initially
    git pull origin master
    sudo python3 setup.py install



.. _kite: https://pyrocko.org/docs/kite/current/
.. _pyrocko: https://pyrocko.org/docs/current/install/
