Installation Instructions
=========================

``grond`` depends on `pyrocko`_, which needs to be installed first.

Then install ``grond``:

.. code-block :: sh

    git clone https://gitext.gfz-potsdam.de/heimann/grond.git
    cd grond
    sudo python3 setup.py install


For updating an existing installation of ``grond``:

.. code-block :: sh

    cd grond  # change to the directory to where you cloned grond initially
    git pull origin master
    sudo python3 setup.py install


For a support of InSAR data modelling you further need to install the module `kite`_.

.. _kite: https://pyrocko.org/docs/kite/current/
.. _pyrocko: https://pyrocko.org/docs/current/install/