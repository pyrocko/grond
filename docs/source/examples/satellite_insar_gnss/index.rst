Rectangular source plane from InSAR and GNSS observations
=========================================================

This example will demonstrate a joint inversion of a finite rectangular fault plane for the 2019 Ridgecrest Earthquake. We will jointly invert InSAR and GNSS data using Grond. Unwrapped InSAR surface displacement data from the `Advanced Rapid Imaging and Analysis (ARIA) Project for Natural Hazards. <https://aria.jpl.nasa.gov/>`_ by JPL/NASA. The coseismic GNSS displacements are delivered by `Nevada Geodetic Laboratory <http://geodesy.unr.edu/index.php>`_.

This tutorial will guide through the conversion of ARIA data products and preparation of the unwrapped InSAR surface displacements in Kite. The GNSS data will be imported and pre-selected.

This is an advanced example to demonstrate a geodetic joint inversion. If you haven't worked with Grond earthquake inversions before, we recommend to exercise a single dataset example first.

Setup
-----

To repeat this exercise on your machine, you should first `install Pyrocko
<https://pyrocko.org/docs/current/install/>`_ and Grond (see :doc:`/install/index`), if you have not already done so. Then create the example project with:

.. code-block :: sh

    grond init example_gnss_insar grond-joint-geodetic/

The project folder
------------------

The project folder now contains a configuration file for Grond and some utility scripts to download pre-calculated static Green's functions and the InSAR and GNSS data:

.. code-block :: sh
    
    grond-joint-geodetic        # project folder
    ├── bin                        # directory with scripts
    │   ├── download_gf_stores.sh  # download pre-calculated Green's functions
    │   ├── download_data.sh      # a simple event-based waveform downloader
    └── config                     # directory for configuration files
        └── insar_rectangular.gronf       # Grond configuration file for this exercise

Green's function download
-------------------------

To download the pre-calculated Green's functions needed in this exercise, run

.. code-block :: sh
    
    bin/download_gf_stores.sh

When the command succeeds, you should have a new subdirectory :file:`gf_stores` in your project folder:

.. code-block :: sh

    gf_stores
    └── crust2_ib_static/... # Green's function store

It contains a Pyrocko Green's function store, named ``crust2_ib_static``, which has been created using the `Fomosto <https://pyrocko.org/docs/current/apps/fomosto/index.html>`_ tool of `Pyrocko <http://pyrocko.org/>`_ and the modelling code `PSGRN/PSCMP <https://pyrocko.org/docs/current/apps/fomosto/backends.html#the-psgrn-pscmp-backend>`_. The Green's functions in this store have been calculated for a regional `CRUST2 <https://igppweb.ucsd.edu/~gabi/crust2.html>`_ earth model for source depths between 0 and 30 km in 500 m steps, and horizontal extent from 0 - 300 km in 500 m steps.

InSAR and GNSS data download
----------------------------

The example includes a script to download unwrapped InSAR and GNSS data from Pyrocko's servers. The unwrapped InSAR Sentinel-1 surface displacement data is provided by `ARIA (NASA) <https://aria.jpl.nasa.gov/>`_. An example how to download and convert ARIA data in Kite can be followed in the `documentation <https://pyrocko.org/docs/kite/>`_.

.. code-block :: sh
    
    bin/download_data.sh

This will the InSAR scenes and GNSS data to :file:`data/events/2019-ridgecrest`. Surface displacement is held in Kite container format and the Pyrocko GNSS containers.

Unwrapped InSAR displacement preparation with Kite
--------------------------------------------------

The downloaded InSAR data has to be prepared for the inversion with the Kite tool. To install the software, follow the `install instructions <https://pyrocko.org/docs/kite/current/installation.html>`_.

Once the software is installed we need to parametrize the two scenes:

    1. The data sub-sampling quadtree: This efficiently reduces the resolution of the scene, yet conserves the important data information. A reduced number of samples will benefit the forward-modelling computing cost.

    2. Estimate the spatial data covariance: By looking at the spatial noise of the scene we can estimate the data covariance. Kite enables us to calculate a covariance matrix for the quadtree, which will be used as a weight matrix in our Grond inversion.


.. note ::
    The scenes come pre-configured. The following steps of defining the quadtree and calculating the covariance matrix are optional.


We start by parametrizing the quadtree: find a good parameters for the sub-sampling quadtree by tuning four parameters:

    1. ``epsilon``, the variance threshold in each quadtree's tile.
    2. ``nan_fraction``, percentage of allowed NaN pixels per tile.
    3. ``tile_size_min``, minimum size of the tiles.
    4. ``tile_size_max``, maximum size of the tiles.

.. code-block :: sh

    spool data/events/2019-ridgecrest/insar/ascending

    spool data/events/2019-ridgecrest/insar/descending

Now we can parametrize the quadtree visually:

.. figure:: ../../images/example_spool-ridgecrest-quadtree.png
    :name: Fig. 1 Example InSAR quadtree
    :width: 100%
    :align: center
    
    **Figure 1**: Parametrizing the quadtree. This efficiently sub-samples the high-resolution Sentinel-1 surface displacement data. (command :command:`spool`; `Kite <https://pyrocko.org/docs/kite/>`_ toolbox).

.. note ::
    
    Delete unnecessary tiles of the quadtree by right-clicking, and delete with :kbd:`Del`.

Once you are done, click on the next tab :guilabel:`Scene.covariance`. Now we will define a window for the data's noise. The window's data will be use for calculating the spatial covariance of the scene (see `details <https://pyrocko.org/kite/docs/current/examples/covariance.html>`_).

Use a spatial window far away from the earthquake signal to capture only the noise, yet the bigger the window is, the better the data covariance estimation.

On the left hand side of the GUI you find parameters to tune the spatial covariance analysis. We now can fit an analytical model to the empirical covariance: :math:`\exp(d)` and :math:`\exp + \sin`. For more details on the method, see `Kite's documentation <https://pyrocko.org/docs/kite/current>`_.

.. figure:: ../../images/example_spool-ridgecrest-covariance.png
    :name: Fig. 2 Example InSAR covariance
    :width: 100%
    :align: center
    
    **Figure 2**: Spatial data covariance inspection and definition of the noise window.

Once we finished parametrisation of the quadtree and covariance, we have to calculate the full covariance and weight matrix from the complete scene resolution:

    1. Calculate the full covariance: :menuselection:`Tools --> Calculate Full Matrix`
       Depending on the scene's resolution this process can take time.
    2. Save the parametrized scene: :menuselection:`File --> Save Scene`.


Grond configuration
-------------------

The project folder already contains a configuration file for rectangular source optimisation with Grond, so let's have a look at it.

It's a `YAML`_ file: This file format has been chosen for the Grond configuration because it can represent arbitrarily nested data structures built from mappings, lists, and scalar values. It also provides an excellent balance between human and machine readability. When working with YAML files, it is good to know that the **indentation is part of the syntax** and that comments can be introduced with the ``#`` symbol. The type markers, like ``!grond.RectangularProblemConfig``, select the Grond object type of the following mapping and it's documentation can likely be found in the :doc:`/library/index`.


.. literalinclude :: ../../../../examples/example_insar_gnss/config/insar_rectangular.gronf
    :language: yaml
    :caption: config/insar_rectangular.gronf (in project folder)


Checking the optimisation setup
-------------------------------

Before running the actual optimisation, we can now use the command

.. code-block :: sh
    
    grond check config/insar_rectangular.gronf

to run some sanity checks. In particular, Grond will try to run a few forward models to see if the modelling works and if it can read the input data. If only one event is available, we can also neglect the event name argument in this and other Grond commands.


Starting the optimisation
-------------------------

Now we are set to start the optimisation with:

.. code-block :: sh

    grond go config/insar_rectangular.gronf


During the optimisation a status monitor will show the optimisation's progress.

.. raw:: html
    
    <script id="asciicast-1wi554jdNaO8Pn2HNx3hJaw9g" src="https://asciinema.org/a/1wi554jdNaO8Pn2HNx3hJaw9g.js" async></script>

Depending on the configured number of iterations and the computer's hardware the optimisation will run several minutes to hours.


Optimisation report
-------------------

Once the optimisation is finished we can generate and open the final report with:

.. code-block :: sh

    grond report -so runs/rectangular_2019ridgecrest.grun


Example report
~~~~~~~~~~~~~~

Explore the `online example reports <https://pyrocko.org/grond/reports>`_ to see what information the inversion reveals.


.. _Kite: https://pyrocko.org/docs/kite/current/
.. _YAML: https://en.wikipedia.org/wiki/YAML
