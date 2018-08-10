Centroid moment tensor from regional surface wave observations
==============================================================

This step-by-step guide explains how to obtain a probabilistic centroid moment tensor (CMT) solution from surface waves for an `Mw 5.9 aftershock <https://geofon.gfz-potsdam.de/eqinfo/event.php?id=gfz2018pmjk>`_ of the Mw 6.9 2018 Indonesia earthquake using Grond.

Setup
-----

To repeat this exercise on your machine, you should first `install Pyrocko <https://pyrocko.org/docs/current/install/>`_ and Grond (see :doc:`/install/index`), if you have not already done so. Then, copy the exercise project directory from Grond's git repos to a place of your choice:

.. code-block :: sh

    # git clone https://gitext.gfz-potsdam.de/heimann/grond.git  # <- done during installation
    cp -r grond/examples/grond-playground-regional ~/grond-playground-regional
    cd ~/grond-playground-regional


The project folder
------------------

The project folder now contains a configuration file for Grond, some utility scripts to download precalculated Green's functions and seismic waveforms from public datacenters.

.. code-block :: sh
    
    grond-playground-regional  # project folder
    ├── bin                         # directory with scripts
    │   ├── download_gf_stores.sh   # download precalculated Green's functions
    │   ├── grondown                # a simple event-based waveform downloader
    │   └── grondown_regional.sh     # downloader configured for this exercise
    └── config                      # directory for configuration files
        └── surface_cmt.gronf       # Grond configuration file for this exercise

Green's function download
-------------------------

To download the precalculated Green's functions needed in this exercise, run

.. code-block :: sh
    
    bin/download_gf_stores.sh

When the command succeeds, you should have a new subdirectory :file:`gf_stores` in your project folder:

.. code-block :: sh

    gf_stores
    └── crust2_j3/... # Green's function store

It contains a Pyrocko Green's function store, named ``crust2_j3``, which has been created using the `Fomosto <https://pyrocko.org/docs/current/apps/fomosto/index.html>`_ tool of `Pyrocko <http://pyrocko.org/>`_ and the modelling code `QSEIS <https://pyrocko.org/docs/current/apps/fomosto/backends.html#the-qseis-backend>`_. The Green's functions in this store have been calculated for a regional `Crust2x2 <https://igppweb.ucsd.edu/~gabi/crust2.html>`_ earth model for a source depths between 0 and 30 km in 1 km steps. It is sampled at 2 Hz, which is sufficient for our target frequency range of 0.01 - 0.1 Hz.

Seismic waveform data download
------------------------------

A preconfigured script is provided to download seismic waveform recordings via FDSN web services from the `IRIS <http://service.iris.edu/fdsnws/>`_ and `GEOFON <https://geofon.gfz-potsdam.de/waveform/webservices.php>`_ datacenters. Just run it with the GEOFON event ID of the study earthquake. The GEOFON event ID of the Mw 5.9 aftershock is ``gfz2018pmjk`` (you can find the ID in the `GEOFON catalog <https://geofon.gfz-potsdam.de/eqinfo/list.php>`_ event links).

Now run:

.. code-block :: sh
    
    bin/grondown_regional.sh gfz2015pmjk

This shell script calls the data downloader :file:`bin/grondown` with parameters appropriate to get a dataset of broadband seismometer recordings, sufficient for a surface wave CMT optimisation. It performs the following steps for us:

* Querry the `GEOFON catalog <https://geofon.gfz-potsdam.de/eqinfo/list.php>`_ for event information about ``gfz2018pmjk``.
* Select time windows based on event origin and time, considering that we want to analyse the signals at very low frequencies (0.01 - 0.1 Hz).
* Querry datacenters for seismic stations with epicentral distance between 0 and 1000 km.
* From the available recorder channels select appropriate ones for a target sampling rate of 2 Hz.
* Download raw waveform data for the selected stations and channels.
* Download instrument transfer function meta-information for all successfully downloaded waveform data.
* Calculate displacement seismograms for quality check (Grond will use the raw data). If all went well, the displacement seismograms should be valid in the frequency range 0.01 - 0.05 Hz, sampled at 1 Hz and rotated to radial, transverse, and vertical components.

After running the download script, the playground directory should contain a new :file:`data` directory with the following content:

.. code-block :: sh

    data
    └── events
        └── gfz2015pmjk
            ├── event.txt                 # catalog information about the event
            └── waveforms
                ├── grondown.command
                ├── prepared/...          # rotated, displacement waveforms
                ├── raw/...               # raw Mini-SEED waveforms
                ├── rest/...
                ├── stations.geofon.xml   # instrument response information
                ├── stations.iris.xml
                ├── stations.orfeus.xml
                ├── stations.prepared.txt # stations files for Snuffler
                └── stations.raw.txt

Because of various data problems, like missing instrument response information, gappy traces, data inconsistencies and what not, only about half of the initially requested stations will be useful in the optimisation. Some problems are not detected by the downloader, so we will have to look at the seismograms.

Data screening
--------------

For a quick visual inspection of the dataset, we can use the `Snuffler <https://pyrocko.org/docs/current/apps/snuffler/index.html>`_ program contained in Pyrocko.

.. code-block :: sh

    cd data/events/gfz2015pmjk/waveforms
    snuffler --event=../event.txt --stations=stations.prepared.txt prepared
    cd -  # change to previous folder

Figure 1 shows our view after some interactive adjustments in Snuffler. In particular, we we may want to

* sort the traces according to epicentral distance (Menu → check *Sort by Distance*).
* configure display style (Menu → uncheck *Show Boxes*, check *Common Scale per Station*, uncheck *Clip Traces*).
* filter between 0.01 and 0.05 Hz.
* add markers for expected P and S phase arrivals, (Menu → *Panels* → *Cake Phase (builtin)*).
* show only vertical components: Command ‣ :command:`c *z`.

.. figure:: ../../images/example_snuffler-gfz2015pmjk.svg
    :name: Fig. 1 Example surface wave CMT inversion
    :width: 100%
    :align: center
    
    **Figure 1**: Displacement seismograms for surface wave CMT optimisation as viewed in the waveform browser Snuffler.

Grond configuration
-------------------

The project folder already contains a configuration file for surface wave CMT optimisation with Grond, so let's have a look at it.

It's a `YAML`_ file: This file format has been choosen for the Grond configuration because it can represent arbitrarily nested data structures built from mappings, lists, and scalar values. It also provides an excellent balance between human and machine readability. When working with YAML files, it is good to know that the **indentation is part of the syntax** and that comments can be introduced with the ``#`` symbol. The type markers, like ``!grond.CMTProblemConfig``, select the Grond object type of the following mapping and it's documentation can likely be found in the :doc:`/library/index`.

.. literalinclude :: ../../../../examples/grond-playground-regional/config/regional_cmt.gronf
    :language: yaml
    :caption: config/regional_cmt.gronf (in project folder)

.. _YAML: https://en.wikipedia.org/wiki/YAML

Checking the optimisation setup
-------------------------------

Before running the actual optimisation, we can now use the command

.. code-block :: sh
    
    grond check config/regional_cmt.gronf gfz2015pmjk

to run some sanity checks. In particular, Grond will try to run a few forward models to see if the modelling works and if it can read the input data. If only one event is available, we can also neglect the event name argument in this and other Grond commands.

To get some more insight into the setup, we can now run

.. code-block :: sh

    grond report -so config/regional_cmt.gronf gfz2015pmjk

This will plot some diagnostic figures, create web pages in a new directory :file:`reports`, and finally open these in a web browser.


Starting the optimisation
-------------------------

Let's start the optimisation with:

.. code-block :: sh

    grond go config/regional_cmt.gronf

During the optimisation a status monitor will show the optimisation's progress.

.. figure:: ../../images/example_grond-run-regional.png
    :width: 100%
    :align: center

    **Figure 3**: Runtime information given by :command:`grond`.

Depending on the configured number of iterations and the computer's hardware the optimisation will run several minutes to hours.


Optimisation report
-------------------

Once the optimisation is finished we can generate and open the final report with:

.. code-block :: sh

    grond report -so config/regional_cmt.gronf
