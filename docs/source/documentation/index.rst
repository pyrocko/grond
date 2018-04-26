Basic Setup 
===========

Grond can be run as a command line tool or by calling Grond's library functions from a Python script. To get a brief description on available options of Grond's command line tool, run:

::
   
	grond --help 
	grond <subcommand> --help


Structuring your ``grond`` Projects 
-----------------------------------

To use ``grond`` with your data and the medium model of your choice, we suggest the following folder stucture. Single files listed here are explained below.

.. code-block :: sh

    ├── +
    ├── config
    │   ├── laquila.gronf
    │   ├── ...
    │   :
    │
    ├── data
    │   └── events  # several events could be set up here
    │       ├── laquila2009   # example for the 2009 L'Aquila event
    │       │   ├── event.txt  
    │       │   ├── insar      # contains kite-prepared InSAR-data
    │       │   │   ├── dsc_LAquila2009_Envisat.npz
    │       │   │   └── dsc_LAquila2009_Envisat.yml
    │       │   ├── waveforms  # contains seismic waveforms and meta data
    │       │   │   ├── raw    # raw downloaded data                           
    │       │   │   │   └── 2009-04-06_MW6.3_Central_Italy.421793.seed    
    │       │   │   └── stations.xml
    │       │   └── gnss
    │       │       └── gnss.yml   
    │       :
    │
    └── gf_stores  # contains Green's function stores for near-field data       
        ├── Abruzzo_Ameri_nearfield
        :   └── ...

Preparing your ``grond`` Data
--------------------------------

**Seismic Waveforms Data:**

Required input files are:
	- waveform data (.seed or .mseed format)
	- response functions (.xml format) 

To download waveform data and instrument response information, you may follow this link: ``downloadwave``
For significant earthquake, you can assemble your dataset from the global seismic network using ``pyrocko`` functions.

**InSAR Data:**

For preparing and incorporating InSAR data see the  `kite documentation`_. ``kite`` is an interactive tool for inspection and transport of static displacement maps. It can be used for data noise estimations, easy quadtree data subsampling and calculation of data error variance-covariance matrices for proper data weighting. 

Required input files are:

:: 

	kite_scene.yml kite_scene.npz


**GNSS Data:**

Requiered input file is a simple ``YAML`` file containing GNSS station positions, displacement values and measurement uncertainties. A gnss.yml-file should look like as follows:

.. code-block :: sh

    --- !pf.gnss.GNSSCampaign
    name: northridge
    stations:
    - !pf.gnss.GNSSStation
    lat: 18.30
    lon: -74.22
    elevation: 0.0
    depth: 0.0
    code: ANGL
    style: static
    north: !pf.gnss.GNSSComponent
        unit: m
        shift: 5.7e-06
        sigma: 5.7e-06
    east: !pf.gnss.GNSSComponent
        unit: m
        shift: -7e-07
        sigma: 8.8e-06
    up: !pf.gnss.GNSSComponent
        unit: m
        shift: 2.93e-05
        sigma: 1.11e-05
    - !pf.gnss.GNSSStation
    lat: 18.28
    ...

(add more station information in the same manner) 

Setup of Green's Functions Databases defining the Medium
--------------------------------------------------------

A Green's functions (GF) database is needed. You can either download from the online repository (`online GF databases`_) or compute them with the `fomosto`_ module of ``pyrocko``. Depending on the data sets, different GF databases are suitable:

.. _fomosto: https://pyrocko.org/docs/current/apps/fomosto/index.html


**GF's for global teleseismic waveform data:**

For you general point-source analysis a global store of Green's functions with a sampling frequency of 2 Hz may suffice. 

::

	fomosto download kinherd global_2s store 

You can browse for more GF stores to `download here`_.

**GF's for regional and local seismic waveform data:**

Regional analyses may require smaller and individual medium GF stores. Suitable are GF stores built with the method ``qseis``.

**GF's for near-field static displacements (e.g. InSAR, GNSS):**

Near-field static displacements require GF stores with high spatial sampling and mostly only little temporal sampling. With the PSGRN/PSCMP GF method you can build for any given local 1d-layered velocity model your own GF store ``psgrn``.


Preparing your ``grond`` Configuration File
-------------------------------------------

You can inititiate a ``grond`` configuration file for a centroid moment tensor optimization based on  global seismic waveforms with: 

.. code-block :: sh

    grond init > <project>.gronf
    
Identically, for static near-field displacement (e.g. InSAR and/or GNSS data sets) and finite source optimisation setup, ``grond`` configuration file can be initialise with: 

.. code-block :: sh

    grond init --insar > <project>.gronf
    grond init --gnss > <project>.gronf
    grond init --gnss --insar > <project>.gronf   
 
The ``targets`` (data and misfit setups for seimsic waveforms, InSAR and or GNSS data) can be combined and sources types can be exchanged. A ``grond`` configuration file showing all possible options with their default values is given using: 

.. code-block :: sh

    grond init --full > <project>.gronf`

Commented snippets of ``grond`` configuration files explaining all options can be found here for 
    * point-source optimizations based on waveforms: :download:`config_example_waveforms.yaml </../../examples/config_example_waveforms.yaml>`
    * finite source optimizations based on InSAR data: :download:`config_example_static.yaml </../../examples/config_example_static.yaml>`
    #* full configuration documentation 

    
.. literalinclude :: /../../examples/config_example_static.yaml
    :language: yaml


Optimisation
------------

You may want to check your dataset and configuration file (see suggestions above) and debug
it if needed with the command:

::

	grond check <configfile> <eventname>

Now, you may start the optimization for a given event using:

::
	
	grond go <configfile> <eventname>

During the optimization, results are aggregated in an output directory, referred to `<rundir>`  in the configuration. 


Results plots, exports and reports
----------------------------------

To visualize the results run:

::

	grond plot <plotnames> <rundir> 

The results can be exported in various ways by running the subcommand:

::

	grond export <what> <rundir>

Finally, you may run:

::
	
	grond report <rundir>
	grond report-index reports 

to aggregate all results to a browsable summary, (by default) under the directory `reports`.


.. _kite documentation: https://pyrocko.org/docs/kite/current/
.. _downloadwave: https://pyrocko.org/docs/current/library/examples/fdsn_download.html
.. _qseis: https://pyrocko.org/docs/current/apps/fomosto/tutorial.html#creating-a-new-green-s-function-store
.. _psgrn: https://pyrocko.org/docs/current/apps/fomosto/tutorial.html#creating-a-new-green-s-function-store
.. _online GF databases: http://kinherd.org:8080/gfws/static/stores/
.. _download here: http://kinherd.org:8080/gfws/



