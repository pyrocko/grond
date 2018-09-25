Dataset
=======

The :term:`dataset` configuration section is the core of data input and data configuration. Here we define where all observed data is stored, define response functions paths for waveforms, blacklist stations or load phase arrival picks. Finally the :class:`~grond.dataset.Dataset` object is managing all available data.

.. literalinclude :: dataset.gronf
    :language: yaml
    :caption: A basic dataset section of a Grond configuration file (``gronf``).

General configuration and templating
------------------------------------

All folder and file paths in the dataset support templating and prefixing:

.. glossary ::

    ``events_path``
        File with hypocenter information and possibly reference solution.

    ``path_prefix``
        defines a prefix which is prepended to all paths in the configuration.

    ``${event_name}``
       will be substituted with the event name defined in your ``events_path`` file.


Waveform data
-------------
Usually raw, unrestituted waveforms are loaded into Grond, together with StationXML data describing the station location and response function - Grond will take care of proper restitution.

.. glossary ::

    ``waveform_paths``
        List of directories with raw waveform data.

    ``stations_stationxml_paths``
        List of files with station coordinates in StationXML format.

    ``stations_path``
        List of files with station coordinates in Pyrocko format.

    ``extend_incomplete``
        Extend incomplete seismic traces: ``true``/``false``.

    ``clippings_path``
        Pyrocko marker file indicating where a seismic trace is masked.

    ``responses_stationxml_paths``
        List of StationXML response files for restitution of the raw waveform data.

    ``responses_sacpz_path``
        List of SACPZ response files for restitution of the raw waveform data.

    ``station_corrections_path``
        File containing station correction informations. See :download:`example station corrections <station_corrections.yaml>`.

    ``apply_correction_factors``
        Apply the correction factors from station corrections: ``true``/``false``.

    ``apply_correction_delays``
        Apply the correction delays from station corrections: ``true``/``false``.

    ``picks_paths``
        List of phase picks in `Pyrocko format <https://pyrocko.org/docs/current/apps/snuffler/tutorial.html#markers>`_.

    ``blacklist``
        List of stations/components to be **excluded** according to their STA, NET.STA, NET.STA.LOC, or NET.STA.LOC.CHA codes

    ``blacklist_paths``
        List of text files with blacklisted stations in NSLC pattern.

    ``whitelist``
        List of stations/components to be **included** according to their STA, NET.STA, NET.STA.LOC, or NET.STA.LOC.CHA codes

        Note: when whitelisting on channel level, both, the raw and the processed channel codes have to be listed.

    ``whitelist_paths``
        List of text files with whitelisted stations in NSLC pattern.

    ``synthetic_test``
        Run a synthetic test: ``true``/``false``


Satellite data
--------------

Unwrapped static surface displacements have to be prepared in Kite format.

.. glossary ::
    
    ``kite_scene_paths``
        List of folders where pre-processed `Kite <https://pyrocko.org/>`_ surface displacement scenes are stored.

GNSS campaign data
------------------

Single measurements of surface displacement data from GNSS campaigns can be loaded from YAML text files.

.. glossary ::

    ``gnss_campaign_paths``
        List of folders where `GNSS data <https://pyrocko.org/docs/current/library/examples/gnss_data.html>`_  of static surface displacements are stored.


