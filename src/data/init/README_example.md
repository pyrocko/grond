Grond Project Layout
====================

For detailed instructions see the documentation https://pyrocko.org/grond/current.

Folder Layout
-------------

```
├── config
│   ├── laquila2009_joint.gronf
│   ├── ...
│   :
│
├── data
│   └── events  # several events could be set up here
│       ├── laquila2009
│       │   ├── event.txt
│       │   ├── insar
│       │   │   ├── dsc_insar.npz
│       │   │   ├── dsc_insar.yml
│       │   │   :
│       │   │
│       │   ├── waveforms
│       │   │   ├── raw    # contains Mini-SEED files
│       │   │   │   ├── trace_BK-CMB--BHE_2009-04-06_00-38-31.mseed
│       │   │   │   ├── trace_BK-CMB--BHN_2009-04-06_00-38-31.mseed
│       │   │   │   :
│       │   │   └── stations.xml
│       │   │
│       │   └── gnss
│       │       └── gnss.yml
│       :
│
├── gf_stores  # contains Green's functions
│   ├── Abruzzo_Ameri_nearfield # static near-field GF store
│   │   └── ...
│   ├── global_2s_25km  # dynamic far-field GF store
│   │   └── ...
│   :
│
├── runs  # created at runtime, contains individual optimisation results
│   └── ...
└── reports
    └── ...
```
