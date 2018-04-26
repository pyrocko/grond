# Grond

## Introduction

is a battering ram to non-linearly explore earthquake source parameters based on geophysical data sets. Grond uses as foundation the pyrocko toolbox and is for some functionalities supported by pyrocko modules. 

Enabled earthquake source types so far are:
    * centroid moment tensors
    * double-couple point-sources
    * finite rectangular sources
    
Enabled data sets to be used for the inverse modelling are 
    * seismic waveforms (all standard and some strange formats)
    * InSAR coseismic near-field displacement maps
    * GNSS coseismic displacement vectors
    
These data sets can used also combined. For using coseismic InSAR displacement maps, a [Kite] InSAR data container needs to be prepared. [Kite] is an interactive InSAR-data post-processor, which imports a lot of standard InSAR processor output formats. 
For using GNSS campaign coseismic displacement vectors, a simple data file needs to be prepared, containing the GNSS station positions, displacement values and measurement uncertainties.

For the forward modelling precalculated and stored Green's functions are used. These stores can be obtained from an online repository or are user-defined and user-calculated using pyrockos fomosto module.

## Installation

First, install [Pyrocko](http://pyrocko.org/docs/current/install/),
then install Grond:

```bash
git clone https://gitext.gfz-potsdam.de/heimann/grond.git
cd grond
sudo python setup.py install [Pyrocko](http://pyrocko.org/docs/current/install/),
```

For support of InSAR data modelling install also [Kite](https://github.com/pyrocko/kite)

## Updating an existing installation

```bash
cd grond  # change to the directory to where you cloned grond initially
git pull origin master
sudo python setup.py install
```

# Documentation

Find the documentation at https://pyrocko.github.io/grond/.



## Basic preparation and setup of Green's functions and data 

### Green's function preparation defining the medium 

To use grond with your data and the medium model of your choice, we suggest the following folder stucture below. 
Single files listed here are explained below.


### data preparation

add pyrocko examples for data download

For incorporating InSAR data see the  [Kite documentation and examples][https://pyrocko.org/docs/kite/current/]

### Configuration file

You can inititiate a configuration file for a centroid moment tensor optimization based on  global seismic waveforms with `grond init > <project>.gronf`.
For specifically static near-field displacement, e.g. measured with InSAR you initiate the 
 grond configuration file with `grond init --insar > <project>.gronf`.
An example configuration file with all possible option is given using `grond init --full > <project>.gronf`. 

Commented snippets of configuration files explaining all options can be found here for 
    * point-source optimizations based on waveforms (links)
    * finite source optimizations based on InSAR data (links)
    * full configuration documentation (links)
    
### Structuring grond projects 

```
├── +
├── config
│   └── laquila.gronf
├── data
│   └── events  # several events could be set up here
│       ├── laquila2009   # example for the 2009 │L'Aquila event 
│           ├── event.txt  
│           ├── insar      # contains kite-prepared InSAR-data
│           │   ├── dsc_LAquila2009_Envisat.npz│
│           │   └── dsc_LAquila2009_Envisat.yml
│           └── waveforms  # contains seismic waveforms and meta data
│               ├── raw    # raw downloaded data 
│               │   └── 2009-04-06_MW6.3_Central_Italy.421793.seed
│               └── stations.xml 
└── gf_stores  # contains Green's function stores for near-field data
    └── Abruzzo_Ameri_nearfield
         └── ...
```
    

## Basic usage 

Grond can be run as a command line tool or by calling Grond's library functions
from a Python script. To get a brief description on available options of
Grond's command line tool, run `grond --help` or `grond <subcommand> --help`.

You may want to check your dataset and configuration (see suggestions above)and debug
it if needed, with the command `grond check <configfile> <eventname>`. 

Now, you may start the optimization for a given event using `grond go <configfile> <eventname>`.

During the optimization, results are aggregated in
a directory, referred to in the configuration as `<rundir>`. To visualize the
results run `grond plot <plotnames> <rundir>`. The results can be exported in
various ways by running the subcommand `grond export <what> <rundir>`. Finally,
you may run `grond report <rundir>` to aggregate results to a browsable
summary, (by default) under the directory `reports`.


## Example configuration file

```yaml
%YAML 1.1
--- !grond.Config
# Path, where to store output (run-directories)
path_prefix: '..'
rundir_template: 'runs/${problem_name}.run'
dataset_config: !grond.DatasetConfig
  # File with hypocenter information and possibly reference solution
  events_path: 'data/events/${event_name}/event.txt'
  extend_incomplete: false

  # List of paths to InSAR data sets prepared with `kite`
  kite_scene_paths: ['data/events/${event_name}/insar']
  # List of paths to GNSS displacement vectors 
  gnss_campaign_paths: ['data/events/${event_name}/gnss']
```
  
# -----------------------------------------------------------------------------
# Configuration section selecting data to be included in the data optimization. In this example these are InSAR displacements and GNSS coseismic displacements. In the  target groups also the objective function for the optimization is defined per group. The targets can be composed of one or more contributions, each represented by a !grond.TargetConfig section.
# -----------------------------------------------------------------------------

```
target_groups:
# setup for InSAR
- !grond.SatelliteTargetGroup

  # misfits are normalized within each normalization_family separately
  normalisation_family: 'statics'

  # Name and group reference of this contribution
  path: 'insar'
  
  # How to weight the target group within in the global misfit
  weight: 1.0

  # Subsection defining the misfit
  misfit_config: !grond.SatelliteMisfitConfig
    # account for a remaining planar ramp in the InSAR data.
    optimize_orbital_ramp: true
    # define ranges for offset and gradients in m/m 
    ranges:
      offset: '-0.05 .. 0.05'
      ramp_east: '-1e-7 .. 1e-7'
      ramp_north: '-1e-7 .. 1e-7'
  
  # Green's function store to be used (see fomosto)
  store_id: 'Abruzzo_Ameri_nearfield'
  
  # list of individual InSAR scenes with scene identities.
  # Can be like here a '*all' wildcard. 
  kite_scenes: ['*all']
  
# setup for coseismic GNSS displacements
- !grond.GNSSCampaignTargetGroup

  # here set to be in the same family as InSAR targets
  normalisation_family: 'statics'
  
  path: 'gnss'
  weight: 1.0
  misfit_config: !grond.GNSSCampaignMisfitConfig {}
  store_id: 'Abruzzo_Ameri_nearfield'
  
  # list of individual GNSS campaign data with campaign identities.
  gnss_campaigns: ['*all']
```
# -----------------------------------------------------------------------------
# Definition of the problem to be solved - source model, parameter space, and
# global misfit configuration settings.
# -----------------------------------------------------------------------------

  
  
```  
  problem_config: !grond.RectangularProblemConfig
  
  # Name used when creating output directory
  name_template: 'rect_source'
  
  # Here a L2-norm is defined. Put to 1 for l1-norm
  norm_exponent: 2
  
  # Definition of model parameter space to be searched in the optimization
  ranges:
    # relative positions [m] to reference location in 'event.txt' 
    north_shift: '-2000 .. 20000'
    east_shift: '-2000 .. 20000'
    
    # depth [m] of the upper fault edge (not centroid!)
    depth: '5000 .. 30000'
    
    # fault dimensions and fault slip [m]
    length: '12000 .. 18000'
    width: '4000 .. 14000'
    slip: '0.2 .. 2.'
    
    # mechanism parameters [degree]
    strike: '80 .. 330'
    dip: '0 .. 60'
    rake: '60 .. 90'
    
  # decimation factor for point sources building the finite source 
  decimation_factor: 8
  
  # Define the norm for the global misfits combining the members of all defined 'normalization_family` 
  norm_exponent: 1

# -----------------------------------------------------------------------------
# Configuration of the optimization procedure. The following example setup will
# run a Bayesian bootstrap optimization (BABO).
# -----------------------------------------------------------------------------
 
optimizer_config: !grond.HighScoreOptimizerConfig
  nbootstrap: 0
  sampler_phases:
  - !grond.UniformSamplerPhase
      niterations: 5000
  - !grond.DirectedSamplerPhase
      niterations: 30000
      scatter_scale_begin: 1.6
      scatter_scale_end: 0.8

# -----------------------------------------------------------------------------
# The engine is the forward modelling machine and is here configured.
# -----------------------------------------------------------------------------

engine_config: !grond.EngineConfig
  # The Green's function stores have been given with the targets?
  gf_stores_from_pyrocko_config: true
  # there is a common folder containing gf stores 
  gf_store_superdirs: ['gf_stores']
  # ...
  gf_store_dirs: []
  
