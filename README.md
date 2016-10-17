# Grond

A bootstrap-based probabilistic battering ram to explore solution spaces in
earthquake source parameter estimation problems.


## Installation

First, install [Pyrocko](http://pyrocko.org/current/install.html),
then install Grond:

```bash
git clone https://gitext.gfz-potsdam.de/heimann/grond.git
cd grond
sudo python setup.py install
```


## Usage

Grond can be run as a command line tool or by calling Grond's library functions
from a Python script. To get a brief description on available options of
Grond's command line tool, run `grond --help` or `grond <subcommand> --help`.
Once dataset and configuration are ready, the command
`grond go <configfile> <eventname>` starts the optimization algorithm for a
selected event. During the optimization, results are aggregated in a directory,
referred to in the configuration as `<rundir>`. To visualize the results run
`grond plot <plotnames> <rundir>`. The results can be exported in various way
by running the subcommand `grond export <what> <rundir>`.


## Example configuration file

```yaml
%YAML 1.1
--- !grond.Config

# Path, where to store output (run-directories)
rundir_template: 'gruns/${problem_name}.run'


# -----------------------------------------------------------------------------
# Configuration section for dataset (input data)
# -----------------------------------------------------------------------------

dataset_config: !grond.DatasetConfig

  # List of files with station coordinates
  stations_stationxml_paths:
  - 'events/${event_name}/responses.stationxml'

  # File with hypocenter information and possibly reference solution
  events_path: 'events/${event_name}/event.txt'

  # List of directories with raw waveform data
  waveform_paths:
  - 'data/${event_name}/raw'

  # List of files with instrument response information
  responses_stationxml_paths:
  - 'meta/${event_name}/responses.stationxml'

  # List with station/channel codes to exclude
  blacklist: [OTAV, RCBR, PTGA, AXAS2, SNAA, PAYG, RAR, SDV, VAL, XMAS, ODZ,
              Z51A, MSVF, SHEL, SUR, ROSA, 'IU.PTCN.00.T', 'C1.VA02..T', RPN]


# -----------------------------------------------------------------------------
# Configuration section for synthetic seismogram engine (configures where
# to look for GF stores)
# -----------------------------------------------------------------------------

engine_config: !grond.EngineConfig

  # Whether to use GF store directories listed in ~/.pyrocko/config.pf
  gf_stores_from_pyrocko_config: false

  # List of directories with GF stores
  gf_store_superdirs:
  - 'gf_stores'


# -----------------------------------------------------------------------------
# Configuration section selecting data to be included in the data fitting
# procedure. This defines the objective function to be minimized in the
# optimization procedure. It can be composed of one or more contributions, each
# represented by a !grond.TargetConfig section.
# -----------------------------------------------------------------------------

target_configs:

- !grond.TargetConfig

  # Name of the super-group to which this contribution belongs
  super_group: 'time_domain'

  # Name of the group to which this contribution belongs
  group: 'rayleigh'

  # Minimum distance of stations to be considered
  distance_min: 1000e3

  # Maximum distance of stations to be considered
  distance_max: 10000e3

  # List with names of channels to be considered
  channels: ['Z']

  # How to weight stations from this contribution in the global misfit
  weight: 1.0

  # Subsection on how to fit the traces
  inner_misfit_config: !grond.InnerMisfitConfig

    # Frequency band [Hz] of acausal filter (flat part of frequency taper)
    fmin: 0.002
    fmax: 0.008

    # Factor defining fall-off of frequency taper
    # (zero at fmin/ffactor, fmax*ffactor)
    ffactor: 1.5

    # Time window to include in the data fitting. Times can be defined offset
    # to given phase arrivals. E.g. 'begin-100' would mean 100 s before arrival
    # of the phase named 'begin', which must be defined in the travel time
    # tables in the GF store.
    tmin: '{stored:begin}-100'
    tmax: '{stored:end}+100'

    # How to fit the data (available choices: 'time_domain',
    # 'frequency_domain', 'absolute', 'envelope', 'cc_max_norm')
    domain: 'time_domain'

  # How to interpolate the Green's functions (available choices:
  # 'nearest_neighbor', 'multilinear')
  interpolation: 'nearest_neighbor'

  # Name of GF store to use
  store_id: 'global_20s_shallow'


# A second contribution to the misfit function (for descriptions, see above)
- !grond.TargetConfig
  super_group: 'time_domain'
  group: 'love'
  distance_min: 1000e3
  distance_max: 10000e3
  channels: ['T']
  weight: 1.0
  inner_misfit_config: !grond.InnerMisfitConfig
    fmin: 0.002
    fmax: 0.008
    ffactor: 1.5
    tmin: '{stored:begin}-100'
    tmax: '{stored:end}+100'
    domain: 'time_domain'
  interpolation: 'nearest_neighbor'
  store_id: 'global_20s_shallow'


# -----------------------------------------------------------------------------
# Definition of the problem to be solved - source model, parameter space, and
# global misfit configuration settings.
# -----------------------------------------------------------------------------

problem_config: !grond.CMTProblemConfig

  # Name used when creating output directory
  name_template: 'cmt_surface_wave_${event_name}'

  # Definition of model parameter space to be searched in the optimization
  ranges:
    # Time relative to hypocenter origin time [s]
    time: '-100 .. 100 | add'

    # Centroid location with respect to hypocenter origin [m]
    north_shift: '-200e3 .. 200e3'
    east_shift: '-200e3 .. 200e3'
    depth: '0 .. 50e3'

    # Range of magnitudes to allow
    magnitude: '7.0 .. 9.0'

    # Relative moment tensor component ranges (don't touch)
    rmnn: '-1.41421 .. 1.41421'
    rmee: '-1.41421 .. 1.41421'
    rmdd: '-1.41421 .. 1.41421'
    rmne: '-1 .. 1'
    rmnd: '-1 .. 1'
    rmed: '-1 .. 1'

    # Source duration range [s]
    duration: '30. .. 120.'

  # Clearance distance around stations (no models with origin closer than this
  # distance to any station are produced by the sampler)
  distance_min: 0.

  # Number of bootstrap realizations to be tracked simultaneously in the
  # optimization
  nbootstrap: 100

  # Type of moment tensor to restrict to (choices: 'full', 'deviatoric')
  mt_type: 'deviatoric'

  # Whether to apply automatic weighting to balance the effects of geometric
  # spreading etc.
  apply_balancing_weights: true


# -----------------------------------------------------------------------------
# Configuration of the optimization procedure
# -----------------------------------------------------------------------------

solver_config: !grond.SolverConfig

  # Distribution used when drawing new candidate models (choices: 'normal',
  # 'multivariate_normal') (used in 'transition', 'explorative', and
  # 'non_explorative' phase)
  sampler_distribution: 'normal'

  # Number of iterations to operate in 'uniform' mode
  niter_uniform: 1000

  # Number of iterations to operate in 'transition' mode
  niter_transition: 2000

  # Number of iterations to operate in 'explorative' mode
  niter_explorative: 2000

  # Number of iterations to operate in 'non_explorative' mode
  niter_non_explorative: 0

  # Multiplicator for width of sampler distribution in 'explorative' and
  # 'non-explorative' phases
  scatter_scale: 0.25

  # Multiplicator for width of sampler distribution at start of 'transition'
  # phase. (From there, it exponentially decreases to the value defined in
  # 'scatter_scale' during the 'transition' phase).
  scatter_scale_transition: 2.


# -----------------------------------------------------------------------------
# Configuration of pre-optimization analysis phase. E.g. balancing weights are
# determined during this phase.
# -----------------------------------------------------------------------------

analyser_config: !grond.AnalyserConfig

  # Number of iterations (number of models to forward model in the analysis,
  # larger number -> better statistics)
  niter: 1000
```
