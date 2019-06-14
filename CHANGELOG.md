# Changelog

All notable changes to Grond are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.3.1] 2019-06-08

### Added
- Allow controlling number of threads in `grond report` with `--threads`

### Changed
- Default number of threads used in `grond go` and `grond report` is now 1.
  Setting both together, `--parallel` and `--threads,` to values higher than 1
  may currently cause hangs.
- Improved control on threading utilization.

### Fixed
- Repaired `grond report` and `grond plot location_mt` which were broken in
  v1.3.0.


## [1.3.0] 2019-06-04

### Added
- Added Covariance weighting from `SatelliteTarget` and `GNSSTarget`.
- Added new Volume source (VLVD).
- Added SVG as `grond plot` export format.

### Changed
- Consistent utilisation of `ModelHistory` for plots and results.
- `MTLocationPlot` supports Gamma scaling (misfit^gamma).
- `SatelliteTarget` speeds up bootstrapping by multi-threading.
- `Envronment` can be initialised from `ProblemConfig`.
- Improved: `grondown` downloading seismic wave forms.


## [1.2.0] 2019-02-19

### Added
- Waveform targets: switch to change to acceleration / velocity fitting.
- Waveform targets: include / exclude stations by pattern on target-group
  level.
- Option to export list of stations used in one or more setups (`grond check
  --save-stations-used=<filename>`).
- Can now handle GNSS stations lacking complete set of component orientations.
- Noise analyser: added possibility to except stations with high S/N ratio from
  removal due to high pre-event noise.
- Improved unit tests.
- Report archive generation can now be skipped via command line flag or report
  configuration setting.
- CMT problem: can now switch between different source time functions.
- Added workaround switch for `"cannot project traces with displaced sampling"`
  issues.

### Changed
- Transparent event loading and checking.
- Noise analyser: target groups are now handled independently. Each group now
  uses its own threshold in weeding mode.
- Improved error handling (`grond check`, instrument responses, 
- Only exclude waveform targets when `distance_min` constraint is given in
  `problem_config`.
- Improved method chapter in documentation.

### Fixed
- Waveform fit plots: fix crashes while plotting results from joint inversions.
- Satellite fit plots: fix bug in source outline drawing
- Waveform targets: fixed handling of channel epochs from StationXML for
  channel orientations.
- Station plots: fixed problems with empty target groups.
- No more MPL warnings 'too many open figures' when creating sequence plots;
  the plots are now created one by one.
- Report: fix dsiplay issue with inaccessible elements in left navigation.
- Fixed crash when `starting_point` setting in highscore optimiser is set to
  `mean`.

## [1.1.1] 2019-02-05

### Fixed
- Bug in volume point source plot causing crashes.

## [1.1.0] 2019-01-22

### Added
- New VolumePointProblem to optimise magmatic and volcanic processes.
- New problem config sections in `grond init`.

### Changed
- Documentation of problems configurations are now centralised at 
  `src/data/snippets`.
- Output of `grond init list`.

### Fixed
- Bug in GNSSMisfitTarget.
- Plotting in GNSS plotting functions.
- Highscore optimiser logging output
- Satellite plot: Setting geographical aspect ratio for LatLon data
- GNSS Plotting function

## [1.0.0] 2019-01-07

### Added
- Cluster analysis on result ensemble with `grond cluster`.
- Plotting of cluster analysis results in `jointpar`, `mt_decompostion`, and
  `mt_fuzzy` plots.
- User can now add self-defined labels to a run dir (`grond tag`). The labels
  are shown in the report list.
- The optimiser's acceptance and choice history is now dumped and plotted.
- Added station distribution plot for seismic and GNSS stations.
- Installation instructions for Anaconda and pip.

### Changed
- Parallelized `grond report`.
- Rewritten `grond init` to deliver examples and commented snippets.
- Plots appearing in report can now be customized in a configuration file.
- Whether reference solutions are shown in plots can be configured in the
  report's plot config.
- `grond plot` now has a `--show` option to display MPL plots interactively.
- Plot `location_mt` improved.
- `ProblemConfig`: added `nthreads` argument.
- Plot `contributions`: show cumulative contributions for targets providing
  multiple misfits.
- Optimiser can now be configured to yield exactly reproducible results by
  providing seed values for all random number generators involved.
- Plots `sequence` and `fits_waveform`: layout as single plot figures by 
  default.

### Fixed
- Plot descriptions and titles improved.
- Reports can now be viewed in IE.
- Improved consistency of log messages.
- Fix display issues in waveform plots (partially hidden labels).
- Various improvements on the documentation.
- Improved robustness of report generation during `grund run`.
