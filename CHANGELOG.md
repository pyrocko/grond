# Changelog

All notable changes to Grond will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- Waveform targets: switch to change to acceleration / velocity fitting
- Waveform targets: include / exclude stations by pattern on target-group level
- Option to export list of stations used in one or more setups
  (`grond check --save-stations-used=<filename>`)

### Changed
- Transparent event loading and checking.

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
