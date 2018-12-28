# Changelog
All notable changes to the Grond project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### New Features
- Cluster analysis on result ensemble with `grond cluster`.
- Plotting of cluster analysis results in `jointpar`, `mt_decompostion`, and
  `mt_fuzzy` plots.
- User can now add self-defined labels to a run dir (`grond tag`). The labels
  are shown in the report list.
- The optimiser's acceptance and choice history is now dumped and plotted.
- Added station distrubtion plot for seismic and GNSS stations.
- Added `AcceptancePlot` to follow acceptance history.

### Bugfixes and Improvements
- Parallelized `grond report`.
- Rewrite `grond init` to deliver examples and commented snippets.
- Plot descriptions and titles improved.
- Reports can now be viewed in IE.
- Plots appearing in report can now be customized in a configuration file.
- Whether reference solutions are shown in plots can be configured in the report's plot config.
- `grond plot` now has a `--show` option to display MPL plots interactively.
- Plot `location_mt` improved.
- Contributions Plot: cummulating `nmisfits` for `ScatteredTargets` (e.g. `SatelliteMisfitTarget`).
- `ProblemConfig` added `nthreads` argument.
- Improved consistency of log messages.
- Plots for waveforms showing window length.
- Optimiser now has a consistent and fully reproducable `RandomState`.
- Various improvements on the documentation.
