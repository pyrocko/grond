# Changelog
All notable changes to the Grond project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]
### Added
- The optimiser's acceptance and choice history is now dumped and plotted.
- Plots appearing in report can now be customized in a configuration file.
- Cluster analysis on result ensemble with `grond cluster`.
- Plotting of cluster analysis results in `jointpar`, `mt_decompostion`, and
  `mt_fuzzy` plots.
- Whether reference solutions are shown in plots can be configured in the
  report's plot config.
- Plot descriptions improved
- Grond plot now has a `--show` option to display MPL plots interactively.
- Grond reports can now be viewed in IE.
- User can now add self-defined labels to a run dir (`grond tag`). The labels
  are shown in the report list.
- Parallelized 'grond report'
- Improved `location_mt` plot
- Improved Plot short descriptions `grond plot list` 
- Contributions Plot: cummulating nmisfits for `ScatteredTargets` (e.g. `SatelliteMisfitTarget`)
- `ProblemConfig` added `nthreads` argument
