# Changelog

All notable changes to the caver-pymol-plugin will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [4.2.1] — 2026-07-01

### Fixed
- Restore `org-openide-util-lookup-8.3.1.jar` — required at runtime by
  `delaunay-cell-discovery` via ServiceLoader (`Lookup.getDefault()`).
  Removed in 4.2.0 by mistake after static import analysis missed the
  runtime class loading path (#29).
- Use `os.pathsep` for Java classpath separator instead of hard-coded `:`
  for Windows compatibility (#29).

### CI
- Add Windows (`windows-latest`) and macOS (`macos-15`) test jobs to the
  Bare Tests matrix (#29).
- Temporarily drop `pymol-bundle` (latest) job due to missing
  `libCatch2.so` shared library in the Schrodinger conda channel (#29).

## [4.2.0] — 2026-07-01

### Changed
- **Standalone CAVER 3.1.0 engine.** Replaced the stripped-down 2012-era Java
  fork (148 source files, simplified Voronoi diagram, no MD trajectory support)
  with the recovered standalone CAVER 3.1.0 engine (186 source files). The new
  engine adds additively weighted (AW) Voronoi diagrams via the awvoronoi
  library, Pulsar straight-path search mode, AllPaths alternative search, and
  cross-snapshot tunnel clustering via `TunnelSetClustering`. The Java
  invocation switches from `-jar caver.jar` to `-cp caver.jar:lib/*
  caver.ui.Launcher` (#28).

- **Progress tracking.** Updated the per-snapshot progress counter from
  `.pdb.obj` to `.obj` suffix to match the new engine's `tunnels_N.obj`
  output format (#28).

- **Dependency refresh.** Replaced the 4 old `lib/` JARs (`AverageLinkClustering`,
  `kd`, `ml`, `vecmath`) with 9 lean standalone JARs: AW Voronoi core
  (`awvoronoi-1.0.2`, `awvoronoi-diagram-1.0.1`), Delaunay triangulation
  (`caver-delaunay-factory-1.0.3`, `caver-qhull-3.0_beta4`,
  `delaunay-cell-discovery-1.0.0`), alternative search (`AllPaths.jar`),
  clustering (`AverageLinkClustering.jar`), KD-tree (`kd.jar`), and vector
  math (`vecmath.jar`). 5 distribution-only JARs removed as dead weight:
  `awvoronoi-run-1.0.0` (standalone CLI), `schema-pdbx-4.1.0` (2.2 MB PDBx/mmCIF),
  `schema-atoms-1.1.0`, `jaxb2-basics-runtime-0.6.4`, and
  `org-openide-util-lookup-8.3.1` (#28).

### Fixed
- **Apple Silicon Qt bus error.** Fixed a `SIGBUS` (`EXC_ARM_DA_ALIGN`) crash on
  macOS ARM64 test teardown. The crash occurred when `QtCore.abi3.so`'s atexit
  handler called `sip_api_visit_wrappers` during `Py_FinalizeEx` after sip
  wrappers had already been partially freed. Calling `sip.setdestroyonexit(False)`
  at session startup prevents the dangling-pointer walk. Makefile's
  `test-skip-ci-segfault` target now also tolerates exit code 138 (SIGBUS
  on macOS, signal 10) alongside the existing 139 (SIGSEGV) (#28).

### Removed
- **`ml.jar` (Weka).** The Weka machine-learning library used by the old
  engine's `Clustering.java` is replaced by `AverageLinkClustering.jar` in
  the standalone engine (#28).

## [4.1.3] — 2025-07-27

### Fixed
- Corrected documentation references and reduced movie file size (#27)

## [4.1.2] — 2025-07-20

### Fixed
- Run ID now correctly resolved from the PyMOL path context, fixing edge
  cases where path-dependent state lookups failed (#23)

## [4.1.1] — 2025-07-13

### Added
- **Tunnel plotter.** Interactive matplotlib-based plotting of tunnel
  profiles (radius vs. length) from within the analysis panel, with
  configurable image size, DPI, colormap selection, and aspect ratio
  locking (#21)

## [4.1.0] — 2025-07-06

### Added
- **MD trajectory analysis.** The analyst panel can now load and analyze
  tunnels computed from multi-state MD trajectories, with per-frame tunnel
  preview and spectrum rendering (#18)

### Fixed
- Radio button toggle behavior for custom starting point selection (#17)
- Failed output directory IDs when directory listing was incomplete (#16)
- Empty output directory edge case (#15)

## [4.0.2] — 2025-06-29

### Added
- **`caver_set` autocomplete.** The `caver_set` PyMOL command now provides
  tab-completion for configuration keys and values via `CaverShortcut` (#13)

### Fixed
- Miscellaneous bug fixes and stability improvements (#14)

## [4.0.1] — 2025-06-22

### Added
- **MD trajectory support.** The plugin can now save multi-state PyMOL
  objects as per-state PDB files for CAVER MD analysis, with configurable
  state range and MD-specific UI controls (#10)

### Changed
- Reorganized monolithic source into separate modules:
  `caver_config.py`, `caver_java.py`, `caver_analysis.py`, `utils/` (#12)

## [4.0.0] — 2025-06-15

### Changed
- **PyQt5 rewrite.** Complete rewrite of the GUI layer from tkinter to PyQt5
  for compatibility with modern PyMOL versions. New UI includes tabbed
  configuration, analysis panel with tunnel spectrum visualization, and
  real-time progress tracking (#8)

## [3.0.3] — 2024-01-01  *(pre-history)*

> **Note:** Versions prior to 4.0.0 predate this changelog. The entry below
> is reconstructed from the original release notes bundled with the plugin.

### Added
- Support for PyMOL >= 2.0

### Changed
- Simplified installation procedure (patch provided by Thomas Holder,
  Schrodinger Inc.)
