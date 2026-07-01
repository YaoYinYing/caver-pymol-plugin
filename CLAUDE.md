# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run full test suite
make prepare-test && make test

# Run tests matching a keyword (e.g. "config", "plotter", "tunnel_jump")
make kw-test PYTEST_KW='keyword'

# Format code (isort + black via pre-commit)
make black

# Install plugin to PyMOL startup directory
make install

# CI test run (tolerates segfault exit codes on macOS ARM64)
make test-skip-ci-segfault
```

All test commands use pytest with coverage (`--cov-config=.coveragerc`). Tests require a PyMOL + Qt environment — CI uses conda + Xvfb for headless display.

## Architecture

This is a **PyMOL plugin** (not a pip-installable package) distributed as a ZIP and installed via PyMOL's Plugin Manager. The plugin entry point is `Caver4/__init__.py` which registers a "Caver NG" menu item.

### Key modules in `Caver4/`

| Module | Responsibility |
|--------|---------------|
| `caver_pymol.py` | Main plugin class (~1320 lines) — GUI dialog, PyMOL integration, command registration (`caver_set`, `caver_tunnel_jump`, `caver_tunnel_jump_to`), workflow orchestration |
| `caver_config.py` | Config serialization/deserialization — reads/writes CAVER text config files, manages type casting ("yes"/"no" ↔ bool) |
| `caver_analysis.py` | Tunnel analysis results — `CaverAnalyst` renders tunnels in PyMOL, `CaverAnalystPlotter` makes matplotlib plots, `CaverAnalystPreviewer` provides per-frame tunnel navigation |
| `caver_java.py` | Java bridge — locates JDK, builds classpath, runs CAVER 3.02 JAR as subprocess, parses progress |
| `utils/ui_tape.py` | Qt abstraction over PyMOL's bundled PyQt5 (`pymol.Qt`) — widget value get/set, signals, threads, colormap icons |
| `utils/live_run.py` | Subprocess runner with real-time stdout/stderr capture via threading |
| `utils/caver_utils.py` | Atom/residue filtering helpers |

### Data Flow

```
User Input (GUI) → CaverPyMOL → CaverConfig.to_txt() → config file
                                    ↓
                              PyJava.run() → CAVER 3.02 JAR (subprocess)
                                    ↓
                              Output dir (tunnels CSV/PDB)
                                    ↓
                              CaverAnalyst → render tunnels in PyMOL
                              CaverAnalystPlotter → matplotlib plots
                              CaverAnalystPreviewer → per-frame tunnel view
```

### UI Files (`Caver4/ui/`)
- `Ui_caver.py`, `Ui_caver_analysis.py`, `Ui_caver_config.py` — auto-generated PyQt5 UI classes from `.ui` designer files

### CAVER Algorithm Core (`Caver4/bin/`)
- `caver.py` — sphere computation, tunnel detection, clustering (standalone Python port)
- `view.py`, `view_plugin.py`, `view_timeless.py` — tunnel visualization in PyMOL
- `zones.py`, `rgb.py`, `void_template.py` — supporting definitions

### Test Structure (`tests/`)

```
tests/
  conftest.py              # CaverPluginWorker (full E2E driver), notify_box_spy fixture, cleanup_pymol
  test_version.py          # Version string validation
  data/test_data.py        # STATIC_WORKFLOW / DYNAMIC_WORKFLOW scenarios
  config/test_caver_config.py
  java/test_caver_java.py
  pymol/test_caver_set.py, test_caver_tunnel_jump.py, test_tunnel_movie.py
  gui/test_caver_pymol.py  # Minimal placeholder
  gui/analyst/test_plotter.py, test_previewer.py
  analysis/test_analysis_utils.py
  utils/test_colormap.py
  utils/gui/test_matplotlib_color_map.py
  cache/                   # Pre-computed CAVER output for regression testing
```

Test infrastructure: `CaverPluginWorker` in conftest.py initializes PyMOL, loads structures, runs analysis, captures screenshots. `notify_box_spy` replaces modal dialogs with a recorder. `cleanup_pymol` prevents macOS ARM64 segfaults.

## Dependencies

| Dependency | Role |
|------------|------|
| **PyMOL** (≥ 2.5, conda-forge) | Molecular visualization engine, bundled PyQt5 (`pymol.Qt`), `pymol.cmd` API |
| **OpenJDK** (8+) | Runs CAVER 3.1.0 engine JAR — CI uses JDK 21, install docs recommend 17+ |
| **matplotlib** | Publication-quality tunnel plots and colormap icons |
| **pytest**, **pytest-qt**, **pytest-cov** (dev) | Test suite |
| **pre-commit**, **isort**, **black** (dev) | Code formatting |

**PyMOL is imported with Qt via `pymol.Qt`** — never import PyQt5 directly. All Qt widgets use `from pymol.Qt import ...`.

## CI

GitHub Actions matrix: Ubuntu/Windows/macOS × Python 3.10–3.13 × multiple PyMOL variants. Conda-based environment setup. Headless display via Xvfb (Linux) or setup-headless-display-action (Windows/macOS).

## Version

Defined in `Caver4/__init__.py` as `VERSION = "4.2.1"`. License: GNU GPL v3.