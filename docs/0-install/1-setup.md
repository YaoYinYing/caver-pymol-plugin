# Installation

Follow this guide the first time you set up the Caver PyMOL plugin or when you need to rebuild the environment from scratch.

## Prerequisites

### JDK

Caver is a Java application, so install a recent OpenJDK (17 or later works well).

- **macOS** – `brew install openjdk` and follow Homebrew's post-install message to add it to `PATH`.
- **Windows** – install an OpenJDK distribution such as Temurin or Zulu and add `java` to the `PATH` environment variable.
- **Linux** – use your distribution package manager (`apt install openjdk-17-jre`, `dnf install java-17-openjdk`, etc.).

### PyMOL environment

Install PyMOL with Qt and matplotlib so the plugin UI renders correctly. Conda keeps the dependency chain simple:

```bash
conda create -n Caver -y
conda install -n Caver -c conda-forge 'pymol-open-source' pyqt=5 matplotlib
# or, for the Schrödinger build
conda install -n Caver -c schrodinger -c conda-forge 'pymol-bundle' pyqt=5 matplotlib
```

Activate the environment whenever you run PyMOL: `conda activate Caver`.

## Install the plugin

1. Download the latest [Caver PyMOL plugin](https://github.com/YaoYinYing/caver-pymol-plugin/archive/refs/heads/master.zip).
2. Start PyMOL inside the `Caver` conda environment.
3. Open **Plugin ▸ Plugin Manager**.
4. Choose **Install New Plugin ▸ Choose File…** and select the downloaded ZIP file.
5. Restart PyMOL so the plugin menu entries load.

## Upgrade to a new release

Open the plugin window and click **Upgrade** to check new release and lead user to GitHub release page to download and install the most recent release directly from PyMOL.
