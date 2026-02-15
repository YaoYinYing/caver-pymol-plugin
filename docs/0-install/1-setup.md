# Installation

## For the very first time

### JDK

- macOS:  run `brew install openjdk` in terminal, then add it to PATH according to the instructions output with HomeBrew.
- Windows:
- Linux: according to the distributions and package management tool

### PyMOL

A proper installation of PyMOL is required, whether from incensive or open-sourced.

Conda installation is personally recommended.

`conda install -c conda-forge 'pymol-open-source' pyqt=5 matplotlib -n Caver -y`
or 
`conda install -c schrodinger -c conda-forge 'pymol-bundle' pyqt=5 matplotlib -n Caver -y`

### Caver

1. Download [plugin](https://github.com/YaoYinYing/caver-pymol-plugin/archive/refs/heads/master.zip)
2. Start PyMOL: `conda activate Caver; pymol`
3. Open Plugin manager (Plugin -> Plugin manager)
4. Install plugin (Install new plugin -> Choose file... and select downloaded zip file)
5. Restart PyMOL


## Upgrade to a new release

In Caver PyMOL plugin, you can use `Update` button to check the release