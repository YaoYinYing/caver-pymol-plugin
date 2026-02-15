# CAVER PyMOL plugin

The CAVER PyMOL plugin enables calculation and visualization of tunnels in
PyMOL. For the calculation of tunnels, the plugin utilizes CAVER 3.02 software
package.

Original Website: <https://www.caver.cz>

## Features of this fork
- Support modern PyMOL driven by PyQt
- Support MD analysis
- Flexibler starting point and residue inclusions
- Results playback
- Extended PyMOL commandline prompt `caver_set`
- Extended PyMOL commandline prompt `caver_tunnel_jump` and `caver_tunnel_jump_to`
- Caver Tunnel Analyst: analyst, visualize, plot (requires `matplotlib`)
- Tunnel Movie

## Image Preview

![](https://github-image-cache.yaoyy.moe/2026/02/da13013716666134c9963de705f59eb1.png)
(A). Refactored pluggin interface. (B). Configuration window. (C). Starting point selection and refines. (D). Starting point from residues, atons or coordinates. (E). Tunnel MD analysis. (F, G). `caver_set` auto-completion for (E) options and (G) values.

![](https://github-image-cache.yaoyy.moe/2026/02/4c7c05ea846b7df446300a2904c2cf93.png)
(A) Analyst window for per-tunnel analysis. (B). Tunnels rendered as lines. (C). Timeline view of tunnel evolution. (D). Per-frame view of tunnel. (E). Tunnel plot settings. (F). Publication-quality tunnel Plot.

## Movie Preview

![](https://github-image-cache.yaoyy.moe/2026/02/5d7683460ae932e9ee695072eac24f87.gif)

## Contributors:
* Human Computer Interaction Laboratory Faculty of Informatics, Masaryk University
* Loschmidt Laboratories, Department of Experimental Biology and Research Centre for Toxic Compounds in the Environment Faculty of Science, Masaryk University
* Thomas Holder, Schrodinger Inc.
* Yinying Yao, Huazhong Agricultural University


## Requirements
* PyMOL  >=  2.5
* OpenJDK

## Installation instructions

See the [installation instruction](docs/0-install/1-setup.md) for more detailed information.

1. Install OpenJDK
  - For macOS:  run `brew install openjdk` in terminal, then add it to PATH according to the instructions output with HomeBrew.
2. Download [plugin](https://github.com/YaoYinYing/caver-pymol-plugin/archive/refs/heads/master.zip)
3. Start PyMOL
4. Open Plugin manager (Plugin -> Plugin manager)
5. Install plugin (Install new plugin -> Choose file... and select downloaded zip file)
6. Restart PyMOL

## License

GNU General Public License, version 3 (GPL-3.0)
