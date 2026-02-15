# Static Analysis of tunnels using Caver



## Output directory

Output directory is essential to any part of Caver analysis.

## Load structure

Load a normal structure model into PyMOL session and launch Caver Plugin. Ensure the correct model is selected in `Input Model`

## Tunnel starting point

The Caver PyMOL plugin provides two ways to input the starting point.

1. (Tab `Select` and `Refine`) Convert and  refine xyz coords from a PyMOL selection object
2. (Tab `Custom`) Direct input as atom ids/residue ids/xyz coordinates

It's recommend to use selection-coordinates for static structure analysis.

The priority of selection-convert is higher than direct inputs. To apply customized overriding, use `Clear` button at `Refine` tab to set all of xyz coordinates to `0`.

### XYZ coordinates From Selection

The COM (center of mass) coordinates of a selection will be converted to `starting_point_coordinates` options.

### Starting point groups

All value list must written as space-separated format. e.g: `111 222 333 444 555` or `A:11 B:22` or `1.1 2.2 3.3`


#### Atoms

refer to `starting_point_atom`.

#### Residues

refer to `starting_point_residue`


#### Coordinates

refer to `starting_point_coordinates`


## Residue Type to take into account

Caver requires customized residue type to get included/excluded from tunnel calculations.

This plugin gives use free options. When a new model is detected, all residue names available will be listed, including the 20 canonical amino acids and possible ligands (water, ion, glycan, small molecule, etc). 

Normally, one can click `Protein` to check all protein code. However this may contains some exceptions, such as non canonical protein codes (`HIE`, `HID`, for example). Also, one can optionally check some ligand residue names if the analysis required (`HEM`, `FAD`, for example).

## Run Caver

After all is set, click `Run` to process the task. The window will be disabled for a while during the processing, and the result shall be loaded into PyMOL also.