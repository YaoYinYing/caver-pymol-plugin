# Dynamic Analysis of tunnels using Caver

## Output directory

Output directory is essential to any part of Caver analysis.

## Prepare reduced trajectory file

Caver analysis is a kin visualization tasks in protein structure. Thus, a proper trajectory file that has PME  (center) as well as rot+trans (fit) removed is necessary.

As large md data may crash PyMOL, a long md simulation trajectory must be reduced before get loaded.

### Gromacs

For example, this command takes one frame per 200 ps from the fitted trajectory, and produce a 500-frame version.

```bash
echo "Protein_HEM_LIG"  | gmx trjconv -s ${traj_dir}/${variant_name}/md_0_100.tpr -f md_0_100_fit.xtc -o ${output_dir}/${variant_name}_md_fit_stride.xtc -dt 200 -n ../index.ndx  # 200 ps per frame

# also copy the structure
cp ${traj_dir}/${variant_name}/md_0_100.gro ${output_dir}/${variant_name}.gro
```

## Load structure

1. Launch PyMOL and load the static structure:  `cmd.load(f'{variant}.gro',discrete=1)`
2. Remove unnecessary atoms (SOL, ions): `cmd.remove('r. SOL+NA+CL')`
3. Load stride trajectory: `cmd.load_traj(f'{variant}_md_fit_stride.xtc',)`
4. Center the protein: `cmd.center(polymer.protein)`
5. Alter chain ids for downstream analyses: `cmd.alter('polymer.protein', "chain='A'");cmd.alter('r. HEM', "chain='B'");cmd.alter('r. LIG', "chain='C'")`

## Tunnel Starting point groups

Although the static selection-convert and customized xyz-coordinates is convinient, for dynamic simulation, the real starting point may fluctuate, causing the failure on detecting tunnel in some situation.

All value list must written as space-separated format. e.g: `111 222 333 444 555` or `A:11 B:22` or `1.1 2.2 3.3`

#### Atoms

refer to `starting_point_atom`. e.g: `111 222 333 444 555`

#### Residues

refer to `starting_point_residue` e.g: `A:11 B:22`

## Residue Type to take into account

Caver requires customized residue type to get included/excluded from tunnel calculations.

This plugin gives use free options. When a new model is detected, all residue names available will be listed, including the 20 canonical amino acids and possible ligands (water, ion, glycan, small molecule, etc). 

Normally, one can click `Protein` to check all protein code. However this may contains some exceptions, such as non canonical protein codes (`HIE`, `HID`, for example). Also, one can optionally check some ligand residue names if the analysis required (`HEM`, `FAD`, for example).

## MD state range

Input the start and end frame id (state number in PyMOL). 


## Speedup

Check `Trim` if you need to clean up PyMOL session to leave more memory to caver


## Cleanup

Check `Prune` if you need to clean up `<run-id>/input` directory to save disk space. The input will not be used again.