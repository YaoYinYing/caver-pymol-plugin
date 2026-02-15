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

Although the static selection-convert and customized xyz-coordinates is convinient, 

#### Atoms

#### Residues


#### Coordinates


## Residue Type to take into account


## MD state range


## Speedup


## Cleanup