# Dynamic analysis with Caver

Dynamic (trajectory) runs follow the same structure as static ones but require extra preparation to keep PyMOL fast and memory-safe.

## 1. Choose the output directory

Every run writes into `output/<run-id>`. Set this field before loading data so the plugin can discover results later for playback, plotting, or movies.

## 2. Prepare a reduced trajectory

PyMOL renders best with pre-centered, fitted, and decimated trajectories. Remove global rotation/translation and keep only the atoms needed for the analysis. Reduce very long MD simulations to a manageable frame count before loading them.

### Example (GROMACS)

The snippet below extracts one frame every 200 ps from a fitted trajectory and copies the matching structure file:

```bash
echo "Protein_HEM_LIG" | gmx trjconv \
  -s ${traj_dir}/${variant_name}/md_0_100.tpr \
  -f md_0_100_fit.xtc \
  -o ${output_dir}/${variant_name}_md_fit_stride.xtc \
  -dt 200 -n ../index.ndx

cp ${traj_dir}/${variant_name}/md_0_100.gro ${output_dir}/${variant_name}.gro
```

## 3. Load the structure and trajectory

Inside PyMOL:

1. `cmd.load(f"{variant}.gro", discrete=1)`
2. Remove solvent/ions you do not need: `cmd.remove('r. SOL+NA+CL')`
3. `cmd.load_traj(f"{variant}_md_fit_stride.xtc")`
4. Center the protein: `cmd.center('polymer.protein')`
5. Normalize chain IDs if desired: `cmd.alter('polymer.protein', "chain='A'")`, `cmd.alter('r. HEM', "chain='B'")`, etc.

## 4. Define starting-point groups

Frame-by-frame fluctuations make automatic selection conversion less reliable than in static runs. Provide redundant inputs (atoms plus residues) whenever possible. Use space-separated lists such as `111 222 333` or `A:11 A:22`.

- `starting_point_atom` – e.g., `111 222 333 444 555`
- `starting_point_residue` – e.g., `A:11 B:22`
- `starting_point_coordinates` – absolute XYZ values if you have a reference point

## 5. Pick residue types to include

Same rules as the static workflow: click **Protein** to grab canonical residues, then add non-standard codes (`HIE`, `HID`, …) or ligands (`HEM`, `FAD`, …) that form tunnel walls.

## 6. Set the MD state range

Provide the starting and ending PyMOL state numbers you want Caver to process. Limiting the range speeds up runs and helps when you only need part of the trajectory.

## 7. Optional speed/cleanup toggles

- **Trim** – removes temporary PyMOL objects before launching Caver to free memory.
- **Prune** – deletes `<run-id>/input` after the run finishes to save disk space (the raw inputs are not reused once processing completes).
