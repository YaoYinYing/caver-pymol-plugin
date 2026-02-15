# Static analysis with Caver

Use this guide when you only need tunnel profiles for a single structure (PDB/CIF). The workflow is linear: set an output directory, choose the model, define the starting point, and launch the calculation.

## 1. Pick the output directory

All results live in the output directory. Set it first so every subsequent step saves into the same run folder.

## 2. Load the structure

Load the structure inside PyMOL and open the Caver plugin. The structure you want to analyse must be selected in **Input model**.

## 3. Define the tunnel starting point

Two inputs are available:

1. **Select ▸ Refine** – convert a PyMOL selection into XYZ coordinates (center of mass). This is the most reproducible option for static structures.
2. **Custom** – type atom IDs, residue IDs, or XYZ values manually if you need complete control.

Selection-derived coordinates override manual entries. Click **Clear** on the **Refine** tab to zero out XYZ values before entering custom coordinates.

### Input formats

Lists are space-separated. Examples: `111 222 333`, `A:11 B:22`, or `1.2 2.4 3.6`.

- **Atoms** – populates `starting_point_atom`.
- **Residues** – populates `starting_point_residue`.
- **Coordinates** – populates `starting_point_coordinates`.

## 4. Choose residue types to include

Caver needs to know which residue types participate in tunnel calculations. When you load a model the plugin lists every residue code it finds: standard amino acids plus ligands (waters, ions, cofactors, glycans, small molecules, etc.).

- Click **Protein** to select canonical residues, then add any uncommon ones (`HIE`, `HID`, …) manually.
- Include ligands that should form part of the tunnel wall (`HEM`, `FAD`, …) if needed.

## 5. Run Caver

Click **Run** and wait for PyMOL to re-enable the window. The plugin automatically loads the generated tunnels back into the session for inspection or playback.
