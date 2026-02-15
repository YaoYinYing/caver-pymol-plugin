# Preview tunnels

Use the Analyst timeline to inspect tunnel geometry across MD frames, step through results, and recolor surfaces on the fly.

## Launch Analyst

> Instead of clicking **Load**, choose **Analysis** in the playback section to open the Analyst window.

1. Set the correct **Output directory**.
2. Enable **Results playback** and click **Refresh**.
3. Select the run ID corresponding to your dataset.
4. Click **Analysis**.

## Timeline setup

1. Switch to the **Timeline** tab inside Analyst.
2. In **Timeline setup**, pick a color palette.
3. Click **Refresh** to populate the tunnel list.
4. Select the tunnel you want to study.
5. Configure the diameter range. Water channels often work well with `[1.5, 2.0]`; ligand channels depend on the ligand size.
6. Choose the representation (mesh, surface, spheres, …).
7. Set the spectrum source. Using `b` or `vdw` encodes radius/bottleneck variation across the tunnel.
8. Click **1. Run** to cache data.
9. Click **2. Render** to build the tunnel objects for every frame.

## Preview timeline

1. Switch to the **Preview** group.
2. Hit **Refresh** to initialize the controls.
3. Use the playback bar (first/previous/play/pause/next/last) or drag the slider to move between frames. Typing a frame number in the spin box jumps directly to that state.
4. The `?` button opens a detail dialog for the current tunnel frame.

## Spectrum vs. ramp coloring

- **Spectrum** – maps arbitrary model properties (`resi`, `b`, `vdw`, etc.) to colors.
- **Ramp** – colors by distance from a reference object.

### Quick ramp example

1. `create startpoint, sele` – define the origin.
2. `cmd.ramp_new('r', 'startpoint', [0, 30], 'rainbow')`
3. Enter `r` into the **Spect.** combo box and click **2. Render**.
4. Re-running `ramp_new` updates the colors instantly without rebuilding the tunnels.
