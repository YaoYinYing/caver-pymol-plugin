# Plot Tunnel for bottlenecks

## Concept

Plotter is a re-implement of Caver's `profile_heat_maps`.

### Tunnel Runner

1. Ensure the `output directory` is set.
2. Enable `results playback` and click `Refresh` to discover results that already processed. It does not matter whether the input model is correct or not. Tunnel Analyst and Previewer will not use it.
3. Pick the correct run id that **matches** the model and trajectory.
4. Click `Analysis` to open Analyst window.

### Tunnel MD Analyst

1. The Analyst window has to tabs. For tunnel preview and timeline tasks, use `Timeline` tab.
2. Navigate to the `Timeline setup` group. 
3. Select color palette you want to use.
4. Click `Refresh` to read the results under run-id directory.
5. Pick the desired tunnel id that you wish to visualize.
6. Setup tunnel color range. For water channel, `[1.5, 2]` is okay. For ligand channel, it depends on what the diameter you need to adjust to.
7. Pick a representation method. I personally use mesh.
8. Pick what way the spectrum will run against. `b` or `vdw` is sufficient for normal way to show the diameters and bottlenecks.
9. Click `1. Run` to load the tunnel. This steps takes no effects on GUI or PyMOL.
10. Click `2. Render` to render the tunnel one-after-another. This step takes a while, which depends on the frame number of trajectory. Now the Previewer is activated.
11. Now switch to `plot` tab.

### Tunnel Plotter

> [!WARNING] Tunnel plotter requires `matplotlib`, which is not a part of PyMOL.
> You should have got message like `"Please install matplotlib before using this feature"` in `Color` combobox if `matplotlib` is not ready.

1. Reset tunnel range to meet the minimum and maximum of tunnel length. Adjust them if you need.
2. Select a proper color for tunnel diameters. default is red - green to fit the original palette.
3. If you need to define colors mapping to tunnel diameters, back to tab `Timeline` and adjust `Range` in group `Timeline setup`.
4. Adjust picture size and DPI if you need to do that. 
5. Optionally, set x/y transformations, labels and ranges(after transformation).
   1. Transformation: Inspired by XMGrace. `$t` refers to the original value. `$t/5` to divide the value by 5. 
   2. Transformation is helpful if you need to convert `frame index` (from caver trajectory) back to `ns(nanosecond)` (used by MD program).
   3. Quick example: x transformation `$t/5`, x label set to `MD Simulation Time (ns)`.
6. Click `Plot` to open matplotlib GUI. if no window opens, use the following and plot again:

    ```python
    # force matplotlib to use Qt backend
    python

    import matplotlib

    matplotlib.use('Qt5Agg')

    python end
    ```

7. Adjust some details from the new figure window. Save the final plot w/ proper filenames.