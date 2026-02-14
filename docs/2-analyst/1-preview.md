# Preview tunnel

## Load structure


## Enable Playback and discover the run ids

> [!WARNING] Instead of clicking `Load` for raw result from Caver, click `Analysis` to launch Analyst.


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
11. Now navigate to the `Preview` group.


### Tunnel MD Previewer

1. Click `Refresh` to get previewer ready.
2. Use play buttons (the first, the previous, auto-play, pause auto-play, the next, the last) as well as adjust auto-play interval to preview the time evolution of selected tunnel. Also, drag the slider or input frame id into the spinbox next to the slider will take similar effects. 
3. To check the detail of a tunnel frame, click `?` to take a look.
   

### Ramp or spectrum

- Spectrum: use any model properties (resi, b, vdw, etc.) to color
- Ramp: use only distance proximals.
  

#### Quick example

To render tunnel by the distance to some startpoint:
1. Use `create startpoint, sele` to create an origin start point object
2. Create an ramp object: `cmd.ramp_new('r', 'startpoint', [0, 30], 'rainbow')`
3. Input `r` into the `Spect.` combobox.
4. Click `2. Render to render`.
5. Post adjustment: rerun `ramp_new` will automaticcally change the colors. no need to re-render.
