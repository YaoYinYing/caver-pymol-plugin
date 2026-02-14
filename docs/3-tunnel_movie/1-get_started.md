# Rendering tunnel animations into a movie


## Note

This doc assumes you have already load an MD trajectory into PyMOL and get a MD trajectory analysis processed.

Also, this doc assumes you have already known the run id and tunnel id you wish to process in this tutorial.

## UI Operations

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

## Quick example for moview scripting

```python
# tell PyMOL command line prompts to inter a python script
python

# set the total frame numbers
num_of_frames=500

# create the moview w/ the frame number
cmd.mset(f'1 x{num_of_frames}')

# loop these frames and set key frames
for frame_id in range(1,1+num_of_frames):
  # set the state as the frame
  cmd.mset(frame_id, frame_id)
  # in the frame, jump to the tunnel frame
  cmd.mdo(frame_id,f'caver_tunnel_jump_to {frame_id}; ')

  # optionaly, to apply any after-jump operations, for example turn camera view, remember 
  # that after state switch, the view will **always** get reset to the initial, that's 
  # why you should use `turn y <frame_id>` instead of `turn y, 1`.
  
  # cmd.mdo(frame_id,f'caver_tunnel_jump_to {frame_id}; turn y, {frame_id}; ')

  # store it
  cmd.mview('store')

# that is it!

# tell PyMOL command line prompts to exit a python script
python end
```

Now the movie sequence is ready to be rendered.

## Dump the movie

Render the movie like what you have done to the official movie tutorials.