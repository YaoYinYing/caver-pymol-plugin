# Render tunnel animations

This walkthrough assumes you already processed a dynamic run and know the run/tunnel IDs you want to animate.

## Prepare the tunnel objects

Follow the usual Analyst workflow:

1. Set the **Output directory**, enable **Results playback**, click **Refresh**, and choose the run ID.
2. Click **Analysis** instead of **Load**.
3. In Analyst ▸ **Timeline setup**, pick your palette, click **Refresh**, select the tunnel, and configure diameter range/representation/spectrum.
4. Click **1. Run** then **2. Render** to build the tunnel objects.
5. Open the **Preview** group and click **Refresh** so you can scrub through the trajectory.

## Build the movie timeline

Use PyMOL's movie commands to keyframe the rendered tunnel. The script below synchronizes the movie with tunnel frames:

```python
python
num_of_frames = 500
cmd.mset(f'1 x{num_of_frames}')
for frame_id in range(1, num_of_frames + 1):
    cmd.mset(frame_id, frame_id)  # sync PyMOL state
    cmd.mdo(frame_id, f'caver_tunnel_jump_to {frame_id}; ')
    # Example camera adjustment:
    # cmd.mdo(frame_id, f'caver_tunnel_jump_to {frame_id}; turn y, {frame_id}; ')
    cmd.mview('store')
python end
```

Adjust `num_of_frames` to match your trajectory length. If you rotate/translate the view, always apply the operation inside `mdo` so the transformation is tied to the frame.

## Render the movie

Use standard PyMOL rendering commands (`mpng`, `movie.produce`, etc.) to export the animation once the timeline looks correct.
