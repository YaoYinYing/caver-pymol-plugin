# Plot tunnel bottlenecks

The Plot tab recreates Caver Analyst's `profile_heat_maps` inside PyMOL, letting you export tunnel-diameter timelines as publication-ready figures.

## Launch Analyst

Follow the same steps listed in the [preview guide](../2-analyst/1-preview.md) to open Analyst and populate the **Timeline** tab. Clicking **1. Run** ensures the tunnel data are cached for plotting.

## Configure the plot

> The Plot tab requires `matplotlib`. If PyMOL displays “Please install matplotlib before using this feature” in the color selector, install it in your environment.

1. Switch to the **Plot** tab.
2. Set the tunnel length range to match the minimum/maximum path distance you care about.
3. Choose a color palette for diameter values. The default red→green scheme mimics the original Analyst output, but any matplotlib colormap works.
4. Adjust the figure size and DPI for your target medium (slides vs. print).
5. Optionally add axis transformations, labels, and ranges:
   - Use `$t` for the raw value, e.g., `$t/5` converts frame indices back to nanoseconds if each frame equals 0.2 ns.
   - Set axis labels to describe the transformed quantity (`MD simulation time (ns)`, `Tunnel length (Å)`, etc.).
6. Click **Plot** to open the matplotlib window. If nothing appears, force Qt usage:

   ```python
   python
   import matplotlib
   matplotlib.use('Qt5Agg')
   python end
   ```

7. Tweak fonts or layout in the matplotlib GUI, then save the figure to PNG, PDF, or SVG.
