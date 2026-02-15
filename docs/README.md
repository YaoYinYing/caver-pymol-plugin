# Documentation Overview

This directory describes the complete workflow for running the Caver PyMOL plugin, verifying the results, and exporting publication-ready media. Each guide focuses on a discrete stage so you can jump directly to the topic you need or work through them sequentially for a full walkthrough.

## Quick navigation

1. `docs/0-install/1-setup.md` – Install Java, PyMOL, and the Caver PyMOL plugin. Includes the recommended conda environment and the in-app update flow.
2. `docs/1-caver/1-static.md` – Configure Caver for static structures, launch calculations, and understand tunnel tables.
3. `docs/1-caver/2-dynamic.md` – Run ensemble/dynamic analyses, including how to import trajectories and switch between frames efficiently.
4. `docs/1-caver/3-playback.md` – Replay calculated tunnels in PyMOL, adjust visualization parameters, and capture screenshots.
5. `docs/2-analyst/1-preview.md` – Open results inside Caver Analyst for rapid previewing, filtering, and clustering.
6. `docs/2-analyst/2-plot.md` – Produce tunnel plots, annotate bottlenecks, and export vector graphics.
7. `docs/3-tunnel_movie/1-get_started.md` – Build tunnel movies from PyMOL sessions, add labels, and render final videos.

## Recommended workflow

1. Follow the installation guide to prepare a working environment.
2. Choose between the static or dynamic tutorials depending on your data. The static guide assumes a single PDB/CIF structure, while the dynamic guide walks through multi-frame trajectories.
3. Use the playback tutorial to inspect tunnels interactively inside PyMOL once calculations finish.
4. Switch to the Analyst tutorials when you need advanced filtering, clustering, or publication graphics.
5. Finish with the tunnel movie guide to generate narrated animations for presentations.

Each tutorial contains prerequisites at the top, screenshots when needed, and command snippets you can copy directly into your shell or the PyMOL console. Keep the README open as a roadmap while you work through the detailed guides.
