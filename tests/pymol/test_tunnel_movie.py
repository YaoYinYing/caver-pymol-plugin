from __future__ import annotations

import pytest

pytest.importorskip("pymol")
import shutil
import subprocess

from Caver4.utils.ui_tape import set_widget_value
from tests.gui.analyst.test_previewer import _ensure_cached_run


def _apply_basic_style(cmd) -> None:
    cmd.set("cartoon_color", "gray70")
    cmd.set("cartoon_transparency", 0.9)
    cmd.hide("sticks", "hydrogens")
    cmd.bg_color("white")
    cmd.origin("all")


def test_tunnel_movie_generation(caver_worker, test_data_dir, results_root):
    worker = caver_worker
    plugin = worker.plugin
    cmd = worker.cmd

    session_path = test_data_dir / "md_snapshots" / "caver_md.snapshots.pze"
    assert session_path.is_file()
    cmd.load(str(session_path))
    worker.process_events()
    _apply_basic_style(cmd)

    cache_root = _ensure_cached_run(test_data_dir)
    set_widget_value(plugin.ui.lineEdit_outputDir, str(cache_root))
    plugin.config.output_dir = str(cache_root)
    plugin.ui.checkBox_EnablePlayBack.setChecked(True)
    plugin.ui.pushButton_RefreshRunID.click()
    worker.process_events()
    set_widget_value(plugin.ui.comboBox_RunID, "1")
    worker.process_events()

    plugin.ui.pushButton_analysis.click()
    worker.process_events()
    analysis_ui = plugin.ui_analyst
    analysis_ui.tabWidget_analysis.setCurrentWidget(analysis_ui.tabTimeline)
    worker.process_events()

    analysis_ui.doubleSpinBox_spectrumMin.setValue(1.0)
    analysis_ui.doubleSpinBox_spectrumMax.setValue(3.0)
    set_widget_value(analysis_ui.comboBox_representation, "mesh")
    set_widget_value(analysis_ui.comboBox_spectrumPalette, "green_red")

    analysis_ui.pushButton_refreshTunnels.click()
    worker.process_events()
    set_widget_value(analysis_ui.comboBox_tunnel, "1")
    worker.process_events()

    analysis_ui.pushButton_runTunnelsSpectrum.click()
    worker.process_events()
    analysis_ui.pushButton_renderTunnelsSpectrum.click()
    worker.process_events()
    analysis_ui.pushButton_refreshTunnelPreview.click()
    worker.process_events()

    previewer = plugin.analyst_previewer
    assert previewer is not None
    previewer.head()
    worker.process_events()

    total_frames = 24
    fps = 24
    assert previewer._max_frame_id - previewer._min_frame_id + 1 >= total_frames

    cmd.set("ray_trace_frames", 0)
    cmd.set("ray_trace_mode", 0)
    cmd.set("cache_frames", 1)

    cmd.orient("cl_000*")

    movie_dir = results_root / "movie"
    frames_dir = movie_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    for leftover in frames_dir.glob("*.png"):
        leftover.unlink()
    movie_path = movie_dir / "tunnel_movie.mp4"
    if movie_path.exists():
        movie_path.unlink()

    for frame in range(1, total_frames + 1):
        if frame > 1:
            cmd.do("caver_tunnel_jump 1")
            worker.process_events()
            cmd.refresh()
        cmd.turn("y", 1)
        cmd.refresh()
        frame_file = frames_dir / f"tunnel_movie{frame:04d}.png"
        cmd.png(str(frame_file), ray=0, quiet=0, dpi=150, width=400, height=300)

    assert previewer._current_frame_id == previewer._min_frame_id + total_frames - 1

    frames = sorted(frames_dir.glob("tunnel_movie*.png"))
    assert len(frames) == total_frames
    assert frames[0].is_file()
    assert frames[-1].is_file()

    ffmpeg_bin = shutil.which("ffmpeg")
    if not ffmpeg_bin:
        pytest.skip("ffmpeg is required to encode MP4 output")

    subprocess.run(
        [
            ffmpeg_bin,
            "-y",
            "-framerate",
            str(fps),
            "-i",
            str(frames_dir / "tunnel_movie%04d.png"),
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            str(movie_path),
        ],
        check=True,
        capture_output=True,
    )
    assert movie_path.is_file()
    assert movie_path.stat().st_size > 0
