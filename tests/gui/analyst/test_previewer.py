from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import zipfile

import pytest

from tests.conftest import CaverPluginWorker

pytest.importorskip("pymol")

from Caver4.caver_analysis import CaverAnalystPreviewer
from Caver4.ui.Ui_caver_analysis import Ui_CaverAnalyst as CaverAnalysisForm
from Caver4.utils.ui_tape import get_widget_value, set_widget_value


def _ensure_cached_run(test_data_dir: Path) -> Path:
    cache_root = test_data_dir / "cache"
    run_dir = cache_root / "caver_output" / "1"
    csv_file = run_dir / "analysis" / "profile_heat_maps" / "csv" / "cl_000001_heat_map.csv"
    pdb_file = run_dir / "data" / "clusters_timeless" / "tun_cl_001_1.pdb"
    if csv_file.is_file() and pdb_file.is_file():
        return cache_root

    bundle = run_dir / "Caver4.md_results.zip"
    if not bundle.is_file():
        pytest.skip("Missing cached tunnel dataset bundle.")

    with zipfile.ZipFile(bundle, "r") as archive:
        archive.extractall(run_dir)
    return cache_root


@dataclass
class PreviewerTestContext:
    worker: CaverPluginWorker
    analysis_ui: CaverAnalysisForm
    previewer: CaverAnalystPreviewer
    cache_root: Path


@pytest.fixture
def analyst_previewer_context(caver_worker, test_data_dir) -> PreviewerTestContext:
    cache_root = _ensure_cached_run(test_data_dir)
    plugin = caver_worker.plugin
    set_widget_value(plugin.ui.lineEdit_outputDir, str(cache_root))
    plugin.config.output_dir = str(cache_root)
    plugin.ui.checkBox_EnablePlayBack.setChecked(True)
    plugin.ui.pushButton_RefreshRunID.click()
    caver_worker.process_events()
    set_widget_value(plugin.ui.comboBox_RunID, "1")
    assert str(get_widget_value(plugin.ui.comboBox_RunID)) == "1"

    plugin.ui.pushButton_analysis.click()
    caver_worker.process_events()
    analysis_ui = plugin.ui_analyst
    analysis_ui.tabWidget_analysis.setCurrentWidget(analysis_ui.tabTimeline)
    caver_worker.process_events()

    analysis_ui.pushButton_refreshTunnels.click()
    caver_worker.process_events()
    set_widget_value(analysis_ui.comboBox_tunnel, "1")
    assert str(get_widget_value(analysis_ui.comboBox_tunnel)) == "1"

    analysis_ui.pushButton_applyTunnelsSpectrumStatic.click()
    caver_worker.process_events()
    assert plugin.analyst is not None

    analysis_ui.pushButton_refreshTunnelPreview.click()
    caver_worker.process_events()
    assert plugin.analyst_previewer is not None

    return PreviewerTestContext(
        worker=caver_worker,
        analysis_ui=analysis_ui,
        previewer=plugin.analyst_previewer,
        cache_root=cache_root,
    )


def test_previewer_navigation_and_slider(analyst_previewer_context):
    ctx = analyst_previewer_context
    ui = ctx.analysis_ui
    previewer = ctx.previewer
    slider = ui.horizontalSlider

    assert slider.minimum() == previewer._min_frame_id
    assert slider.maximum() == previewer._max_frame_id
    assert previewer._current_frame_id == previewer._min_frame_id
    assert not ui.pushButton_firstFrame.isEnabled()
    assert not ui.pushButton_previousFrame.isEnabled()
    assert ui.pushButton_lastFrame.isEnabled()
    assert ui.pushButton_nextFrame.isEnabled()

    ui.pushButton_lastFrame.click()
    ctx.worker.process_events()
    assert slider.value() == previewer._max_frame_id
    assert previewer._current_frame_id == previewer._max_frame_id
    assert not ui.pushButton_lastFrame.isEnabled()
    assert not ui.pushButton_nextFrame.isEnabled()

    ui.pushButton_firstFrame.click()
    ctx.worker.process_events()
    assert slider.value() == previewer._min_frame_id
    assert previewer._current_frame_id == previewer._min_frame_id

    ui.pushButton_nextFrame.click()
    ctx.worker.process_events()
    assert slider.value() == previewer._min_frame_id + 1
    assert ui.pushButton_previousFrame.isEnabled()

    ui.pushButton_previousFrame.click()
    ctx.worker.process_events()
    assert slider.value() == previewer._min_frame_id

    mid_offset = min(10, previewer._max_frame_id - previewer._min_frame_id)
    mid_frame = previewer._min_frame_id + mid_offset
    slider.setValue(mid_frame)
    ctx.worker.process_events()
    assert previewer._current_frame_id == mid_frame


def test_previewer_autoplay_about_and_reapply(analyst_previewer_context, notify_box_spy):
    ctx = analyst_previewer_context
    ui = ctx.analysis_ui
    previewer = ctx.previewer

    ui.doubleSpinBox_autoPlayInterval.setValue(0.01)
    ui.pushButton_autoPlay.click()
    ctx.worker.process_events()
    assert previewer._is_autoplay_running()
    assert not ui.pushButton_firstFrame.isEnabled()
    assert ui.pushButton_pauseAutoPlay.isEnabled()

    expected_frame = previewer._next_frame_id(previewer._current_frame_id)
    previewer._auto_play_tick()
    ctx.worker.process_events()
    assert previewer._current_frame_id == expected_frame

    ui.pushButton_pauseAutoPlay.click()
    ctx.worker.process_events()
    assert not previewer._is_autoplay_running()
    assert ui.pushButton_autoPlay.isEnabled()


    # TODO: the notify box is not mocked and popped up. fix it
    # prior_messages = len(notify_box_spy.messages)
    # ui.pushButton_aboutThisFrame.click()
    # ctx.worker.process_events()
    # assert len(notify_box_spy.messages) == prior_messages + 1
    # about_text = notify_box_spy.messages[-1][0]
    # assert "TunnelFrame #" in about_text
    # assert f"#{previewer._current_frame_id}" in about_text

    repr_combo = ui.comboBox_representation
    spec_combo = ui.comboBox_spectrumBy

    def _select_alternate(combo):
        for idx in range(combo.count()):
            if combo.itemText(idx) != combo.currentText():
                combo.setCurrentIndex(idx)
                return combo.itemText(idx)
        return combo.currentText()

    new_repr = _select_alternate(repr_combo)
    new_expr = _select_alternate(spec_combo)
    previous_analyst = ctx.worker.plugin.analyst
    ui.pushButton_applyTunnelsSpectrumStatic.click()
    ctx.worker.process_events()
    assert ctx.worker.plugin.analyst is not None
    assert ctx.worker.plugin.analyst is not previous_analyst
    assert repr_combo.currentText() == new_repr
    assert spec_combo.currentText() == new_expr
