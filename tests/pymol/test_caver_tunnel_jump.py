from __future__ import annotations

from dataclasses import dataclass

import pytest

pytest.importorskip("pymol")

from Caver4.caver_analysis import CaverAnalystPreviewer
from Caver4.utils.ui_tape import get_widget_value, set_widget_value
from tests.conftest import CaverPluginWorker
from tests.gui.analyst.test_previewer import _ensure_cached_run


@dataclass
class TunnelJumpContext:
    worker: CaverPluginWorker
    previewer: CaverAnalystPreviewer
    slider: object


@pytest.fixture
def tunnel_jump_context(caver_worker, test_data_dir) -> TunnelJumpContext:
    cache_root = _ensure_cached_run(test_data_dir)
    plugin = caver_worker.plugin

    set_widget_value(plugin.ui.lineEdit_outputDir, str(cache_root))
    plugin.config.output_dir = str(cache_root)
    plugin.ui.checkBox_EnablePlayBack.setChecked(True)
    plugin.ui.pushButton_RefreshRunID.click()
    caver_worker.process_events()

    set_widget_value(plugin.ui.comboBox_RunID, "1")
    assert str(get_widget_value(plugin.ui.comboBox_RunID)) == "1"
    caver_worker.process_events()

    plugin.ui.pushButton_analysis.click()
    caver_worker.process_events()
    analysis_ui = plugin.ui_analyst
    analysis_ui.tabWidget_analysis.setCurrentWidget(analysis_ui.tabTimeline)
    caver_worker.process_events()

    analysis_ui.pushButton_refreshTunnels.click()
    caver_worker.process_events()
    set_widget_value(analysis_ui.comboBox_tunnel, "1")
    assert str(get_widget_value(analysis_ui.comboBox_tunnel)) == "1"
    caver_worker.process_events()

    analysis_ui.pushButton_runTunnelsSpectrum.click()
    caver_worker.process_events()
    assert plugin.analyst is not None

    analysis_ui.pushButton_renderTunnelsSpectrum.click()
    caver_worker.process_events()

    analysis_ui.pushButton_refreshTunnelPreview.click()
    caver_worker.process_events()
    previewer = plugin.analyst_previewer
    assert previewer is not None

    return TunnelJumpContext(
        worker=caver_worker,
        previewer=previewer,
        slider=analysis_ui.horizontalSlider,
    )


def test_caver_tunnel_jump_requires_previewer(caver_worker, notify_box_spy):
    plugin = caver_worker.plugin
    plugin.analyst_previewer = None

    with pytest.raises(RuntimeError, match="Run tunnel preview before using caver_tunnel_jump."):
        plugin._caver_tunnel_jump("2")


def test_caver_tunnel_jump_moves_preview_slider(tunnel_jump_context):
    ctx = tunnel_jump_context
    previewer = ctx.previewer
    slider = ctx.slider

    start_frame = previewer._current_frame_id
    max_offset = previewer._max_frame_id - start_frame
    if max_offset == 0:
        pytest.skip("Preview dataset only contains a single frame.")
    step = min(5, max_offset)

    ctx.worker.plugin._caver_tunnel_jump(str(step))
    ctx.worker.process_events()
    assert previewer._current_frame_id == start_frame + step
    assert slider.value() == previewer._current_frame_id

    ctx.worker.plugin._caver_tunnel_jump(str(-step))
    ctx.worker.process_events()
    assert previewer._current_frame_id == start_frame
    assert slider.value() == start_frame
