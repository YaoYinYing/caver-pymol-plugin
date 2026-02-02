from __future__ import annotations

import pytest

pytest.importorskip("pymol")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from Caver4.utils.ui_tape import list_color_map


def test_combobox_matches_matplotlib_colormaps(caver_worker):
    plugin = caver_worker.plugin
    plugin.ui.pushButton_analysis.click()
    caver_worker.process_events()

    analysis_ui = plugin.ui_analyst
    analysis_ui.tabWidget_analysis.setCurrentWidget(analysis_ui.tabPlot)
    caver_worker.process_events()

    combo = analysis_ui.comboBox_plotColormap
    expected = list_color_map()

    assert combo.count() == len(expected)
    assert [combo.itemText(idx) for idx in range(combo.count())] == list(expected.keys())
