from __future__ import annotations

import zipfile
from dataclasses import dataclass
from pathlib import Path

import pytest

from tests.conftest import CaverPluginWorker

pytest.importorskip("pymol")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from Caver4.caver_analysis import CaverAnalystPlotter
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
class PlotterTestContext:
    worker: CaverPluginWorker
    plotter: CaverAnalystPlotter
    analysis_ui: CaverAnalysisForm
    cache_root: Path
    snapshot_path: Path


@pytest.fixture
def analyst_plotter_context(caver_worker, test_data_dir) -> PlotterTestContext:
    """
    Prepare the analyst UI so the plotter can operate on cached tunnel data.
    """

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

    snapshot_path = caver_worker.capture_gui_snapshot("analyst_plotter_ready.png")
    assert snapshot_path.is_file()

    analysis_ui.pushButton_applyTunnelsSpectrumStatic.click()
    caver_worker.process_events()
    assert plugin.analyst is not None
    assert plugin.analyst_plotter is not None

    return PlotterTestContext(
        worker=caver_worker,
        plotter=plugin.analyst_plotter,
        analysis_ui=analysis_ui,
        cache_root=cache_root,
        snapshot_path=snapshot_path,
    )


def test_plotter_generates_heatmap(analyst_plotter_context, monkeypatch):
    ctx = analyst_plotter_context
    analysis_ui = ctx.analysis_ui
    analysis_ui.tabWidget_analysis.setCurrentWidget(analysis_ui.tabPlot)
    ctx.worker.process_events()

    analysis_ui.spinBox_tunnelStart.setValue(5)
    analysis_ui.spinBox_tunnelEnd.setValue(5)
    analysis_ui.pushButton_resetPlotTunnelRange.click()
    ctx.worker.process_events()
    assert get_widget_value(analysis_ui.spinBox_tunnelStart) == 1
    assert get_widget_value(analysis_ui.spinBox_tunnelEnd) == ctx.plotter._max_section

    from matplotlib.figure import Figure

    shown_sizes: list[tuple[float, float]] = []

    def fake_show(self: Figure) -> None:  # type: ignore[override]
        shown_sizes.append(tuple(self.get_size_inches()))

    monkeypatch.setattr(Figure, "show", fake_show, raising=False)
    analysis_ui.pushButton_tunnelPlot.click()
    ctx.worker.process_events()
    assert shown_sizes, "Plotter never displayed the heatmap figure"


def test_plotter_reset_controls(analyst_plotter_context):
    ctx = analyst_plotter_context
    analysis_ui = ctx.analysis_ui
    analysis_ui.tabWidget_analysis.setCurrentWidget(analysis_ui.tabPlot)
    ctx.worker.process_events()

    analysis_ui.spinBox_imageSizeWidthCm.setValue(3)
    analysis_ui.spinBox_imageSizeHightCm.setValue(3)
    analysis_ui.spinBox_imageSizeWidthPx.setValue(64)
    analysis_ui.spinBox_imageSizeHightPx.setValue(32)

    cmap_combo = analysis_ui.comboBox_plotColormap
    for idx in range(cmap_combo.count()):
        text = cmap_combo.itemText(idx)
        if text != CaverAnalystPlotter._DEFAULT_CMAP:
            cmap_combo.setCurrentIndex(idx)
            break
    assert cmap_combo.currentText() != CaverAnalystPlotter._DEFAULT_CMAP

    analysis_ui.pushButton_resetPlotTunnelImageSize.click()
    analysis_ui.pushButton_resetPlotTunnelCmap.click()
    ctx.worker.process_events()

    assert get_widget_value(analysis_ui.spinBox_imageSizeWidthCm) == CaverAnalystPlotter._DEFAULT_IMAGE_WIDTH_CM
    assert get_widget_value(analysis_ui.spinBox_imageSizeHightCm) == CaverAnalystPlotter._DEFAULT_IMAGE_HEIGHT_CM
    width_px, height_px = ctx.plotter._default_pixel_dimensions()
    assert get_widget_value(analysis_ui.spinBox_imageSizeWidthPx) == width_px
    assert get_widget_value(analysis_ui.spinBox_imageSizeHightPx) == height_px
    assert cmap_combo.currentText() == CaverAnalystPlotter._DEFAULT_CMAP
