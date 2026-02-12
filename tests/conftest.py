from __future__ import annotations

import warnings
from collections.abc import Iterable
from contextlib import contextmanager
from pathlib import Path
from typing import cast

import pytest

from tests.data.test_data import DYNAMIC_WORKFLOW, STATIC_WORKFLOW, WorkflowScenario

TESTS_DIR = Path(__file__).resolve().parent
RESULTS_DIR = TESTS_DIR / "results"
CACHE_DIR = RESULTS_DIR / "caver_cache"
SNAPSHOT_DIR = RESULTS_DIR / "snapshots"
PYMOL_IMAGE_DIR = RESULTS_DIR / "pymol_images"


def _ensure_directories() -> None:
    for directory in (RESULTS_DIR, CACHE_DIR, SNAPSHOT_DIR, PYMOL_IMAGE_DIR):
        directory.mkdir(parents=True, exist_ok=True)


class _NotifyRecorder:
    """
    Captures calls to notify_box without showing real modal dialogs.
    """

    def __init__(self) -> None:
        self.messages: list[tuple[str, type[Exception] | None, str | None]] = []

    def __call__(
        self, message: str = "", error_type: type[Exception] | None = None, details: str | None = None
    ) -> None:
        self.messages.append((message, error_type, details))
        if error_type is None:
            return
        if issubclass(error_type, Warning):
            warnings.warn(error_type(message))
            return
        raise error_type(message)


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """
    Root directory for bundled test data files.
    """
    return TESTS_DIR / "data"


@pytest.fixture(scope="session")
def results_root() -> Path:
    """
    Ensure each run shares a deterministic tree for cached analysis data and snapshots.
    """
    _ensure_directories()
    return RESULTS_DIR


@pytest.fixture
def notify_box_spy(monkeypatch: pytest.MonkeyPatch) -> _NotifyRecorder:
    """
    Replace notify_box with a recorder so tests do not open modal dialogs.
    """

    pytest.importorskip("pymol")
    from Caver4 import caver_analysis, caver_pymol
    from Caver4.utils import ui_tape

    recorder = _NotifyRecorder()
    monkeypatch.setattr(ui_tape, "notify_box", recorder)
    monkeypatch.setattr(caver_pymol, "notify_box", recorder)
    monkeypatch.setattr(caver_analysis, "notify_box", recorder)
    return recorder


class CaverPluginWorker:
    """
    Helper that drives the PyMOL GUI plugin the same way a user would.
    """

    STATIC_SCENARIO: WorkflowScenario = STATIC_WORKFLOW
    DYNAMIC_SCENARIO: WorkflowScenario = DYNAMIC_WORKFLOW
    STATIC_COORDS = cast(tuple[float, float, float], STATIC_SCENARIO.start_value)
    STATIC_OBJECT = STATIC_SCENARIO.object_name
    DYNAMIC_OBJECT = DYNAMIC_SCENARIO.object_name
    DYNAMIC_ATOMS = cast(str, DYNAMIC_SCENARIO.start_value)
    MD_STATE_RANGE = DYNAMIC_SCENARIO.md_state_range or (1, 50)

    def __init__(self, qtbot, data_dir: Path, results_root: Path) -> None:
        pytest.importorskip("pymol")
        import pymol

        from Caver4.caver_pymol import CaverPyMOL, QtWidgets
        from Caver4.utils.ui_tape import get_widget_value, set_widget_value

        self._qtbot = qtbot
        self._data_dir = data_dir
        self._cache_dir = results_root / "caver_cache"
        self.snapshot_dir = results_root / "snapshots"
        self.pymol_image_dir = results_root / "pymol_images"
        for directory in (self._cache_dir, self.snapshot_dir, self.pymol_image_dir):
            directory.mkdir(parents=True, exist_ok=True)

        self._qt = QtWidgets
        self._set_widget_value = set_widget_value
        self._get_widget_value = get_widget_value
        self._app = self._qt.QApplication.instance() or self._qt.QApplication([])

        self._pymol = pymol
        self.cmd = pymol.cmd
        self.cmd.reinitialize()
        self.plugin = CaverPyMOL()
        self.plugin.run_plugin_gui()
        self._qtbot.addWidget(self.plugin.dialog)
        self._configure_output_dir()
        self.process_events()

    def process_events(self) -> None:
        """
        Flush pending Qt events to keep widgets responsive during the tests.
        """
        self._app.processEvents()

    def _configure_output_dir(self) -> None:
        self._set_widget_value(self.plugin.ui.lineEdit_outputDir, str(self._cache_dir))
        self.plugin.config.output_dir = str(self._cache_dir)

    @staticmethod
    @contextmanager
    def _block_signals(widgets: Iterable) -> Iterable[bool]:
        states = []
        for widget in widgets:
            states.append(widget.blockSignals(True))
        try:
            yield states
        finally:
            for widget, state in zip(widgets, states):
                widget.blockSignals(state)

    def _clear_custom_startpoint(self) -> None:
        widgets = [
            self.plugin.ui.textEdit_startpoint,
            self.plugin.ui.radioButton_startAsAtoms,
            self.plugin.ui.radioButton_startAsResidues,
            self.plugin.ui.radioButton_startAsCoords,
        ]
        with self._block_signals(widgets):
            self.plugin.ui.textEdit_startpoint.clear()
            self.plugin.ui.radioButton_startAsAtoms.setChecked(False)
            self.plugin.ui.radioButton_startAsResidues.setChecked(False)
            self.plugin.ui.radioButton_startAsCoords.setChecked(False)
        for key in ("starting_point_atom", "starting_point_residue", "starting_point_coordinates"):
            if self.plugin.config.has(key):
                self.plugin.config.set(key, "???")

    def reset(self) -> None:
        """
        Reset PyMOL, clear input widgets, and bring the plugin back to a neutral state.
        """
        self.cmd.reinitialize()
        self.plugin.register_pymol_commands()
        self.plugin._clear_pymol_sel_and_coords()
        self._clear_custom_startpoint()
        self.plugin.ui.checkBox_MD.setChecked(False)
        self.plugin.update_model_list()
        self.process_events()

    def load_structure(self, file_path: Path, object_name: str) -> str:
        """
        Load a structure file and refresh the model combo box.
        """
        self.cmd.load(str(file_path), object_name)
        self.plugin.update_model_list()
        self._set_widget_value(self.plugin.ui.comboBox_inputModel, object_name)
        self.process_events()

        self.plugin.ui.pushButton_allAA.click()
        self.process_events()
        return object_name

    @staticmethod
    def _format_coordinate_string(coords: tuple[float, float, float]) -> str:
        return " ".join(f"{value:.3f}" for value in coords)

    def set_start_coordinates(self, coords: tuple[float, float, float]) -> None:
        self._clear_custom_startpoint()
        for axis, value in zip("xyz", coords):
            spin = getattr(self.plugin.ui, f"doubleSpinBox_{axis}")
            spin.setValue(float(value))
            setattr(self.plugin.config, f"start_point_{axis}", float(value))
        self.plugin.config.set("starting_point_coordinates", self._format_coordinate_string(coords))
        self.process_events()

    def set_custom_startpoint(self, mode: str, start_value: str) -> None:
        radio_map = {
            "atoms": self.plugin.ui.radioButton_startAsAtoms,
            "residues": self.plugin.ui.radioButton_startAsResidues,
            "coords": self.plugin.ui.radioButton_startAsCoords,
        }
        widgets = [self.plugin.ui.textEdit_startpoint, *radio_map.values()]
        with self._block_signals(widgets):
            self.plugin.ui.textEdit_startpoint.setPlainText(start_value)
            for key, radio in radio_map.items():
                radio.setChecked(key == mode)
        self.plugin._use_custom_startpoint()
        self.process_events()

    def _configure_md_states(self, md_range: tuple[int, int] | None) -> None:
        if not md_range:
            self.plugin.ui.checkBox_MD.setChecked(False)
            return
        min_state, max_state = md_range
        self.plugin.ui.checkBox_MD.setChecked(True)
        self.plugin.ui.spinBox_MD_StateMin.setValue(min_state)
        self.plugin.ui.spinBox_MD_StateMax.setValue(max_state)

    def _run_workflow(self, scenario: WorkflowScenario) -> Path:
        self.reset()
        structure = scenario.structure_path(self._data_dir)
        self.load_structure(structure, scenario.object_name)
        self._configure_md_states(scenario.md_state_range)
        if scenario.start_mode == "coords":
            coords = cast(tuple[float, float, float], scenario.start_value)
            self.set_start_coordinates(coords)
        else:
            start_value = cast(str, scenario.start_value)
            self.set_custom_startpoint(scenario.start_mode, start_value)
        snapshot_prefix = scenario.snapshot_prefix
        self.capture_gui_snapshot(f"{snapshot_prefix}_ui.png")
        self.process_events()
        self.capture_pymol_scene(f"{snapshot_prefix}_scene.png")
        self.process_events()
        self.plugin.execute()
        self.process_events()
        return Path(self.plugin.out_dir)

    def capture_gui_snapshot(self, filename: str) -> Path:
        self.process_events()
        target = self.snapshot_dir / filename
        pixmap = self.plugin.dialog.grab()
        pixmap.save(str(target))
        return target

    def capture_pymol_scene(self, filename: str) -> Path:
        target = self.pymol_image_dir / filename
        # PyMOL segfaults on macOS when OpenGL screenshots are requested without a live viewer.
        # Force the CPU ray-tracer so we can capture scenes headlessly during tests.
        self.cmd.png(str(target), ray=1)
        return target

    def run_static_analysis(self) -> Path:
        """
        Execute the static analysis workflow end-to-end.
        """
        return self._run_workflow(self.STATIC_SCENARIO)

    def run_dynamic_analysis(self) -> Path:
        """
        Execute the MD workflow end-to-end.
        """
        return self._run_workflow(self.DYNAMIC_SCENARIO)

    def shutdown(self) -> None:
        """
        Clean up between tests.
        """
        self.cmd.reinitialize()
        self.process_events()


@pytest.fixture
def caver_worker(qtbot, test_data_dir: Path, results_root: Path, notify_box_spy: _NotifyRecorder) -> CaverPluginWorker:
    """
    Provide a fully-initialized plugin driver for GUI tests.
    """
    worker = CaverPluginWorker(qtbot, test_data_dir, results_root)
    try:
        yield worker
    finally:
        worker.shutdown()
