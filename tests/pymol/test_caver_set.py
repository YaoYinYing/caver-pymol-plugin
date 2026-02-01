from __future__ import annotations

import pytest
from types import SimpleNamespace

pytest.importorskip("pymol")

from Caver4.caver_config import CaverConfig, CaverShortcut
from Caver4.caver_pymol import CaverPyMOL


class DummyWidget:
    def __init__(self) -> None:
        self.value = None


@pytest.fixture
def caver_set_driver(monkeypatch: pytest.MonkeyPatch):
    calls: list[tuple[object, object]] = []

    def fake_set_widget_value(widget, value):
        widget.value = value
        calls.append((widget, value))

    monkeypatch.setattr("Caver4.caver_pymol.set_widget_value", fake_set_widget_value)

    driver = SimpleNamespace(
        config=CaverConfig(),
        ui=SimpleNamespace(lineEdit_outputDir=DummyWidget()),
        ui_config=SimpleNamespace(spinBox_maxJavaHeapSize=DummyWidget()),
        config_bindings_main_rev=CaverPyMOL.config_bindings_main_rev,
        config_bindings_config_rev=CaverPyMOL.config_bindings_config_rev,
    )
    return SimpleNamespace(driver=driver, calls=calls)


@pytest.mark.parametrize(
    ("key", "value", "owner", "widget"),
    (
        ("output_dir", "/tmp/caver", "ui", "lineEdit_outputDir"),
        ("customized_java_heap", 8192, "ui_config", "spinBox_maxJavaHeapSize"),
    ),
)
def test_caver_set_updates_bound_widgets(caver_set_driver, key, value, owner, widget):
    widget_obj = getattr(getattr(caver_set_driver.driver, owner), widget)

    CaverPyMOL.caver_set(caver_set_driver.driver, key, value)

    assert widget_obj.value == value
    assert caver_set_driver.calls == [(widget_obj, value)]


def test_caver_set_updates_raw_config_when_widget_missing(caver_set_driver):
    caver_set_driver.driver.config.set("starting_point_coordinates", "???")

    CaverPyMOL.caver_set(caver_set_driver.driver, "starting_point_coordinates", "1 2 3")

    assert caver_set_driver.driver.config.get("starting_point_coordinates") == "1 2 3"
    assert caver_set_driver.calls == []


def test_caver_set_warns_and_sets_unknown_key(caver_set_driver):
    with pytest.warns(UserWarning, match="not a valid Caver setting"):
        CaverPyMOL.caver_set(caver_set_driver.driver, "experimental_mode", "enabled")

    assert caver_set_driver.driver.config.get("experimental_mode") == "enabled"
    assert caver_set_driver.calls == []


def test_caver_shortcut_tracks_previous_completion(monkeypatch: pytest.MonkeyPatch):
    shortcut = CaverShortcut(config=CaverConfig(), keywords=[])

    def fake_interpret(self, keyword, mode=False):
        outcomes = {
            "valid": "valid",
            "ambiguous": ["valid"],
        }
        return outcomes.get(keyword)

    monkeypatch.setattr(CaverShortcut, "_interpret", fake_interpret, raising=False)

    assert shortcut.interpret("valid") == "valid"
    assert shortcut.config._complete_temp == "valid"

    assert shortcut.interpret("ambiguous") == ["valid"]
    assert shortcut.config._complete_temp == ""
