from __future__ import annotations

import builtins
import importlib
import logging
import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


class DummyQIcon:
    def __init__(self, payload):
        self.payload = payload


class _DynamicQtNamespace:
    def __init__(self, **defaults):
        self._values = dict(defaults)

    def __getattr__(self, name):
        if name not in self._values:
            self._values[name] = type(name, (), {})
        return self._values[name]


@pytest.fixture
def ui_tape_module(monkeypatch):
    """Load Caver4.utils.ui_tape with stubbed PyMOL/Qt dependencies."""

    qt_module = types.ModuleType("pymol.Qt")
    qt_module.QtCore = _DynamicQtNamespace(Qt=type("Qt", (), {}), pyqtSignal=lambda *_, **__: object())
    qt_module.QtWidgets = _DynamicQtNamespace()
    qt_module.QtGui = _DynamicQtNamespace(QIcon=DummyQIcon)

    pymol_module = types.ModuleType("pymol")
    pymol_module.Qt = qt_module
    monkeypatch.setitem(sys.modules, "pymol", pymol_module)
    monkeypatch.setitem(sys.modules, "pymol.Qt", qt_module)

    caver_pkg = types.ModuleType("Caver4")
    caver_pkg.__path__ = [str(REPO_ROOT / "Caver4")]
    monkeypatch.setitem(sys.modules, "Caver4", caver_pkg)

    importlib.import_module("Caver4.utils")

    caver_pymol_stub = types.ModuleType("Caver4.caver_pymol")
    caver_pymol_stub.ROOT_LOGGER = logging.getLogger("Caver4TestRoot")
    monkeypatch.setitem(sys.modules, "Caver4.caver_pymol", caver_pymol_stub)

    module = sys.modules.get("Caver4.utils.ui_tape")
    if module is None:
        module = importlib.import_module("Caver4.utils.ui_tape")
    else:
        module = importlib.reload(module)
    return module


def test_list_color_map_returns_icons_when_matplotlib_available(ui_tape_module, monkeypatch):
    fake_matplotlib = types.ModuleType("matplotlib")
    fake_matplotlib.colormaps = lambda: ["Alpha", "Beta"]
    monkeypatch.setitem(sys.modules, "matplotlib", fake_matplotlib)

    created = []

    def fake_create_cmap_icon(*, cmap):
        created.append(cmap)
        return f"pixmap-{cmap}"

    monkeypatch.setattr(ui_tape_module, "create_cmap_icon", fake_create_cmap_icon)

    result = ui_tape_module.list_color_map()

    assert list(result.keys()) == ["Alpha", "Beta"]
    assert created == ["Alpha", "Beta"]
    for cmap_name, icon in result.items():
        assert isinstance(icon, DummyQIcon)
        assert icon.payload == f"pixmap-{cmap_name}"


def test_list_color_map_prompts_install_without_matplotlib(ui_tape_module, monkeypatch):
    monkeypatch.delitem(sys.modules, "matplotlib", raising=False)
    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "matplotlib":
            raise ImportError("missing matplotlib")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    assert ui_tape_module.list_color_map() == "Please install matplotlib before using this feature".split()
