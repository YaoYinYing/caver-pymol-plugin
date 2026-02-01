import importlib.util
import logging
import subprocess
import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_caver_java():
    package_name = "Caver4"
    if package_name not in sys.modules:
        package = types.ModuleType(package_name)
        package.__path__ = [str(REPO_ROOT / package_name)]
        sys.modules[package_name] = package

    if "Caver4.caver_pymol" not in sys.modules:
        caver_pymol_stub = types.ModuleType("Caver4.caver_pymol")
        caver_pymol_stub.ROOT_LOGGER = logging.getLogger("Caver4TestRoot")
        sys.modules["Caver4.caver_pymol"] = caver_pymol_stub

    spec = importlib.util.spec_from_file_location(
        "Caver4.caver_java",
        REPO_ROOT / "Caver4" / "caver_java.py",
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules["Caver4.caver_java"] = module
    spec.loader.exec_module(module)
    return module


caver_java = _load_caver_java()


def _install_run_command_stub(monkeypatch, success_heaps, final_returncode=0):
    """
    Patch Caver4.caver_java.run_command while still invoking the real Java binary
    for `java -version`. Heap probing commands are short-circuited to keep tests
    fast deterministically.
    """

    recorded = {"heaps": [], "final_cmd": []}
    real_run_command = caver_java.run_command

    def fake_run_command(cmd, verbose=False, env=None):
        if "-version" in cmd:
            return real_run_command(cmd, verbose=verbose, env=env)

        if "do_nothing" in cmd:
            for token in cmd:
                if token.startswith("-Xmx") and token.endswith("m"):
                    heap_value = int(token[4:-1])
                    recorded["heaps"].append(heap_value)
                    break
            rc = 0 if recorded["heaps"] and recorded["heaps"][-1] in success_heaps else 1
            return subprocess.CompletedProcess(cmd, rc, "", "")

        recorded["final_cmd"].append(cmd)
        return subprocess.CompletedProcess(cmd, final_returncode, "caver run", "")

    monkeypatch.setattr(caver_java, "run_command", fake_run_command)
    return recorded


def _build_pyjava(
    monkeypatch,
    tmp_path,
    customized_heap,
    success_heaps,
    final_returncode=0,
):
    records = _install_run_command_stub(monkeypatch, success_heaps, final_returncode)
    repo_root = Path(__file__).resolve().parents[2]
    caver_folder = repo_root / "Caver4"
    jar_path = tmp_path / "dummy.jar"
    jar_path.write_text("not a real jar")
    pdb_dir = tmp_path / "pdb"
    config_file = tmp_path / "config.txt"
    output_dir = tmp_path / "output"
    pdb_dir.mkdir()
    config_file.write_text("CONFIG")
    output_dir.mkdir()

    instance = caver_java.PyJava(
        customized_heap,
        str(caver_folder),
        str(jar_path),
        str(pdb_dir),
        str(config_file),
        str(output_dir),
    )
    return instance, records, {
        "caver_folder": caver_folder,
        "jar_path": jar_path,
        "pdb_dir": pdb_dir,
        "config_file": config_file,
        "output_dir": output_dir,
    }


@pytest.mark.parametrize(
    ("customized_heap", "success_heaps", "expected_heap", "expected_calls"),
    [
        (
            1150,
            {500, 800, 900, 1000, 1050},
            1050,
            [500, 800, 900, 950, 1000, 1050, 1100, 1150, 1150],
        ),
        (
            300,
            {300},
            300,
            [300],
        ),
    ],
)
def test_optimize_memory_tracks_highest_success(monkeypatch, tmp_path, customized_heap, success_heaps, expected_heap, expected_calls):
    pyjava, records, _ = _build_pyjava(monkeypatch, tmp_path, customized_heap, success_heaps)
    assert pyjava.memory_heap_level == expected_heap
    assert records["heaps"] == expected_calls


def test_run_caver_uses_constructed_command(monkeypatch, tmp_path):
    success_heaps = {500, 800}
    pyjava, records, paths = _build_pyjava(monkeypatch, tmp_path, 800, success_heaps)

    expected_cmd = [
        pyjava.java_bin,
        f"-Xmx{pyjava.memory_heap_level}m",
        "-cp",
        str(paths["caver_folder"] / "lib"),
        "-jar",
        str(paths["jar_path"]),
        "-home",
        str(paths["caver_folder"]),
        "-pdb",
        str(paths["pdb_dir"]),
        "-conf",
        str(paths["config_file"]),
        "-out",
        str(paths["output_dir"]),
    ]
    assert pyjava.cmd == expected_cmd

    result = pyjava.run_caver()
    assert records["final_cmd"] == [expected_cmd]
    assert result.stdout == "caver run"
