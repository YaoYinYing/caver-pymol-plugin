import shutil

import pytest

pytest.importorskip("pymol")

if shutil.which("java") is None:
    pytest.skip("Java runtime is required for the Caver4 PyMOL integration tests.", allow_module_level=True)


def test_static_analysis_workflow(caver_worker) -> None:
    """
    Run the plugin against a single PDB and verify that artifacts are produced.
    """
    run_dir = caver_worker.run_static_analysis()

    assert run_dir.exists(), f"Missing output for static analysis: {run_dir}"
    assert (run_dir / "input").is_dir(), "Input directory not created for static run"
    assert any(run_dir.joinpath("input").glob("config_*.txt")), "Config file was not exported"
    assert (run_dir / "pymol" / "view_plugin.py").is_file(), "PyMOL visualization script missing"

    start_coords = caver_worker.plugin.config.get("starting_point_coordinates")
    assert isinstance(start_coords, str) and "17.012" in start_coords

    static_ui = caver_worker.snapshot_dir / "static_analysis_ui.png"
    static_scene = caver_worker.pymol_image_dir / "static_analysis_scene.png"
    assert static_ui.is_file(), "UI snapshot for static run was not captured"
    assert static_scene.is_file(), "PyMOL scene for static run was not captured"


def test_dynamic_analysis_workflow(caver_worker) -> None:
    """
    Execute the dynamic MD workflow and ensure expected caches are generated.
    """
    run_dir = caver_worker.run_dynamic_analysis()

    assert run_dir.exists(), "Dynamic analysis did not create an output directory"
    md_state_file = run_dir / "md_state_number.txt"
    assert md_state_file.is_file(), "MD state tracking file missing"
    contents = md_state_file.read_text().strip().splitlines()
    min_state, max_state = caver_worker.MD_STATE_RANGE
    assert (min_state, max_state) == (1, 5), "Dynamic tests should cover MD states 1-5"
    assert contents[0] == str(min_state) and contents[-1] == str(max_state)

    input_pdbs = list(run_dir.joinpath("input").glob("*.pdb"))
    assert input_pdbs, "MD frames were not exported as PDBs"

    assert caver_worker.plugin.ui.checkBox_MD.isChecked()
    assert caver_worker.plugin.config.get("starting_point_atom") == caver_worker.DYNAMIC_ATOMS

    dynamic_ui = caver_worker.snapshot_dir / "dynamic_analysis_ui.png"
    dynamic_scene = caver_worker.pymol_image_dir / "dynamic_analysis_scene.png"
    assert dynamic_ui.is_file(), "UI snapshot for MD run missing"
    assert dynamic_scene.is_file(), "PyMOL scene for MD run missing"
