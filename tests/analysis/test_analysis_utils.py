from pathlib import Path
import shutil
import sys

import pytest
from pytest import approx

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Caver4.caver_analysis import (
    TunnelDynamic,
    read_tunnel_csv,
    read_tunnel_pdb,
)


DATA_FILE = (
    Path(__file__).resolve().parents[1] / "data" / "csv" / "cl_000001_heat_map.csv"
)
PDB_FILE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "tunnel_pdb"
    / "tun_cl_001_1.pdb"
)


def test_read_tunnel_csv_returns_column_major_structure():
    columns = read_tunnel_csv(str(DATA_FILE))

    assert len(columns) == 50
    assert all(len(column) == 27 for column in columns)


def test_read_tunnel_csv_preserves_column_values():
    columns = read_tunnel_csv(str(DATA_FILE))

    assert columns[0][:5] == approx([-1.0, -1.0, -1.0, -1.0, -1.0])
    assert columns[0][-5:] == approx(
        [
            2.1201782975193733,
            2.0120324202156965,
            1.939933117327831,
            1.890524522545811,
            1.9191565337526026,
        ]
    )
    assert columns[49][-5:] == approx(
        [
            2.0226632238628115,
            1.894042750044157,
            1.843130431605548,
            1.8737900234769984,
            1.983665716291471,
        ]
    )


def test_read_tunnel_pdb_splits_frames_and_preserves_content():
    frames = read_tunnel_pdb(str(PDB_FILE))

    assert len(frames) == 48
    assert all(frame.endswith("\n") for frame in frames)
    assert frames[0].splitlines()[0] == (
        "ATOM      1  H   FIL T   1       8.794  -0.208   0.610        1.90              "
    )
    assert frames[0].splitlines()[-1] == "CONECT   18   19"
    assert frames[-1].splitlines()[0] == (
        "ATOM    741  H   FIL T   1       9.138  -0.032   0.881        1.98              "
    )
    assert frames[-1].splitlines()[-1] == "CONECT  754  755"


def test_read_tunnel_pdb_handles_empty_file(tmp_path):
    empty_file = tmp_path / "empty.pdb"
    empty_file.write_text("")

    assert read_tunnel_pdb(str(empty_file)) == []


def test_tunnel_dynamic_from_result_dir_builds_frames(tmp_path):
    res_dir = tmp_path / "results"
    csv_dir = res_dir / "1" / "analysis" / "profile_heat_maps" / "csv"
    pdb_dir = res_dir / "1" / "data" / "clusters_timeless"
    csv_dir.mkdir(parents=True)
    pdb_dir.mkdir(parents=True)

    csv_target = csv_dir / "cl_000001_heat_map.csv"
    pdb_target = pdb_dir / "tun_cl_001_1.pdb"
    shutil.copy(DATA_FILE, csv_target)
    shutil.copy(PDB_FILE, pdb_target)

    dynamic = TunnelDynamic.from_result_dir(str(res_dir), 1, 1)

    assert dynamic.name == "cl_000001"
    assert len(dynamic.frames) == 49

    tunnel_frames = [frame for frame in dynamic.frames if not frame.is_empty]
    null_frames = [frame for frame in dynamic.frames if frame.is_empty]

    assert len(tunnel_frames) == 1
    assert tunnel_frames[0].frame_id == 36
    assert len(null_frames) == 48
    assert null_frames[0].frame_id == 1


def test_tunnel_dynamic_from_result_dir_requires_csv(tmp_path):
    res_dir = tmp_path / "results"
    pdb_dir = res_dir / "1" / "data" / "clusters_timeless"
    pdb_dir.mkdir(parents=True)
    shutil.copy(PDB_FILE, pdb_dir / "tun_cl_001_1.pdb")

    with pytest.raises(FileNotFoundError) as excinfo:
        TunnelDynamic.from_result_dir(str(res_dir), 1, 1)

    assert "cl_000001_heat_map.csv" in str(excinfo.value)


def test_tunnel_dynamic_from_result_dir_requires_pdb(tmp_path):
    res_dir = tmp_path / "results"
    csv_dir = res_dir / "1" / "analysis" / "profile_heat_maps" / "csv"
    csv_dir.mkdir(parents=True)
    shutil.copy(DATA_FILE, csv_dir / "cl_000001_heat_map.csv")

    with pytest.raises(FileNotFoundError) as excinfo:
        TunnelDynamic.from_result_dir(str(res_dir), 1, 1)

    assert "tun_cl_001_1.pdb" in str(excinfo.value)
