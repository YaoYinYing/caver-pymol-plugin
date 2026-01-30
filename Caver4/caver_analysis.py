from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Optional, Union
from pymol import cmd

# pandas is not supposed to be installed with PyMOL

from .caver_pymol import ROOT_LOGGER

logging=ROOT_LOGGER.getChild('Analysis')

def exists(name: str):
    return name in cmd.get_names("all")

def read_tunnel_csv(file: str)-> list[list[float]]:
    columns: list[list[float]]=[]
    with open(file, "r", newline="") as f:
        for line in f:
            stripped=line.strip()
            if not stripped:
                continue
            values=[float(value) for value in stripped.split(",") if value]
            for idx, value in enumerate(values):
                if idx==len(columns):
                    columns.append([])
                columns[idx].append(value)
    return columns

def read_tunnel_pdb(file_path: str)-> list[str]:
    with open(file_path, "r", newline="") as f:
        data=f.read()
    if not data:
        return []
    frames=[]
    for chunk in data.replace("\r\n", "\n").split("\n\n"):
        stripped=chunk.strip("\n")
        if stripped:
            if not stripped.endswith("\n"):
                stripped+= "\n"
            frames.append(stripped)
    return frames
    

@dataclass
class NullFrame:
    """
    Null frame
    """
    frame_id: int


@dataclass
class TunnelFrame:
    
    pdb_strings: str
    frame_id: int

    def load(self, name: str, group: str, apply_to_vdw: bool=True) -> str:
        obj_name=f'{name}.{self.frame_id}'

        if exists(obj_name):
            logging.warning(f'Object {obj_name} already exists. Deleted.')
            cmd.delete(obj_name)
        
        cmd.load_raw(self.pdb_strings, obj_name)
        cmd.group(group, obj_name)
        if apply_to_vdw:
            # backpropagate vdw radius from b-factor
            cmd.alter(obj_name, "vdw=b")
        return obj_name
    
    def render(self, obj_name: str, minimum: float=1.5, maximum: float=3.0, palette: str='red_green') -> None:
        try:
            cmd.spectrum('vdw', palette, obj_name, minimum=minimum, maximum=maximum)
        except Exception as e:
            logging.error(f"Error rendering {obj_name} due to:\n {e}")
        



@dataclass
class TunnelDynamic:
    name: str
    frames: list[Union[TunnelFrame, NullFrame]]
    
    @classmethod
    def from_result_dir(cls, res_dir: str, run_id: int, tunnel_id) -> TunnelDynamic:
        res_dir=os.path.abspath(res_dir)

        # <res_dir>/<run_id>/analysis/profile_heat_maps/csv/cl_000001_heat_map.csv
        csv_file = os.path.join(res_dir, str(run_id), "analysis", "profile_heat_maps", "csv", f"cl_{tunnel_id:06d}_heat_map.csv")

        # <res_dir>/<run_id>/data/clusters_timeless/tun_cl_001_1.pdb
        pdb_file = os.path.join(res_dir, str(run_id), "data", "clusters_timeless", f"tun_cl_{tunnel_id:03d}_1.pdb")

        if not os.path.exists(csv_file):
            raise FileNotFoundError(f"{csv_file} does not exist")
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"{pdb_file} does not exist")
        
        csv_data = read_tunnel_csv(csv_file)
        pdb_data = read_tunnel_pdb(pdb_file)

        frames: list[Union[TunnelFrame, NullFrame]] = []
        pdb_idx = 0
        pdb_len = len(pdb_data)
        append_frame = frames.append

        for frame_id, column in enumerate(csv_data[1:], start=1):
            if column and all(value != -1.0 for value in column):
                if pdb_idx >= pdb_len:
                    logging.warning(
                        "Not enough PDB frames for tunnel %s (frame %s)", tunnel_id, frame_id
                    )
                    append_frame(NullFrame(frame_id))
                    continue
                append_frame(TunnelFrame(pdb_data[pdb_idx], frame_id))
                pdb_idx += 1
            else:
                append_frame(NullFrame(frame_id))

        return cls(name=f"cl_{tunnel_id:06d}", frames=frames)



        
