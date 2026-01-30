from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Optional, Union
from pymol import cmd

# pandas is not supposed to be installed with PyMOL

from .caver_pymol import ROOT_LOGGER
from .utils.ui_tape import get_widget_value, set_widget_value
from .ui.Ui_caver_analysis import Ui_CaverAnalysis as CaverAnalysisForm

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
class TunnelFrame:
    
    
    frame_id: int
    pdb_strings: str = ''

    @property
    def is_empty(self) -> bool:
        return not self.pdb_strings

    def load(self, name: str, group: str, apply_to_vdw: bool=True) -> str:
        obj_name=f'{name}.{self.frame_id}'

        if exists(obj_name):
            logging.warning(f'Object {obj_name} already exists. Deleted.')
            cmd.delete(obj_name)
        
        cmd.load_raw(self.pdb_strings, format='pdb', object=obj_name)
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
    frames: list[TunnelFrame]
    
    @classmethod
    def from_result_dir(cls, res_dir: str, run_id: int, tunnel_id: int) -> TunnelDynamic:
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

        frames: list[TunnelFrame] = []
        pdb_idx = 0
        append_frame = frames.append

        for frame_id, column in enumerate(csv_data[1:], start=1):
            # check if frame is missing
            missing_frame= not all(value != -1.0 for value in column)
            # if the frame is missing, use empty string for pdb_datum
            pdb_datum= pdb_data[pdb_idx] if not missing_frame else ''

            # if not missing frame, count the pdb index for next iteration
            if not missing_frame:
                pdb_idx+=1
            
            # assemble and append frame
            append_frame(TunnelFrame(frame_id, pdb_datum))

        return cls(name=f"cl_{tunnel_id:06d}", frames=frames)


class CaverAnalyst:
    # salute to the original caver analyst package
    
    def __init__(self, res_dir: str, run_id: Union[int, str], tunnel_id: int, pallete: str='red_green'):
        self.res_dir = res_dir
        self.run_id = run_id
        self.tunnel_id = tunnel_id
        self.palette = pallete
        self.tunnels: TunnelDynamic = TunnelDynamic.from_result_dir(res_dir, run_id, tunnel_id)
    
    def render(self, minimum: float, maximum:float, palette: Optional[str]='red_green') -> None:
        for frame in self.tunnels.frames:
            logging.info(f"Rendering frame {frame.name}")
            frame_name=frame.load(self.tunnels.name, group=f'{self.tunnels.name}_{self.run_id}_t{self.tunnel_id:03d}')
            frame.render(frame_name, minimum=minimum, maximum=maximum, palette=palette or self.palette)
        logging.info(f"Rendered tunnel {self.tunnels.name}")
            



def run_analysis(form: CaverAnalysisForm, run_id: Union[str, int], res_dir: str):
    palette=get_widget_value(form.comboBox_spectrumPalette)
    run_id=str(run_id)
    tunnel_id=get_widget_value(form.comboBox_tunnel)
    spectrum_min=get_widget_value(form.doubleSpinBox_min)
    spectrum_max=get_widget_value(form.doubleSpinBox_max)

    analyst=CaverAnalyst(res_dir=res_dir, run_id=run_id, tunnel_id=tunnel_id, pallete=palette)
    analyst.render(minimum=spectrum_min, maximum=spectrum_max)


    