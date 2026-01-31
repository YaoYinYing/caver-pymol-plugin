from __future__ import annotations

from dataclasses import dataclass,field
import os
from typing import Optional, Union
from pymol import cmd
from pymol.constants_palette import palette_dict

# pandas is not supposed to be installed with PyMOL

from .caver_pymol import ROOT_LOGGER
from .utils.ui_tape import get_widget_value
from .ui.Ui_caver_analysis import Ui_CaverAnalysis as CaverAnalysisForm

logging=ROOT_LOGGER.getChild('Analyst')

palette_tuple = tuple(palette_dict.keys())



def list_palettes() -> tuple[str, ...]:
    return palette_tuple



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
    diameters: list[float] = field(default_factory=list)

    @property
    def is_empty(self) -> bool:
        return not self.pdb_strings
    
    @property
    def node_number(self) -> int:
        return len(x for x in self.pdb_strings.split('\n') if x.startswith('ATOM'))
    
    @property
    def diameter_records_number(self) -> int:
        return len(self.diameters)

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
            missing_frame= all(value == -1.0 for value in column)
            # if the frame is missing, use empty string for pdb_datum
            pdb_datum= pdb_data[pdb_idx] if not missing_frame else ''
            atom_count=[x for x in pdb_datum.split("\n") if x.startswith("ATOM")]
            logging.debug(f'frame_id: {frame_id} ({len(atom_count)}), missing_frame: {missing_frame}')

            # if not missing frame, count the pdb index for next iteration
            if not missing_frame:
                pdb_idx+=1
            
            # assemble and append frame
            
            append_frame(TunnelFrame(frame_id, pdb_datum, column))

        return cls(name=f"cl_{tunnel_id:06d}", frames=frames)


class CaverAnalyst:
    # salute to the original caver analyst package
    
    def __init__(self, res_dir: str, run_id: int, tunnel_id: int, palette: str='red_green'):
        self.res_dir = res_dir
        self.run_id = run_id
        self.tunnel_id = tunnel_id
        self.palette = palette

        self.tunnels: TunnelDynamic = TunnelDynamic.from_result_dir(res_dir, run_id, tunnel_id)
    
    def render(self, minimum: float, maximum:float, palette: Optional[str]='red_green', show_as: str='lines') -> None:
        for frame in self.tunnels.frames:
            logging.info(f"Rendering frame {frame.frame_id} ({len(frame.pdb_strings)}) from {self.tunnels.name} ...")
            frame_name=frame.load(self.tunnels.name, group=f'{self.tunnels.name}_{self.run_id}_t{self.tunnel_id:03d}')
            cmd.show(show_as, frame_name)
            frame.render(frame_name, minimum=minimum, maximum=maximum, palette=palette or self.palette)
        logging.info(f"Rendered tunnel {self.tunnels.name}")
        

def run_analysis(form: CaverAnalysisForm, run_id: Union[str, int], res_dir: str) -> CaverAnalyst:
    palette=get_widget_value(form.comboBox_spectrumPalette)
    run_id=int(run_id)
    tunnel_id=int(get_widget_value(form.comboBox_tunnel))
    spectrum_min=get_widget_value(form.doubleSpinBox_spectrumMin)
    spectrum_max=get_widget_value(form.doubleSpinBox_spectrumMax)

    repre=get_widget_value(form.comboBox_representation)

    analyst=CaverAnalyst(res_dir=res_dir, run_id=run_id, tunnel_id=tunnel_id, palette=palette)
    analyst.render(minimum=spectrum_min, maximum=spectrum_max, palette=palette, show_as=repre)

    return analyst

class CaverAnalystPreviewer:
    def __init__(self, form: CaverAnalysisForm,analyst:CaverAnalyst, res_dir: str, run_id: int):
        
        self.form=form
        self.analyst=analyst
        if not analyst:
            raise ValueError("No analyst provided")

        self.res_dir=analyst.res_dir
        self.run_id=analyst.run_id
        self.tunnel_id=analyst.tunnel_id
        self.tunnel_name=analyst.tunnels.name

        self.slider=form.horizontalSlider
        self.autoplay_interval=get_widget_value(form.doubleSpinBox_autoPlayInterval)   
        
        md_state_file=os.path.join(res_dir, str(run_id), "md_state_number.txt")
        with open(md_state_file, "r") as f:
            self.frame_ids=[int(line.strip()) for line in f.readlines()]
        
        self.num_frames=len(self.frame_ids)
        self.init_slider_range()
        self._current_frame_id=min(self.frame_ids)
    def init_slider_range(self):
        self.slider.setRange(min(self.frame_ids), max(self.frame_ids))
        # only released signal is emitted when the slider is released,
        # so we can use it to trigger the frame switch and skip the middle frames
        self.slider.valueChanged.connect(self._switch_frame)
    
    @property
    def tunnel_objects_to_hide(self, ) -> str:
        return ' or '.join(f'{self.tunnel_name}.{frame_id}' for frame_id in self.frame_ids if frame_id != self._current_frame_id)

    @property
    def tunnel_objects_to_show(self) -> str:
        return f'{self.tunnel_name}.{self._current_frame_id}'
    

    def _update_button_status(self):
        self.form.pushButton_firstFrame.setEnabled(self._current_frame_id != min(self.frame_ids))
        self.form.pushButton_lastFrame.setEnabled(self._current_frame_id != max(self.frame_ids))

        self.form.pushButton_nextFrame.setEnabled(self._current_frame_id != max(self.frame_ids))
        self.form.pushButton_previousFrame.setEnabled(self._current_frame_id != min(self.frame_ids))

    # the real work is done here
    def _switch_frame(self):
        logging.debug(f"Switching frame to {self._current_frame_id}")
        cmd.frame(self._current_frame_id)
        cmd.refresh()
        # logging.debug(f"Hiding frames {self.tunnel_objects_to_hide}")
        cmd.disable(f'( {self.tunnel_objects_to_hide} )')
        cmd.refresh()
        # logging.debug(f"Showing frame {self.tunnel_objects_to_show}")
        cmd.enable(self.tunnel_objects_to_show)

        cmd.refresh()
        self._update_button_status()

    def _update_index_to_slider(self):
        # logging.debug(f"Syncing slider to frame {self._current_frame_id}")
        self.slider.setValue(self._current_frame_id)
    
    def forward(self):
        
        if self._current_frame_id >= self.num_frames:
            raise IndexError("Reached the end of the tunnel")
        # logging.debug(f"Moving forward to frame {self._current_frame_id+1}")
        self._current_frame_id+=1
        self._update_index_to_slider()

    def backward(self):
        
        if self._current_frame_id <= 0:
            raise IndexError("Reached the beginning of the tunnel")
        # logging.debug(f"Moving backward to frame {self._current_frame_id-1}")
        self._current_frame_id-=1
        self._update_index_to_slider()
    def head(self):
        logging.debug(f"Moving to head of tunnel")
        self._current_frame_id=min(self.frame_ids)
        self._update_index_to_slider()
    def tail(self):
        logging.debug(f"Moving to tail of tunnel")
        self._current_frame_id=max(self.frame_ids)
        self._update_index_to_slider()

    