from __future__ import annotations

from dataclasses import dataclass,field
import os
import threading
from typing import Optional, Union
from PyQt5 import QtCore
from pymol import cmd
from pymol.constants_palette import palette_dict

# pandas is not supposed to be installed with PyMOL

from .caver_pymol import ROOT_LOGGER
from .utils.ui_tape import get_widget_value, notify_box
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
    


@dataclass(frozen=True)
class TunnelFrame:
    
    
    frame_id: int
    pdb_strings: str = ''
    diameters: list[float] = field(default_factory=list)


    def __repr__(self):
        return f'''TunnelFrame #{self.frame_id}
-=-=-=-=
Atoms: {self.node_number}
Bonds: {self.num_bonds}
Tunnel Length: {self.diameter_records_number}
Is empty: {self.is_empty}
-=-=-=-=
Valid Diameters: {[round(x, 2) for x in self.valid_diameters]}
Full Diameters: {[round(x, 2) for x in self.diameters]}
'''

    @property
    def is_empty(self) -> bool:
        return not self.pdb_strings
    
    @property
    def num_bonds(self) -> int:
        return len([x for x in self.pdb_strings.split("\n") if x.startswith("CONECT")])
    
    @property
    def valid_diameters(self) -> list[float]:
        return [x for x in self.diameters if x > 0]
    
    @property
    def node_number(self) -> int:
        return len([x for x in self.pdb_strings.split('\n') if x.startswith('ATOM')])
    
    @property
    def diameter_records_number(self) -> int:
        return len(self.valid_diameters)

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

    def which_frame(self, frame_id: int) -> TunnelFrame:
        """
        Returns the frame with the given frame_id.
        """
        for frame in self.frames:
            if frame.frame_id == frame_id:
                return frame
        raise IndexError(f"Frame {frame_id} not found in tunnel {self.name}")

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
            logging.info(f"Rendering frame {frame.frame_id} ({repr(frame)}) ...")
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
        self.autoplay_interval=float(get_widget_value(form.doubleSpinBox_autoPlayInterval))   
        
        md_state_file=os.path.join(res_dir, str(run_id), "md_state_number.txt")
        with open(md_state_file, "r") as f:
            self.frame_ids=[int(line.strip()) for line in f.readlines()]

        self._min_frame_id=min(self.frame_ids)
        self._max_frame_id=max(self.frame_ids)
        self._current_frame_id=self._min_frame_id

        self.init_slider_range()
        self._setup_autoplay()
        
    def init_slider_range(self):
        self.slider.setRange(self._min_frame_id, self._max_frame_id)
        # only released signal is emitted when the slider is released,
        # so we can use it to trigger the frame switch and skip the middle frames
        self.slider.valueChanged.connect(self._switch_frame)
        self._update_button_status()
        # update the frame by the initial preview
        self._switch_frame()

    def about_this_frame(self):
        the_frame=self.analyst.tunnels.which_frame(self._current_frame_id)
        notify_box(
            message=f'{repr(the_frame)}',
        )
    
    @property
    def tunnel_objects_to_hide(self, ) -> str:
        return ' or '.join(f'{self.tunnel_name}.{frame_id}' for frame_id in self.frame_ids if frame_id != self._current_frame_id)

    @property
    def tunnel_objects_to_show(self) -> str:
        return f'{self.tunnel_name}.{self._current_frame_id}'
    

    def _update_button_status(self):
        self.form.pushButton_firstFrame.setEnabled(self._current_frame_id != self._min_frame_id)
        self.form.pushButton_lastFrame.setEnabled(self._current_frame_id != self._max_frame_id)

        self.form.pushButton_nextFrame.setEnabled(self._current_frame_id != self._max_frame_id)
        self.form.pushButton_previousFrame.setEnabled(self._current_frame_id != self._min_frame_id)

    # the real work is done here
    def _switch_frame(self):
        # force to sync the slider index 
        self._current_frame_id = get_widget_value(self.form.horizontalSlider)
        logging.debug(f"Switching frame to {self._current_frame_id}")
        cmd.frame(self._current_frame_id)
        # cmd.refresh()
        cmd.disable(f'( {self.tunnel_objects_to_hide} )')
        # cmd.refresh()
        cmd.enable(self.tunnel_objects_to_show)

        cmd.refresh()
        if not self._is_autoplay_running():
            self._update_button_status()

    def _update_index_to_slider(self):

        self.slider.setValue(self._current_frame_id)
    
    def forward(self):
        self._current_frame_id+=1
        self._update_index_to_slider()

    def backward(self):
        
        self._current_frame_id-=1
        self._update_index_to_slider()
    def head(self):
        logging.debug(f"Moving to head of tunnel")
        self._current_frame_id=self._min_frame_id
        self._update_index_to_slider()
    def tail(self):
        logging.debug(f"Moving to tail of tunnel")
        self._current_frame_id=self._max_frame_id
        self._update_index_to_slider()

    def _setup_autoplay(self) -> None:
        self._autoplay_thread: Optional[threading.Thread] = None
        self._autoplay_stop_event = threading.Event()
        self._autoplay_nav_buttons = (
            self.form.pushButton_firstFrame,
            self.form.pushButton_lastFrame,
            self.form.pushButton_nextFrame,
            self.form.pushButton_previousFrame,
        )
        self._autoplay_other_buttons = (
            self.form.pushButton_aboutThisFrame,
            self.form.pushButton_refreshTunnelPreview,
            self.form.pushButton_applyTunnelsSpectrumStatic,
            self.form.pushButton_clearTunnelsSpectrumStatic,
        )

        self.form.pushButton_autoPlay.clicked.connect(self.start_auto_play)
        self.form.pushButton_pauseAutoPlay.clicked.connect(self.pause_auto_play)
        self.form.doubleSpinBox_autoPlayInterval.valueChanged.connect(self._autoplay_interval_changed)
        self.form.pushButton_pauseAutoPlay.setEnabled(False)

    def _is_autoplay_running(self) -> bool:
        thread = getattr(self, "_autoplay_thread", None)
        return bool(thread and thread.is_alive())

    def _autoplay_interval_changed(self, value: float) -> None:
        try:
            self.autoplay_interval = float(value)
        except (TypeError, ValueError):
            self.autoplay_interval = float(get_widget_value(self.form.doubleSpinBox_autoPlayInterval))

    def _set_autoplay_running(self, running: bool) -> None:
        if running:
            for button in (*self._autoplay_nav_buttons, *self._autoplay_other_buttons):
                button.setEnabled(False)
            self.form.pushButton_autoPlay.setEnabled(False)
            self.form.pushButton_pauseAutoPlay.setEnabled(True)
        else:
            for button in self._autoplay_other_buttons:
                button.setEnabled(True)
            self._update_button_status()
            self.form.pushButton_autoPlay.setEnabled(True)
            self.form.pushButton_pauseAutoPlay.setEnabled(False)

    def _next_frame_id(self, current: int) -> int:
        return self._min_frame_id if current >= self._max_frame_id else current + 1

    def _auto_play_worker(self) -> None:
        next_frame = self._current_frame_id
        while not self._autoplay_stop_event.is_set():
            next_frame = self._next_frame_id(next_frame)
            QtCore.QMetaObject.invokeMethod(
                self.slider,
                "setValue",
                QtCore.Qt.QueuedConnection,
                QtCore.Q_ARG(int, next_frame),
            )
            interval = max(float(self.autoplay_interval), 0.01)
            if self._autoplay_stop_event.wait(interval):
                break

        self._autoplay_thread = None

    def start_auto_play(self) -> None:
        if self._autoplay_thread and self._autoplay_thread.is_alive():
            return
        self.autoplay_interval = float(get_widget_value(self.form.doubleSpinBox_autoPlayInterval))
        self._autoplay_stop_event.clear()
        self._set_autoplay_running(True)
        self._autoplay_thread = threading.Thread(target=self._auto_play_worker, daemon=True)
        self._autoplay_thread.start()

    def pause_auto_play(self) -> None:
        thread = self._autoplay_thread
        if not thread:
            return
        self._autoplay_stop_event.set()
        thread.join()
        self._autoplay_thread = None
        self._autoplay_stop_event.clear()
        self._set_autoplay_running(False)
