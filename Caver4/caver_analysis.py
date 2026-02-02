from __future__ import annotations

from dataclasses import dataclass,field
import os
from typing import Optional, Union

from pymol import cmd
from pymol.constants_palette import palette_dict

# pandas is not supposed to be installed with PyMOL

from .caver_pymol import ROOT_LOGGER
from .utils.ui_tape import get_widget_value, notify_box, QtCore
from .ui.Ui_caver_analysis import Ui_CaverAnalyst as CaverAnalysisForm

logging=ROOT_LOGGER.getChild('Analyst')

palette_tuple = tuple(palette_dict.keys())

TUNNEL_REPRE=(
    'lines',
    'sticks',
    'spheres',
    'mesh',
    'surface',

)

TUNNEL_SPECTRUM_EXPRE=(
    'b',
    'vdw',
    'resi',
    'index'
)


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
Node: {self.node_number}
Connections: {self.num_bonds}
Length: {self.diameter_records_number}
Is empty: {self.is_empty}
-=-=-=-=
Valid Diameters: {[round(x, 2) for x in self.valid_diameters]}
Full Diameters: {[round(x, 2) for x in self.diameters]}
-=-=-=-=
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
            # delete the object for session safety
            cmd.delete(obj_name)
        
        cmd.load_raw(self.pdb_strings, format='pdb', object=obj_name)
        cmd.group(group, obj_name)
        if apply_to_vdw:
            # backpropagate vdw radius from b-factor
            cmd.alter(obj_name, "vdw=b")
        return obj_name
    def render(self, *args, **kwargs) -> None:
        expression=kwargs.get('expression')
        # spectrum as normal way of PyMOL
        # using vdw, b, index, etc.
        if expression in TUNNEL_SPECTRUM_EXPRE:
            logging.debug(f'rendering spectrum {expression}')
            return self.spectrum(*args, **kwargs)
        
        # apply existing ramp to the object
        # eg: 
        # ```pymol
        # create startpoint, sele
        # cmd.ramp_new('r', 'startpoint', [0, 10], 'rainbow')
        # ```
        # use `r` as ramp name will force to use the proximity from the `startpoint`
        # 
        if expression in cmd.get_names_of_type('object:ramp'):
            logging.debug(f'rendering existing ramp {expression}')
            return cmd.color(expression, kwargs.get('obj_name'))

        # build new ramp based on expression
        # eg:
        # ```pymol
        # create startpoint, sele
        # ```
        # use `startpoint` as ramp starting point will force to use the proximity from the `startpoint`
        logging.debug(f'rendering new ramp {expression}')
        return self.ramp(*args, **kwargs)
        
    
    def ramp(self, obj_name: str, expression: str,minimum: float=1.5, maximum: float=3.0,palette: str='red_green' ) -> None:
        if palette.startswith('rainbow'):
            # fix to rainbow if palette is rainbow_*
            palette='rainbow'
        else:
            # otherwise split it as a list
            palette=palette.split('_')
            
        
        # ramp for each frame and get it colored.
        try:
            new_ramp=f'r_{obj_name}'
            
            cmd.ramp_new(
                name=new_ramp, 
                map_name=expression, 
                color=palette, 
                selection=obj_name,
                range=[minimum, maximum],
                state=self.frame_id
                )
            # hide all ramps
            cmd.group('all_ramps', new_ramp)
            
            cmd.color(new_ramp, obj_name)
        except Exception as e:
            logging.error(e)


    def spectrum(self, obj_name: str, minimum: float=1.5, maximum: float=3.0, palette: str='red_green', expression: str='vdw') -> None:
        try:
            cmd.spectrum(expression, palette, obj_name, minimum=minimum, maximum=maximum)
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
    
    def render(self, minimum: float, maximum:float, palette: Optional[str]='red_green', show_as: str='lines', expression: str='vdw') -> None:
        for frame in self.tunnels.frames:
            logging.info(f"Rendering frame {frame.frame_id} ({repr(frame)}) ...")
            frame_name=frame.load(self.tunnels.name, group=f'{self.tunnels.name}_{self.run_id}_t{self.tunnel_id:03d}')
            cmd.show(show_as, frame_name)
            frame.render(
                obj_name=frame_name, 
                minimum=minimum, 
                maximum=maximum, 
                palette=palette or self.palette, 
                expression=expression
            )
        logging.info(f"Rendered tunnel {self.tunnels.name}")
        

def run_analysis(form: CaverAnalysisForm, run_id: Union[str, int], res_dir: str) -> CaverAnalyst:
    palette=get_widget_value(form.comboBox_spectrumPalette)
    run_id=int(run_id)
    tunnel_id=int(get_widget_value(form.comboBox_tunnel))
    spectrum_min=get_widget_value(form.doubleSpinBox_spectrumMin)
    spectrum_max=get_widget_value(form.doubleSpinBox_spectrumMax)

    spectrum_expression=get_widget_value(form.comboBox_spectrumBy) or 'vdw'

    repre=get_widget_value(form.comboBox_representation)

    analyst=CaverAnalyst(res_dir=res_dir, run_id=run_id, tunnel_id=tunnel_id, palette=palette)
    analyst.render(
        minimum=spectrum_min, 
        maximum=spectrum_max, 
        palette=palette, 
        show_as=repre, 
        expression=spectrum_expression
        )

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
        self._autoplay_timer: Optional[QtCore.QTimer] = None
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
        self._autoplay_timer = QtCore.QTimer(self.slider)
        self._autoplay_timer.setSingleShot(False)
        self._autoplay_timer.timeout.connect(self._auto_play_tick)

    def _is_autoplay_running(self) -> bool:
        timer = getattr(self, "_autoplay_timer", None)
        return bool(timer and timer.isActive())

    def _autoplay_interval_changed(self, value: float) -> None:
        try:
            self.autoplay_interval = float(value)
        except (TypeError, ValueError):
            self.autoplay_interval = float(get_widget_value(self.form.doubleSpinBox_autoPlayInterval))
        if self._is_autoplay_running():
            self._start_autoplay_timer(restart=True)

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

    def _auto_play_tick(self) -> None:
        self.slider.setValue(self._next_frame_id(self._current_frame_id))

    def _start_autoplay_timer(self, restart: bool = False) -> None:
        if not self._autoplay_timer:
            return
        interval_ms = max(int(float(self.autoplay_interval) * 1000), 10)
        if restart and self._autoplay_timer.isActive():
            self._autoplay_timer.stop()
        self._autoplay_timer.start(interval_ms)

    def start_auto_play(self) -> None:
        if self._is_autoplay_running():
            return
        self.autoplay_interval = float(get_widget_value(self.form.doubleSpinBox_autoPlayInterval))
        self._set_autoplay_running(True)
        self._start_autoplay_timer()

    def pause_auto_play(self) -> None:
        timer = self._autoplay_timer
        if timer and timer.isActive():
            timer.stop()
        self._set_autoplay_running(False)

# TODO: class analyst plotter (CaverAnalystPlotter)
# function: plot tunnel data from a given analyst object.
# input: 
#  - analyst object, which contains run id, tunnel id, tunnel data and Form
#  - bool crop empty frames: crop empty frames or not. if true, empty frames will be cropped.
# Widgets:
#  - tunnel start/end (spinBox_tunnelStart, spinBox_tunnelEnd)
#  - colormap (comboBox_plotColormap, default bwr_r)
#  - size of the plot (spinBox_imageSizeWidthCm, spinBox_imageSizeHightCm,spinBox_imageSizeWidthPx , spinBox_imageSizeHightPx,)
#  - dpi (comboBox_DPI, default 150)
#  - aspect ratio lock (checkBox_lockAspectRatio)
#  - save path (lineEdit_imageSavePath) and save path button (pushButton_openSaveImage)
#  - plot button (pushButton_tunnelPlot)
# mechanism:
#  - when clicking the plot button, read all the widget data, and call matplotlib to plot the tunnel data, preview the plot, and save the plot to the save path
# methods:
#  - initialize with:
#   - analyst object, read the min and max tunnel start and end of a tunnel frame (min is basically zero, while max is the length of self.diameters); fill them to `spinBox_tunnelStart` and `spinBox_tunnelEnd`
#   - CaverAnalysisForm: for reading widget data
#  - plot tunnel data
#   - retrieve tunnel start/end, colormap, etc
#   call matplotlib to plot tunnel data for time-series
#   - 
#  - save plot in preview window
#   - open file dialog to select save path and save the plot if given
