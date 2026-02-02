from __future__ import annotations

from dataclasses import dataclass,field
import os
from typing import Optional, Union

from pymol import cmd
from pymol.constants_palette import palette_dict

# pandas is not supposed to be installed with PyMOL

from .caver_pymol import ROOT_LOGGER
from .utils.ui_tape import get_widget_value, notify_box, set_widget_value, QtCore, QtWidgets
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
    def reversed_diameters(self) -> list[float]:
        return list(reversed(self.diameters))
    
    @property
    def reversed_valid_diameters(self) -> list[float]:
        return list(reversed(self.valid_diameters))

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

class CaverAnalystPlotter:
    """
    Plot time-series tunnel diameter heat maps from an Analyst instance.

    The plotter wires UI widgets to Matplotlib so users can preview and export plots
    without leaving the plugin window.
    """

    _DEFAULT_CMAP = "bwr_r"
    _DEFAULT_DPI = 150
    _DPI_CHOICES = ("72", "96", "150", "200", "300", "600")
    _FILE_FILTER = (
        "PNG Image (*.png);;"
        "JPEG Image (*.jpg *.jpeg);;"
        "TIFF Image (*.tif *.tiff);;"
        "PDF Document (*.pdf);;"
        "SVG Image (*.svg);;"
        "All Files (*)"
    )

    def __init__(self, form: CaverAnalysisForm, analyst: CaverAnalyst, crop_empty_frames: bool = False):
        if analyst is None:
            raise ValueError("Analyst instance is required for plotting.")

        self.form = form
        self.analyst = analyst
        self.crop_empty_frames = crop_empty_frames

        self._save_path_widget: QtWidgets.QLineEdit = self._locate_save_path_widget()
        self._aspect_checkbox: Optional[QtWidgets.QCheckBox] = getattr(
            form, "checkBox_lockAspectRatio", getattr(form, "checkBox", None)
        )

        self._frames = self._prepare_frames()
        self._max_section = self._compute_max_section()
        if not self._max_section:
            raise ValueError("No tunnel diameter data available for plotting.")

        self._init_widgets()

        self.form.pushButton_tunnelPlot.clicked.connect(self.plot)  # type: ignore[attr-defined]
        self.form.pushButton_openSaveImage.clicked.connect(self._select_save_path)  # type: ignore[attr-defined]

    def _locate_save_path_widget(self) -> QtWidgets.QLineEdit:
        widget = getattr(self.form, "lineEdit_imageSavePath", None)
        if widget is None:
            widget = getattr(self.form, "lineEdit", None)
        if widget is None:
            raise AttributeError("Save path line edit is missing from the form.")
        return widget

    def _prepare_frames(self) -> list[TunnelFrame]:
        frames = list(self.analyst.tunnels.frames)
        if self.crop_empty_frames:
            filtered = [frame for frame in frames if not frame.is_empty]
            if filtered:
                frames = filtered
        return frames

    def _compute_max_section(self) -> int:
        return max((len(frame.diameters) for frame in self.analyst.tunnels.frames if frame.diameters), default=0)

    def _init_widgets(self) -> None:
        self._init_range_inputs()
        self._init_colormap_combo()
        self._init_dpi_combo()
        self._ensure_default_save_path()
        self._ensure_default_sizes()

    def _init_range_inputs(self) -> None:
        start_spin: QtWidgets.QSpinBox = self.form.spinBox_tunnelStart  # type: ignore[attr-defined]
        end_spin: QtWidgets.QSpinBox = self.form.spinBox_tunnelEnd  # type: ignore[attr-defined]
        start_spin.setRange(1, self._max_section)
        end_spin.setRange(1, self._max_section)
        if start_spin.value() <= 0:
            start_spin.setValue(1)
        if end_spin.value() <= 0 or end_spin.value() < start_spin.value():
            end_spin.setValue(self._max_section)

    # no need to init the colormap combo, it has been done in the parent class
    def _init_colormap_combo(self) -> None:
        combo: QtWidgets.QComboBox = self.form.comboBox_plotColormap  # type: ignore[attr-defined]
        if combo.count() == 0:
            try:
                import matplotlib

                cmap_names = sorted(matplotlib.colormaps(), key=str.lower)
            except Exception:
                cmap_names = ["viridis", "plasma", "inferno", "magma", "cividis", "coolwarm", self._DEFAULT_CMAP]
            combo.addItems(cmap_names)
        idx = combo.findText(self._DEFAULT_CMAP)
        if idx >= 0:
            combo.setCurrentIndex(idx)

    def _init_dpi_combo(self) -> None:
        combo: QtWidgets.QComboBox = self.form.comboBox_DPI  # type: ignore[attr-defined]
        current_items = [combo.itemText(i) for i in range(combo.count())]
        if current_items != list(self._DPI_CHOICES):
            set_widget_value(combo, self._DPI_CHOICES)
        set_widget_value(combo, str(self._DEFAULT_DPI))

    def _ensure_default_save_path(self) -> None:
        if not get_widget_value(self._save_path_widget).strip():
            self._save_path_widget.setText(self._default_filename())

    def _ensure_default_sizes(self) -> None:
        width_cm: QtWidgets.QSpinBox = self.form.spinBox_imageSizeWidthCm  # type: ignore[attr-defined]
        height_cm: QtWidgets.QSpinBox = self.form.spinBox_imageSizeHightCm  # type: ignore[attr-defined]
        width_px: QtWidgets.QSpinBox = self.form.spinBox_imageSizeWidthPx  # type: ignore[attr-defined]
        height_px: QtWidgets.QSpinBox = self.form.spinBox_imageSizeHightPx  # type: ignore[attr-defined]

        if width_cm.value() == 0:
            width_cm.setValue(15)
        if height_cm.value() == 0:
            height_cm.setValue(10)
        estimated_pixels = max(self._max_section * 4, 200)
        if width_px.value() == 0:
            width_px.setValue(estimated_pixels)
        if height_px.value() == 0:
            height_px.setValue(estimated_pixels // 2)

    def _default_filename(self) -> str:
        base = f"run_{self.analyst.run_id}_tunnel_{self.analyst.tunnel_id}.png"
        return os.path.join(self.analyst.res_dir, base)

    def _select_save_path(self) -> None:
        default = self._save_path_widget.text() or self._default_filename()
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self.form.tabPlot, "Save tunnel plot", default, self._FILE_FILTER  # type: ignore[attr-defined]
        )
        if filename:
            self._save_path_widget.setText(filename)

    def _read_range(self) -> tuple[int, int]:
        start = max(1, min(self._max_section, int(get_widget_value(self.form.spinBox_tunnelStart))))  # type: ignore[attr-defined]
        end = max(1, min(self._max_section, int(get_widget_value(self.form.spinBox_tunnelEnd))))  # type: ignore[attr-defined]
        if start > end:
            start, end = end, start
        return start, end

    def _gather_heatmap_data(self, start: int, end: int) -> tuple[list[list[float]], list[int]]:
        width = end - start + 1
        matrix: list[list[float]] = []
        ids: list[int] = []
        lower = start - 1
        upper = min(end, self._max_section)
        for frame in self._frames:
            slice_data = frame.diameters[lower:upper]
            sanitized = [value if value >= 0 else float("nan") for value in slice_data]
            if not sanitized:
                sanitized = [float("nan")] * width
            if len(sanitized) < width:
                sanitized.extend(float("nan") for _ in range(width - len(sanitized)))
            matrix.append(sanitized)
            ids.append(frame.frame_id)
        return matrix, ids

    def _figure_size_inches(self, dpi: int) -> tuple[float, float]:
        width_cm: int = get_widget_value(self.form.spinBox_imageSizeWidthCm)  # type: ignore[attr-defined]
        height_cm: int = get_widget_value(self.form.spinBox_imageSizeHightCm)  # type: ignore[attr-defined]
        if width_cm > 0 and height_cm > 0:
            return width_cm / 2.54, height_cm / 2.54
        width_px: int = get_widget_value(self.form.spinBox_imageSizeWidthPx)  # type: ignore[attr-defined]
        height_px: int = get_widget_value(self.form.spinBox_imageSizeHightPx)  # type: ignore[attr-defined]
        if width_px > 0 and height_px > 0 and dpi > 0:
            return width_px / dpi, height_px / dpi
        return 8.0, 5.0

    def _get_selected_dpi(self) -> int:
        text = str(get_widget_value(self.form.comboBox_DPI))  # type: ignore[attr-defined]
        try:
            return max(1, int(float(text)))
        except ValueError:
            return self._DEFAULT_DPI

    def _get_colormap(self):
        cmap_name = get_widget_value(self.form.comboBox_plotColormap).strip()  # type: ignore[attr-defined]
        if not cmap_name:
            cmap_name = self._DEFAULT_CMAP
        import matplotlib.pyplot as plt

        try:
            return plt.get_cmap(cmap_name)
        except ValueError:
            notify_box(f"Colormap '{cmap_name}' not found, falling back to {self._DEFAULT_CMAP}.", Warning)
            return plt.get_cmap(self._DEFAULT_CMAP)

    def _color_limits(self) -> tuple[Optional[float], Optional[float]]:
        try:
            vmin = float(get_widget_value(self.form.doubleSpinBox_spectrumMin))  # type: ignore[attr-defined]
            vmax = float(get_widget_value(self.form.doubleSpinBox_spectrumMax))  # type: ignore[attr-defined]
        except (AttributeError, TypeError, ValueError):
            return None, None
        if vmin > vmax:
            vmin, vmax = vmax, vmin
        if vmin == vmax:
            vmax = vmin + 1e-9
        return vmin, vmax

    # TODO: the start and end of tunnels are currently reversed, use `reversed_diameters` property to reverse it
    def plot(self) -> None:
        if not self._frames:
            notify_box("No tunnel frames available for plotting.", RuntimeError)
            return

        start, end = self._read_range()
        if start == end:
            notify_box("Tunnel start and end must be different positions.", Warning)
            return

        data, frame_ids = self._gather_heatmap_data(start, end)
        if not data:
            notify_box("No tunnel data exists in the selected range.", Warning)
            return
        plot_data = [list(column) for column in zip(*data)]

        dpi = self._get_selected_dpi()
        figsize = self._figure_size_inches(dpi)
        cmap = self._get_colormap()
        vmin, vmax = self._color_limits()

        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:
            notify_box(
                "Matplotlib is required to plot tunnel data. Install matplotlib if you need plotting.",
                Warning,
                details=str(exc),
            )
            return

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        aspect = "equal" if self._aspect_checkbox and self._aspect_checkbox.isChecked() else "auto"
        first_id = frame_ids[0]
        last_id = frame_ids[-1]
        x_min = first_id - 0.5
        x_max = last_id + 0.5 if last_id != first_id else first_id + 0.5
        y_min = start - 0.5
        y_max = end + 0.5
        im = ax.imshow(
            plot_data,
            aspect=aspect,
            cmap=cmap,
            origin="lower",
            interpolation="nearest",
            extent=[x_min, x_max, y_min, y_max],
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_xlabel("Frame ID")
        ax.set_ylabel("Tunnel position (index)")
        ax.set_title(f"Tunnel {self.analyst.tunnel_id} · Run {self.analyst.run_id}")
        from matplotlib.ticker import FuncFormatter, MaxNLocator

        def _format_tick(value: float, _pos: int) -> str:
            try:
                return str(int(round(value)))
            except (TypeError, ValueError):
                return ""

        max_x_ticks = 12
        max_y_ticks = 15
        ax.xaxis.set_major_locator(MaxNLocator(nbins=max_x_ticks, integer=True, min_n_ticks=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=max_y_ticks, integer=True, min_n_ticks=4))
        formatter = FuncFormatter(_format_tick)
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
        rotation = 45 if len(frame_ids) > max_x_ticks else 0
        ax.tick_params(axis="x", rotation=rotation)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Diameter (Å)")
        fig.tight_layout()
        manager = getattr(fig.canvas, "manager", None)
        if manager is not None:
            try:
                manager.set_window_title(f"Caver Tunnel Plot - Run {self.analyst.run_id}")
            except Exception:
                pass

        save_path = self._save_path_widget.text().strip()
        if save_path:
            save_dir = os.path.dirname(save_path)
            if save_dir and not os.path.exists(save_dir):
                os.makedirs(save_dir, exist_ok=True)
            fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        fig.show()
