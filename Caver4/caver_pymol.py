# CAVER Copyright Notice
# ============================
#

'''
"THE BEERWARE LICENSE" (Revision 42):

Yinying wrote the refactored code.
As long as you retain this notice, you can do whatever you want with this stuff.
If we meet someday, and you think this stuff is worth it, you can buy me a beer
in return.
-- Yinying Yao

'''

import logging
import math
import os
import re
from contextlib import contextmanager
from functools import partial
from typing import TYPE_CHECKING, Any, Dict, List, Optional
import warnings
import webbrowser

from pymol import cmd

if TYPE_CHECKING:
    from PyQt5 import QtWidgets
else:
    from pymol.Qt import QtWidgets

import time

# TODO: deprecated
from pymol import stored
from pymol.shortcut import Shortcut
from pymol.cgo import BEGIN, END, LINE_STRIP, LINEWIDTH, VERTEX
from pymol.Qt.utils import getSaveFileNameWithExt

from .ui.Ui_caver import Ui_CaverUI as CaverUI
from .ui.Ui_caver_config import Ui_CaverConfigForm as CaverConfigForm
from .utils.ui_tape import (CheckableListView, get_widget_value,
                            getExistingDirectory, getOpenFileNameWithExt,
                            hold_trigger_button, notify_box,
                            run_worker_thread_with_progress, set_widget_value,
                            widget_signal_tape)


from .caver_config import CaverConfig, THIS_DIR,CONFIG_TXT, CaverShortcut
from .caver_java import PyJava


VERSION = '4.0.2'
url = "https://www.caver.cz/index.php?sid=123"


THE_20s=['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


class CaverPyMOL(QtWidgets.QWidget):
    # configuration binding from UI to CaverConfig
    config_bindings_main: Dict[str, str] = {
        'lineEdit_outputDir': 'output_dir',
        'comboBox_startPointSele': 'selection_name',
        'doubleSpinBox_x': 'start_point_x',
        'doubleSpinBox_y': 'start_point_y',
        'doubleSpinBox_z': 'start_point_z',
    }

    config_bindings_config: Dict[str, str] = {
        # from config.txt
        'spinBox_maxJavaHeapSize': 'customized_java_heap',
        'doubleSpinBox_maxProRad': 'probe_radius',
        'doubleSpinBox_shellRad': 'shell_radius',
        'spinBox_shellDepth': 'shell_depth',
        'doubleSpinBox_clusterThreshold': 'clustering_threshold',
        'comboBox_numApproxBalls': 'number_of_approximating_balls',

        'doubleSpinBox_maxDist': 'max_distance',
        'doubleSpinBox_desiredDist': 'desired_radius'
    }
    
    # reversed
    config_bindings_main_rev={v:k for k,v in config_bindings_main.items()}
    config_bindings_config_rev={v:k for k,v in config_bindings_config.items()}
    
    

    def bind_config(self):
        """
        Binds UI widgets to their corresponding configuration items.
        
        Iterates through the config_bindings dictionary where keys are widget names 
        and values are the corresponding configuration item names. For each widget, 
        it connects the widget's signal to the _wiget_link method using the widget name 
        as an argument. This ensures that changes in the UI are reflected in the configuration.
        """
        for wn in self.config_bindings_main:
            widget = getattr(self.ui, wn)
            widget_signal_tape(widget, partial(self._wiget_link, wn))

        for wn in self.config_bindings_config:
            widget = getattr(self.ui_config, wn)
            widget_signal_tape(widget, partial(self._wiget_link, wn))

    def _wiget_link(self, widget_name):
        """
        Links a widget's value with the corresponding configuration item.
    
        Retrieves the configuration item associated with the widget from the config_bindings dictionary.
        If the configuration item does not exist in the config object, an AttributeError is raised.
        
        Parameters:
        - widget_name: The name of the widget, used to look up the corresponding configuration item.
        
        Raises:
        - AttributeError: If the configuration item is not found in the config object.
        """
        # Retrieve the configuration item associated with the widget
        if widget_name in self.config_bindings_config: 
            config_item = self.config_bindings_config[widget_name]
            ui=self.ui_config
        else:
            config_item = self.config_bindings_main[widget_name]
            ui=self.ui

        
        # Check if the configuration item exists in the config object
        if not hasattr(self.config, config_item):
            raise AttributeError(f'{config_item} not found in config')
    
        # Get the current value of the configuration item
        pv = getattr(self.config, config_item)
        
        # Get the widget object
        widget = getattr(ui, widget_name)
        
        # Get the current value of the widget
        nv = get_widget_value(widget)
        
        # Update the value of the configuration item with the widget's value
        self.config.set_value(config_item, nv)
    
        # Get the updated value of the configuration item
        uv = getattr(self.config, config_item)
        
        # Log the change in value
        logging.info(f'{widget_name} {pv} -> {nv} -> {uv}')

    def refresh_window_from_cfg(self):
        '''
        Refresh the windows from the configuration file
        '''
        for wn, cn in self.config_bindings_main.items():
            widget = getattr(self.ui, wn)
            set_widget_value(widget, getattr(self.config, cn))
        
        for wn, cn in self.config_bindings_config.items():
            widget = getattr(self.ui_config, wn)
            set_widget_value(widget, getattr(self.config, cn))

    def make_window(self):
        '''
        Make windows and load config
        '''
        main_window = QtWidgets.QWidget()
        self.ui = CaverUI()
        self.ui.setupUi(main_window)

        self.ui_config=CaverConfigForm()
        config_window=QtWidgets.QWidget()
        self.ui_config.setupUi(config_window)

        self.ui.pushButton_help.clicked.connect(self.launchHelp)
        self.bind_config()
        self.ui.pushButton_compute.clicked.connect(self.execute)
        self.ui.pushButton_convertStartPointSele.clicked.connect(self.convert_sele_to_coords)
        self.ui.pushButton_reloadInputModel.clicked.connect(self.update_model_list)

        self.ui.pushButton_loadConfig.clicked.connect(self.configin)
        self.ui.pushButton_saveConfig.clicked.connect(self.configout)

        self.ui.doubleSpinBox_x.valueChanged.connect(self.changeCoords)
        self.ui.doubleSpinBox_y.valueChanged.connect(self.changeCoords)
        self.ui.doubleSpinBox_z.valueChanged.connect(self.changeCoords)
        self.ui.pushButton_openOutputDir.clicked.connect(
            lambda: self.ui.lineEdit_outputDir.setText(
                getExistingDirectory()))
        self.ui.comboBox_startPointSele.currentIndexChanged.connect(self._analysis_model_resn)
        self.ui.textEdit_startpoint.textChanged.connect(self._use_custom_startpoint)

        return main_window, config_window
    
    @contextmanager
    def freeze_window(self):
        """
        Freezes the dialog while the plugin is running.
        """
        self.dialog.setEnabled(False)
        try:
            yield
        except Exception as e:
            logging.error(f"Error occurred: {e}")
        self.dialog.setEnabled(True)

    def run_plugin_gui(self):
        """PyMOL entry for running the plugin"""
        super().__init__()
        # global reference to avoid garbage collection of our dialog
        self.dialog = None
        self.config = CaverConfig()

        # ignore structures which match the follwing regexps
        self.ignoreStructures = [r"^origins$", r"_origins$", r"_v_origins$", r"_t\d\d\d_\d$"]

        if self.dialog is None:
            self.dialog,self.config_dialog = self.make_window()
        self.dialog.show()

        self.ui.pushButton_editConfig.clicked.connect(self.config_dialog.show)


        # aa bias
        self.checktable_aa = CheckableListView(self.ui.listView_residueType, {aa: aa for aa in THE_20s})
        self.ui.pushButton_allAA.clicked.connect(self.checktable_aa.check_all)
        self.ui.pushButton_noneAA.clicked.connect(self.checktable_aa.uncheck_all)
        self.ui.pushButton_reverseAAsel.clicked.connect(self.checktable_aa.reverse_check)
        self.checktable_aa.checkStateChanged.connect(self._update_aa_sel)
        self.ui.pushButton_proteinResn.clicked.connect(lambda: self.checktable_aa.check_these(THE_20s, clear_before_check=True))
        self.ui.pushButton_ligandResn.clicked.connect(
            lambda: self.checktable_aa.check_these([x for x in self.checktable_aa.items.keys() if x not in THE_20s], clear_before_check=False)
            )
        self.ui.pushButton_RefreshSelection.clicked.connect(self._update_pymol_sel)
        self.ui.pushButton_clearStartPointSele.clicked.connect(self._clear_pymol_sel_and_coords)
        self.ui.pushButton_RefreshRunID.clicked.connect(self._update_run_id)
        self.ui.pushButton_LoadRunID.clicked.connect(lambda: self._playback_run_id(get_widget_value(self.ui.comboBox_RunID)))

        self.ui.pushButton_cite.clicked.connect(self.cite_info)
        self.ui.pushButton_doc.clicked.connect(self.open_doc_pdf)

        # register as a pymol command
        cmd.extend('caver_set', self.caver_set)

        # autocompletion for key
        cmd.auto_arg[0]["caver_set"] = [lambda : CaverShortcut(config=self.config,keywords=self.config.all_keys), "Caver setting key", ", "]

        # dynamic autocompletion for value (current)
        cmd.auto_arg[1]["caver_set"] = [
            lambda: Shortcut(
                keywords=[
                    self.config.get(self.config._complete_temp) if self.config.has(self.config._complete_temp) else ' '
                    ]
                ), "Caver setting value", ""]

        self.configin(CONFIG_TXT)

        self.update_model_list()

        self._analysis_model_resn()


    def _update_run_id(self):
        run_ids=self.get_run_ids()[1]
        set_widget_value(self.ui.comboBox_RunID, run_ids or [])

    def _playback_run_id(self, run_id: str):
        '''
        Play Back historical runs according to the run id.

        Parameters:
            - run_id [str]: Run ID of Caver.
        '''
        if not run_id.isdigit():
            notify_box(f"Run ID '{run_id}' is not a valid number.", ValueError)

        
    
        run_id= int(run_id)
        out_home, idxs = self.get_run_ids()

        if not run_id in idxs:
            notify_box(f"Run ID '{run_id}' does not exist in the output directory.", ValueError)

        expected_view_file = os.path.join(out_home, str(run_id), "pymol", "view_plugin.py")
        if not os.path.isfile(expected_view_file):
            notify_box(f"Run ID '{run_id}' does not contain a valid output file ({expected_view_file}).", FileNotFoundError)
        
        with self.freeze_window():
            # Run the PyMOL view plugin to visualize the results
            runview = f"run {expected_view_file}" 
            logging.debug(runview)
            cmd.do(runview)

    def _update_pymol_sel(self):
        selections: list[str]=cmd.get_names(type="selections")
        set_widget_value(self.ui.comboBox_startPointSele, selections)

    def _use_custom_startpoint(self):
        startpoint_val=get_widget_value(self.ui.textEdit_startpoint)
        
        if not startpoint_val:
            notify_box("Please specify starting point like atom ids, residue ids or x y z coordinates.", ValueError)

        for i in ('starting_point_atom', 'starting_point_residue', 'starting_point_coordinates'):
            self.config.set(i, '???')
        
        if self.ui.radioButton_startAsAtoms.isChecked():
            self.config.set('starting_point_atom', startpoint_val)
        elif self.ui.radioButton_startAsResidues.isChecked():
            self.config.set('starting_point_residue', startpoint_val)
        elif self.ui.radioButton_startAsCoords.isChecked():
            self.config.set('starting_point_coordinates', startpoint_val)
        else:
            notify_box('To use custom starting point, please select the corresponding radio button.', ValueError)

    def _clear_pymol_sel_and_coords(self):
        self.ui.comboBox_startPointSele.clear()
        self.ui.doubleSpinBox_x.setValue(0)
        self.ui.doubleSpinBox_y.setValue(0)
        self.ui.doubleSpinBox_z.setValue(0)
        cmd.delete("crisscross")

    def _update_aa_sel(self, aa_sel: Optional[List[str]]):
        if not aa_sel:
            if not self.config.has('include_residue_names'):
                return
            self.config.delete('include_residue_names')
            return
        self.config.set('include_residue_names', ' '.join(aa_sel))
        return

    def showCrisscross(self):
        cmd.delete("crisscross")
        CaverPyMOL.crisscross(
            self.config.start_point_x,
            self.config.start_point_y,
            self.config.start_point_z,
            0.5,
            "crisscross")

    def changeCoords(self, *args):
        self.showCrisscross()
        self.config.set(
            'starting_point_coordinates', 
            f'{get_widget_value(self.ui.doubleSpinBox_x)} {get_widget_value(self.ui.doubleSpinBox_y)} {get_widget_value(self.ui.doubleSpinBox_z)}'
            )
        
        

    def structureIgnored(self, name):
        for key in self.ignoreStructures:
            if re.search(key, name):
                return 1
        return 0

    def update_model_list(self):
        '''
        Update the list of models in the combobox
        and set the current model to the first one

        '''
        models=[str(i) for i in cmd.get_object_list() if not self.structureIgnored(str(i))]
        set_widget_value(self.ui.comboBox_inputModel, models)

        self._analysis_model_resn()

    def launchHelp(self):
        webbrowser.open(url)

    def loadFileContent(self, file):
        handler = open(file)
        lines = handler.readlines()
        wresult = ""
        for line in lines:
            wresult += line
        return wresult


    def get_run_ids(self) -> tuple[str, list[int]]:
        '''
        Get all run IDs in the output directory.
        Returns:
        - out_home (str): The path to the output directory.
        - idxs (list[int]): A list of all run IDs found in the output directory.
        '''
        dir = self.config.output_dir

        if not dir:
            notify_box("Please specify output directory.", ValueError)

        os.makedirs(dir, exist_ok=True)

        dir = dir.replace("\\", "/")
        if (dir.endswith("/")):
            dir = dir[:-1]

        out_home = os.path.join(dir, "caver_output")
        os.makedirs(out_home, exist_ok=True)

        # only successful runs are saved with view_plugin.py file
        idxs = [int(x) for x in os.listdir(out_home) if x.isdigit() if os.path.isfile(os.path.join(out_home, str(x), "pymol", "view_plugin.py"))]
        
        if not idxs:
            notify_box(f"Output directory '{out_home}' does not contain any run files, perhaps you dont have run any analysis yet.", ValueError)
        
        return out_home,idxs

    def initialize_out_dir(self):
        '''
        Initialize the output directory for the current run.
        A new run ID will be generated by adding 1 to the maximum run ID found in the output directory.

        '''

        out_home, idxs = self.get_run_ids()
        max_idx = max(idxs)

        new_dir = os.path.join(out_home, str(max_idx + 1))
        os.makedirs(new_dir)

        self.out_dir = new_dir
        logging.info("Output will be stored in " + self.out_dir)

    @property
    def coordinatesNotSet(self) -> bool:
        '''
        Checks if the coordinates are set.

        Returns
        bool
            True if the coordinates are not set, False otherwise.
        '''
        return all(self.config.get(f'start_point_{i}') == 0 for i in 'xyz')

    @property
    def customStartPointNotSet(self) -> bool:
        '''
        Check if a custom starting point is set in the UI.
        Returns:
            - bool: True if no custom starting point is set, False otherwise.
        '''
        return get_widget_value(self.ui.textEdit_startpoint).strip() == ''

    def execute(self):
        """
        Executes the analysis process, including checking prerequisites, preparing the environment,
        creating configuration files, running Caver through Java, and handling the results.
        """
        # Check if coordinates are set, if not, prompt the user to set them
        if self.coordinatesNotSet and self.customStartPointNotSet:
            notify_box(
                "Please specify starting point - "
                "e.g. by selecting atoms or residues and clicking at the button 'Convert to x, y, z'.",
                ValueError)

        # Display the crisscross structure
        self.showCrisscross()
    
        # Get the selected model's name from the UI list widget
        selected_model = get_widget_value(self.ui.comboBox_inputModel)
    
        # Initialize the output directory
        self.initialize_out_dir()
    
        # Create a subdirectory for inputs
        outdirInputs = os.path.join(self.out_dir, 'input')
        os.makedirs(outdirInputs, exist_ok=True)
    
        if not self.ui.checkBox_MD.isChecked():
            # Save the selected model as a PDB file in the input directory
            input = os.path.join(outdirInputs, f'{selected_model}.pdb')
            cmd.set('retain_order', 1)
            cmd.sort()
            cmd.save(input, selected_model)  # to by ulozilo cely model selected_model.
        else:
            self.prepare_md_pdb_traj(selected_model, outdirInputs)
            
        
        # Get the path to the Caver JAR file
        caverjar = os.path.join(THIS_DIR, "caver.jar")
    
        # Create a new configuration file with a timestamp
        cfgTimestamp = time.strftime("%Y-%m-%d-%H-%M")
        cfgnew = os.path.join(outdirInputs, f"config_{cfgTimestamp}.txt")
        self.configout(cfgnew)
    
        # Initialize and run Caver through the PyJava interface
        try:
            pj = PyJava(self.config.customized_java_heap, THIS_DIR, caverjar, outdirInputs, cfgnew, self.out_dir)
        except RuntimeError as e:
            notify_box(str(e), RuntimeError)

        with hold_trigger_button(self.ui.pushButton_compute), self.freeze_window():
            ret=run_worker_thread_with_progress(pj.run_caver)

        # Check for out of memory errors in Caver's output
        if 'OutOfMemory' in ret.stdout or 'OutOfMemory' in ret.stderr:
            notify_box(
                'Insufficient memory.',
                details=f"Available memory ({pj.memory_heap_level} MB) is not sufficient to analyze this structure. "
                "Try to allocate more memory. 64-bit operating system and Java are needed to get over 1200 MB. "
                "Using smaller 'Number of approximating balls' can also help, but at the cost of decreased accuracy of computation.")
    
        # Store the current working directory
        prevDir = os.getcwd()
        logging.debug(prevDir)

        runview_file=os.path.join(self.out_dir, "pymol","view_plugin.py")
        if not os.path.isfile(runview_file):
            notify_box('No tunnel found!', RuntimeError)
    
        # Run the PyMOL view plugin to visualize the results
        runview = f"run {runview_file}" 
        logging.info(runview)
        cmd.do(runview)

    @staticmethod
    def fixPrecision(numberStr: Any) -> float:
        return math.floor(float(numberStr) * 1000) / 1000
    

    def convert_sele_to_coords(self):
        """
        Converts the current selection to coordinates and displays them.
        
        This method calculates the center of the selected area and displays its coordinates in the UI.
        If the selection does not exist or the session is empty, it will notify the user with a prompt.
        """
    
        # Check if the session is empty
        if len([a for a in cmd.get_model('(all)').atom]) == 0:
            notify_box("Session is empty", ValueError)
    
        # Explicitly prohibit PyMOL selection syntax to avoid confusion
        if self.config.selection_name not in cmd.get_names('selections'):
            notify_box("Selection does not exist. If you are using PyMOL selection syntax, "
                       "please create a new selection in PyMOL, then input the correct selection name", ValueError)
    
        # Calculate the center of the selection
        startpoint = self.compute_center(self.config.selection_name)
        if not startpoint:
            return
        # Update the UI with the coordinates of the center
        self.ui.doubleSpinBox_x.setValue(CaverPyMOL.fixPrecision(startpoint[0]))
        self.ui.doubleSpinBox_y.setValue(CaverPyMOL.fixPrecision(startpoint[1]))
        self.ui.doubleSpinBox_z.setValue(CaverPyMOL.fixPrecision(startpoint[2]))
    
        # Mark the center point in the 3D space
        CaverPyMOL.crisscross(startpoint[0], startpoint[1], startpoint[2], 0.5, "crisscross")
        self.showCrisscross()

    def configin(self, filepath: Optional[str] = None):
        """
        Load configuration from a file. If no filepath is provided, a file dialog will prompt the user to select a file.
        
        Parameters:
        - filepath: Optional[str] - The path to the configuration file. If not provided, the user will be prompted to select a file.
        
        Returns:
        None
        """
        # Prompt the user to select a configuration file if no filepath is provided
        filepath = filepath or getOpenFileNameWithExt(
            self.dialog,
            "Select configuration file",
            filter="JSON ( *.json );;TXT ( *.txt )")
        if not filepath:
            return
    
        # Load the configuration from the selected file, depending on its extension
        self.config = CaverConfig.from_json(filepath) if filepath.endswith(".json") else CaverConfig.from_txt(filepath)
        # refresh window wiget from input config
        self.refresh_window_from_cfg()
        # Update the start point based on the loaded configuration if it exists
        self.refresh_start_point_from_cfg()
        # Update the configuration status label with the loaded filename
        set_widget_value(self.ui.label_configStatus, f'Loaded from {os.path.basename(filepath)}')
        # Perform post-processing on the configuration
        self.config_post_process()

    def refresh_start_point_from_cfg(self):

        if not self.config.has('starting_point_coordinates') or self.config.get('starting_point_coordinates') == '???':
            return

        coords_from_config = tuple(map(float, self.config.get('starting_point_coordinates').split(' ')))
        if not len(coords_from_config) == 3:
            notify_box("Invalid starting point coordinates in configuration file", ValueError)

        for (i, axis), (j, coord) in zip(enumerate('xyz'), enumerate(coords_from_config)):
            set_widget_value(getattr(self.ui, f'doubleSpinBox_{axis}'), coord)
        notify_box(f"Starting point coordinates loaded from configuration file: {coords_from_config}")

    def configout(self, filepath: Optional[str] = None):
        filepath = filepath or getSaveFileNameWithExt(
            self.dialog,
            "Select configuration file",
            filter="JSON ( *.json );;TXT ( *.txt )")
        if not filepath:
            return
        self.config.to_json(filepath) if filepath.endswith(".json") else self.config.to_txt(filepath)
        set_widget_value(self.ui.label_configStatus, f'Saved as {os.path.basename(filepath)}')

    def config_post_process(self):

        conflict_flags = {
            'starting_point_coordinates': ['starting_point_atom', 'starting_point_residue']
        }
        for primary_flag, conflict_flag_list in conflict_flags.items():
            if not self.config.has(primary_flag):
                continue
            for flag in conflict_flag_list:
                if self.config.has(flag) and self.config.get(flag) != '???':
                    notify_box(
                        f'Conflict between {primary_flag} and {flag}',
                        details=f'Simultaneous usage of {primary_flag} parameter with {flag} parameters is not supported by plugin. '
                        f'{flag} is ignored.')
                    self.config.delete(flag)


        self.ensure_residue_names_to_checktable()

    def ensure_residue_names_to_checktable(self):
        if self.config.has("include_residue_names"):
            aa_from_config: List[str] = self.config.get("include_residue_names").split(" ")
            if aa_from_config:
                self.checktable_aa.check_these(aa_from_config)

    def _analysis_model_resn(self):
        sel = cmd.get_model('(all)')

        for a in sel.atom:
            if a.resn in self.checktable_aa.items:
                continue
            if a.resn not in THE_20s:
                self.checktable_aa.update({a.resn: a.resn})

    def getAtoms(self, selection="(all)"):
        return cmd.identify(selection, 0)

    def getResids(self, selection="(all)"):
        stored.list = []
        cmd.iterate(selection, "stored.list.append((resi,chain))")
        return set(stored.list)

    def getObjectName(self, selection="(all)"):
        pairs = cmd.identify(selection, 1)
        name = None
        names = set()
        for p in pairs:
            names.add(p[0])
        if 0 == len(names):
            notify_box("Selection is empty.")
        elif 1 == len(names):
            name = names.pop()
        else:
            s = "Starting point selection need to be limited to one object. Currently, it includes these objects: "
            for n in names:
                s += n + ' '
            notify_box(s)
        return name

    def compute_center(self, selection="(all)"):
        if selection not in cmd.get_names("selections") and selection not in cmd.get_names("objects"):
            notify_box(f"Selection '{selection}' does not exist, using all atoms.")
            selection = "all"
        object = self.getObjectName(selection)
        if None is object:
            return None
        Ts = []
        residues = self.getResids(selection)  # SET1
        atoms = self.getAtoms(selection)     # SET2
        for r in residues:
            r_sel = 'resi ' + str(r[0]) + ' and chain ' + r[1] + ' and object ' + object
            residue_atoms = self.getAtoms(r_sel)
            all = []
            for a in residue_atoms:
                if a in atoms:
                    all = all + [a]
            if len(all) == len(residue_atoms):
                Ts = Ts + [self.computecenterRA(r_sel)]
            else:
                for a in all:
                    Ts = Ts + [self.computecenterRA('id ' + str(a) + ' and object ' + object)]

        logging.info('Centers: %s' % ', '.join(map(str, Ts)))
        sumx = 0
        sumy = 0
        sumz = 0
        if len(Ts) == 0:
            return (0, 0, 0)
        l = len(Ts)
        for center in Ts:
            sumx += center[0]
            sumy += center[1]
            sumz += center[2]
        logging.info('Starting point: ' + str(sumx) + " " + str(sumy) + " " + str(sumz) + " " + str(l))
        return (sumx / l, sumy / l, sumz / l)

    # compute center for given selection
    def computecenterRA(self, selection="(all)"):
        stored.xyz = []
        cmd.iterate_state(1, selection, "stored.xyz.append([x,y,z])")
        centx = 0
        centy = 0
        centz = 0
        cnt = 0
        for a in stored.xyz:
            centx += a[0]
            centy += a[1]
            centz += a[2]
            cnt += 1
        centx /= cnt
        centy /= cnt
        centz /= cnt
        return (centx, centy, centz)

    # TODO: unused?
    def computecenter(self, selection="(all)"):
        gcentx = 0
        gcenty = 0
        gcentz = 0
        gcnt = 0
        for selstr in selection.split():
            sel = cmd.get_model(selstr)

            centx = 0
            centy = 0
            centz = 0
            cnt = len(sel.atom)
            if (cnt == 0):
                warnings.warn(UserWarning('selection used to compute starting point is empty'))
                return (0, 0, 0)
            for a in sel.atom:
                centx += a.coord[0]
                centy += a.coord[1]
                centz += a.coord[2]
            centx /= cnt
            centy /= cnt
            centz /= cnt

            gcentx += centx
            gcenty += centy
            gcentz += centz
            gcnt += 1

        gcentx /= gcnt
        gcenty /= gcnt
        gcentz /= gcnt
        return (gcentx, gcenty, gcentz)

    @staticmethod
    def crisscross(x, y, z, d, name="crisscross"):

        obj = [
            LINEWIDTH, 3,

            BEGIN, LINE_STRIP,
            VERTEX, float(x - d), float(y), float(z),
            VERTEX, float(x + d), float(y), float(z),
            END,

            BEGIN, LINE_STRIP,
            VERTEX, float(x), float(y - d), float(z),
            VERTEX, float(x), float(y + d), float(z),
            END,

            BEGIN, LINE_STRIP,
            VERTEX, float(x), float(y), float(z - d),
            VERTEX, float(x), float(y), float(z + d),
            END

        ]
        view = cmd.get_view()
        cmd.load_cgo(obj, name)
        cmd.set_view(view)


    def prepare_md_pdb_traj(self, selected_model: str, input_dir: str):
        '''
        Generate pdbs from multi-state PyMOL session for selected model

        Parameters
            selected_model: str
            input_dir: str

        '''
        state_start=get_widget_value(self.ui.spinBox_MD_StateMin)
        state_stop =get_widget_value(self.ui.spinBox_MD_StateMax)

        if state_start>state_stop:
            notify_box("Starting state is greater than the stoping state", ValueError)
        
        for state in range(state_start+1, state_stop+1):
            cmd.save(os.path.join(input_dir, f'{state}.pdb'), state=state, selection=selected_model)

    def cite_info(self):
        '''
        Cite info
        '''
        citation_file=os.path.join(THIS_DIR, 'bin', 'citation.txt')
        with open(citation_file, 'r') as f:
            citation_text=f.read()
        
        notify_box('Thank you for using Caver, please cite the following paper.', details=citation_text)

    def open_doc_pdf(self):
        '''
        Open doc pdf
        '''
        doc_file=os.path.join(THIS_DIR, 'config', 'caver_userguide.pdf')

        webbrowser.open(f'file://{doc_file}')

    def caver_set(self, key, value):
        '''
        Setting caver configurations

        Parameters
        key : str
            key of the configuration item
        value : str
            value of the configuration item
        '''
        if not self.config.has(key):
            warnings.warn(UserWarning(f'{key} is not a valid Caver setting. Will added it to the config.'))
        
        d={}
        ui=None
        if key in self.config_bindings_main_rev:
            ui=self.ui
            d=self.config_bindings_main_rev
        
        elif key in self.config_bindings_config_rev:
            ui=self.ui_config
            d=self.config_bindings_config_rev

        if d:
            wn=d.get(key)
            logging.debug(f'Setting {wn} -> {value}')
            set_widget_value(getattr(ui, wn), value)
        else:
            logging.debug(f'Setting config {key} -> {value}')
            self.config.set(key, value)
        

