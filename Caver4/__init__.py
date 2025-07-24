#!/usr/bin/env python

# CAVER Copyright Notice
# ============================
#


from contextlib import contextmanager
import json
import logging
import math
import os
import re
from dataclasses import dataclass
from functools import cached_property, partial
from typing import TYPE_CHECKING, Any, Dict, List, Optional

from pymol import cmd

if TYPE_CHECKING:
    from PyQt5 import QtWidgets
else:
    from pymol.Qt import QtWidgets

import time

from pymol import stored
from pymol.cgo import BEGIN, END, LINE_STRIP, LINEWIDTH, VERTEX
from pymol.plugins import addmenuitemqt
from pymol.Qt.utils import getSaveFileNameWithExt

from .ui.Ui_caver import Ui_CaverUI as CaverUI
from .utils.live_run import run_command
from .utils.ui_tape import (CheckableListView, get_widget_value,
                            getExistingDirectory, getOpenFileNameWithExt,
                            notify_box, set_widget_value, widget_signal_tape, hold_trigger_button, run_worker_thread_with_progress)

THIS_DIR = os.path.dirname(__file__)
CONFIG_TXT = os.path.join(THIS_DIR, "config", "config.txt")


VERSION = '4.0.0'


THE_20s=['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


@dataclass
class CaverConfig:

    # connect to wigets
    output_dir: str = ""

    customized_java_heap: int = 6000

    shell_radius: float = 3.0
    shell_depth: int = 2
    probe_radius: float = 0.7
    clustering_threshold: float = 1.5

    number_of_approximating_balls: int = 4
    ignore_water: bool = False

    selection_name: str = ''

    start_point_x: float = 0.0
    start_point_y: float = 0.0
    start_point_z: float = 0.0

    max_distance: float = 4.0
    desired_radius: float = 1.8

    def has(self, key: str) -> bool:
        return hasattr(self, key)

    def delete(self, key: str):
        delattr(self, key)

    def get(self, key: str) -> Any:
        return getattr(self, key)

    def set(self, key: str, value: Any):
        setattr(self, key, value)

    @classmethod
    def from_json(cls, json_file: str) -> 'CaverConfig':
        """
        Load configuration information from a JSON file and create a CaverConfig object.
    
        This method is a class method, which means it can be called directly on the class rather than an instance of the class.
        It takes a path to a JSON file as input and returns an instance of CaverConfig initialized with the data from the JSON file.
    
        Parameters:
        - json_file: str - The path to the JSON file containing the configuration information.
    
        Returns:
        - CaverConfig: An instance of the CaverConfig class initialized with the data from the JSON file.
        """
        # Load JSON data from the specified file path
        json_data = json.load(open(json_file))
        # filter dataclass keys with fixed typing
        # `__match_args__` has all dataclass keys
        # `__dataclass_fields__` has all dataclass key fields
        # Initialize a new instance of the class using the data from the JSON, converting types as necessary
        new_self = cls(**{k: cls.__dict__['__dataclass_fields__'][k].type(v)
                       for k, v in json_data.items() if k in cls.__dict__['__match_args__']})
        # add extra keys that come from json
        # Iterate through the JSON data again, adding any extra keys not included in the initial dataclass keys
        for k, v in json_data.items():
            if k in cls.__dict__['__match_args__']:
                continue
            new_self.set(k, v)
        return new_self

    def to_json(self, json_file: str):
        json.dump(self.__dict__, open(json_file, 'w'))

    def set_value(self, key: str, value: Any):
        setattr(self, key, value)

    @classmethod
    def from_txt(cls, txt_file: Optional[str] = None):
        '''
        Load configuration from a specified text file or the default configuration file if none is provided.
        
        This method reads the configuration from a text file, where each line represents a configuration item.
        The format of each line is "key value", with the following rules:
        - Lines starting with "#" are comments and are ignored.
        - If "#" appears in the middle of a line, the part after "#" is treated as a comment and is ignored.
        - Boolean values: "yes" represents True, "no" represents False.
        - String values: If they contain spaces, they must be enclosed in double quotes "".
        - List values: Space-separated.
        
        Parameters:
        - txt_file (Optional[str]): The path to the text file containing the configuration. If not provided, the default configuration file is used.
        
        Returns:
        - An instance of the class initialized with the configuration read from the file.
        '''
        # Initialize a new instance of the class
        new_self = cls()
        
        # Open the specified configuration file or the default configuration file
        with open(txt_file or CONFIG_TXT) as f:
            # Read each line of the file
            for l in f.readlines():
                # Remove leading and trailing whitespace
                l = l.strip()
                
                # If the line contains "#" and does not start with "#", remove the comment part
                if '#' in l and not l.startswith('#'):
                    l = l[0:l.rfind("#") - 1]
                    # remove trailing whitespaces
                    l = l.rstrip(' ')
                
                # If the line starts with "#" or is empty, skip it
                if l.startswith("#") or not l:
                    continue
                
                # Split the line into key and value parts
                parsed = l.split(' ', 1)
                key = parsed[0]
                # If the key part is empty or there is no value part, skip it
                if len(parsed) <= 1:
                    print('skipping ' + key)
                    continue
                
                # Set the value for the key in the new instance
                new_self._set_value(key, parsed[1])
        
        # Return the initialized new instance
        return new_self

    def _set_value(self, key: str, new_value: str):
        if hasattr(self, key):
            v_type = type(getattr(self, key))
            # align to the self variable type
            setattr(self, key, v_type(new_value) if v_type != bool else CaverConfig._true_or_false(new_value))
            return
        # unknown variable, just set it as str
        setattr(self, key, new_value)

    def _get_value(self, key: str) -> str:
        # if the variable is a defined slot in this dataclass, return with the correct type
        value = getattr(self, key)
        v_type = type(value)

        # if it's a boolean, return yes or no according to the config.txt
        if v_type == bool:
            return CaverConfig._yes_or_no(value)
        return str(value)

    @staticmethod
    def _yes_or_no(v: bool) -> str:
        return 'yes' if v else 'no'

    @staticmethod
    def _true_or_false(v: str) -> bool:
        return True if v == 'yes' else False

    def to_txt(self, txt_file: str):
        '''
        Export the current configuration to a TXT file for the JAR file to read.
        The format of the TXT file is as follows:
        - Boolean values are represented by "True" or "False"
        - List values are separated by spaces
    
        Parameters:
        - txt_file (str): The path to the TXT file to be exported
    
        This function does not return any value.
        '''
        # Set the starting point coordinates if any of the start_point_x, start_point_y, or start_point_z attributes are not 0
        if any(getattr(self, f'start_point_{a}') != 0 for a in 'xyz'):
            self.set('starting_point_coordinates', " ".join(str(getattr(self, f'start_point_{a}')) for a in 'xyz'))
            logging.debug(f"starting point: {self.get('starting_point_coordinates')}")
    
        # Initialize a new list to store the updated configuration content
        new_txt_contents = []
        # Read the configuration template file
        with open(CONFIG_TXT) as template_f:
            template = template_f.readlines()
        for l in template:
            l = l.strip()
            # Directly append comment lines and empty lines to the new configuration content
            if l.startswith('#') or not l:
                new_txt_contents.append(l)
                continue
            # Remove everything after the last occurrence of the # character and trailing whitespaces for lines containing inline comments
            if '#' in l and not l.startswith('#'):
                l = l[0:l.rfind("#") - 1]
                l = l.rstrip(' ')
    
            # Split each line into a key-value pair
            parsed = l.split(' ', 1)
            key = parsed[0]
            # Skip lines with no value set
            if len(parsed) <= 1:
                print('skipping ' + key)
            # Update keys that exist in the TXT file but not in the current configuration
            elif hasattr(self, key):
                self_val = getattr(self, key)
                _val_type = type(self_val)
                # Handle boolean values
                if _val_type == bool:
                    new_txt_contents.append(f'{key} {CaverConfig._yes_or_no(self_val)}')
                # Skip unset values
                elif self_val == '???':
                    continue
                # Handle float values
                elif _val_type == float:
                    new_txt_contents.append(f'{key} {self_val:.1f}')
                # Handle other types of values
                else:
                    new_txt_contents.append(f'{key} {str(self_val)}')
    
        # Create the directory for the TXT file if it does not exist
        os.makedirs(os.path.dirname(txt_file), exist_ok=True)
        # Write the updated configuration content to the TXT file
        with open(txt_file, 'w') as f:
            f.write('\n'.join(new_txt_contents))



url = "http://www.caver.cz/index.php?sid=123"


class PyJava:
    def __init__(self, customized_memory_heap, caverfolder, caverjar, outdirInputs, cfgnew, out_dir):
        self.jar = caverjar

        print("\n*** Testing if Java is installed ***")
        if not self.java_present:
            notify_box("Java is not installed. Please install Java and try again.", RuntimeError)

        print("\n*** Optimizing memory allocation for Java ***")
        self.optimize_memory(customized_memory_heap)
        self.cmd = [
            "java",
            "-Xmx%dm" % self.memory_heap_level,
            "-cp", os.path.join(caverfolder, "lib"),
            "-jar", caverjar,
            "-home", caverfolder,
            "-pdb", outdirInputs,
            "-conf", cfgnew,
            "-out", out_dir,
        ]
        print("*** Caver will be called using command ***")
        print(" ".join(['"%s"' % t if t != "java" and t[0] != "-" else t for t in self.cmd]))
        print("******************************************")

    @cached_property
    def java_present(self):
        return run_command(["java", "-version"], verbose=True).returncode == 0

    def run_caver(self):
        return run_command(self.cmd, verbose=True)

    def optimize_memory(self, customized_memory_heap):
        customized_memory_heap = int(customized_memory_heap)
        memory_allocate_options = [500, 800, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1400, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 14000, 16000, 20000, 32000, 48000, 64000]
        memory_allocate_options.append(customized_memory_heap)
        memory_allocate_options.sort()
        # sorted(values)
        self.memory_heap_level = memory_allocate_options[0]
        for heap_level in memory_allocate_options:
            if int(heap_level) <= customized_memory_heap:
                cmd = ["java", "-Xmx%dm" % heap_level, "-jar", self.jar, "do_nothing"]
                heap_level_test = run_command(cmd)

                if heap_level_test.returncode == 0:
                    self.memory_heap_level = heap_level
                    print("Memory heap level: " + str(self.memory_heap_level))
        print("*** Memory for Java: " + str(self.memory_heap_level) + " MB ***")


class AnBeKoM(QtWidgets.QWidget):
    config_bindings: Dict[str, str] = {
        'lineEdit_outputDir': 'output_dir',
        'spinBox_maxJavaHeapSize': 'customized_java_heap',
        'doubleSpinBox_maxProRad': 'probe_radius',
        'doubleSpinBox_shellRad': 'shell_radius',
        'spinBox_shellDepth': 'shell_depth',
        'doubleSpinBox_clusterThreshold': 'clustering_threshold',
        'comboBox_numApproxBalls': 'number_of_approximating_balls',
        'checkBox_ignoreWater': 'ignore_water',
        'lineEdit_startPointSele': 'selection_name',
        'doubleSpinBox_x': 'start_point_x',
        'doubleSpinBox_y': 'start_point_y',
        'doubleSpinBox_z': 'start_point_z',
        'doubleSpinBox_maxDist': 'max_distance',
        'doubleSpinBox_desiredDist': 'desired_radius'
    }
    

    def bind_config(self):
        """
        Binds UI widgets to their corresponding configuration items.
        
        Iterates through the config_bindings dictionary where keys are widget names 
        and values are the corresponding configuration item names. For each widget, 
        it connects the widget's signal to the _wiget_link method using the widget name 
        as an argument. This ensures that changes in the UI are reflected in the configuration.
        """
        for wn in self.config_bindings:
            widget = getattr(self.ui, wn)
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
        config_item = self.config_bindings.get(widget_name)
        
        # Check if the configuration item exists in the config object
        if not hasattr(self.config, config_item):
            raise AttributeError(f'{config_item} not found in config')
    
        # Get the current value of the configuration item
        pv = getattr(self.config, config_item)
        
        # Get the widget object
        widget = getattr(self.ui, widget_name)
        
        # Get the current value of the widget
        nv = get_widget_value(widget)
        
        # Update the value of the configuration item with the widget's value
        self.config.set_value(config_item, nv)
    
        # Get the updated value of the configuration item
        uv = getattr(self.config, config_item)
        
        # Log the change in value
        logging.info(f'{widget_name} {pv} -> {nv} -> {uv}')

    def refresh_window_from_cfg(self):
        for wn, cn in self.config_bindings.items():
            widget = getattr(self.ui, wn)
            set_widget_value(widget, getattr(self.config, cn))

    def make_window(self):
        main_window = QtWidgets.QWidget()
        self.ui = CaverUI()
        self.ui.setupUi(main_window)

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
        self.ui.lineEdit_startPointSele.textChanged.connect(self._analysis_model_resn)

        return main_window
    
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
            self.dialog = self.make_window()
        self.dialog.show()

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

        self.configin(CONFIG_TXT)

        self.update_model_list()

        self._analysis_model_resn()

    def _update_pymol_sel(self, selection: str):
        set_widget_value(self.ui.lineEdit_startPointSele, selection)

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
        AnBeKoM.crisscross(
            self.config.start_point_x,
            self.config.start_point_y,
            self.config.start_point_z,
            0.5,
            "crisscross")

    def changeCoords(self, *args):
        self.showCrisscross()

    def structureIgnored(self, name):
        for key in self.ignoreStructures:
            if re.search(key, name):
                return 1
        return 0

    def update_model_list(self):
        self.ui.listWidget_inputModel.clear()
        self.ui.listWidget_inputModel.addItems(
            [str(i) for i in cmd.get_object_list() if not self.structureIgnored(str(i))])
        self.ui.listWidget_inputModel.setCurrentItem(self.ui.listWidget_inputModel.item(0))

        self._analysis_model_resn()

    def launchHelp(self):
        import webbrowser
        webbrowser.open(url)

    def loadFileContent(self, file):
        handler = open(file)
        lines = handler.readlines()
        wresult = ""
        for line in lines:
            wresult += line
        return wresult

    def initialize_out_dir(self):

        dir = self.config.output_dir

        os.makedirs(dir, exist_ok=True)

        dir = dir.replace("\\", "/")
        if (dir.endswith("/")):
            dir = dir[:-1]

        out_home = os.path.join(dir, "caver_output")
        os.makedirs(out_home, exist_ok=True)

        idxs = map(int, [x for x in os.listdir(out_home) if x.isdigit()] or [0])
        max_idx = max(idxs)

        new_dir = os.path.join(out_home, str(max_idx + 1))
        os.makedirs(new_dir)

        self.out_dir = new_dir
        print("Output will be stored in " + self.out_dir)

    @property
    def coordinatesNotSet(self) -> bool:
        return all(self.config.get(f'start_point_{i}') == 0 for i in 'xyz')

    def execute(self):
        # Check if coordinates are set, if not, prompt the user to set them
        if self.coordinatesNotSet:
            notify_box(
                "Please specify starting point - "
                "e.g. by selecting atoms or residues and clicking at the button 'Convert to x, y, z'.",
                ValueError)

        with hold_trigger_button(self.ui.pushButton_compute), self.freeze_window():
            ret=run_worker_thread_with_progress(self._execute)

        # Check for out of memory errors in Caver's output
        if 'OutOfMemory' in ret.stdout or 'OutOfMemory' in ret.stderr:
            notify_box(
                'Insufficient memory.',
                details=f"Available memory ({pj.memory_heap_level} MB) is not sufficient to analyze this structure. "
                "Try to allocate more memory. 64-bit operating system and Java are needed to get over 1200 MB. "
                "Using smaller 'Number of approximating balls' can also help, but at the cost of decreased accuracy of computation.")
    
        # Store the current working directory
        prevDir = os.getcwd()
        print(prevDir)
    
        # Run the PyMOL view plugin to visualize the results
        runview = "run " + self.out_dir + "/pymol/view_plugin.py"
        print(runview)
        cmd.do(runview)
    def _execute(self):
        """
        Executes the analysis process, including checking prerequisites, preparing the environment,
        creating configuration files, running Caver through Java, and handling the results.
        """

        # Display the crisscross structure
        self.showCrisscross()
    
        # Get the selected model's name from the UI list widget
        selected_model = self.ui.listWidget_inputModel.currentItem().text()
    
        # Initialize the output directory
        self.initialize_out_dir()
    
        # Create a subdirectory for inputs
        outdirInputs = os.path.join(self.out_dir, 'input')
        os.makedirs(outdirInputs, exist_ok=True)
    
        # Save the selected model as a PDB file in the input directory
        input = os.path.join(outdirInputs, f'{selected_model}.pdb')
        cmd.set('retain_order', 1)
        cmd.sort()
        cmd.save(input, selected_model)  # to by ulozilo cely model selected_model.
    
        # Get the path to the Caver JAR file
        caverjar = os.path.join(THIS_DIR, "caver.jar")
    
        # Create a new configuration file with a timestamp
        cfgTimestamp = time.strftime("%Y-%m-%d-%H-%M")
        cfgnew = os.path.join(outdirInputs, f"config_{cfgTimestamp}.txt")
        self.configout(cfgnew)
    
        # Initialize and run Caver through the PyJava interface
        pj = PyJava(self.config.customized_java_heap, THIS_DIR, caverjar, outdirInputs, cfgnew, self.out_dir)
        return pj.run_caver()

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
        self.ui.doubleSpinBox_x.setValue(AnBeKoM.fixPrecision(startpoint[0]))
        self.ui.doubleSpinBox_y.setValue(AnBeKoM.fixPrecision(startpoint[1]))
        self.ui.doubleSpinBox_z.setValue(AnBeKoM.fixPrecision(startpoint[2]))
    
        # Mark the center point in the 3D space
        AnBeKoM.crisscross(startpoint[0], startpoint[1], startpoint[2], 0.5, "crisscross")
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
        # Update the start point based on the loaded configuration
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

        print('Centers: %s' % ', '.join(map(str, Ts)))
        # print('Centers: ' + Ts)
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
        print('Starting point: ' + str(sumx) + " " + str(sumy) + " " + str(sumz) + " " + str(l))
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
                print('warning: selection used to compute starting point is empty')
                return (0, 0, 0)
            for a in sel.atom:
                centx += a.coord[0]
                centy += a.coord[1]
                centz += a.coord[2]
            centx /= cnt
            centy /= cnt
            centz /= cnt
        #       fmttext="%lf\t%lf\t%lf\n" % (centx,centy,centz)
#               print(centx,centy,centz)
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


def __init_plugin__(app=None):
    """
    Add an entry to the PyMOL "Plugin" menu
    """

    plugin = AnBeKoM()
    addmenuitemqt("Caver NG", plugin.run_plugin_gui)
