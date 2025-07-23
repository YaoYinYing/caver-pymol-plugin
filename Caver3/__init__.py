#!/usr/bin/env python

# CAVER Copyright Notice
# ============================
#


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
                            notify_box, set_widget_value, widget_signal_tape)

THIS_DIR = os.path.dirname(__file__)
CONFIG_TXT = os.path.join(THIS_DIR, "config", "config.txt")


VERSION = '4.0.0'

CAVER3_LOCATION = os.path.dirname(__file__)

OUTPUT_LOCATION = os.path.abspath(".")

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
        json_data = json.load(open(json_file))
        # filter dataclass keys with fixed typing
        # `__match_args__` has all dataclass keys
        # `__dataclass_fields__` has all dataclass key fields
        new_self = cls(**{k: cls.__dict__['__dataclass_fields__'][k].type(v)
                       for k, v in json_data.items() if k in cls.__dict__['__match_args__']})
        # add extra keys that come from json
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
        comment:
            - started by #: ignored
        bool:
            - True: yes
            - False: no
        str:
            - quoted by "" if has space

        list:
            - space separated

        '''
        new_self = cls()
        with open(txt_file or CONFIG_TXT) as f:
            for l in f.readlines():
                l = l.strip()
                # 'option yes # comment' -> 'option yes'
                if '#' in l and not l.startswith('#'):
                    l = l[0:l.rfind("#") - 1]
                    # remove trailing whitespaces
                    l = l.rstrip(' ')

                if l.startswith("#") or not l:
                    continue

                parsed = l.split(' ', 1)
                key = parsed[0]
                if len(parsed) <= 1:
                    print('skipping ' + key)
                    continue

                new_self._set_value(key, parsed[1])

        return new_self

    @property
    def has_include_exclude(self) -> bool:
        not_allowed = [
            "include_residue_names",
            "include_residue_ids",
            "include_atom_numbers",
            "exclude_residue_names",
            "exclude_residue_ids",
            "exclude_atom_numbers"]
        return any(self.has(key=k) for k in not_allowed)

    def _set_value(self, key: str, new_value: str):
        if hasattr(self, key):
            v_type = type(getattr(self, key))
            # align to the self variable type
            setattr(self, key, v_type(new_value) if v_type != bool else CaverConfig._true_or_false(new_value))
            return
        # unknown variable, just set it as str
        setattr(self, key, new_value)

    def _get_value(self, key: str) -> str:
        value = getattr(self, key)
        v_type = type(value)
        if v_type == bool:
            return CaverConfig._yes_or_no(value)
        return CaverConfig._need_quote(str(value))

    @staticmethod
    def _yes_or_no(v: bool) -> str:
        return 'yes' if v else 'no'

    @staticmethod
    def _true_or_false(v: str) -> bool:
        return True if v == 'yes' else False

    @staticmethod
    def _need_quote(v: str) -> str:
        if ' ' not in v:
            return v
        # simple digit int array
        if all(_v.isdigit() for _v in v.split(' ') if v):
            return v
        # complex float array
        if all(re.match(r'^(-)?\d+(\.?\d+)?$', _v) for _v in v.split(' ') if v):
            return v

        # normal string by lacking quotes
        return '"' + v + '"' if ' ' in v and (v[0] != '"' and v[-1] != '"') else v

    def to_txt(self, txt_file: str):
        '''
        export to txt file so jar can read it.
        bool:
            - True: yes
            - False: no
        str:
            - quoted by "" if has space

        list:
            - space separated

        '''
        if any(getattr(self, f'start_point_{a}') != 0 for a in 'xyz'):
            self.set('starting_point_coordinates', " ".join(str(getattr(self, f'start_point_{a}')) for a in 'xyz'))
            logging.debug(f"starting point: {self.get('starting_point_coordinates')}")

        new_txt_contents = []
        with open(CONFIG_TXT) as template_f:
            template = template_f.readlines()
        for l in template:
            l = l.strip()
            if l.startswith('#') or not l:
                # direct append comment and empty lines
                new_txt_contents.append(l)
                continue
            if '#' in l and not l.startswith('#'):
                # remove everything after last occurence of # char
                l = l[0:l.rfind("#") - 1]
                # remove trailing whitespaces
                l = l.rstrip(' ')

            parsed = l.split(' ', 1)
            key = parsed[0]
            if len(parsed) <= 1:
                print('skipping ' + key)

            # update keys that only exists in the txt file
            elif hasattr(self, key):
                self_val = getattr(self, key)
                _val_type = type(self_val)
                if _val_type == bool:
                    new_txt_contents.append(f'{key} {CaverConfig._yes_or_no(self_val)}')
                elif self_val == '???':  # drop unset values
                    continue
                elif _val_type == str:
                    new_txt_contents.append(f'{key} {CaverConfig._need_quote(self_val)}')
                elif _val_type == float:
                    new_txt_contents.append(f'{key} {self_val:.1f}')
                else:
                    new_txt_contents.append(f'{key} {str(self_val)}')

        os.makedirs(os.path.dirname(txt_file), exist_ok=True)
        with open(txt_file, 'w') as f:
            f.write('\n'.join(new_txt_contents))


defaults = CaverConfig()

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
    }

    def bind_config(self):
        for wn in self.config_bindings:
            widget = getattr(self.ui, wn)
            widget_signal_tape(widget, partial(self._wiget_link, wn))

    def _wiget_link(self, widget_name):
        config_item = self.config_bindings.get(widget_name)
        if not hasattr(self.config, config_item):
            raise AttributeError(f'{config_item} not found in config')

        pv = getattr(self.config, config_item)
        widget = getattr(self.ui, widget_name)
        nv = get_widget_value(widget)
        self.config.set_value(config_item, nv)

        uv = getattr(self.config, config_item)
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
        self.ui.pushButton_reloadInputModel.clicked.connect(self.updateList)

        self.ui.pushButton_loadConfig.clicked.connect(self.configin)
        self.ui.pushButton_saveConfig.clicked.connect(self.configout)

        self.ui.doubleSpinBox_x.valueChanged.connect(self.changeCoords)
        self.ui.doubleSpinBox_y.valueChanged.connect(self.changeCoords)
        self.ui.doubleSpinBox_z.valueChanged.connect(self.changeCoords)
        self.ui.pushButton_openOutputDir.clicked.connect(
            lambda: self.ui.lineEdit_outputDir.setText(
                getExistingDirectory()))
        self.ui.lineEdit_startPointSele.textChanged.connect(self._analysis_sel_resn)

        return main_window

    def run_plugin_gui(self):
        """PyMOL entry for running the plugin"""
        super().__init__()
        # global reference to avoid garbage collection of our dialog
        self.window = None
        self.config = CaverConfig()

        # ignore structures which match the follwing regexps
        self.ignoreStructures = [r"^origins$", r"_origins$", r"_v_origins$", r"_t\d\d\d_\d$"]

        if self.window is None:
            self.window = self.make_window()
        self.window.show()

        # aa bias
        self.checktable_aa = CheckableListView(self.ui.listView_residueType, {aa: aa for aa in THE_20s})
        self.ui.pushButton_allAA.clicked.connect(self.checktable_aa.check_all)
        self.ui.pushButton_noneAA.clicked.connect(self.checktable_aa.uncheck_all)
        self.ui.pushButton_reverseAAsel.clicked.connect(self.checktable_aa.reverse_check)
        self.checktable_aa.checkStateChanged.connect(self._update_aa_sel)

        self.configin(CONFIG_TXT)

        self.updateList()

        self.config_post_process()

        self._analysis_sel_resn()

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

    def updateList(self):
        self.ui.listWidget_inputModel.clear()
        self.ui.listWidget_inputModel.addItems(
            [str(i) for i in cmd.get_object_list() if not self.structureIgnored(str(i))])
        self.ui.listWidget_inputModel.setCurrentItem(self.ui.listWidget_inputModel.item(0))

        self._analysis_sel_resn()

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

        if self.coordinatesNotSet:
            notify_box(
                "Please specify starting point - "
                "e.g. by selecting atoms or residues and clicking at the button 'Convert to x, y, z'.",
                ValueError)

        self.showCrisscross()
        selected_model = self.ui.listWidget_inputModel.currentItem().text()

        self.initialize_out_dir()

        # create subdirectory for inputs

        outdirInputs = os.path.join(self.out_dir, 'input')
        os.makedirs(outdirInputs, exist_ok=True)

        input = os.path.join(outdirInputs, f'{selected_model}.pdb')
        cmd.set('retain_order', 1)
        cmd.sort()
        cmd.save(input, selected_model)  # to by ulozilo cely model selected_model.

        caverjar = os.path.join(THIS_DIR, "caver.jar")

        # create new config
        cfgTimestamp = time.strftime("%Y-%m-%d-%H-%M")
        cfgnew = os.path.join(outdirInputs, f"config_{cfgTimestamp}.txt")
        self.configout(cfgnew)

        # set correct java options
        # javaOpts = JOPTS.replace("@", self.javaHeap.getvalue())

        pj = PyJava(self.config.customized_java_heap, THIS_DIR, caverjar, outdirInputs, cfgnew, self.out_dir)

        ret = pj.run_caver()

        if 'OutOfMemory' in ret.stdout or 'OutOfMemory' in ret.stderr:
            notify_box(
                'Insufficient memory.',
                details=f"Available memory ({pj.memory_heap_level} MB) is not sufficient to analyze this structure. "
                "Try to allocate more memory. 64-bit operating system and Java are needed to get over 1200 MB. "
                "Using smaller 'Number of approximating balls' can also help, but at the cost of decreased accuracy of computation.")

        prevDir = os.getcwd()
        print(prevDir)

        runview = "run " + self.out_dir + "/pymol/view_plugin.py"
        print(runview)
        cmd.do(runview)

    @staticmethod
    def fixPrecision(numberStr: Any) -> float:
        return math.floor(float(numberStr) * 1000) / 1000

    def convert_sele_to_coords(self):

        if len([a for a in cmd.get_model('(all)').atom]) == 0:
            notify_box("Session is empty", ValueError)

        # prohibit pymol selection syntax, explicitly.
        if self.config.selection_name not in cmd.get_names('selections'):
            notify_box("Selection does not exist. If you are using PyMOL selection syntax, "
                       "please create a new selection in PyMOL, then refresh the list", ValueError)

        startpoint = self.compute_center(self.config.selection_name)
        if not startpoint:
            return
        self.ui.doubleSpinBox_x.setValue(AnBeKoM.fixPrecision(startpoint[0]))
        self.ui.doubleSpinBox_y.setValue(AnBeKoM.fixPrecision(startpoint[1]))
        self.ui.doubleSpinBox_z.setValue(AnBeKoM.fixPrecision(startpoint[2]))

        AnBeKoM.crisscross(startpoint[0], startpoint[1], startpoint[2], 0.5, "crisscross")
        self.showCrisscross()

    def configin(self, filepath: Optional[str] = None):
        filepath = filepath or getOpenFileNameWithExt(
            self.window,
            "Select configuration file",
            filter="JSON ( *.json );;TXT ( *.txt )")
        if not filepath:
            return

        self.config = CaverConfig.from_json(filepath) if filepath.endswith(".json") else CaverConfig.from_txt(filepath)
        # refresh window wiget from input config
        self.refresh_window_from_cfg()
        self.refresh_start_point_from_cfg()
        set_widget_value(self.ui.label_configStatus, f'Loaded from {os.path.basename(filepath)}')

    def refresh_start_point_from_cfg(self):

        if not self.config.has('starting_point_coordinates') or self.config.get('starting_point_coordinates') == '???':
            return

        coords_from_config = tuple(map(float, self.config.get('starting_point_coordinates').split(' ')))
        if not len(coords_from_config) == 3:
            notify_box("Invalid starting point coordinates in configuration file")
            return

        for i, axis, j, coord in zip(enumerate('xyz'), enumerate(coords_from_config)):
            set_widget_value(getattr(self.ui, f'doubleSpinBox_{axis}'), coord)
        notify_box(f"Starting point coordinates loaded from configuration file: {coords_from_config}")

    def configout(self, filepath: Optional[str] = None):
        filepath = filepath or getSaveFileNameWithExt(
            self.window,
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
                        details=f'Simultaneous usage of {primary_flag} parameter with {flag} parameters is not supported by plugin. Now ignoring atom.')
                    delattr(self.config, flag)

        # test include/exclude
        if self.config.has_include_exclude:
            notify_box(
                'include_ and exclude_ parameters are not supported by plugin. '
                'Please use the plugin to specify residues to be analyzed.')

        self.ensure_residue_names_to_checktable()

    def ensure_residue_names_to_checktable(self):
        if self.config.has("include_residue_names"):
            aa_from_config: List[str] = self.config.get("include_residue_names").split(" ")
            if aa_from_config:
                self.checktable_aa.check_these(aa_from_config)

    def _analysis_sel_resn(self):
        if not self.config.selection_name:
            return

        sel = cmd.get_model(self.config.selection_name)

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
