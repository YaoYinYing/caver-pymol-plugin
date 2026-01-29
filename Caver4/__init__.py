from pymol.plugins import addmenuitemqt

from .caver_pymol import CaverPyMOL


def __init_plugin__(app=None):
    """
    Add an entry to the PyMOL "Plugin" menu
    """

    plugin = CaverPyMOL()
    addmenuitemqt("Caver NG", plugin.run_plugin_gui)
