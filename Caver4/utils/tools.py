

import os
import webbrowser

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

from .ui_tape import notify_box
def open_doc_pdf():
    """
    Open doc pdf
    """
    doc_file = os.path.join(THIS_DIR, '..',"config", "caver_userguide.pdf")

    webbrowser.open(f"file://{doc_file}")


def cite_info():
    """
    Cite info
    """
    citation_file = os.path.join(THIS_DIR, '..',"bin", "citation.txt")
    with open(citation_file) as f:
        citation_text = f.read()

    notify_box("Thank you for using Caver, please cite the following paper.", details=citation_text)

