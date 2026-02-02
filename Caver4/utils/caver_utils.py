import os

from ..caver_pymol import ROOT_LOGGER

logging = ROOT_LOGGER.getChild("Utils")


IGNORED_STRUCTURES = [r"^origins$", r"_origins$", r"_v_origins$", r"_t\d\d\d_\d$"]


THE_20s = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]


def find_centrial_pdb(out_home: str, run_id: int) -> str:
    """
    Find the central PDB file in the specified output directory for a given run ID.

    Args:
        out_home (str): The root output directory path
        run_id (int): The run identifier used to construct the subdirectory path

    Returns:
        str: The full path to the central PDB file

    Raises:
        RuntimeError: If no PDB files are found or more than one PDB file is found in the directory
    """

    # Construct the path to the data directory containing PDB files
    pdb_path = os.path.join(out_home, str(run_id), "data")

    # Get all .pdb files in the directory excluding v_origins.pdb and origins.pdb
    pdb_files = [x for x in os.listdir(pdb_path) if x.endswith(".pdb") and x != "v_origins.pdb" and x != "origins.pdb"]

    logging.debug(f"pdb_files: {pdb_files}")

    if not pdb_files:
        raise RuntimeError(f"No PDB files found in the output directory: {pdb_path}")

    # Verify that only one PDB file exists
    if len(pdb_files) > 1:
        raise RuntimeError(f"More than one PDB file found in the output directory: {pdb_path}: {pdb_files}")
    input_pdb = pdb_files[0]
    logging.debug(f"input_pdb: {input_pdb}")
    return os.path.join(pdb_path, input_pdb)
