from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

StartMode = Literal["atoms", "residues", "coords"]


@dataclass(frozen=True)
class WorkflowScenario:
    """
    Describes a PyMOL workflow that can be executed by the GUI test worker.
    """

    slug: str
    structure_relpath: str
    object_name: str
    start_mode: StartMode
    start_value: tuple[float, float, float] | str
    md_state_range: tuple[int, int] | None = None

    def structure_path(self, data_root: Path) -> Path:
        """
        Resolve the workflow's structure file relative to the bundled test data.
        """
        return data_root / self.structure_relpath

    @property
    def snapshot_prefix(self) -> str:
        return self.slug


STATIC_WORKFLOW = WorkflowScenario(
    slug="static_analysis",
    structure_relpath="pdb/1AKD.pdb",
    object_name="caver_static",
    start_mode="coords",
    start_value=(17.012, 24.139, 7.790),
)


DYNAMIC_WORKFLOW = WorkflowScenario(
    slug="dynamic_analysis",
    structure_relpath="md_snapshots/caver_md.snapshots.pze",
    object_name="caver_test_md",
    start_mode="atoms",
    start_value="578 1609 3258",
    md_state_range=(1, 50),
)


__all__ = ["WorkflowScenario", "STATIC_WORKFLOW", "DYNAMIC_WORKFLOW"]
