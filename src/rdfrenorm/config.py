from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class ProjectPaths:
    config_json: Path = Path("data/lrp1_clusters.json")


@dataclass
class ClusterConfig:
    name: str
    topology: Path
    trajectory: Path
    pdb_to_gro_offset: int
    rdf_example: Path
    reference_atoms: dict[str, str]
    residues_pdb: dict[str, list[int]]
    vdw_radii_angstrom: dict[str, float]
