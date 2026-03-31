from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import MDAnalysis as mda
import numpy as np

from .config import ClusterConfig, ProjectPaths


def load_json(path: str | Path) -> dict:
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def load_cluster_config(cluster_name: str, paths: ProjectPaths = ProjectPaths()) -> ClusterConfig:
    data = load_json(paths.config_json)

    common = data["common"]
    clusters = data["clusters"]

    if cluster_name not in clusters:
        raise KeyError(f"Cluster '{cluster_name}' not found in {paths.config_json}")

    cluster = clusters[cluster_name]

    reference_atoms = common["default_reference_atoms"].copy()
    reference_atoms.update(cluster.get("reference_atoms", {}))

    return ClusterConfig(
        name=cluster_name,
        topology=Path(cluster["topology"]),
        trajectory=Path(cluster["trajectory"]),
        pdb_to_gro_offset=int(cluster["pdb_to_gro_offset"]),
        rdf_example=Path(cluster["rdf_example"]),
        reference_atoms=reference_atoms,
        residues_pdb=cluster["residues_pdb"],
        vdw_radii_angstrom=common["vdw_radii_angstrom"],
    )


def load_universe(cluster_cfg: ClusterConfig) -> mda.Universe:
    return mda.Universe(str(cluster_cfg.topology), str(cluster_cfg.trajectory))


def pdb_to_gro_resids(pdb_resids: Iterable[int], offset: int) -> np.ndarray:
    return np.asarray(list(pdb_resids), dtype=int) - int(offset)


def build_residue_maps(cluster_cfg: ClusterConfig) -> dict[str, np.ndarray]:
    gro_map: dict[str, np.ndarray] = {}
    for residue_name, pdb_resids in cluster_cfg.residues_pdb.items():
        gro_map[residue_name] = pdb_to_gro_resids(
            pdb_resids,
            offset=cluster_cfg.pdb_to_gro_offset,
        )
    return gro_map


def get_reference_atom_indices(
    universe: mda.Universe,
    gro_resids: Iterable[int],
    atom_name: str,
) -> np.ndarray:
    indices: list[int] = []

    for resid in gro_resids:
        selection = universe.select_atoms(f"name {atom_name} and resid {int(resid)}")
        if len(selection) == 0:
            raise ValueError(
                f"No atom found for selection: name {atom_name} and resid {int(resid)}"
            )
        if len(selection) > 1:
            raise ValueError(
                f"More than one atom found for selection: name {atom_name} and resid {int(resid)}"
            )
        indices.append(int(selection[0].index))

    return np.asarray(indices, dtype=int)


def build_reference_index_map(
    universe: mda.Universe,
    cluster_cfg: ClusterConfig,
) -> dict[str, np.ndarray]:
    gro_residue_map = build_residue_maps(cluster_cfg)
    index_map: dict[str, np.ndarray] = {}

    for residue_name, gro_resids in gro_residue_map.items():
        atom_name = cluster_cfg.reference_atoms[residue_name]
        index_map[residue_name] = get_reference_atom_indices(
            universe=universe,
            gro_resids=gro_resids,
            atom_name=atom_name,
        )

    return index_map


def read_xvg(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    x, y = np.loadtxt(path, comments=("@", "#"), unpack=True)
    return x, y


def load_r_grid_from_example(example_xvg: str | Path) -> np.ndarray:
    r, _ = read_xvg(example_xvg)
    return r


def atom_radii_from_types(atom_types: Iterable[str], vdw_radii: dict[str, float]) -> np.ndarray:
    radii = []
    for atom_type in atom_types:
        if atom_type not in vdw_radii:
            raise KeyError(f"Missing van der Waals radius for atom type: {atom_type}")
        radii.append(vdw_radii[atom_type])
    return np.asarray(radii, dtype=float)
