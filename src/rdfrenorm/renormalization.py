from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from tqdm.auto import tqdm

from .data import atom_radii_from_types, read_xvg


@dataclass
class RenormConfig:
    dr_angstrom: float = 0.02
    frame_stride: int = 250
    contact_margin_angstrom: float = 1.8
    n_phi: int = 50
    n_theta: int = 25

    accurate_indices: tuple[int, ...] = tuple(range(100, 751, 25))
    approximate_indices: tuple[int, ...] = tuple(range(800, 2000, 100))

    input_pattern: str = "comparison_{residue}/comp_traj5_loopCL2_CB{cb}_{residue}{idx}.xvg"
    sampled_output_pattern: str = "comparison_{residue}/full_calc_comp_traj5_loopCL2_CB{cb}_{residue}{idx}.xvg"
    corrected_output_pattern: str = "comparison_{residue}/for_corr_full_comp_traj5_loopCL2_CB{cb}_{residue}{idx}.xvg"
    figure_output_pattern: str = "comparison_{residue}/CB{cb}_{residue}{idx}_renorm.png"

    excluded_volume_output: str = "{cluster}/renorm_excluded_volume_{residue}_full.npy"
    accessible_volume_output: str = "{cluster}/renorm_accessible_volume_{residue}_full.npy"


def shell_volume_angstrom3(radius_angstrom: float, dr_angstrom: float) -> float:
    return (4.0 / 3.0) * np.pi * ((radius_angstrom + dr_angstrom) ** 3 - radius_angstrom**3)


def make_surface_grid(n_phi: int, n_theta: int) -> tuple[np.ndarray, np.ndarray, float, float]:
    phi_values = np.linspace(-np.pi, np.pi, n_phi)
    theta_values = np.linspace(0.0, np.pi, n_theta)
    dphi = phi_values[1] - phi_values[0]
    dtheta = theta_values[1] - theta_values[0]
    return phi_values, theta_values, dphi, dtheta


def make_unit_sphere_points(phi_values: np.ndarray, theta_values: np.ndarray) -> np.ndarray:
    points = []
    for phi in phi_values:
        for theta in theta_values:
            points.append(
                [
                    np.sin(theta) * np.cos(phi),
                    np.sin(theta) * np.sin(phi),
                    np.cos(theta),
                ]
            )
    return np.asarray(points, dtype=float)


def build_shell_selection(reference_index: int, radius_angstrom: float, margin_angstrom: float) -> str:
    r_min = max(radius_angstrom - margin_angstrom, 0.0)
    r_max = radius_angstrom + margin_angstrom
    return f"sphlayer {r_min:.6f} {r_max:.6f} (index {reference_index})"


def accessible_volume_surface_method(
    reference_position: np.ndarray,
    contact_positions: np.ndarray,
    contact_radii: np.ndarray,
    shell_radius_angstrom: float,
    dr_angstrom: float,
    unit_sphere_points: np.ndarray,
    theta_values: np.ndarray,
    dphi: float,
    dtheta: float,
) -> float:
    if contact_positions.size == 0:
        return shell_volume_angstrom3(shell_radius_angstrom, dr_angstrom)

    accessible_volume = 0.0
    n_theta = len(theta_values)

    for idx, direction in enumerate(unit_sphere_points):
        theta = theta_values[idx % n_theta]
        surface_center = reference_position + shell_radius_angstrom * direction
        occupied = np.any(np.linalg.norm(contact_positions - surface_center, axis=1) <= contact_radii)

        if not occupied:
            accessible_volume += (
                (shell_radius_angstrom**2)
                * np.sin(theta)
                * dphi
                * dtheta
                * dr_angstrom
            )

    return accessible_volume


def excluded_volume_analytical_method(
    reference_position: np.ndarray,
    contact_positions: np.ndarray,
    contact_radii: np.ndarray,
    shell_radius_angstrom: float,
    dr_angstrom: float,
) -> float:
    if contact_positions.size == 0:
        return 0.0

    excluded = 0.0
    distances = np.linalg.norm(contact_positions - reference_position, axis=1)

    for d, radius in zip(distances, contact_radii, strict=True):
        if shell_radius_angstrom < d < shell_radius_angstrom + dr_angstrom:
            excluded += np.pi * dr_angstrom * radius**2 + (np.pi / 6.0) * dr_angstrom**3
        elif d - radius < shell_radius_angstrom and shell_radius_angstrom + dr_angstrom < d:
            term = radius**2 - (d - shell_radius_angstrom - dr_angstrom) ** 2
            if term > 0.0:
                excluded += np.pi * dr_angstrom * term + (np.pi / 6.0) * dr_angstrom**3

    return excluded


def compute_shell_volumes(
    universe: mda.Universe,
    residue_reference_indices: list[int],
    r_grid_nm: np.ndarray,
    radii_by_type: dict[str, float],
    config: RenormConfig,
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    sampled_indices = list(config.accurate_indices) + list(config.approximate_indices)
    n_residues = len(residue_reference_indices)
    n_radii = len(sampled_indices)

    excluded_volume = np.zeros((n_residues, n_radii), dtype=float)
    accessible_volume = np.zeros((n_residues, n_radii), dtype=float)

    phi_values, theta_values, dphi, dtheta = make_surface_grid(config.n_phi, config.n_theta)
    unit_sphere_points = make_unit_sphere_points(phi_values, theta_values)

    n_frames = len(range(0, len(universe.trajectory), config.frame_stride))
    total_steps = n_residues * n_radii * n_frames

    with tqdm(total=total_steps, desc="Computing shell volumes") as pbar:
        for residue_idx, reference_index in enumerate(residue_reference_indices):
            for sampled_col, rdf_idx in enumerate(sampled_indices):
                shell_radius_angstrom = r_grid_nm[rdf_idx] * 10.0
                shell_total = shell_volume_angstrom3(shell_radius_angstrom, config.dr_angstrom)

                excluded_accumulator = 0.0
                accessible_accumulator = 0.0

                for _ts in universe.trajectory[:: config.frame_stride]:
                    reference_position = universe.atoms[reference_index].position.copy()

                    selection = build_shell_selection(
                        reference_index=reference_index,
                        radius_angstrom=shell_radius_angstrom,
                        margin_angstrom=config.contact_margin_angstrom,
                    )
                    contacts = universe.select_atoms(selection)

                    if len(contacts) > 0:
                        mask = contacts.indices != reference_index
                        contacts = contacts[mask]

                    if len(contacts) == 0:
                        frame_excluded = 0.0
                        frame_accessible = shell_total
                    else:
                        contact_positions = contacts.positions.copy()
                        contact_radii = atom_radii_from_types(contacts.types, radii_by_type)

                        if rdf_idx in config.accurate_indices:
                            frame_accessible = accessible_volume_surface_method(
                                reference_position=reference_position,
                                contact_positions=contact_positions,
                                contact_radii=contact_radii,
                                shell_radius_angstrom=shell_radius_angstrom,
                                dr_angstrom=config.dr_angstrom,
                                unit_sphere_points=unit_sphere_points,
                                theta_values=theta_values,
                                dphi=dphi,
                                dtheta=dtheta,
                            )
                            frame_excluded = max(shell_total - frame_accessible, 0.0)
                        else:
                            frame_excluded = excluded_volume_analytical_method(
                                reference_position=reference_position,
                                contact_positions=contact_positions,
                                contact_radii=contact_radii,
                                shell_radius_angstrom=shell_radius_angstrom,
                                dr_angstrom=config.dr_angstrom,
                            )
                            frame_excluded = min(max(frame_excluded, 0.0), shell_total)
                            frame_accessible = max(shell_total - frame_excluded, 0.0)

                    excluded_accumulator += frame_excluded
                    accessible_accumulator += frame_accessible
                    pbar.update(1)

                excluded_volume[residue_idx, sampled_col] = excluded_accumulator / n_frames
                accessible_volume[residue_idx, sampled_col] = accessible_accumulator / n_frames

    return excluded_volume, accessible_volume, sampled_indices


def correct_sampled_rdf(
    raw_rdf: np.ndarray,
    sampled_indices: list[int],
    r_grid_nm: np.ndarray,
    excluded_volume_row: np.ndarray,
    dr_angstrom: float,
    eps: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray]:
    sampled_r_nm = []
    sampled_corrected = []

    for col, rdf_idx in enumerate(sampled_indices):
        radius_angstrom = r_grid_nm[rdf_idx] * 10.0
        shell_total = shell_volume_angstrom3(radius_angstrom, dr_angstrom)
        accessible = max(shell_total - excluded_volume_row[col], eps)

        corrected_value = raw_rdf[rdf_idx] * shell_total / accessible
        sampled_r_nm.append(r_grid_nm[rdf_idx])
        sampled_corrected.append(corrected_value)

    return np.asarray(sampled_r_nm), np.asarray(sampled_corrected)


def interpolate_full_correction(
    x_nm: np.ndarray,
    raw_rdf: np.ndarray,
    sampled_r_nm: np.ndarray,
    sampled_corrected: np.ndarray,
    sampled_indices: list[int],
) -> tuple[np.ndarray, np.ndarray]:
    raw_sampled = np.asarray([raw_rdf[idx] for idx in sampled_indices], dtype=float)

    correction_factor = np.divide(
        sampled_corrected,
        raw_sampled,
        out=np.ones_like(sampled_corrected),
        where=raw_sampled != 0.0,
    )

    full_factor = np.interp(
        x_nm,
        sampled_r_nm,
        correction_factor,
        left=correction_factor[0],
        right=correction_factor[-1],
    )
    corrected_full = raw_rdf * full_factor
    return full_factor, corrected_full


def save_array(path: str | Path, array: np.ndarray) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, array)


def save_text(path: str | Path, x: np.ndarray, y: np.ndarray) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.column_stack([x, y]))


def plot_correction(
    x_nm: np.ndarray,
    raw_rdf: np.ndarray,
    sampled_r_nm: np.ndarray,
    sampled_corrected: np.ndarray,
    corrected_full: np.ndarray,
    title: str,
    output_path: str | Path | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.plot(x_nm, raw_rdf, color="darkcyan", alpha=0.5, label="raw g(r)")
    ax.plot(sampled_r_nm, sampled_corrected, "o", color="green", markersize=7, label="sampled corrected g'(r)")
    ax.plot(x_nm, corrected_full, color="black", linewidth=2, label="interpolated corrected g'(r)")

    ax.set_xlabel(r"r (nm)", fontsize=14)
    ax.set_ylabel("RDF", fontsize=14)
    ax.set_xlim(0.0, 3.0)
    ax.set_title(title, fontsize=16)
    ax.legend()
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200)

    plt.show()


def default_file_patterns(
    cluster_name: str,
    residue_name: str,
) -> dict[str, str]:
    input_pattern = f"data/rdf_raw/{cluster_name}/{residue_name}/comp_{{cb}}_{{idx}}.xvg"
    sampled_pattern = f"outputs/rdf_corrected/{cluster_name}/{residue_name}/sampled_{{cb}}_{{idx}}.xvg"
    corrected_pattern = f"outputs/rdf_corrected/{cluster_name}/{residue_name}/corrected_{{cb}}_{{idx}}.xvg"
    figure_pattern = f"outputs/figures/{cluster_name}/{residue_name}/plot_{{cb}}_{{idx}}.png"

    return {
        "input_pattern": input_pattern,
        "sampled_output_pattern": sampled_pattern,
        "corrected_output_pattern": corrected_pattern,
        "figure_output_pattern": figure_pattern,
        "excluded_volume_output": f"outputs/volumes/{cluster_name}/excluded_{residue_name}.npy",
        "accessible_volume_output": f"outputs/volumes/{cluster_name}/accessible_{residue_name}.npy",
    }


def process_cluster_residue(
    universe: mda.Universe,
    cluster_name: str,
    residue_name: str,
    residue_reference_indices: list[int],
    r_grid_nm: np.ndarray,
    radii_by_type: dict[str, float],
    config: RenormConfig,
) -> None:
    excluded_volume, accessible_volume, sampled_indices = compute_shell_volumes(
        universe=universe,
        residue_reference_indices=residue_reference_indices,
        r_grid_nm=r_grid_nm,
        radii_by_type=radii_by_type,
        config=config,
    )

    save_array(config.excluded_volume_output, excluded_volume)
    save_array(config.accessible_volume_output, accessible_volume)

    residue_counter = 0
    max_idx = len(residue_reference_indices)

    # Historical CB ranges differ slightly by cluster family in your notebooks,
    # but the generic script keeps the same 3..10 loop unless the user edits patterns.
    for cb in range(3, 11):
        for idx in range(1, max_idx + 1):
            input_path = Path(config.input_pattern.format(cb=cb, idx=idx, residue=residue_name))
            if not input_path.is_file():
                continue

            if residue_counter >= len(residue_reference_indices):
                raise IndexError(
                    f"More RDF files were found than reference indices for {cluster_name} {residue_name}."
                )

            x_nm, raw_rdf = read_xvg(input_path)

            sampled_r_nm, sampled_corrected = correct_sampled_rdf(
                raw_rdf=raw_rdf,
                sampled_indices=sampled_indices,
                r_grid_nm=r_grid_nm,
                excluded_volume_row=excluded_volume[residue_counter],
                dr_angstrom=config.dr_angstrom,
            )

            _full_factor, corrected_full = interpolate_full_correction(
                x_nm=x_nm,
                raw_rdf=raw_rdf,
                sampled_r_nm=sampled_r_nm,
                sampled_corrected=sampled_corrected,
                sampled_indices=sampled_indices,
            )

            sampled_output = Path(config.sampled_output_pattern.format(cb=cb, idx=idx, residue=residue_name))
            corrected_output = Path(config.corrected_output_pattern.format(cb=cb, idx=idx, residue=residue_name))
            figure_output = Path(config.figure_output_pattern.format(cb=cb, idx=idx, residue=residue_name))

            save_text(sampled_output, sampled_r_nm, sampled_corrected)
            save_text(corrected_output, x_nm, corrected_full)

            plot_correction(
                x_nm=x_nm,
                raw_rdf=raw_rdf,
                sampled_r_nm=sampled_r_nm,
                sampled_corrected=sampled_corrected,
                corrected_full=corrected_full,
                title=f"{cluster_name} CB{cb} {residue_name}{idx}",
                output_path=figure_output,
            )

            residue_counter += 1
