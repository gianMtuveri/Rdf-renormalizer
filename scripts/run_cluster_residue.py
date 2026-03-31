import argparse
from pathlib import Path

from rdfrenorm import (
    ProjectPaths,
    build_reference_index_map,
    load_cluster_config,
    load_r_grid_from_example,
    load_universe,
)
from rdfrenorm.renormalization import RenormConfig, process_cluster_residue, default_file_patterns


def parse_args():
    parser = argparse.ArgumentParser(description="Run RDF renormalization")

    parser.add_argument("--cluster", required=True, help="Cluster name (cl2, cl3, cl4)")
    parser.add_argument("--residue", required=True, help="Residue family (ASP, CYS, ...)")

    parser.add_argument("--light", action="store_true", help="Run lightweight test")

    return parser.parse_args()


def main():
    args = parse_args()

    project_paths = ProjectPaths(config_json=Path("data/lrp1_clusters.json"))
    cluster_cfg = load_cluster_config(args.cluster, project_paths)

    print(f"Loading universe for {args.cluster}...")
    u = load_universe(cluster_cfg)

    print("Building residue index map...")
    ref_map = build_reference_index_map(u, cluster_cfg)

    print("Loading RDF grid...")
    r_grid = load_r_grid_from_example(cluster_cfg.rdf_example)

    patterns = default_file_patterns(args.cluster, args.residue)

    # -------- LIGHT TEST MODE --------
    if args.light:
        print("Running LIGHT configuration")

        renorm_cfg = RenormConfig(
            dr_angstrom=0.02,
            frame_stride=1000,
            contact_margin_angstrom=1.8,
            n_phi=20,
            n_theta=10,
            accurate_indices=(100, 125),
            approximate_indices=(800,),
            input_pattern=patterns["input_pattern"],
            sampled_output_pattern=patterns["sampled_output_pattern"],
            corrected_output_pattern=patterns["corrected_output_pattern"],
            figure_output_pattern=patterns["figure_output_pattern"],
            excluded_volume_output=patterns["excluded_volume_output"],
            accessible_volume_output=patterns["accessible_volume_output"],
        )

    # -------- FULL RUN --------
    else:
        print("Running FULL configuration")

        renorm_cfg = RenormConfig(
            dr_angstrom=0.02,
            frame_stride=250,
            contact_margin_angstrom=1.8,
            n_phi=50,
            n_theta=25,
            accurate_indices=tuple(range(100, 751, 25)),
            approximate_indices=tuple(range(800, 2000, 100)),
            input_pattern=patterns["input_pattern"],
            sampled_output_pattern=patterns["sampled_output_pattern"],
            corrected_output_pattern=patterns["corrected_output_pattern"],
            figure_output_pattern=patterns["figure_output_pattern"],
            excluded_volume_output=patterns["excluded_volume_output"],
            accessible_volume_output=patterns["accessible_volume_output"],
        )

    process_cluster_residue(
        universe=u,
        cluster_name=args.cluster,
        residue_name=args.residue,
        residue_reference_indices=ref_map[args.residue].tolist(),
        r_grid_nm=r_grid,
        radii_by_type=cluster_cfg.vdw_radii_angstrom,
        config=renorm_cfg,
    )

    print("Done.")


if __name__ == "__main__":
    main()
