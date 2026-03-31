# RDF Volume Renormalization for Protein–Solvent Systems

This repository implements a geometrical renormalization of radial distribution functions (RDFs) to correct for excluded volume effects in heterogeneous environments such as proteins in solvent.

The method was developed in the context of the LRP1 receptor analysis and is designed to be reusable and extensible for similar molecular systems.

## Motivation

Standard RDFs assume a homogeneous, fully accessible volume.

However, near macromolecules (e.g., proteins), a fraction of space is physically inaccessible due to steric exclusion.

This leads to systematic bias:

- Underestimation of local densities
- Incorrect normalization
- Distortion of structural interpretation

This project corrects RDFs by explicitly accounting for excluded volume.

## Method Overview

The corrected RDF is defined as:

g'(r) = g(r) * V_shell(r) / V_accessible(r)

where:

- V_shell(r) = 4πr²Δr
- V_accessible(r) = V_shell(r) − V_excluded(r)

### Two complementary approaches

#### 1. Surface sampling (short range)

- Sample directions around the reference residue
- Detect intersections with protein atoms (vdW radii)
- Estimate excluded solid angle

Accurate but computationally expensive.

#### 2. Analytical approximation (long range)

- Sum excluded contributions from nearby atoms
- Based on atom distances and vdW radii

Efficient but slightly overestimates volume because overlaps are neglected.

### Hybrid strategy

- Short distances → surface sampling
- Long distances → analytical method
- Smooth interpolation between regimes

This provides a good balance between accuracy and computational cost.

## Repository Structure

rdf-volume-renormalization/
├── data/
│   ├── config/
│   │   └── lrp1_clusters.json
│   └── rdf_raw/
│       └── cl2/
│           └── ASP/
├── systems/
│   └── cl2/
├── outputs/
├── scripts/
│   └── run_cluster_residue.py
└── src/
    └── rdfrenorm/

## Installation

Create a virtual environment and install:

pip install -r requirements.txt
pip install -e .

## Usage

Run the renormalization for a given cluster and residue family:

python scripts/run_cluster_residue.py --cluster cl2 --residue ASP

For a lightweight test:

python scripts/run_cluster_residue.py --cluster cl2 --residue ASP --light

## Input Requirements

To run the pipeline you need:

1. MD system
- .gro topology file
- .xtc trajectory file

2. RDF data
- Raw RDF .xvg files for each residue instance

3. Configuration file
- data/config/lrp1_clusters.json

This file contains:
- residue lists (PDB numbering)
- PDB to GRO offset
- reference atoms
- system paths

## Outputs

The code generates:

- Corrected RDFs
- Sampled RDFs
- Plots
- Excluded and accessible volume arrays

All outputs are written to:

outputs/

## Important Note on Provided Data

This repository includes:

- Example trajectory (cluster 2)
- Example RDF files for ASP residues only

These are provided solely to verify that the code runs correctly.

## Access to Full Dataset

The complete dataset (all clusters and residue types) is not included in this repository.

To access the full data, please contact:

Giuseppe Battaglia, Institute for Bioengineering of Catalunya

## Limitations

- Analytical method neglects overlap between atoms
- Slight overestimation of excluded volume at large distances
- Assumes spherical shell geometry
- Performance depends on trajectory size

## Validation

A correct implementation should satisfy:

g'(r) → 1 as r → ∞

This is a key sanity check for the renormalization.

## Scientific Context

This method relates to:

- Confined fluids
- Protein hydration layers
- Density normalization in heterogeneous systems

## Future Improvements

Possible extensions:

- Overlap-corrected excluded volume
- Grid-based masking approaches
- Parallelization
- Automatic RDF discovery
- Integration with MDAnalysis pipelines

## Author

Gian Marco Tuveri

PhD in Computational Biophysics
Universitat de Barcelona

## License

MIT license

## Notes for Future Users

This code was developed as part of a research workflow.

It is intended as:
- a reproducible implementation
- a starting point for further development
- a reference for RDF correction methods

If you extend this work, please document changes clearly.
