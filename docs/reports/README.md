# Thesis Report

This folder contains the final compiled report (PDF) that documents the theory,
implementation and validation of the custom OpenFOAM Lagrangian spray solver
hosted in this repository.

## Master's thesis

**Title.** A Diffusion-Based Lagrangian Spray Solver in OpenFOAM for Liquid Jet
Atomization in Subsonic Crossflow

**Author.** Lorenzo Proli — Politecnico di Torino, M.Sc. in Aerospace
Engineering, 2026.

**File.** [`MSc_Thesis_Proli_2026.pdf`](MSc_Thesis_Proli_2026.pdf)

### Abstract

A custom OpenFOAM Lagrangian spray solver implementing the Madabhushi–Lambert
breakup model for liquid jets in crossflow (LJICF) and gas-film (LJIGF/GLOP)
configurations. The solver couples the breakup model with a Diffusion-Based
Method (DBM) for momentum source-term smoothing, eliminating the velocity-locking
artefact of the standard Particle Centroid Method (PCM) in dense spray regions.
The complete model chain — KH/Wave surface stripping, Madabhushi column breakup,
the Lambert ligament correction, Pilch–Erdman secondary breakup and a
deformation-dependent drag law — is assessed against the Lambert LJICF benchmark
(J = 10 and J = 20) and the Zhang LJIGF/GLOP cases at two gas-film Weber numbers.

### Contents at a glance

1. Introduction — motivation, objectives and contributions
2. Physical background and Lagrangian model formulation
3. OpenFOAM implementation (`sprayFoamDBM`)
4. Lambert LJICF validation: J10 sensitivity and J20 follow-up
5. Zhang LJIGF/GLOP validation and DBM necessity
6. Conclusions, limitations and future developments

## Citation

If you use the solver or refer to this work, please cite the archived software
record (Zenodo DOI
[10.5281/zenodo.20666067](https://doi.org/10.5281/zenodo.20666067)); see
[`CITATION.cff`](../../CITATION.cff) in the repository root.

## Related material

- Solver source code and documentation: top-level
  [`README.md`](../../README.md).
- Curated bibliography of the supporting literature:
  [`docs/papers/README.md`](../papers/README.md).
