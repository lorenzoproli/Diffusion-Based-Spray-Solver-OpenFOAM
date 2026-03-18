# Thesis Lagrangian Spray Solver – Crossflow Atomization

This repository contains the development and validation of a custom Lagrangian spray solver based on OpenFOAM, carried out within a Master's thesis in Aerospace Engineering.

The focus of the work is the modelling of liquid atomization in subsonic crossflow conditions, with particular attention to breakup dynamics and droplet aerodynamic behaviour.

---

## Project objectives

The main objectives of the project are:

- Develop custom breakup and drag models for Lagrangian spray simulations
- Improve the physical representation of crossflow atomization
- Validate the numerical framework against literature and experimental data
- Perform detailed post-processing of spray characteristics (e.g. D32, dispersion)

---

## Key features

- Custom multi-stage breakup model:
  - Liquid column behaviour
  - Kelvin–Helmholtz wave stripping
  - Primary breakup transition
  - Secondary breakup (Pilch–Erdman regime)

- Custom drag model:
  - Blob/column regime
  - Spherical droplet regime
  - Deformation-based drag during breakup

- Crossflow atomization setup:
  - Gas–liquid interaction in subsonic regime
  - Validation at LMR = 1

- Structured repository for:
  - solver development
  - case management
  - post-processing

---

## Repository structure
