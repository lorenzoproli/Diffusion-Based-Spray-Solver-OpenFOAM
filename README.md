# Thesis Lagrangian Spray Solver – Crossflow Atomization

This repository contains the development and validation of a custom Lagrangian spray solver based on OpenFOAM, carried out within a Master's thesis in Aerospace Engineering.

The project focuses on the numerical modelling of liquid atomization in gas crossflow conditions, with particular attention to breakup dynamics, droplet aerodynamic behaviour, and validation against reference literature and experimental data.

---

## Project objectives

The main objectives of the project are:

- develop custom breakup and drag models for Lagrangian spray simulations
- improve the physical representation of crossflow atomization
- validate the numerical framework against literature and experimental references
- perform detailed post-processing of spray characteristics such as droplet size distribution, D32, and dispersion behaviour

---

## Key features

- Custom multi-stage breakup model:
  - liquid column behaviour
  - Kelvin–Helmholtz wave stripping
  - primary breakup transition
  - secondary breakup based on Pilch–Erdman-type behaviour

- Custom drag model:
  - blob/column regime
  - spherical droplet regime
  - deformation-based drag during breakup

- Crossflow atomization framework:
  - gas–liquid interaction in crossflow conditions
  - reproducible OpenFOAM case organization
  - dedicated setup, source-code and post-processing structure
  - validation-oriented workflow against literature and experimental references

- Structured repository for:
  - solver development
  - case management
  - post-processing
  - documentation and reference material

---

## Repository structure

```text
Thesis-Lagrangian-Solver/
│
├── src/
│   ├── breakup/           # Custom breakup models
│   ├── drag/              # Custom drag models
│   └── evaporation/       # Future extensions
│
├── cases/
│   ├── subsonic/
│   │   └── Wu_case/       # Validation case
│   ├── supersonic/
│   └── evaporation/
│
├── postProcessing/
│   ├── common/
│   ├── subsonic/
│   ├── supersonic/
│   └── evaporation/
│
├── docs/
│   ├── papers/
│   ├── notes/
│   └── figures/
│
└── README.md
