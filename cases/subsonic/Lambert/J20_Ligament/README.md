# Zhang subsonic crossflow case

This folder contains the subsonic validation case derived from the reference configuration reported by Wu et al., restricted to the condition with liquid-to-gas momentum ratio:

LMR = 1

The purpose of this case is not to reproduce the full experimental and numerical campaign, but to use the selected LMR = 1 condition as a validation benchmark for the thesis spray solver.

## Validation objective

This case is used to compare:

- the numerical results obtained with the present OpenFOAM implementation
- the numerical results reported in the reference study
- the available experimental measurements

The comparison is focused on the validation of the Lagrangian spray modelling strategy in subsonic crossflow conditions, with particular attention to atomization behaviour and droplet statistics.

## Physical configuration

- Gas phase: Nitrogen (N₂)
- Liquid phase: Water
- Injection in subsonic crossflow
- Validation condition: LMR = 1

The case is intended as a controlled benchmark for assessing the predictive capability of the implemented custom breakup and drag models.

## Folder structure

- `setup/`: complete OpenFOAM case
  - `0/`: initial and boundary conditions
  - `constant/`: thermophysical properties, spray setup, geometry, and mesh
  - `system/`: numerical setup and control dictionaries
- `results/`: reserved for lightweight tracked outputs when needed

## Custom models

This case relies on the custom user libraries implemented in:

- `../../../src/breakup/Madabhushi/`
- `../../../src/drag/Madabhushi/`

These models are used to improve the representation of spray breakup and droplet drag under crossflow conditions.

## Mesh and reproducibility

The final mesh is included in:

`setup/constant/polyMesh/`

This choice ensures reproducibility of the validated configuration currently used in the thesis work.

Although a `blockMeshDict` is available, the full mesh generation pipeline is not yet completely reconstructed in this repository through all intermediate dictionaries and utilities. For this reason, the final mesh is retained in the setup.

## Post-processing

Case-specific post-processing scripts are stored in:

`../../../postProcessing/subsonic/Zhang_case/`

These scripts are used to extract validation quantities such as:

- droplet size distribution
- Sauter Mean Diameter (D32)
- velocity-related spray statistics

## Purpose of the case

This case is used to:

- validate the custom Lagrangian spray solver at LMR = 1
- compare the present numerical results with literature numerical data
- compare the present numerical results with experimental reference data
- assess the behaviour of the implemented breakup and drag models in a controlled subsonic crossflow benchmark
