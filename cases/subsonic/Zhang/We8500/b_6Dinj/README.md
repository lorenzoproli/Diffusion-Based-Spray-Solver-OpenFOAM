# Zhang – We8500 – b_6Dinj

Subsonic crossflow validation case from the Zhang configuration (We8500), used as a kernel-study point for the Diffusion-Based smoothing of the Lagrangian source terms.

## Kernel setting

- Bandwidth label: `b_6Dinj` (6 × D_inj)
- `constant/smoothingProperties` → `smoothBandwidth = 0.0096` [m]
- Interpretation: Diffusion-Based smoothing with bandwidth equal to six injector diameters


## Folder structure

- `setup/0/`: initial and boundary conditions (the `U` field is taken from the
  original `U_backup/U` of the source case, before any overwrite by
  `potentialFoam`)
- `setup/constant/`: thermophysical properties, spray cloud setup, smoothing
  configuration, gravity, turbulence, radiation, species
- `setup/chemkin/`: CHEMKIN-format thermodynamic/transport data
- `setup/system/`: numerical setup, controlDict, fvSchemes, fvSolution,
  blockMeshDict, decomposeParDict and case-specific dictionaries
- `runScripts/`: `Allrun_automated` and `monitor_breakup.sh` for cluster execution

## Mesh

The `polyMesh/` directory is not stored in the repository. The mesh is
regenerated through `blockMesh` (and any case-specific topology utilities
referenced in `system/`) as part of the `Allrun_automated` workflow.

## Results

Simulation results (numeric time directories, processor folders, `postProcessing`
and log files) are intentionally excluded. They should be regenerated locally or
on the cluster from the provided setup.
