# Lambert – J20_kernelStudy – b_0_PCM

Lambert (ILASS 2019) reference case with momentum ratio J20, used as a kernel-study point for the Diffusion-Based smoothing of the Lagrangian source terms.

## Kernel setting

- Bandwidth label: `b_0_PCM` (0)
- `constant/smoothingProperties` → `smoothBandwidth = 0` [m]
- Interpretation: PCM reference (no DBM smoothing; cellPoint projection)

## Source folder

Imported from:

`~/runs/subsonic/Lambert_case_J20_kernelStudy_cluster/b_0_PCM_cell`

The original time directories, processor decomposition, `postProcessing` outputs
and run logs are intentionally excluded. Only the clean setup needed to rerun
the case is retained.

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
