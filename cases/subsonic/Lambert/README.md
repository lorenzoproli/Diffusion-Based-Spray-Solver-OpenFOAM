# Lambert subsonic crossflow cases

This directory groups the Lambert (ILASS 2019) reference cases used as
validation points for the custom Lagrangian spray solver.

```
Lambert/
├── J10/                  Reference case, momentum ratio J = 10
├── J20/                  Reference case, momentum ratio J = 20
├── J10_Ligament/         J10 with updated ligament-aware libraries
├── J20_Ligament/         J20 with updated ligament-aware libraries
├── J10_kernelStudy/      DBM kernel study at J = 10 (b_0_PCM, b_Dinj, b_3Dinj, b_6Dinj)
└── J20_kernelStudy/      DBM kernel study at J = 20 (b_0_PCM, b_Dinj, b_3Dinj, b_6Dinj)
```

The `*_kernelStudy/` folders contain the OpenFOAM setups used for the
parametric study of the Diffusion-Based smoothing bandwidth. Each kernel-study
group shares the same physical configuration and differs only in
`constant/smoothingProperties → smoothBandwidth`.

These cases are setup-only. Time directories, processor decomposition,
`postProcessing` outputs and run logs are intentionally excluded.
