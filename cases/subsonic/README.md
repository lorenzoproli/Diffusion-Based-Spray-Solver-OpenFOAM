# Subsonic crossflow cases

This directory contains the OpenFOAM setups used for the subsonic-crossflow
validation and kernel-study activities of the thesis.

```
subsonic/
├── Zhang_case/             Original Zhang subsonic validation case (LMR = 1)
├── Zhang/                  Zhang kernel study (We500, We3500, We8500 × 4 bandwidths)
└── Lambert/
    ├── J10/, J20/                       Lambert J = 10, J = 20 reference cases
    ├── J10_Ligament/, J20_Ligament/     Lambert cases with ligament-aware libraries
    └── J10_kernelStudy/, J20_kernelStudy/  DBM bandwidth kernel study at J = 10, 20
```

All cases are setup-only. Numeric time directories, processor folders,
`postProcessing` outputs, log files and other heavy runtime artefacts are
intentionally excluded; results must be regenerated locally or on the cluster
from the provided setup.

The DBM bandwidth is configured per-case in `constant/smoothingProperties`
through the parameter `smoothBandwidth`. See the individual case READMEs for
the exact values.
