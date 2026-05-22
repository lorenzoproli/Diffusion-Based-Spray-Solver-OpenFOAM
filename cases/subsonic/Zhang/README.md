# Zhang subsonic crossflow – kernel study

This directory contains the OpenFOAM setups for the kernel study performed on
the Zhang subsonic crossflow validation configuration. The study spans three
Weber-number conditions and four DBM bandwidth values:

```
Zhang/
├── We500/
│   ├── b_0_PCM/   smoothBandwidth = 0           (PCM reference)
│   ├── b_Dinj/    smoothBandwidth = 1 × D_inj
│   ├── b_3Dinj/   smoothBandwidth = 3 × D_inj
│   └── b_6Dinj/   smoothBandwidth = 6 × D_inj
├── We3500/        (same four bandwidths)
└── We8500/        (same four bandwidths)
```

The injector diameter for these cases is `D_inj = 1.6 mm`, so:

- `b_Dinj`  ↔ `smoothBandwidth = 0.0016 m`
- `b_3Dinj` ↔ `smoothBandwidth = 0.0048 m`
- `b_6Dinj` ↔ `smoothBandwidth = 0.0096 m`

These cases are setup-only. Time directories, processor decomposition,
`postProcessing` outputs and run logs are intentionally excluded; they must be
regenerated locally or on the cluster from the provided setup using the
`runScripts/Allrun_automated` workflow.

See the parent `cases/subsonic/Zhang_case/` for the original LMR = 1
validation case from which the kernel study was derived.
