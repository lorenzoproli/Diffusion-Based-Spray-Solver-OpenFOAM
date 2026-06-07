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

These cases are setup-only. Numeric time directories, processor decomposition,
function-object outputs, parcel-level CSV dumps and run logs are intentionally
excluded; the solver fields must be regenerated locally or on the cluster from
the provided setup using the `runScripts/Allrun_automated` workflow.

The folder layout is:

```
<We>/
├── postProcessing/        Cross-bandwidth (kernel) comparison assets for this We
│   ├── figures/            Combined / pair / single-location comparison plots
│   ├── tables/             Kernel-comparison summary + reference profiles
│   ├── paraview/           Placeholder for ParaView assets
│   ├── scripts/            Python scripts used for the comparison
│   └── README.md
├── b_0_PCM/
│   ├── setup/              OpenFOAM input deck (0, constant, system, chemkin)
│   ├── runScripts/         Allrun_automated and monitor scripts
│   ├── postProcessing/     Curated per-case figures, tables, ParaView assets
│   │   ├── figures/
│   │   ├── tables/
│   │   ├── paraview/
│   │   ├── scripts/
│   │   └── README.md
│   └── README.md
├── b_Dinj/   (same layout)
├── b_3Dinj/  (same layout)
└── b_6Dinj/  (same layout)
```

The `postProcessing/` subfolders are reserved for curated outputs only — final
plots (PNG), lightweight processed summaries (CSV, JSON), ParaView
colormaps/state files/screenshots and the Python scripts that produced them.
Raw OpenFOAM runtime outputs remain excluded. See each `postProcessing/README.md`
for the per-case and per-group inventories.
