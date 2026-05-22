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

Each case folder follows the layout:

```
<We>/<band>/
├── setup/        OpenFOAM input deck (0, constant, system, chemkin)
├── runScripts/   Allrun_automated and monitor scripts
├── postProcessing/   Curated lightweight figures, tables and ParaView assets
└── README.md
```

The `postProcessing/` subfolder is reserved for curated outputs only — final
plots (PNG), lightweight processed summaries (CSV, JSON) and ParaView
colormaps/state files/screenshots. Raw OpenFOAM runtime outputs remain
excluded. See each `<band>/postProcessing/README.md` for the per-case
inventory.

See the parent `cases/subsonic/Zhang_case/` for the original LMR = 1
validation case from which the kernel study was derived.
