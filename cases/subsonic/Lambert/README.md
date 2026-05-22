# Lambert subsonic crossflow cases

This directory groups the Lambert (ILASS 2019) Liquid Jet In Crossflow (LJICF)
reference cases used as validation points for the custom Lagrangian spray
solver. The repository snapshot only tracks the DBM kernel-study setups:

```
Lambert/
├── J10_kernelStudy/      DBM kernel study at J = 10 (b_0_PCM, b_Dinj, b_3Dinj, b_6Dinj)
└── J20_kernelStudy/      DBM kernel study at J = 20 (b_0_PCM, b_Dinj, b_3Dinj, b_6Dinj)
```

The standalone reference folders that previously existed in this directory
(`J10/`, `J20/`, `J10_Ligament/`, `J20_Ligament/`) have intentionally been
removed from this snapshot and are not part of the currently tracked case set;
the kernel-study folders supersede them for the validation of the
Diffusion-Based smoothing.

The `*_kernelStudy/` folders contain the OpenFOAM setups used for the
parametric study of the Diffusion-Based Method (DBM) smoothing bandwidth.
Each kernel-study group shares the same physical configuration and differs
only in `constant/smoothingProperties → smoothBandwidth`. The four bandwidth
variants are:

| Case label | `smoothBandwidth` | Meaning                                |
|------------|-------------------|----------------------------------------|
| `b_0_PCM`  | `0`               | Particle Centroid Method (PCM) reference, no DBM smoothing |
| `b_Dinj`   | `1 × D_inj`       | DBM smoothing bandwidth = 1 injector diameter |
| `b_3Dinj`  | `3 × D_inj`       | DBM smoothing bandwidth = 3 injector diameters |
| `b_6Dinj`  | `6 × D_inj`       | DBM smoothing bandwidth = 6 injector diameters |

The `b_0_PCM` folder is the normalized in-repository name of the PCM/cellPoint
source case and provides the PCM baseline against which the DBM bandwidths are
compared.

All cases are setup-only. Numeric time directories, processor decomposition,
function-object outputs, parcel-level CSV dumps and run logs are intentionally
excluded; the `polyMesh/` directory is regenerated locally through `blockMesh`
(plus any case-specific topology utilities) via the `Allrun_automated` workflow
under `runScripts/`.

The folder layout is:

```
<J*_kernelStudy>/
├── postProcessing/        Cross-bandwidth (kernel) comparison assets for this J
│   ├── figures/            Combined / pair / single-location comparison plots
│   ├── tables/             Kernel-comparison summary + exp/Fluent reference profiles
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
