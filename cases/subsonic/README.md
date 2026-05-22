# Subsonic crossflow cases

This directory contains the OpenFOAM setups used for the subsonic-crossflow
validation and kernel-study activities of the thesis. Two reference families
are tracked:

- **Zhang** – Liquid Jet In subsonic Gas Film (LJIGF) configuration derived
  from the gas-liquid orifice-type pintle (GLOP) experiments of Zhang et al.,
  used as a benchmark for the Diffusion-Based smoothing (DBM) at three Weber
  numbers (We500, We3500, We8500);
- **Lambert** – Liquid Jet In Crossflow (LJICF) reference case from Lambert
  et al. (ILASS-Europe 2019), used at two momentum ratios (J = 10, J = 20)
  to assess the kernel-study cases of the Diffusion-Based smoothing.

```
subsonic/
├── Zhang_case/             Original Zhang subsonic validation case (LMR = 1)
├── Zhang/                  Zhang kernel study (We500, We3500, We8500 × 4 bandwidths)
└── Lambert/
    ├── J10_kernelStudy/    DBM bandwidth kernel study at J = 10
    └── J20_kernelStudy/    DBM bandwidth kernel study at J = 20
```

All cases under this directory are **setup-only**. Numeric time directories,
processor folders, function-object outputs, log files and other heavy runtime
artefacts are intentionally excluded; the solver fields must be regenerated
locally or on the cluster from the provided setup. The `polyMesh/` directory
is regenerated through `blockMesh` (plus any case-specific topology utilities
referenced in `system/`) as part of the `Allrun_automated` workflow stored
under each case's `runScripts/`.

Each kernel-study case follows the layout:

```
<case>/
├── setup/        OpenFOAM input deck (0, constant, system, chemkin)
├── runScripts/   Allrun_automated and monitor scripts
├── postProcessing/   Curated lightweight figures, tables and ParaView assets
└── README.md
```

The `postProcessing/` subfolder, when present, holds only curated lightweight
assets (final plots, processed CSV/JSON summaries, ParaView colormaps and
state files); raw OpenFOAM runtime outputs remain excluded.

The DBM bandwidth is configured per-case in `constant/smoothingProperties`
through the parameter `smoothBandwidth`. The four bandwidth variants used
throughout the kernel-study cases are:

| Case label | `smoothBandwidth` | Meaning                                |
|------------|-------------------|----------------------------------------|
| `b_0_PCM`  | `0`               | PCM reference, no DBM smoothing        |
| `b_Dinj`   | `1 × D_inj`       | DBM smoothing at one injector diameter |
| `b_3Dinj`  | `3 × D_inj`       | DBM smoothing at three D_inj           |
| `b_6Dinj`  | `6 × D_inj`       | DBM smoothing at six D_inj             |

See the per-family `README.md` (`Zhang/README.md`, `Lambert/README.md`) and
the per-case `README.md` for the exact bandwidth values in metres.
