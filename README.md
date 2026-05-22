# Thesis Lagrangian Spray Solver – Crossflow Atomization

This repository contains the development and validation of a custom Lagrangian spray solver based on OpenFOAM, carried out within a Master's thesis in Aerospace Engineering at Politecnico di Torino.

The project focuses on the numerical modelling of liquid atomization in gas crossflow conditions, with particular attention to breakup dynamics, droplet aerodynamic behaviour, Lagrangian–Eulerian momentum coupling, and validation against reference literature and experimental data.

---

## Project objectives

- Develop custom breakup and drag models for Lagrangian spray simulations
- Improve the physical representation of liquid jet atomization in gas crossflow
- Correct the velocity locking artefact of OpenFOAM's Particle Centroid Method (PCM) via Gaussian kernel source-term smoothing
- Validate the numerical framework against literature and experimental references
- Perform detailed post-processing of spray characteristics such as droplet size distribution, SMD/D32, droplet velocity, and dispersion behaviour
- Provide reproducible OpenFOAM case setups for subsonic crossflow atomization benchmarks

---

## Key features

**Custom multi-stage breakup model (Madabhushi + Lambert):**
- liquid column/blob behaviour at injection
- Kelvin–Helmholtz wave stripping
- primary breakup transition with ligament mechanism based on Lambert (2019)
- secondary breakup based on Pilch–Erdman-type catastrophic breakup
- state tracking for core, Wave-generated children, and post-catastrophic children
- support for multiple child parcels and parent-fragment reclassification

**Custom drag model (MadabhushiDragForce):**
- blob/column regime drag
- spherical droplet regime drag
- deformation-based drag during Pilch–Erdman breakup transition
- sphere-to-disc transition with frontal-area correction for PE-active droplets

**Modified spray parcel infrastructure:**
- extended `SprayParcel.C` for parent/child state handling
- extended `BreakupModel.H` interface for child state initialization, pending children, and parent user-state update
- structured debug support for breakup diagnostics

**Parcel plane sampling:**
- custom `ParcelPlaneSampler` cloud function
- crossing-event sampling on user-defined measurement planes
- output of parcel position, velocity, diameter, `nParticle`, `user`, `ms`, `tc`, `KHindex`, `yState`, and `yDotState`
- designed for Lambert-style post-processing of SMD/D32 and droplet velocity profiles at fixed downstream locations

**Gaussian kernel source-term smoothing (sprayFoamDBM):**
- replaces OpenFOAM's PCM point-source projection with a spatially distributed Gaussian kernel
- mathematically equivalent to ANSYS Fluent's "Enable Node Based Averaging + Gaussian kernel" option
- eliminates velocity locking artefact that suppresses jet deflection in crossflow
- implemented via implicit screened-Poisson solve on `UTrans` and `UCoeff` fields
- tunable via `constant/smoothingProperties` through the bandwidth parameter `b`

**Crossflow atomization framework:**
- gas–liquid interaction in subsonic crossflow conditions
- reproducible OpenFOAM case organization
- dedicated setup, source-code, run-script, and post-processing structure
- validation-oriented workflow against Zhang-type subsonic crossflow benchmarks and Lambert J = 10 / J = 20 kernel-study cases, with experimental comparisons from the literature

---

## Repository structure

```text
Thesis-Lagrangian-Solver/
│
├── src/
│   ├── breakup/                         # Custom breakup models
│   │   └── Madabhushi/                  # Madabhushi + Lambert breakup model
│   │       ├── Madabhushi.H/.C
│   │       └── Make/
│   │
│   ├── drag/                            # Custom drag models
│   │   └── Madabhushi/                  # MadabhushiDragForce model
│   │       ├── MadabhushiDragForce.H/.C
│   │       └── Make/
│   │
│   ├── sprayParcel/                     # Modified OpenFOAM spray parcel infrastructure
│   │   ├── SprayParcel.C                # Parent/child state handling and debug logging
│   │   └── BreakupModel.H               # Extended breakup interface
│   │
│   ├── sampling/                        # Custom cloud function objects
│   │   └── ParcelPlaneSampler/          # Crossing-plane parcel sampler
│   │       ├── ParcelPlaneSampler.H/.C
│   │       ├── makeParcelPlaneSampler.C
│   │       └── Make/
│   │
│   ├── sprayFoamDBM/                    # Custom solver: sprayFoam + DBM smoothing
│   │   ├── sprayFoamDBM.C               # Main solver
│   │   ├── createDBMFields.H            # Reads smoothingProperties and creates DBM fields
│   │   ├── smoothSourceTerms.H          # Core DBM algorithm
│   │   ├── Make/
│   │   │   ├── files
│   │   │   └── options
│   │   └── README_sprayFoamDBM.txt      # Compilation and case setup instructions
│   │
│   └── evaporation/                     # Future extensions
│
├── cases/
│   ├── subsonic/
│   │   ├── Zhang_case/                  # Original Zhang subsonic validation case (LMR = 1)
│   │   ├── Zhang/                       # Zhang DBM kernel study (setup-only)
│   │   │   ├── We500/
│   │   │   │   ├── b_0_PCM/             # PCM reference (smoothBandwidth = 0)
│   │   │   │   ├── b_Dinj/              # DBM, b = 1 D_inj
│   │   │   │   ├── b_3Dinj/             # DBM, b = 3 D_inj
│   │   │   │   └── b_6Dinj/             # DBM, b = 6 D_inj
│   │   │   ├── We3500/                  # same four bandwidths
│   │   │   └── We8500/                  # same four bandwidths
│   │   │
│   │   └── Lambert/
│   │       ├── J10_kernelStudy/         # Lambert DBM kernel study at J = 10
│   │       └── J20_kernelStudy/         # Lambert DBM kernel study at J = 20
│   │
│   ├── supersonic/
│   └── evaporation/
│
├── postProcessing/
│   ├── common/
│   ├── subsonic/
│   ├── supersonic/
│   └── evaporation/
│
├── docs/
│   ├── papers/
│   ├── notes/
│   └── figures/
│
└── README.md
```

---

## Validation cases

The `cases/` directory provides reproducible **setup-only** OpenFOAM cases for
the subsonic-crossflow validation and the Diffusion-Based smoothing kernel
study. The repository deliberately does **not** track simulation results,
processor decomposition, `postProcessing` outputs or run logs: only the clean
setup needed to rerun the cases is committed.

For each imported case the following are tracked:

- `setup/0/` – initial and boundary conditions (with the original `U` field
  restored from the source case's `U_backup/` before any overwrite by
  `potentialFoam`);
- `setup/constant/` – thermophysical properties, spray cloud setup, gravity,
  turbulence, radiation and the `smoothingProperties` dictionary controlling
  the DBM bandwidth;
- `setup/system/` – `controlDict`, `fvSchemes`, `fvSolution`, `blockMeshDict`
  and the other case-level dictionaries;
- `setup/chemkin/` – CHEMKIN-format thermodynamic and transport data, when
  required by the case;
- `runScripts/` – `Allrun_automated` and `monitor_*` helpers used for cluster
  execution;
- `postProcessing/` – when present, holds **only** curated lightweight assets
  for the case (final figures in `figures/`, processed CSV/JSON summaries in
  `tables/`, and ParaView colormaps/state files/screenshots in `paraview/`).
  Raw OpenFOAM runtime outputs (numeric time directories, function-object
  outputs, processor folders, parcel-level CSV dumps, logs) are excluded.

The `polyMesh/` directory is regenerated locally via `blockMesh` (plus any
case-specific topology utilities referenced in `system/`) as part of the
`Allrun_automated` workflow; it is not committed in the kernel-study cases.

The top-level `postProcessing/` directory at the repository root collects
reusable post-processing scripts and shared assets; case-specific curated
figures and tables live in each `cases/**/postProcessing/` folder.

### Subsonic kernel-study cases

The `cases/subsonic/Zhang/` and `cases/subsonic/Lambert/{J10,J20}_kernelStudy/`
folders contain the OpenFOAM setups used for the parametric study of the
Diffusion-Based smoothing bandwidth `b`. Each kernel-study group shares the
same physical configuration and differs only in the value of
`smoothBandwidth` in `constant/smoothingProperties`:

| Case label | `smoothBandwidth` | Meaning                          |
|------------|-------------------|----------------------------------|
| `b_0_PCM`  | `0`               | PCM reference, no DBM smoothing  |
| `b_Dinj`   | `1 × D_inj`       | DBM smoothing at one D_inj       |
| `b_3Dinj`  | `3 × D_inj`       | DBM smoothing at three D_inj     |
| `b_6Dinj`  | `6 × D_inj`       | DBM smoothing at six D_inj       |

For the Zhang configuration the injector diameter is `D_inj = 1.6 mm`, so the
recorded values are `0`, `0.0016`, `0.0048` and `0.0096` metres. For the
Lambert configuration the injector diameter gives `0`, `0.000457`, `0.001371`
and `0.002742` metres. The Lambert `b_0_PCM` folder corresponds to the source
case named `b_0_PCM_cell` in the kernel-study run directory and provides the
PCM baseline used for the kernel comparison.

Results for these cases are intentionally excluded from the repository and
should be obtained externally or regenerated from the provided setup.

### Custom libraries used by the cases

All validation cases rely on the updated custom libraries in `src/`:

- `src/breakup/Madabhushi/` – Madabhushi + Lambert multi-stage breakup model;
- `src/drag/Madabhushi/` – `MadabhushiDragForce` drag model with sphere/disc
  deformation correction;
- `src/sprayParcel/` – modified `SprayParcel.C` and `BreakupModel.H` providing
  parent/child state handling for the ligament mechanism;
- `src/sampling/ParcelPlaneSampler/` – crossing-plane parcel sampler used by
  the post-processing pipeline;
- `src/sprayFoamDBM/` – custom solver implementing the Diffusion-Based
  smoothing of the Lagrangian source terms, controlled per-case via
  `constant/smoothingProperties`.

---

## sprayFoamDBM – Gaussian kernel source-term smoothing

### Background

OpenFOAM's standard Lagrangian coupling uses the **Particle Centroid Method (PCM)**: the entire momentum exchange between a parcel and the gas phase is accumulated into the single cell containing the parcel centroid. In dense spray regions near the injector, this creates an artificial velocity deficit in the gas – the so-called **velocity locking** artefact – which suppresses the crossflow momentum transfer to the liquid jet and produces insufficient deflection compared to experiments and ANSYS Fluent results.

ANSYS Fluent addresses this with an optional **Gaussian kernel averaging** that distributes each parcel's source contribution over a neighbourhood of cells weighted by a Gaussian function.

### Method

`sprayFoamDBM` replicates this Gaussian smoothing via the **Diffusion-Based Method (DBM)** of Sun & Xiao (2015). After `parcels.evolve()` accumulates the PCM source terms, and before the PIMPLE loop uses them, the solver applies the following implicit screened-Poisson solve to both `UTrans` (explicit momentum source) and `UCoeff` (implicit drag coefficient):

```text
φ_smooth − σ² ∇²φ_smooth = φ_PCM
```

with zero-gradient (Neumann) boundary conditions on all patches to guarantee global momentum conservation:

```text
∫φ_smooth dV = ∫φ_PCM dV
```

This single implicit step is mathematically equivalent to convolving the PCM field with a 3D Gaussian kernel of bandwidth:

```text
b = √(4σ²)
```

### Tuning

The bandwidth is set in `constant/smoothingProperties`:

```c++
smoothBandwidth   0.0048;   // [m] example: 3 × Dinj for Dinj = 1.6 mm
```

| Value | Effect |
|-------|--------|
| `0` | PCM pure, equivalent to standard sprayFoam coupling |
| `1 × Dinj` | Minimal smoothing |
| `3 × Dinj` | Recommended starting point for Zhang-type subsonic cases |
| `6 × Dinj` | Aggressive smoothing |

---

## ParcelPlaneSampler – crossing-plane post-processing

The `ParcelPlaneSampler` cloud function records parcel crossing events on user-defined downstream planes. It is intended for comparison with reference measurements and simulations where SMD/D32 and droplet velocity profiles are evaluated at fixed locations such as `x/d = 30` and `x/d = 60`.

The sampler writes one event per parcel crossing with the following fields:

```text
time x y z Ux Uy Uz d nParticle user ms tc KHindex yState yDotState faceID
```

These quantities allow the post-processing scripts to compute:

```text
D32(y) = Σ(n_i d_i³) / Σ(n_i d_i²)
Ux(y)  = Σ(n_i Ux_i) / Σ(n_i)
```

and to separate the contributions of different parcel families:

| `user` | Meaning |
|--------|---------|
| `0` | core / parent parcel |
| `1` | Wave-generated child |
| `2` | post-catastrophic child |

The additional state variables `ms`, `tc`, `KHindex`, `yState`, and `yDotState` are diagnostic fields used to identify the breakup state of the sampled parcels.

---

## Compilation

### Local compilation – ESI OpenFOAM v2406, Int32

The modified spray parcel infrastructure changes OpenFOAM template/base files. Therefore, the OpenFOAM spray library and the custom solver must both be recompiled after modifying `SprayParcel.C` or `BreakupModel.H`.

```bash
# Modified OpenFOAM spray library
cd $FOAM_SRC/lagrangian/spray
wclean libso
wmake libso
```

```bash
# Custom breakup model
cd src/breakup/Madabhushi
wclean libso
wmake libso
```

```bash
# Custom drag model
cd src/drag/Madabhushi
wclean libso
wmake libso
```

```bash
# Parcel plane sampler
cd src/sampling/ParcelPlaneSampler
wclean libso
wmake libso
```

```bash
# Custom solver
cd src/sprayFoamDBM
wclean
wmake
```

### Cluster Legion – Politecnico di Torino

```bash
module load openfoam/v2312
cd src/sprayFoamDBM
wclean
wmake
```

The solver binary is installed to `$FOAM_USER_APPBIN`. The custom libraries must be listed in the case dictionary where the Lagrangian cloud is created, typically in `system/controlDict` or in the relevant cloud dictionary.

```c++
libs
(
    "liblagrangianSpray.so"          // modified SprayParcel/BreakupModel; must be loaded first
    "libMadabhushiSprayModels.so"
    "libMadabhushiDragModels.so"
    "libParcelPlaneSampler.so"       // required when parcelPlaneSampler is active
);
```

---

## Case workflow

A typical subsonic crossflow workflow consists of:

1. preparing the initial carrier-gas field;
2. running the gas-only or warm-up phase if required;
3. enabling the Lagrangian spray injection;
4. running `sprayFoamDBM` with the selected breakup, drag and DBM settings;
5. sampling parcel crossings with `ParcelPlaneSampler`;
6. post-processing D32, droplet velocity and parcel-family contributions.

Case-specific run scripts are stored under each case folder in `runScripts/`. These scripts should be treated as case utilities and may assume a specific run-from-zero or restart workflow.

---

## References

The repository keeps a curated bibliography of the publications that support
the implementation and validation of the solver in
[`docs/papers/README.md`](docs/papers/README.md). The entries are grouped by
topic (LJICF/LJIGF experiments, breakup models, drag modelling, coupling and
smoothing methods, evaporation and supersonic references) and document the
papers that are physically tracked under `docs/papers/` together with any
external citations used in the thesis.
