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
- validation-oriented workflow against Zhang-type subsonic crossflow benchmarks and Lambert J10/J20 reference cases, with experimental comparisons from the literature

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
│   │   ├── Zhang_case/                  # Subsonic crossflow validation case
│   │   │   ├── setup/                   # OpenFOAM setup files
│   │   │   ├── runScripts/              # Run and monitoring scripts
│   │   │   └── results/                 # Optional lightweight reference outputs
│   │   │
│   │   └── Lambert/
│   │       ├── J10/                     # Lambert momentum-ratio J=10 case
│   │       ├── J20/                     # Lambert momentum-ratio J=20 case
│   │       ├── J10_Ligament/            # J10 with updated ligament-aware libraries
│   │       └── J20_Ligament/            # J20 with updated ligament-aware libraries
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

- Wu, P.K. et al. (1997). *Breakup processes of liquid jets in subsonic crossflows*. Atomization and Sprays.
- Zhang, Y. et al. – subsonic crossflow atomization benchmark.
- Madabhushi, R.K. (2003). *A Model for Numerical Simulation of Breakup of a Liquid Jet in Crossflow*. Atomization and Sprays, 13(4).
- Lambert, A. et al. (2019). *Enhancement of the Madabhushi LJICF Breakup Model by a Ligament Breakup Mechanism*. ILASS-Europe.
- Pilch, M. & Erdman, C.A. (1987). *Use of breakup time data and velocity history data to predict the maximum size of stable fragments*. International Journal of Multiphase Flow.
- Sun, D. & Xiao, H. (2015). *Diffusion-based coarse graining in hybrid continuum-discrete solvers: Theoretical formulation and a priori tests*. International Journal of Multiphase Flow, 77, 142–157.
- Esgandari, B. et al. (2025). *Diffusion-based smoothing of Lagrangian source terms in Euler-Lagrange simulations*. Journal of Multiphase Flow, 188, 105223.
- ANSYS Fluent Theory Guide, Ch. 12: *Gaussian kernel for mesh node averaging*.
