# Thesis Lagrangian Spray Solver – Crossflow Atomization

This repository contains the development and validation of a custom Lagrangian spray solver based on OpenFOAM, carried out within a Master's thesis in Aerospace Engineering at Politecnico di Torino.

The project focuses on the numerical modelling of liquid atomization in gas crossflow conditions, with particular attention to breakup dynamics, droplet aerodynamic behaviour, and validation against reference literature and experimental data.

---

## Project objectives

- Develop custom breakup and drag models for Lagrangian spray simulations
- Improve the physical representation of crossflow atomization
- Correct the velocity locking artefact of OpenFOAM's Particle Centroid Method (PCM) via Gaussian kernel source-term smoothing
- Validate the numerical framework against literature and experimental references
- Perform detailed post-processing of spray characteristics such as droplet size distribution, D32, and dispersion behaviour

---

## Key features

**Custom multi-stage breakup model (Madabhushi + Lambert):**
- liquid column/blob behaviour at injection
- Kelvin–Helmholtz wave stripping
- primary breakup transition with ligament mechanism (Lambert 2019)
- secondary breakup based on Pilch–Erdman-type catastrophic breakup

**Custom drag model (MadabhushiDragForce):**
- blob/column regime drag
- spherical droplet regime drag
- deformation-based drag during breakup transition

**Gaussian kernel source-term smoothing (sprayFoamDBM):**
- replaces OpenFOAM's PCM point-source projection with a spatially distributed Gaussian kernel
- mathematically equivalent to ANSYS Fluent's "Enable Node Based Averaging + Gaussian kernel" option
- eliminates velocity locking artefact that suppresses jet deflection in crossflow
- implemented via implicit screened-Poisson solve on `UTrans` and `UCoeff` fields
- tunable via `constant/smoothingProperties` (bandwidth parameter `b`)

**Crossflow atomization framework:**
- gas–liquid interaction in subsonic crossflow conditions
- reproducible OpenFOAM case organization (warm-up + spray phases)
- dedicated setup, source-code and post-processing structure
- validation-oriented workflow against Zhang et al. and Wu et al. experimental references

---

## Repository structure

```text
Thesis-Lagrangian-Solver/
│
├── src/
│   ├── breakup/                    # Custom breakup model (Madabhushi + Lambert)
│   │   ├── Madabhushi.H/.C
│   │   ├── BreakupModel.H          # Extended base class (childMsInit, childUserInit,
│   │   │                           #   hasPendingChildren, nPendingChildren, debugLog)
│   │   └── Make/
│   │
│   ├── drag/                       # Custom drag model (MadabhushiDragForce)
│   │   ├── MadabhushiDragForce.H/.C
│   │   └── Make/
│   │
│   ├── sprayParcel/                # Modified OpenFOAM SprayParcel
│   │   └── SprayParcel.C           # FIX-SP1: parent user_ flag correction
│   │                               # FIX-SP2: structured CSV debug logging
│   │
│   ├── sprayFoamDBM/               # Custom solver: sprayFoam + DBM smoothing
│   │   ├── sprayFoamDBM.C          # Main solver (sprayFoam + smoothSourceTerms call)
│   │   ├── createDBMFields.H       # Reads smoothingProperties, creates σ² parameter
│   │   ├── smoothSourceTerms.H     # Core DBM algorithm (screened-Poisson solve)
│   │   ├── Make/
│   │   │   ├── files
│   │   │   └── options
│   │   └── README_sprayFoamDBM.txt # Compilation and case setup instructions
│   │
│   └── evaporation/                # Future extensions
│
├── cases/
│   └── subsonic/
│       └── Wu_case/                # Validation case (Wu et al. crossflow)
│           ├── 0/                  # Initial conditions
│           ├── constant/           # Mesh + physical properties
│           │   ├── smoothingProperties   # DBM bandwidth parameter
│           │   └── sprayCloudProperties  # Injection + breakup + drag setup
│           └── system/             # controlDict, fvSchemes, fvSolution
│               └── fvSolution      # Includes convergenceSU/Sp GAMG solvers
├── cases/supersonic/
├── cases/evaporation/
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

```
φ_smooth − σ² ∇²φ_smooth = φ_PCM
```

with zero-gradient (Neumann) boundary conditions on all patches to guarantee global momentum conservation:

```
∫φ_smooth dV = ∫φ_PCM dV
```

This single implicit step is mathematically equivalent to convolving the PCM field with a 3D Gaussian kernel of bandwidth `b = √(4σ²)`.

### Tuning

The bandwidth is set in `constant/smoothingProperties`:

```c++
smoothBandwidth   0.0048;   // [m] → 3 × Dinj for Wu case (Dinj = 1.6 mm)
```

| Value | Effect |
|-------|--------|
| `0` | PCM pure (identical to standard sprayFoam) |
| `1 × Dinj` | Minimal smoothing |
| `3 × Dinj` | Recommended starting point |
| `6 × Dinj` | Aggressive smoothing |

### References

1. Sun, D. & Xiao, H. (2015). *Diffusion-based coarse graining in hybrid continuum-discrete solvers: Theoretical formulation and a priori tests*. International Journal of Multiphase Flow, 77, 142–157. [arXiv:1409.0001](https://arxiv.org/abs/1409.0001)
2. Esgandari, B. et al. (2025). *Diffusion-based smoothing of Lagrangian source terms in Euler-Lagrange simulations*. Journal of Multiphase Flow, 188, 105223.
3. ANSYS Fluent Theory Guide, Ch. 12: *Gaussian kernel for mesh node averaging*.
4. Lambert, A. et al. (2019). *Enhancement of the Madabhushi LJICF Breakup Model by a Ligament Breakup Mechanism*. ILASS-Europe.
5. Madabhushi, R.K. (2003). *A Model for Numerical Simulation of Breakup of a Liquid Jet in Crossflow*. Atomization and Sprays, 13(4).

---

## Compilation

### Local (ESI OpenFOAM v2406, Int32)

```bash
# Custom libraries
cd src/breakup && wclean && wmake
cd src/drag    && wclean && wmake

# Custom solver
cd src/sprayFoamDBM && wclean && wmake
```

### Cluster Legion – Politecnico di Torino (v2312, Int64)

```bash
module load openfoam/v2312
cd src/sprayFoamDBM && wclean && wmake
```

The solver binary is installed to `$FOAM_USER_APPBIN`. The custom libraries must be listed in `system/controlDict`:

```c++
libs
(
    "liblagrangianSpray.so"       // modified SprayParcel – must be FIRST
    "libMadabhushiSprayModels.so"
    "libMadabhushiDragModels.so"
);
```

---

## References

- Wu, P.K. et al. (1997). *Breakup processes of liquid jets in subsonic crossflows*. Atomization and Sprays.
- Zhang, Y. et al. – crossflow atomization VOF-to-DPM benchmark (Aerospace Science and Technology).
- Madabhushi, R.K. (2003). Atomization and Sprays, 13(4).
- Lambert, A. et al. (2019). ILASS-Europe.
- Pilch, M. & Erdman, C.A. (1987). *Use of breakup time data and velocity history data to predict the maximum size of stable fragments*. Int. J. Multiphase Flow.
- Sun & Xiao (2015). Int. J. Multiphase Flow, 77, 142–157.
- Esgandari et al. (2025). J. Multiphase Flow, 188, 105223.
