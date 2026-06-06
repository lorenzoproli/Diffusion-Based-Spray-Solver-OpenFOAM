# Madabhushi drag model (custom implementation)

This folder contains a custom drag force model for Lagrangian spray particles, developed for crossflow atomization simulations within the thesis work.

## Purpose

The model provides a physically consistent drag formulation across different stages of the atomization process.

It is tightly coupled with the custom breakup model and is validated within the LMR = 1 crossflow benchmark case.

## Modelling approach

The drag coefficient is computed using a multi-regime formulation depending on the particle state.

### 1. Blob / liquid column regime

- Applied before primary breakup
- Represents intact liquid structures exposed to crossflow
- Constant drag coefficient:

  Cd ≈ 1.48

- Captures the aerodynamic behaviour of elongated liquid columns

### 2. Spherical droplet regime

- Applied to standard droplets after breakup
- Classical drag law based on Reynolds number
- Suitable for fully formed droplets

### 3. Deformation / breakup regime (Pilch–Erdman-based)

- Applied during active breakup
- Drag is computed using a deformation-based formulation:

  Cd = CdDisc · (d_ref / d)

- Accounts for increased drag due to droplet deformation
- Ensures consistency with the breakup dynamics

## Dictionary coefficients

The following parameters are read from `MadabhushiDragForce` (or `MadabhushiLigamentDragForce`) in `constant/sprayCloudProperties`:

| Key | Default | Description |
|-----|---------|-------------|
| `CdBlob` | 1.48 | Drag coefficient for liquid column / blob regime |
| `C0` | 3.44 | Column breakup constant (shared with breakup model) |
| `Dinj` | 0.0016 | Injector diameter [m] |
| `debug` | false | Enable CSV debug logging |

The following parameter is **hardcoded** in the source and is **not** read from the dictionary:

| Parameter | Fixed value | Description |
|-----------|-------------|-------------|
| `CdDisc` | 1.2 | Drag coefficient for deformed droplets (Madabhushi 2003 Eq. 7) |

## Implementation features

- Dynamic selection of drag regime based on particle state
- Coupling with breakup stage information
- Consistent handling of child droplets
- Compatible with OpenFOAM Lagrangian framework

## Files

- `MadabhushiDragForce.C`: main implementation
- `MadabhushiDragForce.H`: class definition
- `makeMadabhushiDrag.C`: runtime selection registration
- `Make/files`, `Make/options`: compilation settings

## Role in the thesis

This model is used to:

- improve prediction of droplet trajectories in crossflow
- correctly represent aerodynamic forces during breakup
- ensure consistency between breakup and drag physics
- support validation against numerical and experimental data at LMR = 1

