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

## Key parameters

- `CdBlob`: drag coefficient for liquid column / blob regime
- `CdDisc`: drag coefficient for deformed droplets
- `UgRef`: reference gas velocity for deformation scaling
- `Dinj`: injector diameter

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

