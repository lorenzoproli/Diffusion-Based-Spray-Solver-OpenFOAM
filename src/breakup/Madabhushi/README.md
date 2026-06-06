# Madabhushi breakup model (custom implementation)

This folder contains a custom breakup model for Lagrangian spray simulations in OpenFOAM, developed within the thesis work for crossflow atomization applications.

## Purpose

The model is designed to improve the physical representation of liquid jet atomization in subsonic crossflow conditions, and is validated against a reference case at LMR = 1.

It is not a direct implementation of a single literature model, but a hybrid formulation combining multiple breakup mechanisms.

## Modelling approach

The breakup process is described through a multi-stage framework:

### 1. Liquid column / wave stage

- The injected liquid initially behaves as a continuous column exposed to crossflow
- Kelvin–Helmholtz (KH) instability governs the early dynamics
- A wave growth model is used to estimate instability development
- Progressive mass stripping occurs from the liquid core

### 2. Wave-shed droplets

- Droplets are generated from KH-driven stripping
- These droplets are tracked separately from standard parcels
- Special internal flags are used to identify their origin

### 3. Primary breakup transition

- Triggered based on local Weber number and instability evolution
- Marks the transition from column-dominated behaviour to dispersed phase

### 4. Secondary breakup (Pilch–Erdman regime)

- Droplet breakup follows Pilch–Erdman-type correlations
- Breakup time and deformation depend on Weber number and flow conditions
- Droplet size evolution is modelled dynamically

### 5. Child droplet generation

- Child parcels are created with corrected diameter and velocity
- Additional corrections are applied to improve physical consistency
- Inspired by literature formulations (e.g. Lambert-type approaches)

## Internal state handling

The model relies on internal variables to manage the breakup stages:

- `ms`: breakup state identifier
- `y`, `yDot`: deformation / instability variables
- `user`: tracking of breakup stage transitions
- `KHindex`: wave growth tracking parameter

This effectively defines a state-machine governing the evolution of each parcel.

## Key features

- Multi-regime breakup modelling
- Explicit treatment of liquid column behaviour
- KH-driven wave stripping
- Pilch–Erdman secondary breakup
- Physically consistent child droplet generation
- Designed for crossflow atomization scenarios

## Files

- `Madabhushi.C`: main implementation
- `Madabhushi.H`: class definition
- `makeMadabhushiBreakup.C`: runtime selection registration
- `Make/files`, `Make/options`: compilation settings

## Dictionary coefficients

The following parameters are read from `MadabhushiCoeffs` (or `MadabhushiLigamentCoeffs`) in `constant/sprayCloudProperties`:

| Key | Default | Description |
|-----|---------|-------------|
| `C0` | 3.44 | Column breakup time constant (Lambert 2019 Eq. 1) |
| `FLig` | 0.4 | Ligament weighting factor (Lambert 2019 Eq. 11) |
| `b0` | 0.61 | Reitz Wave model constant B0 |
| `b1` | 10.0 | Reitz Wave model constant B1 |
| `Dinj` | 0.0016 | Injector diameter [m] |
| `debug` | false | Enable CSV debug logging |

The following parameters are **hardcoded** in the source and are **not** read from the dictionary:

| Parameter | Fixed value | Description |
|-----------|-------------|-------------|
| `nChildren` | 5 | Number of child parcels at catastrophic breakup |
| `minChildDiameter` | 1e-7 m | Numerical lower bound for sampled child diameters |

## Role in the thesis

This model is a core component of the thesis solver and is used to:

- control breakup onset and evolution
- reproduce atomization behaviour in crossflow
- enable comparison with literature and experimental data
- support validation at LMR = 1
