## Automated Run Script (`Allrun_automated`)

This script executes a three-stage transient simulation using `sprayFoam`, designed to ensure a stable and physically consistent gas field before spray injection.

---

### Simulation Workflow

The simulation is divided into three phases:

**Phase 1a — Gas Warm-up**
- Time: `0 → 0.005 s`
- Schemes: Euler (1st order), upwind
- Spray: OFF  
- Includes `potentialFoam` initialization

Purpose: stabilize the flow and damp numerical oscillations from uniform initial conditions.

---

**Phase 1b — Gas Stabilization**
- Time: `0.005 → 0.015 s`
- Schemes: backward (2nd order), linearUpwind
- Spray: OFF

Purpose: improve accuracy once the flow is sufficiently developed.

---

**Phase 2 — Spray Simulation**
- Time: `0.015 → 0.030 s`
- Schemes: backward, linearUpwind
- Spray: ON

Notes:
- Start of Injection (SOI) is relative to Phase 2 start time
- Injection ends at `t = 0.0295 s`
- Uses Madabhushi breakup and Lambert drag models

---

### Usage

Run from the case directory:

```bash
./Allrun_automated           # complete run
./Allrun_automated --gas     # run only Phase 1 (no spray)
./Allrun_automated --1b      # skip Phase 1a, start from latestTime
./Allrun_automated --spray   # run only Phase 2 from latestTime
./Allrun_automated --clean   # remove processor folders and logs
