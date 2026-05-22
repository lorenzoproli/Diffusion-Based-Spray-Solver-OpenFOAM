# Lambert – J20_kernelStudy – b_3Dinj – postProcessing

This folder contains curated lightweight post-processing assets for this specific case. Raw OpenFOAM outputs (numeric time directories, processor decomposition, function-object outputs, parcel-level CSV dumps, logs) are intentionally excluded; only final figures, lightweight processed tables, ParaView visualization assets and the post-processing scripts that produced them are tracked here.

## Structure

- `figures/` – final post-processing plots (PNG).
- `tables/` – lightweight processed summaries (CSV, JSON).
- `paraview/` – ParaView colormaps, state files, presets and screenshots.
- `scripts/` – Python scripts used to generate the curated outputs.

## figures/

- `Lambert_J20_b_3Dinj_Lambert_OpenFOAM_Comparison.png`

## tables/

- `Lambert_J20_b_3Dinj_profile_OpenFOAM_30d.csv`
- `Lambert_J20_b_3Dinj_profile_OpenFOAM_60d.csv`

## paraview/

_Reserved for ParaView colormaps, state files (`.pvsm`), presets and screenshots to be added later._

## scripts/

- `postprocess_lambert_openfoam.py`
