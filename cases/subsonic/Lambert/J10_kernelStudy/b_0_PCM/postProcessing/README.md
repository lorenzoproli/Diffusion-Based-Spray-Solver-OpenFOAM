# Lambert – J10_kernelStudy – b_0_PCM – postProcessing

This folder contains curated lightweight post-processing assets for this specific case. Raw OpenFOAM outputs (numeric time directories, processor decomposition, function-object outputs, parcel-level CSV dumps, logs) are intentionally excluded; only final figures, lightweight processed tables, ParaView visualization assets and the post-processing scripts that produced them are tracked here.

## Structure

- `figures/` – final post-processing plots (PNG).
- `tables/` – lightweight processed summaries (CSV, JSON).
- `paraview/` – ParaView colormaps, state files, presets and screenshots.
- `scripts/` – Python scripts used to generate the curated outputs.

## figures/

- `Lambert_J10_b_0_PCM_Lambert_OpenFOAM_Comparison.png`

## tables/

- `Lambert_J10_b_0_PCM_profile_OpenFOAM_30d.csv`
- `Lambert_J10_b_0_PCM_profile_OpenFOAM_60d.csv`

## paraview/

- `j10_U_PCM.png`
- `j10_Ux_PCM.png`
- `j10_d_PCM.png`
- `j10_p_PCM.png`

## scripts/

- `postprocess_lambert_openfoam.py`
- `postprocess_lambert_openfoam_filter.py`
