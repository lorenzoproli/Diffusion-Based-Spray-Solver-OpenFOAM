# Zhang – We500 – b_6Dinj – postProcessing

This folder contains curated lightweight post-processing assets for this specific case. Raw OpenFOAM outputs (numeric time directories, processor decomposition, function-object outputs, parcel-level CSV dumps, logs) are intentionally excluded; only final figures, lightweight processed tables, ParaView visualization assets and the post-processing scripts that produced them are tracked here.

## Structure

- `figures/` – final post-processing plots (PNG).
- `tables/` – lightweight processed summaries (CSV, JSON).
- `paraview/` – ParaView colormaps, state files, presets and screenshots.
- `scripts/` – Python scripts used to generate the curated outputs.

## figures/

- `profile_Zhang_b_6Dinj_x060mm_D10.png`
- `profile_Zhang_b_6Dinj_x060mm_D32.png`
- `profile_Zhang_b_6Dinj_x060mm_Umean.png`
- `profile_Zhang_b_6Dinj_x060mm_count.png`

## tables/

- `metadata_b_6Dinj.json`
- `profile_Zhang_b_6Dinj_x060mm.csv`

## paraview/

- `We_500_U_6Dinj.png`
- `We_500_Up_6Dinj.png`
- `We_500_d_6Dinj.png`
- `We_500_p_6Dinj.png`

## scripts/

- `calcolo_D32.py`
- `postprocess_zhang_kernel_single.py`
