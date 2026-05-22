# Zhang – We3500 – b_Dinj – postProcessing

This folder contains curated lightweight post-processing assets for this specific case. Raw OpenFOAM outputs (numeric time directories, processor decomposition, function-object outputs, parcel-level CSV dumps, logs) are intentionally excluded; only final figures, lightweight processed tables, ParaView visualization assets and the post-processing scripts that produced them are tracked here.

## Structure

- `figures/` – final post-processing plots (PNG).
- `tables/` – lightweight processed summaries (CSV, JSON).
- `paraview/` – ParaView colormaps, state files, presets and screenshots.
- `scripts/` – Python scripts used to generate the curated outputs.

## figures/

- `Zhang_We3500_b_Dinj_boundary_front_xz.png`
- `Zhang_We3500_b_Dinj_boundary_left_yz.png`
- `Zhang_We3500_b_Dinj_boundary_top_xy.png`
- `Zhang_We3500_b_Dinj_overlay_D10_vs_zhang2025.png`
- `Zhang_We3500_b_Dinj_overlay_D32_vs_zhang2025.png`
- `Zhang_We3500_b_Dinj_overlay_U_vs_zhang2025.png`
- `Zhang_We3500_b_Dinj_overlay_count_vs_zhang2025.png`
- `Zhang_We3500_b_Dinj_projection_front_xz.png`
- `Zhang_We3500_b_Dinj_projection_left_yz.png`
- `Zhang_We3500_b_Dinj_projection_top_xy.png`
- `Zhang_We3500_b_Dinj_x_060mm_downstream_D10.png`
- `Zhang_We3500_b_Dinj_x_060mm_downstream_D32.png`
- `Zhang_We3500_b_Dinj_x_060mm_downstream_Umean.png`
- `Zhang_We3500_b_Dinj_x_060mm_downstream_count.png`
- `Zhang_We3500_b_Dinj_x_060mm_downstream_dropletCount.png`

## tables/

- `Zhang_We3500_b_Dinj_metadata.json`
- `Zhang_We3500_b_Dinj_x_060mm_downstream.csv`

## paraview/

- `Zhang_We3500_b_Dinj_diameter.png`
- `Zhang_We3500_b_Dinj_pressure.png`
- `Zhang_We3500_b_Dinj_velocity.png`
- `Zhang_We3500_b_Dinj_velocity_magnitude.png`
- `Zhang_We3500_b_Dinj_weber.png`

## scripts/

- `calcolo_D32.py`
- `postprocess_zhang_single.py`
