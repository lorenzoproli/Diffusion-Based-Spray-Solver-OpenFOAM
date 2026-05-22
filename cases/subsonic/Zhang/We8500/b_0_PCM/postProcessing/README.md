# Zhang – We8500 – b_0_PCM – postProcessing

This folder contains curated lightweight post-processing assets for this specific case. Raw OpenFOAM outputs (numeric time directories, processor decomposition, function-object outputs, parcel-level CSV dumps, logs) are intentionally excluded; only final figures, lightweight processed tables, ParaView visualization assets and the post-processing scripts that produced them are tracked here.

## Structure

- `figures/` – final post-processing plots (PNG).
- `tables/` – lightweight processed summaries (CSV, JSON).
- `paraview/` – ParaView colormaps, state files, presets and screenshots.
- `scripts/` – Python scripts used to generate the curated outputs.

## figures/

- `Zhang_We8500_b_0_PCM_boundary_front_xz.png`
- `Zhang_We8500_b_0_PCM_boundary_left_yz.png`
- `Zhang_We8500_b_0_PCM_boundary_top_xy.png`
- `Zhang_We8500_b_0_PCM_overlay_D10_vs_zhang2025.png`
- `Zhang_We8500_b_0_PCM_overlay_D32_vs_zhang2025.png`
- `Zhang_We8500_b_0_PCM_overlay_U_vs_zhang2025.png`
- `Zhang_We8500_b_0_PCM_overlay_count_vs_zhang2025.png`
- `Zhang_We8500_b_0_PCM_projection_front_xz.png`
- `Zhang_We8500_b_0_PCM_projection_left_yz.png`
- `Zhang_We8500_b_0_PCM_projection_top_xy.png`
- `Zhang_We8500_b_0_PCM_x_060mm_downstream_D10.png`
- `Zhang_We8500_b_0_PCM_x_060mm_downstream_D32.png`
- `Zhang_We8500_b_0_PCM_x_060mm_downstream_Umean.png`
- `Zhang_We8500_b_0_PCM_x_060mm_downstream_count.png`
- `Zhang_We8500_b_0_PCM_x_060mm_downstream_dropletCount.png`

## tables/

- `Zhang_We8500_b_0_PCM_metadata.json`
- `Zhang_We8500_b_0_PCM_x_060mm_downstream.csv`

## paraview/

_Reserved for ParaView colormaps, state files (`.pvsm`), presets and screenshots to be added later._

## scripts/

- `calcolo_D32.py`
- `postprocess_zhang_single.py`
