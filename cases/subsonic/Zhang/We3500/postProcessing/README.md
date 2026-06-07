# Zhang – We3500 – kernel-comparison postProcessing

This folder contains the curated **cross-bandwidth** post-processing assets that compare the four DBM bandwidth subcases (`b_0_PCM`, `b_Dinj`, `b_3Dinj`, `b_6Dinj`) of this group against each other and against reference data, together with the Python scripts that produced them.

## Structure

- `figures/` – aggregate comparison plots. Where available the figures are tracked at three granularities: the combined multi-panel overview (`*_all.png`), the per-variable comparison across multiple axial locations (`*_<var>.png`), and the single-location comparison (`*_<var>_<location>.png`).
- `tables/` – kernel-comparison error/summary tables and reference experimental/Fluent profiles when used as overlays.
- `paraview/` – placeholder for cross-bandwidth ParaView assets (colormaps, state files, presets).
- `scripts/` – Python scripts used to generate these comparisons.

## figures/

- `Zhang_kernel_We3500_D10_x060mm.png`
- `Zhang_kernel_We3500_D32_x060mm.png`
- `Zhang_kernel_We3500_Umean_x060mm.png`
- `Zhang_kernel_We3500_b_0_PCM_x060mm.png`
- `Zhang_kernel_We3500_b_3Dinj_x060mm.png`
- `Zhang_kernel_We3500_b_6Dinj_x060mm.png`
- `Zhang_kernel_We3500_b_Dinj_x060mm.png`
- `Zhang_kernel_We3500_count_x060mm.png`

## tables/

_Currently a placeholder; no tables have been added yet._

## paraview/

_Reserved for ParaView colormaps, state files (`.pvsm`), presets and screenshots to be added later._

## scripts/

- `postprocess_zhang_kernel_single.py`
- `postprocess_zhang_single.py`

Per-bandwidth curated assets (per-case profiles, ParaView field screenshots and tables) live in each subcase folder:

- `b_0_PCM/postProcessing/`
- `b_Dinj/postProcessing/`
- `b_3Dinj/postProcessing/`
- `b_6Dinj/postProcessing/`
