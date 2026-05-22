# Lambert – J20_kernelStudy – kernel-comparison postProcessing

This folder contains the curated **cross-bandwidth** post-processing assets that compare the four DBM bandwidth subcases (`b_0_PCM`, `b_Dinj`, `b_3Dinj`, `b_6Dinj`) of this group against each other and against reference data, together with the Python scripts that produced them.

## Structure

- `figures/` – aggregate comparison plots. Where available the figures are tracked at three granularities: the combined multi-panel overview (`*_all.png`), the per-variable comparison across multiple axial locations (`*_<var>.png`), and the single-location comparison (`*_<var>_<location>.png`).
- `tables/` – kernel-comparison error/summary tables and reference experimental/Fluent profiles when used as overlays.
- `paraview/` – placeholder for cross-bandwidth ParaView assets (colormaps, state files, presets).
- `scripts/` – Python scripts used to generate these comparisons.

## figures/

- `Lambert_J20_kernel_comparison_D32.png`
- `Lambert_J20_kernel_comparison_D32_30d.png`
- `Lambert_J20_kernel_comparison_D32_60d.png`
- `Lambert_J20_kernel_comparison_Ux.png`
- `Lambert_J20_kernel_comparison_Ux_30d.png`
- `Lambert_J20_kernel_comparison_Ux_60d.png`
- `Lambert_J20_kernel_comparison_all.png`

## tables/

- `Lambert_J20_kernel_error_summary.csv`
- `Lambert_J20_reference_SMD_x_d_30_exp.csv`
- `Lambert_J20_reference_SMD_x_d_30_fluent.csv`
- `Lambert_J20_reference_SMD_x_d_60_exp.csv`
- `Lambert_J20_reference_SMD_x_d_60_fluent.csv`
- `Lambert_J20_reference_Ux_x_d_60_exp.csv`
- `Lambert_J20_reference_Ux_x_d_60_fluent.csv`
- `Lambert_J20_reference_ux_x_d_30_exp.csv`
- `Lambert_J20_reference_ux_x_d_30_fluent.csv`

## paraview/

_Reserved for ParaView colormaps, state files (`.pvsm`), presets and screenshots to be added later._

## scripts/

- `compare_kernels_with_refs.py`
- `postprocess_lambert_openfoam.py`

Per-bandwidth curated assets (per-case profiles, ParaView field screenshots and tables) live in each subcase folder:

- `b_0_PCM/postProcessing/`
- `b_Dinj/postProcessing/`
- `b_3Dinj/postProcessing/`
- `b_6Dinj/postProcessing/`
