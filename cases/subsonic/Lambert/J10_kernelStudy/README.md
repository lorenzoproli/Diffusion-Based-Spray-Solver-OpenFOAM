# Lambert – J10 – kernel study

Kernel-study points for the Lambert (ILASS 2019) reference case at momentum
ratio `J10`. All four cases share the same physical setup and differ only in
the `smoothBandwidth` parameter of the Diffusion-Based smoothing.

| Case        | smoothBandwidth | Description                              |
|-------------|-----------------|------------------------------------------|
| `b_0_PCM`   | 0               | PCM reference (no smoothing)             |
| `b_Dinj`    | 1 × D_inj       | Diffusion-Based smoothing, b = 1 D_inj   |
| `b_3Dinj`   | 3 × D_inj       | Diffusion-Based smoothing, b = 3 D_inj   |
| `b_6Dinj`   | 6 × D_inj       | Diffusion-Based smoothing, b = 6 D_inj   |

The `b_0_PCM` folder corresponds to the source case named `b_0_PCM_cell` in
the kernel-study run directory; it provides the PCM reference with the
corrected cell-point projection used as baseline for the kernel comparison.

The exact bandwidth in metres is reported in each case's
`constant/smoothingProperties` and in the case-level `README.md`.

Only the clean OpenFOAM setup is included. Numeric time directories, processor
folders, `postProcessing` outputs and run logs are intentionally excluded.
