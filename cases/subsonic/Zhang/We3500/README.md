# Zhang – We3500

Kernel-study points at Weber number `We3500` for the Zhang subsonic crossflow
configuration. All four cases share the same physical setup and differ only in
the `smoothBandwidth` parameter of the Diffusion-Based smoothing.

| Case        | smoothBandwidth | Description                              |
|-------------|-----------------|------------------------------------------|
| `b_0_PCM`   | 0               | PCM reference (no smoothing)             |
| `b_Dinj`    | 1 × D_inj       | Diffusion-Based smoothing, b = 1 D_inj   |
| `b_3Dinj`   | 3 × D_inj       | Diffusion-Based smoothing, b = 3 D_inj   |
| `b_6Dinj`   | 6 × D_inj       | Diffusion-Based smoothing, b = 6 D_inj   |

The exact bandwidth in metres is reported in each case's
`constant/smoothingProperties` and in the case-level `README.md`.

Only the clean OpenFOAM setup is included. Numeric time directories, processor
folders, `postProcessing` outputs and run logs are intentionally excluded.
