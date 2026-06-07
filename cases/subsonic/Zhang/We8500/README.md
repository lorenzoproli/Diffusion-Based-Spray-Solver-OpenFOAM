# Zhang – We8500

Kernel-study points at Weber number `We8500` for the Zhang subsonic crossflow
configuration. All four cases share the same physical setup and differ only in
the `smoothBandwidth` parameter of the Diffusion-Based smoothing.

> **⚠️ Convergence note.** This case is published for completeness, but the
> `We8500` simulations **did not reach convergence**: the most energetic Zhang
> condition (Sim. 4, throttle level 100 %) **crashed** before producing usable
> statistics. As discussed in the thesis, the failure is attributed to the fact
> that the Diffusion-Based smoothing (DBM) is currently applied only to the
> **momentum** source term, whereas at this very high Weber number the **energy**
> (and mass) exchange source terms grow large enough to dominate the
> Lagrangian–Eulerian coupling. Extending the kernel smoothing to the remaining
> Lagrangian (DPM) source terms — in particular the **energy** source — is
> required to stabilise the case and is left as future work. `We8500` is
> therefore provided as a setup-only reference, **without validated results**.

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
