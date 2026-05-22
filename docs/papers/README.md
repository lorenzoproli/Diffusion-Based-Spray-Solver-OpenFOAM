# Bibliography – Thesis Lagrangian Spray Solver

This file lists the reference publications that support the implementation
and validation of the custom OpenFOAM Lagrangian spray solver developed in
this thesis. The corresponding PDF files are tracked locally under
`docs/papers/` and can be opened directly from the repository.

References are grouped by topic to make it easy to map each component of the
solver (breakup, drag, source-term smoothing, evaporation) to the literature
on which it is based, and each validation case to the experimental or
numerical work that defines it.

Where author lists, journal names, volumes, pages or DOIs cannot be inferred
reliably from the local PDF filenames, the corresponding fields are marked as
`metadata to verify` rather than guessed.

---

## 1. Liquid jet in subsonic crossflow – experiments and reference data

- **Wu, P.-K., Kirkendall, K.A., Fuller, R.P., Nejad, A.S.** (1997).
  *Breakup processes of liquid jets in subsonic crossflows*.
  Journal of Propulsion and Power, **13**(1), 64–73.
  – Local PDF: `Breakup Processes of Liquid Jets in Subsonic Crossflows.pdf`.
  – Used as a reference for the breakup regimes of liquid jets in subsonic
  gaseous crossflows and as one of the experimental data sources for the
  Lagrangian validation activities.

- **Wu, P.-K., Kirkendall, K.A., Fuller, R.P., Nejad, A.S.** (1998).
  *Spray trajectories of liquid fuel jets in subsonic crossflows*.
  International Journal of Fluid Mechanics Research / AIAA – metadata to verify.
  – Local PDF: `Spray Trajectories of Liquid Fuel Jets in Subsonic Crossflows.pdf`.
  – Companion paper providing penetration and trajectory correlations used to
  benchmark the predicted spray penetration in the subsonic crossflow cases.

- **Lambert, A., Wirth, C., Roussel, S., Sloan, G.** (2019).
  *Enhancement of the Madabhushi LJICF breakup model by a ligament breakup
  mechanism* – metadata to verify.
  ILASS-Europe 2019, 29th Conference on Liquid Atomization and Spray Systems.
  – Local PDF: `Paper_Lambert_ILASS2019_v10.pdf`.
  – Defines the J = 10 and J = 20 Liquid Jet In Crossflow (LJICF) cases used
  in the `cases/subsonic/Lambert/{J10,J20}_kernelStudy/` kernel-study setups
  and introduces the ligament-extended breakup mechanism implemented in
  `src/breakup/Madabhushi/`.

---

## 2. Liquid jet in subsonic gas film – GLOP / Zhang cases

- **Zhang, M. et al.** *Atomization dynamics of gas-liquid orifice type pintle
  injector under different operating conditions* – metadata to verify.
  – Local PDF: `Atomization dynamics of gas-liquid orifice type pintle
  injector under.pdf`.
  – Gas-Liquid Orifice type Pintle (GLOP) injector study, used as a reference
  for the Liquid Jet In subsonic Gas Film (LJIGF) configuration that
  motivates the `cases/subsonic/Zhang_case/` LMR = 1 validation setup.

- **Zhang, M. et al.** *Atomization dynamics of gas-liquid orifice type pintle
  injectors at different local momentum ratios*.
  Aerospace Science and Technology – metadata to verify (year, volume, pages).
  – Local PDF: `Atomization dynamics of gas-liquid orifice type pintle
  injectors at different_zhang_compressed.pdf`.
  – Provides the parametric data set at different Weber numbers used to
  derive the Zhang kernel-study cases stored under
  `cases/subsonic/Zhang/{We500,We3500,We8500}/`.

---

## 3. Primary and secondary breakup modelling

- **Madabhushi, R.K.** (2003).
  *A model for numerical simulation of breakup of a liquid jet in crossflow*.
  Atomization and Sprays, **13**(4). Begell House.
  DOI: `10.1615/AtomizSpr.v13.i4` – metadata to verify (issue page range).
  – Local PDF: `Madabhushi, Ravi K. (author) - A MODEL FOR NUMERICAL
  SIMULATION OF BREAKUP OF A LIQUID JET IN CROSSFLOW (2003, Begell House
  Inc.) [10.1615_atomizspr.v13.i4..pdf`.
  – Base formulation of the LJICF column/Kelvin–Helmholtz breakup model
  implemented in `src/breakup/Madabhushi/` and extended with the Lambert
  ligament mechanism.

- **Pilch, M., Erdman, C.A.** (1987).
  *Use of breakup time data and velocity history data to predict the maximum
  size of stable fragments for acceleration-induced breakup of a liquid drop*.
  International Journal of Multiphase Flow, **13**(6), 741–757.
  – Local PDF: `Use of breakup time data and velocity history data to predict
  the maximum size of stable fragment.pdf`.
  – Source of the Pilch–Erdman (PE) catastrophic secondary-breakup criteria
  used inside the Madabhushi + Lambert model to handle post-column droplets.

- **Reitz, R.D.** (1987).
  *Modeling atomization processes in high-pressure vaporizing sprays*.
  Atomization and Spray Technology, **3**, 309–337 – metadata to verify.
  – Local PDF: `Modeling_atomization_process_in_high_pre.pdf`.
  – Reference for the Kelvin–Helmholtz Wave (KH) atomization model that
  underlies the wave-stripping stage of the Madabhushi-type breakup
  implementation.

---

## 4. Diffusion-based / coarse-graining source-term smoothing

- **Sun, R., Xiao, H.** (2015).
  *Diffusion-based coarse graining in hybrid continuum-discrete solvers:
  Theoretical formulation and a priori tests*.
  International Journal of Multiphase Flow, **77**, 142–157.
  – Local PDF: `Diffusion-based coarse graining in hybrid
  continuum–discrete solvers.pdf`.
  – Theoretical basis of the Diffusion-Based Method (DBM) implemented in
  `src/sprayFoamDBM/`, in which the Lagrangian source terms are smoothed via
  an implicit screened-Poisson solve equivalent to a Gaussian convolution.

---

## 5. Evaporation and spray modelling in combustion chambers

- **Schmehl, R., Klose, G., Maier, G., Wittig, S.** (1997).
  *Efficient numerical calculation of evaporating sprays in combustion
  chamber flows* – metadata to verify (proceedings, page range).
  – Local PDF: `Efficient_Numerical_Calculation_of_Evapo.pdf`.
  – Reference for the evaporation modelling strategy considered for future
  extensions of the solver (see `src/evaporation/`).

---

## 6. Supersonic crossflow references (context)

- **Lin, K.-C., Kennedy, P.J., Jackson, T.A.** (2004).
  *Structures of water jets in a Mach 1.94 supersonic crossflow* – metadata
  to verify (AIAA paper number, proceedings).
  – Local PDF: `STRUCTURES OF WATER JETS IN A MACH 1.94 SUPERSONIC
  CROSSFLOW.pdf`.
  – Used as background reference for the supersonic crossflow extension
  considered in `cases/supersonic/` and not as a current validation target
  of this snapshot.

---

## Cross-references

- The solver implementation that builds on these references is documented in
  the top-level [`README.md`](../../README.md) and in
  [`src/sprayFoamDBM/README_sprayFoamDBM.txt`](../../src/sprayFoamDBM/README_sprayFoamDBM.txt).
- The validation cases that rely on these references are documented in:
  - [`cases/subsonic/README.md`](../../cases/subsonic/README.md)
  - [`cases/subsonic/Zhang/README.md`](../../cases/subsonic/Zhang/README.md)
  - [`cases/subsonic/Lambert/README.md`](../../cases/subsonic/Lambert/README.md)

## Notes on metadata

The bibliographic entries above are reconstructed from the file names of the
PDFs tracked locally in this directory. Some PDFs only carry partial or
truncated identifying information in their filename; the corresponding
journal, volume, page range, year or DOI fields are explicitly flagged with
`metadata to verify` and should be confirmed against the original PDF
metadata or the publisher's record before being used in formal citations of
the thesis.
