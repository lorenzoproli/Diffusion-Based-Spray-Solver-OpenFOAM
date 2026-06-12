# Bibliography — Diffusion-Based Lagrangian Spray Solver

This file lists the reference publications that support the implementation and
validation of the custom OpenFOAM Lagrangian spray solver developed in this
thesis.

> **Copyright notice.** The PDF files of these papers are **not** redistributed
> in this repository. Each entry below provides the full bibliographic
> reference together with a DOI or publisher link; please obtain the papers
> directly from their publishers or official repositories. The complete,
> machine-readable bibliography is maintained in
> [`thesis/bibliografia.bib`](../../thesis/bibliografia.bib).

References are grouped by topic so that each component of the solver (breakup,
drag, source-term smoothing, evaporation) and each validation case can be
mapped to the literature on which it is based. Fields that could not be
confirmed against the publisher record are flagged as `metadata to verify`.

---

## 1. Liquid jet in subsonic crossflow — experiments and reference data

- **Wu, P.-K., Kirkendall, K. A., Fuller, R. P., Nejad, A. S.** (1997).
  *Breakup processes of liquid jets in subsonic crossflows*.
  Journal of Propulsion and Power, **13**(1), 64–73.
  DOI: [10.2514/2.5151](https://doi.org/10.2514/2.5151).
  — Reference for the breakup regimes of liquid jets in subsonic gaseous
  crossflows and for the column drag coefficient used in the LJICF cases.

- **Wu, P.-K., Kirkendall, K. A., Fuller, R. P., Nejad, A. S.** (1998).
  *Spray trajectories of liquid fuel jets in subsonic crossflows*.
  Journal of Propulsion and Power — `metadata to verify` (volume, issue, pages).
  — Companion paper providing penetration and trajectory correlations used as
  background for the predicted spray penetration.

- **Gopala, Y., Zhang, P., Bibik, O., Lubarsky, E., Zinn, B. T.** (2010).
  *Liquid jet in crossflow — trajectory correlations based on the column
  breakup point*. 48th AIAA Aerospace Sciences Meeting, AIAA Paper 2010-214,
  Orlando, Florida, USA.
  DOI: [10.2514/6.2010-214](https://doi.org/10.2514/6.2010-214).
  — One of the experimental sources for the Lambert J10/J20 validation data.

- **Sekar, J., Rao, A., Pillutla, S., Danis, A., Hsieh, S. Y.** (2014).
  *Liquid jet in cross flow modeling*. Proceedings of ASME Turbo Expo 2014,
  Düsseldorf, Germany — `metadata to verify` (paper number GT2014-..., DOI).
  — Numerical reference adopted, together with the experimental data, for the
  Lambert LJICF benchmark.

- **Lambert, M., Esch, T., Braun, M., Elasrag, H.** (2019).
  *Enhancement of the Madabhushi liquid jet in crossflow breakup model by a
  ligament breakup mechanism*. ILASS – Europe 2019, 29th Conference on Liquid
  Atomization and Spray Systems, Paris, France.
  — Defines the J = 10 and J = 20 LJICF cases and introduces the ligament
  correction implemented in the breakup model.

- **Faeth, G. M., Hsiang, L.-P., Wu, P.-K.** (1995).
  *Structure and breakup properties of sprays*.
  International Journal of Multiphase Flow, **21**(Suppl.), 99–127.
  DOI: [10.1016/0301-9322(95)00059-7](https://doi.org/10.1016/0301-9322(95)00059-7).
  — General review of spray structure and droplet-breakup properties.

---

## 2. Liquid jet in subsonic gas film — GLOP / Zhang cases

- **Zhang, M., Zhou, S., Song, L., Zhang, X., Lyu, J.-Y., Dong, F.** (2025).
  *Atomization dynamics of gas-liquid orifice type pintle injectors at
  different throttling levels*. Energy, **335**, 138192.
  DOI: [10.1016/j.energy.2025.138192](https://doi.org/10.1016/j.energy.2025.138192).
  — Primary reference for the LJIGF/GLOP configuration; provides the
  VOF-to-DPM/LES data used for the Zhang gas-film Weber-number validation cases.

- **Zhang, M. et al.**
  *Atomization dynamics of gas-liquid orifice type pintle injector under
  different operating conditions* — `metadata to verify` (journal, year,
  volume, pages, DOI).
  — Companion GLOP study used as additional background for the gas-film
  configuration.

---

## 3. Primary and secondary breakup modelling

- **Madabhushi, R. K.** (2003).
  *A model for numerical simulation of breakup of a liquid jet in crossflow*.
  Atomization and Sprays, **13**(4–5), 413–424. Begell House.
  DOI: [10.1615/AtomizSpr.v13.i45.10](https://doi.org/10.1615/AtomizSpr.v13.i45.10)
  (`metadata to verify`).
  — Base formulation of the LJICF column / Kelvin–Helmholtz breakup model.

- **Reitz, R. D.** (1987).
  *Modeling atomization processes in high-pressure vaporizing sprays*.
  Atomisation and Spray Technology, **3**, 309–337.
  — Kelvin–Helmholtz (KH) Wave atomization model underlying the
  wave-stripping stage.

- **Pilch, M., Erdman, C. A.** (1987).
  *Use of breakup time data and velocity history data to predict the maximum
  size of stable fragments for acceleration-induced breakup of a liquid drop*.
  International Journal of Multiphase Flow, **13**(6), 741–757.
  DOI: [10.1016/0301-9322(87)90063-2](https://doi.org/10.1016/0301-9322(87)90063-2).
  — Pilch–Erdman secondary-breakup criteria for post-column droplets.

- **Liu, A. B., Mather, D., Reitz, R. D.** (1993).
  *Modeling the effects of drop drag and breakup on fuel sprays*.
  SAE Technical Paper 930072.
  DOI: [10.4271/930072](https://doi.org/10.4271/930072).
  — Deformation-dependent drag treatment used for PE-active parcels.

---

## 4. Drag modelling

- **Putnam, A.** (1961).
  *Integrable form of droplet drag coefficient*.
  ARS Journal, **31**(10), 1467–1468.
  DOI: [10.2514/8.5825](https://doi.org/10.2514/8.5825).
  — Sphere-drag correlation used by the standard OpenFOAM drag law.

---

## 5. Source-term coupling and diffusion-based smoothing

- **Sun, R., Xiao, H.** (2015).
  *Diffusion-based coarse graining in hybrid continuum–discrete solvers:
  theoretical formulation and a priori tests*.
  International Journal of Multiphase Flow, **77**, 142–157.
  DOI: [10.1016/j.ijmultiphaseflow.2015.08.014](https://doi.org/10.1016/j.ijmultiphaseflow.2015.08.014).
  — Theoretical basis of the Diffusion-Based Method (DBM) implemented in the
  solver.

---

## 6. Turbulence and two-phase-flow modelling

- **Menter, F. R.** (1994).
  *Two-equation eddy-viscosity turbulence models for engineering applications*.
  AIAA Journal, **32**(8), 1598–1605.
  DOI: [10.2514/3.12149](https://doi.org/10.2514/3.12149).
  — k–ω SST turbulence closure used for the carrier phase.

- **Jiang, X., Siamas, G. A., Jagus, K., Karayiannis, T. G.** (2010).
  *Physical modelling and advanced simulations of gas–liquid two-phase jet
  flows in atomization and sprays*.
  Progress in Energy and Combustion Science, **36**(2), 131–167.
  DOI: [10.1016/j.pecs.2009.09.002](https://doi.org/10.1016/j.pecs.2009.09.002).
  — Review of high-fidelity (VOF/DNS/LES) atomization simulation approaches.

- **Nordin, P. A. N.** (1998).
  *Complex chemistry modeling of diesel spray combustion*.
  PhD thesis, Chalmers University of Technology, Gothenburg, Sweden.
  — Stochastic trajectory collision/coalescence model used in the spray cloud.

---

## 7. Evaporation and spray modelling (future extensions)

- **Schmehl, R., Klose, G., Maier, G., Wittig, S.** (1998).
  *Efficient numerical calculation of evaporating sprays in combustion chamber
  flows* — `metadata to verify` (proceedings, volume, page range, DOI).
  — Reference for the evaporation modelling strategy considered for future
  extensions of the solver.

---

## 8. Supersonic crossflow references (context)

- **Lin, K.-C., Kennedy, P. J., Jackson, T. A.** (2004).
  *Structures of water jets in a Mach 1.94 supersonic crossflow* —
  `metadata to verify` (AIAA paper number, proceedings, DOI).
  — Background reference for the supersonic crossflow extension discussed as
  future work; not a current validation target.

---

## Cross-references

- The solver implementation that builds on these references is documented in
  the top-level [`README.md`](../../README.md).
- The complete thesis bibliography (BibTeX) is in
  [`thesis/bibliografia.bib`](../../thesis/bibliografia.bib).

## Notes on metadata

Some entries are reconstructed from the thesis bibliography and from the
identifying information of the original sources. Fields explicitly flagged with
`metadata to verify` should be confirmed against the publisher record before
being used in formal citations.
