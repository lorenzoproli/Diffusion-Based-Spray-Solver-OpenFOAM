// Copyright (C) 2026 Lorenzo Proli, Politecnico di Torino
// SPDX-License-Identifier: GPL-3.0-or-later
// smoothSourceTerms.H
// ------------------------------------------------------------------
// Diffusion-Based Method (DBM) smoothing of Lagrangian source terms.
//
// Replaces the Particle Centroid Method (PCM) point-source projection
// with a Gaussian-equivalent spatial distribution, eliminating the
// "velocity locking" artefact in dense spray regions.
//
// Algorithm (per accumulated source field phi):
//   1. Build the volumetric density on each cell:
//        phi_PCM = phi_cloud / V
//   2. Solve the screened Poisson equation:
//        phi_s - sigma^2 * laplacian(phi_s) = phi_PCM
//      with zeroGradient on physical patches (Neumann -> conservation)
//      and the patch's own type on coupled patches (processor/cyclic).
//   3. Multiply back by V and store in the cloud.
//
// The single implicit step is mathematically equivalent to convolving
// phi_PCM with a Helmholtz/Yukawa kernel of characteristic length
// b = 2*sigma. For sigma << L_domain this kernel is well-approximated
// by a Gaussian of bandwidth b.
//
// Conservation:
//   With Neumann BCs on all physical patches, the screened Poisson
//   operator preserves the integral, i.e.
//        sum_i phi_s_i V_i = sum_i phi_PCM_i V_i
//   to within the linear-solver tolerance.
//
// Notes on the implementation:
//   * The RHS *must* be passed explicitly via "== phi_PCM"; without it
//     the matrix equation is homogeneous and the solver returns ~0,
//     which silently zeroes the Lagrangian source (one-way coupling).
//   * setReference() must NOT be used here: the operator (I - s2*L)
//     with Neumann BCs is non-singular thanks to the diagonal
//     contribution of fvm::Sp(1, .).
//   * The same wordList of patch types is reused for both the vector
//     field (UTrans) and the scalar field (UCoeff): zeroGradient and
//     coupled-patch type names are identical for both ranks.
//
// References:
//   [1] Sun & Xiao (2015), arXiv:1409.0001
//   [2] Capecelatro & Desjardins (2013), JCP 238, 1-31
//   [3] ANSYS Fluent Theory Guide,
//       "Gaussian kernel for mesh node averaging"
//   [4] Esgandari et al. (2025), J. Multiphase Flow 188, 105223
// ------------------------------------------------------------------

if (dbmActive)
{
    // -------------------------------------------------------------
    // Patch types:
    //   physical patches -> zeroGradient (Neumann, conservation)
    //   coupled patches  -> their own type (processor / cyclic / ...)
    // -------------------------------------------------------------
    wordList convergencePatchTypes(mesh.boundary().size(), "zeroGradient");
    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].coupled())
        {
            convergencePatchTypes[patchi] = mesh.boundary()[patchi].type();
        }
    }

    // =============================================================
    // 1. Smooth explicit momentum source: UTrans  [kg.m/s per cell]
    // =============================================================
    {
        // Volumetric density [kg/(m^2.s)]
        volVectorField convergenceSU
        (
            IOobject
            (
                "convergenceSU",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "zero",
                parcels.UTrans().dimensions()/dimVolume,
                Zero
            ),
            convergencePatchTypes
        );

        // Initialise from the PCM accumulated source
        forAll(convergenceSU, celli)
        {
            convergenceSU[celli] =
                parcels.UTrans()[celli]/mesh.V()[celli];
        }
        convergenceSU.correctBoundaryConditions();

        // RHS of the screened Poisson equation
        const volVectorField convergenceSU_PCM
        (
            "convergenceSU_PCM",
            convergenceSU
        );

        // Solve  (I - sigma^2 * laplacian) phi_s = phi_PCM
        fvVectorMatrix convergenceSUEqn
        (
            fvm::Sp
            (
                dimensionedScalar("one", dimless, 1.0),
                convergenceSU
            )
          - fvm::laplacian(dbmSigma2, convergenceSU)
         ==
            convergenceSU_PCM
        );

        convergenceSUEqn.solve(mesh.solverDict("convergenceSU"));

        // Write back: volumetric density -> per-cell accumulator
        forAll(parcels.UTrans(), celli)
        {
            parcels.UTrans()[celli] =
                convergenceSU[celli]*mesh.V()[celli];
        }
    }

    // =============================================================
    // 2. Smooth implicit drag coefficient: UCoeff  [kg per cell]
    // =============================================================
    {
        volScalarField convergenceSp
        (
            IOobject
            (
                "convergenceSp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "zero",
                parcels.UCoeff().dimensions()/dimVolume,
                0.0
            ),
            convergencePatchTypes
        );

        forAll(convergenceSp, celli)
        {
            convergenceSp[celli] =
                parcels.UCoeff()[celli]/mesh.V()[celli];
        }
        convergenceSp.correctBoundaryConditions();

        const volScalarField convergenceSp_PCM
        (
            "convergenceSp_PCM",
            convergenceSp
        );

        fvScalarMatrix convergenceSpEqn
        (
            fvm::Sp
            (
                dimensionedScalar("one", dimless, 1.0),
                convergenceSp
            )
          - fvm::laplacian(dbmSigma2, convergenceSp)
         ==
            convergenceSp_PCM
        );

        convergenceSpEqn.solve(mesh.solverDict("convergenceSp"));

        // Defensive non-negativity floor: UCoeff is a drag coefficient
        // and must be >= 0 to keep fvm::Sp(UCoeff/V, U) stable in UEqn.
        // The screened Poisson operator with f >= 0 and Neumann BCs
        // preserves positivity (max principle); the clamp guards
        // against tiny negative residuals from the linear solver.
        forAll(parcels.UCoeff(), celli)
        {
            parcels.UCoeff()[celli] =
                max(scalar(0), convergenceSp[celli]*mesh.V()[celli]);
        }
    }

    // -------------------------------------------------------------
    // Optional sanity check (uncomment for debugging):
    //   the L1 norm of the smoothed field should match the PCM one
    //   to within solver tolerance (global momentum conservation).
    // -------------------------------------------------------------
    // Info<< "DBM check: ||UTrans||_1 = "
    //     << gSum(mag(parcels.UTrans())()) << endl;
}

// ************************************************************************* //
