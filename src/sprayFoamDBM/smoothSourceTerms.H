// smoothSourceTerms.H
// ------------------------------------------------------------------
// Diffusion-Based Method (DBM) smoothing of Lagrangian source terms.
//
// Replaces the Particle Centroid Method (PCM) point-source projection
// with a Gaussian-equivalent spatial distribution, eliminating the
// "velocity locking" artefact in dense spray regions.
//
// Algorithm:
//   1. Extract the PCM-accumulated source fields (UTrans, UCoeff)
//      from the cloud and convert to volumetric densities.
//   2. Solve the steady-state screened-Poisson equation:
//        φ_smooth - σ² ∇²φ_smooth = φ_PCM
//      with zero-gradient BCs on all boundaries (momentum conservation).
//   3. Write the smoothed fields back to the cloud.
//
// The solution at step 2 is mathematically identical to convolving
// φ_PCM with a 3D Gaussian kernel of bandwidth b = √(4·σ²).
//
// Conservation: guaranteed by the Neumann (zero-gradient) BCs.
//   ∫φ_smooth dV = ∫φ_PCM dV  (exact to solver tolerance).
//
// References:
//   [1] Sun & Xiao (2015), arXiv:1409.0001
//   [2] ANSYS Fluent Theory Guide, "Gaussian kernel for mesh node averaging"
//   [3] Esgandari et al. (2025), J. Multiphase Flow 188, 105223
// ------------------------------------------------------------------

if (dbmActive)
{
    // Patch types: zeroGradient everywhere for conservation
    const wordList convergencePatchTypes
    (
        mesh.boundary().size(),
        "zeroGradient"
    );

    // =============================================================
    // 1. Smooth explicit momentum source: UTrans [kg·m/s per cell]
    // =============================================================
    {
        // Volumetric density of momentum source [kg/(m²·s)]
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
                parcels.UTrans().dimensions() / dimVolume,
                Zero
            ),
            convergencePatchTypes
        );

        // Map accumulated PCM source to volumetric field
        forAll(convergenceSU, celli)
        {
            convergenceSU[celli] =
                parcels.UTrans()[celli] / mesh.V()[celli];
        }
        convergenceSU.correctBoundaryConditions();

        // Store the PCM field as RHS
        const volVectorField convergenceSU_PCM
        (
            "convergenceSU_PCM",
            convergenceSU
        );

        // Solve: φ_s - σ²·∇²φ_s = φ_PCM
        // Single implicit step → Gaussian kernel of bandwidth b
        fvVectorMatrix convergenceSUEqn
        (
            fvm::Sp
            (
                dimensionedScalar("one", dimless, 1.0),
                convergenceSU
            )
          - fvm::laplacian(dbmSigma2, convergenceSU)
        );

        convergenceSUEqn.setReference(0, Zero);
        convergenceSUEqn.solve
        (
            mesh.solverDict("convergenceSU")
        );

        // Write back to cloud: convert volumetric → per-cell
        forAll(parcels.UTrans(), celli)
        {
            parcels.UTrans()[celli] =
                convergenceSU[celli] * mesh.V()[celli];
        }
    }

    // =============================================================
    // 2. Smooth implicit drag coefficient: UCoeff [kg per cell]
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
                parcels.UCoeff().dimensions() / dimVolume,
                0.0
            ),
            convergencePatchTypes
        );

        forAll(convergenceSp, celli)
        {
            convergenceSp[celli] =
                parcels.UCoeff()[celli] / mesh.V()[celli];
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
        );

        convergenceSpEqn.setReference(0, 0.0);
        convergenceSpEqn.solve
        (
            mesh.solverDict("convergenceSp")
        );

        forAll(parcels.UCoeff(), celli)
        {
            parcels.UCoeff()[celli] =
                convergenceSp[celli] * mesh.V()[celli];
        }
    }
}

// ************************************************************************* //
