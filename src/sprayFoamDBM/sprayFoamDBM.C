/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    sprayFoamDBM — sprayFoam with Diffusion-Based source term smoothing.

    Identical to standard sprayFoam except that after parcels.evolve()
    the accumulated Lagrangian momentum source terms (UTrans, UCoeff)
    are smoothed via an implicit Laplacian solve (DBM), replicating the
    Gaussian kernel smoothing available in ANSYS Fluent DPM.

    This eliminates the "velocity locking" artefact of the Particle
    Centroid Method (PCM), restoring the correct crossflow deflection
    for liquid-jet-in-crossflow simulations.

    [Sun & Xiao (2015), "Diffusion-based coarse graining in hybrid
     continuum-discrete solvers", arXiv:1409.0001]
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "basicSprayCloud.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "fvOptions.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for compressible, laminar or turbulent"
        " reacting and spraying flow with Lagrangian source-term"
        " smoothing (DBM Gaussian kernel)."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    // DBM smoothing: read parameters and create work fields
    #include "createDBMFields.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // ---------------------------------------------------------------
        // Lagrangian cloud evolution (source terms accumulated via PCM)
        // ---------------------------------------------------------------
        parcels.evolve();

        // ---------------------------------------------------------------
        // DBM: smooth source terms BEFORE the PIMPLE loop uses them.
        // This is the Gaussian kernel equivalent.
        // [Sun & Xiao (2015); Fluent: "Gaussian kernel averaging"]
        // ---------------------------------------------------------------
        #include "smoothSourceTerms.H"

        // ---------------------------------------------------------------
        // Eulerian phase solution
        // ---------------------------------------------------------------
        if (pimple.solveFlow())
        {
            #include "rhoEqn.H"

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                #include "UEqn.H"
                #include "YEqn.H"
                #include "EEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            rho = thermo.rho();

            if (runTime.write())
            {
                combustion->Qdot()().write();
            }
        }
        else
        {
            if (runTime.writeTime())
            {
                parcels.write();
            }
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
