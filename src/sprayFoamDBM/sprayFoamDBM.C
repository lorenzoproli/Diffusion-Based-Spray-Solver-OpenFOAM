/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    sprayFoamDBM - sprayFoam with Diffusion-Based source-term smoothing.

    Identical to standard sprayFoam except that, after parcels.evolve(),
    the accumulated Lagrangian momentum source terms (UTrans, UCoeff)
    are smoothed via a single implicit screened-Poisson solve (DBM).
    This replicates the Gaussian kernel smoothing available in ANSYS
    Fluent DPM ("Gaussian Factor" / node-based averaging) and removes
    the "velocity locking" artefact of the Particle Centroid Method
    (PCM) in dense regions. For LJICF (pintle) configurations this
    restores the correct crossflow deflection of the spray.

    NOTE: only the momentum sources (UTrans, UCoeff) are smoothed.
    Mass (rhoTrans), enthalpy (hsTrans, hsCoeff) and species sources
    are left at their PCM (per-cell) values. This is acceptable for
    cold or weakly-evaporating LJICF; for strongly reactive cases the
    same smoothing should be applied to all coupling terms for global
    consistency.

    References:
      [1] Sun & Xiao (2015), "Diffusion-based coarse graining in hybrid
          continuum-discrete solvers", arXiv:1409.0001
      [2] Capecelatro & Desjardins (2013), JCP 238, 1-31
      [3] ANSYS Fluent Theory Guide
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
        " smoothing (DBM Gaussian-equivalent kernel)."
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

    // DBM smoothing parameters (constant/smoothingProperties)
    #include "createDBMFields.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // ---------------------------------------------------------------
        // Lagrangian cloud evolution
        // (PCM source terms accumulated into UTrans, UCoeff, ...)
        // ---------------------------------------------------------------
        parcels.evolve();

        // ---------------------------------------------------------------
        // DBM: smooth UTrans and UCoeff before they enter UEqn.
        // Equivalent to Fluent's Gaussian node-based averaging.
        // ---------------------------------------------------------------
        #include "smoothSourceTerms.H"

        // ---------------------------------------------------------------
        // Eulerian phase (PIMPLE)
        // ---------------------------------------------------------------
        if (pimple.solveFlow())
        {
            #include "rhoEqn.H"

            while (pimple.loop())
            {
                #include "UEqn.H"
                #include "YEqn.H"
                #include "EEqn.H"

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
