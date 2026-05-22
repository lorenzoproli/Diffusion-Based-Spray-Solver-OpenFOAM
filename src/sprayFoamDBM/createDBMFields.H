// createDBMFields.H
// ------------------------------------------------------------------
// Read DBM (Diffusion-Based Method) smoothing parameters.
//
// Dictionary: constant/smoothingProperties
//
// Parameters:
//   smoothBandwidth  [m]   Gaussian-equivalent kernel bandwidth b.
//                          Effective sigma^2 = b^2/4 controls the spread
//                          of the source terms. Typical values:
//                          1-6 x D_inj (injector diameter).
//                          Equivalent to Fluent "Gaussian Factor" x dx.
//                          Set to 0 to disable smoothing (pure PCM).
//
// References:
//   [1] Sun & Xiao (2015), arXiv:1409.0001
//   [2] Capecelatro & Desjardins (2013), JCP 238, 1-31
//   [3] ANSYS Fluent Theory Guide,
//       "Gaussian kernel for mesh node averaging"
// ------------------------------------------------------------------

IOdictionary smoothingProperties
(
    IOobject
    (
        "smoothingProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

// Bandwidth b [m] - default 0 = disabled (fall back to pure PCM)
const dimensionedScalar smoothBandwidth
(
    "smoothBandwidth",
    dimLength,
    smoothingProperties.getOrDefault<scalar>("smoothBandwidth", 0.0)
);

// Diffusion coefficient sigma^2 = b^2/4  [m^2]
const dimensionedScalar dbmSigma2
(
    "dbmSigma2",
    sqr(smoothBandwidth)/scalar(4)
);

// Whether smoothing is active
const bool dbmActive = (smoothBandwidth.value() > SMALL);

if (dbmActive)
{
    Info<< "DBM source-term smoothing ACTIVE" << nl
        << "    bandwidth b  = " << smoothBandwidth.value() << " m" << nl
        << "    sigma^2      = " << dbmSigma2.value()       << " m^2" << nl
        << endl;
}
else
{
    Info<< "DBM source-term smoothing DISABLED (bandwidth = 0)" << nl
        << "    Running with standard PCM (no kernel smoothing)." << nl
        << endl;
}

// ************************************************************************* //
