// createDBMFields.H
// ------------------------------------------------------------------
// Read DBM (Diffusion-Based Method) smoothing parameters.
//
// Dictionary: constant/smoothingProperties
//
// Parameters:
//   smoothBandwidth  [m]   Gaussian kernel bandwidth b.
//                          The effective σ² = b²/4 controls how far
//                          the source terms are spread spatially.
//                          Typical values: 1–6 × Dinj (injector diameter).
//                          Equivalent to Fluent "Gaussian Factor" × Δx.
//                          Set to 0 to disable smoothing (pure PCM).
//
// Reference: Sun & Xiao (2015), arXiv:1409.0001
//            Fluent Theory Guide: "Gaussian kernel for mesh node averaging"
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

// Bandwidth [m] — default 0 = disabled (fallback to pure PCM)
const dimensionedScalar smoothBandwidth
(
    "smoothBandwidth",
    dimLength,
    smoothingProperties.getOrDefault<scalar>("smoothBandwidth", 0.0)
);

// Derived: σ² = b²/4  [m²]
const dimensionedScalar dbmSigma2
(
    "dbmSigma2",
    dimArea,
    sqr(smoothBandwidth.value()) / 4.0
);

// Whether smoothing is active
const bool dbmActive = (smoothBandwidth.value() > SMALL);

if (dbmActive)
{
    Info<< "DBM source-term smoothing ACTIVE" << nl
        << "    bandwidth b  = " << smoothBandwidth.value() << " m" << nl
        << "    sigma^2      = " << dbmSigma2.value() << " m^2" << nl
        << endl;
}
else
{
    Info<< "DBM source-term smoothing DISABLED (bandwidth = 0)" << nl
        << "    Running with standard PCM (no kernel smoothing)." << nl
        << endl;
}

// ************************************************************************* //
