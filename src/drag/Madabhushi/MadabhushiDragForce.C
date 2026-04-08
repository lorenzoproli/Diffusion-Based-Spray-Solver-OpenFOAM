/*------------------------------------------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Madabhushi custom drag force for jet-in-crossflow breakup

    Physical model backbone:
      - Liquid-column / blob drag with constant CdBlob_ = 1.48
        following Madabhushi (2003), jet-regime description
      - Column-breakup time from Madabhushi (2003), Eq. (1)
      - Post-column deformation and breakup timing based on
        Pilch & Erdman correlations
      - Deformed-droplet reference diameter and staged drag
        treatment consistent with Schmehl et al. (1998) and
        Lambert et al. (2019)

    Regimes implemented:
      - user_ = 0:
          * intact liquid column (blob drag)
          * core parcel in PE deformation / disc drag after column breakup
      - user_ = 1 or 2:
          * spherical drag before STD_PE latch
          * PE deformation / disc drag after STD_PE latch

    Notes:
      - Blob regime uses the explicit vector force formulation of
        Lambert et al. (2019), Fluent Theory Guide Eq. (12-452):
          F_D = 0.5 * rhoG * CdBlob * A_ref * |Urel| * Urel
        where A_ref = pi/4 * Dinj^2 is the projected area of the intact
        liquid column.
        The force is split into Su + Sp*(Uc-U) to keep tMom finite:
          Sp  = 0.5 * rhoG * CdBlob * Aref * |Urel|          [scalar, N·s/m]
          Su  = F_blob - Sp*(Uc-U)  = Zero  (identically)
        Because Urel = Uc-U, F_blob = Sp*(Uc-U) exactly — so Su = Zero
        and the full physics is captured in Sp alone, while tMom =
        mass/Sp remains finite and well-conditioned.
      - PE-active parcels use deformation / disc drag law
      - UgRef_ is read from the model coefficients dictionary and used
        consistently for tColumnBreakup evaluation

    Modifications (Thesis — Madabhushi crossflow model):
      FIX-DR1: weCritPE in the PE drag branch now uses the viscosity-corrected
               formula We_crit = 12*(1 + 1.077*Oh^1.6) consistent with
               Schmehl et al. (1998) Eq.(38), Lambert et al. (2019) onset
               criterion, and the identical formula already used in
               Madabhushi.C PATH 1 and PATH 4.
               [Schmehl et al. (1998), Eq. (38);
                Lambert et al. (2019), PE onset criterion;
                Pilch & Erdman (1987), Eq. (5)]

      FIX-DR2: Debug output restructured as CSV (dragDebug.log).
               [All per-tag counter limits removed; disable debug in production]

      FIX-DR3: weExcess for tb/t* correlations uses (We - 12) as per
               Pilch & Erdman (1987) original calibration and Madabhushi
               (2003) Eq. (5).
               [Pilch & Erdman (1987), Eqs. (8)-(12);
                Madabhushi (2003), Eq. (5)]

      FIX-DR4: BLOB regime now implements Lambert et al. (2019) Fluent
               Theory Guide Eq. (12-452) SCALED FOR PARCEL MASS:
                 F_D = 0.5 * rhoG * CdBlob * (pi/4 * Dinj^2) * |Urel| * Urel
               using the factored forceSuSp form:
                 Sp = nParticle_equivalent * 0.5 * rhoG * CdBlob * Aref * |Urel|
                 Su = Zero
               Since Urel = Uc - U, the product Sp*(Uc-U) reproduces F_D
               exactly. This form:
               (a) scales with Dinj^2 as prescribed, not with d_curr^3;
               (b) keeps tMom = mass/Sp finite and numerically stable;
               (c) prevents two-way coupling lock-up by scaling the force
                   to the actual mass fraction of the computational parcel.
               [Lambert et al. (2019), Fluent Theory Guide Eq. (12-452);
                Madabhushi (2003), jet-regime drag description;
                Wu et al. (1997), column aerodynamic force balance, Eq. (1)]

      NOTE-1:  Deformation ramp for We >= 100 uses (1 + 1.9*frac)*Dparent0,
               i.e. Madabhushi (2003) Eq. (6) evaluated at We = 100.
               [Madabhushi (2003), Eq. (6); Lambert et al. (2019), Fluent doc]

      NOTE-2:  tb/t* piecewise correlations use strict-less-than (<) at
               interval boundaries. Physical impact negligible.
               [Pilch & Erdman (1987), Eqs. [8]-[12]]

\*------------------------------------------------------------------------------------------------------------*/

#include "MadabhushiDragForce.H"
#include "mathematicalConstants.H"
#include "OFstream.H"
#include "autoPtr.H"

namespace
{
    // -------------------------------------------------------------------------
    // FIX-DR2: CSV drag log
    // -------------------------------------------------------------------------
    Foam::OFstream& dragLogFile()
    {
        static Foam::autoPtr<Foam::OFstream> logPtr;
        if (!logPtr.valid())
        {
            logPtr.reset(new Foam::OFstream("dragDebug.log"));
            logPtr() << "phase,"
                     << "tc,user,ms,tCb,tSecStart,"
                     << "d,dForceEff,dRefEff,"
                     << "Re,ReEff,Cd,dragScale,"
                     << "WeInitPE,tElapsed,tDef,tb"
                     << Foam::endl;
        }
        return logPtr();
    }

    inline void dragDebugCSV
    (
        const bool debug,
        const Foam::word& phase,
        const Foam::scalar tc,
        const Foam::scalar userFlag,
        const Foam::scalar msState,
        const Foam::scalar tCb,
        const Foam::scalar tSecStart,
        const Foam::scalar d,
        const Foam::scalar dForceEff,
        const Foam::scalar dRefEff,
        const Foam::scalar Re,
        const Foam::scalar ReEff,
        const Foam::scalar Cd,
        const Foam::scalar dragScale,
        const Foam::scalar WeInitPE  = 0.0,
        const Foam::scalar tElapsed  = 0.0,
        const Foam::scalar tDef      = 0.0,
        const Foam::scalar tb        = 0.0
    )
    {
        if (!debug) return;
        dragLogFile()
            << phase      << ","
            << tc         << ","
            << userFlag   << ","
            << msState    << ","
            << tCb        << ","
            << tSecStart  << ","
            << d          << ","
            << dForceEff  << ","
            << dRefEff    << ","
            << Re         << ","
            << ReEff      << ","
            << Cd         << ","
            << dragScale  << ","
            << WeInitPE   << ","
            << tElapsed   << ","
            << tDef       << ","
            << tb
            << Foam::endl;
    }
}


template<class CloudType>
Foam::MadabhushiDragForce<CloudType>::MadabhushiDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    CdBlob_(this->coeffs().getOrDefault("CdBlob", 1.48)),
    CdDisc_(this->coeffs().getOrDefault("CdDisc", 1.2)),
    C0_(this->coeffs().getOrDefault("C0", 3.44)),
    Dinj_(this->coeffs().getOrDefault("Dinj", 0.0016)),
    sigma_(this->coeffs().getOrDefault("sigma", 0.072)),
    UgRef_(this->coeffs().getOrDefault("UgRef", 62.5)),
    rhoRef_(this->coeffs().getOrDefault("rhoRef", 1.186)),
    debug_(this->coeffs().getOrDefault("debug", false))
{}


template<class CloudType>
Foam::MadabhushiDragForce<CloudType>::MadabhushiDragForce
(
    const MadabhushiDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    CdBlob_(df.CdBlob_),
    CdDisc_(df.CdDisc_),
    C0_(df.C0_),
    Dinj_(df.Dinj_),
    sigma_(df.sigma_),
    UgRef_(df.UgRef_),
    rhoRef_(df.rhoRef_),
    debug_(df.debug_)
{}


template<class CloudType>
Foam::MadabhushiDragForce<CloudType>::~MadabhushiDragForce()
{}


template<class CloudType>
Foam::forceSuSp Foam::MadabhushiDragForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    // The full drag force is computed inside calcCoupled.
    // Returning Zero here prevents double-counting.
    // [OpenFOAM ParticleForce API convention]
    return forceSuSp(Zero, 0.0);
}


template<class CloudType>
Foam::forceSuSp Foam::MadabhushiDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    const scalar tc    = p.age();
    const scalar rhoL  = p.rho();
    const scalar rhoG  = td.rhoc();
    const scalar dCurr = max(p.d(), 1.0e-12);

    const scalar sigmaLocal = max(p.sigma(), VSMALL);

    // ------------------------------------------------------------------
    // Persistent parcel state written by the breakup model
    // ------------------------------------------------------------------
    const scalar userFlag  = p.user();
    const scalar msState   = p.ms();
    const scalar tSecStart = p.y();
    const scalar Dparent0  = max(p.yDot(),    1.0e-12);
    const scalar UrelPE0   = max(p.KHindex(), 1.0e-12);

    // ------------------------------------------------------------------
    // Column-breakup time.
    // [Madabhushi (2003), Eq. (1); Lambert et al. (2019), Eq. (1)]
    // ------------------------------------------------------------------
    const scalar tColumnBreakup =
        C0_*(Dinj_/(UgRef_ + VSMALL))*Foam::sqrt(rhoL/(rhoRef_ + VSMALL));

    // ------------------------------------------------------------------
    // PE-active child check.
    // ------------------------------------------------------------------
    const bool isStdPEChild =
    (
        userFlag  > 0.5
     && msState   > 1.5
     && msState   < 2.5
     && tSecStart > 0.0
     && p.KHindex() > 0.0
    );

    // ------------------------------------------------------------------
    // Child parcels not yet in standard PE: spherical drag.
    // ------------------------------------------------------------------
    if (userFlag > 0.5 && !isStdPEChild)
    {
        const scalar CdEffChild = CdSphere(Re);

        dragDebugCSV
        (
            debug_, "CHILD",
            tc, userFlag, msState, tColumnBreakup, tSecStart,
            dCurr, dCurr, dCurr,
            Re, Re, CdEffChild, 1.0
        );

        return forceSuSp
        (
            Zero,
            mass*0.75*CdEffChild*Re*muc/(rhoL*sqr(dCurr))
        );
    }

    // ------------------------------------------------------------------
    // BLOB / intact liquid-column regime (core parcel only).
    //
    // FIX-DR4: implement Lambert et al. (2019) Eq. (12-452):
    //
    //   F_D = 0.5 * rhoG * CdBlob * Aref * |Urel| * Urel
    //
    // with Aref = pi/4 * Dinj^2.
    //
    // Since Urel = Uc - U, this is factored as:
    //
    //   F_D = Sp * (Uc - U)   with   Sp = 0.5 * rhoG * CdBlob * Aref * |Urel|
    //
    // The forceSuSp API then gives: Su = Zero, Sp = scalar above.
    // This is MATHEMATICALLY IDENTICAL to the explicit vector form, and:
    //   (a) scales with Dinj^2 as prescribed by Lambert, not d_curr^3;
    //   (b) keeps tMom = mass/Sp finite — no FPE crash;
    //   (c) prevents two-way coupling lock-up by scaling the force
    //       to the actual mass fraction of the computational parcel.
    //
    // [Lambert et al. (2019), Fluent Theory Guide Eq. (12-452);
    //  Madabhushi (2003), jet-regime drag;
    //  Wu et al. (1997), Eq. (1)]
    // ------------------------------------------------------------------
    if (userFlag < 0.5 && tSecStart <= 0.0 && tc < tColumnBreakup)
    {
        // |Urel| = magnitude of relative velocity (gas - parcel)
        const scalar Urmag_local = mag(td.Uc() - p.U()) + VSMALL;

        // Calcolo il volume e la massa di un singolo blob fisico ideale
        const scalar volume_blob = constant::mathematical::pi / 6.0 * pow3(Dinj_);
        const scalar mass_blob = rhoL * volume_blob;
        
        // Frazione di blob fisico contenuta in questo parcel computazionale
        // (equivale a nParticle se d_parcel fosse forzato a Dinj)
        const scalar nParticle_equivalent = mass / (mass_blob + VSMALL);

        // Projected area of the intact liquid column: A_ref = pi/4 * Dinj^2
        // [Lambert et al. (2019), Fluent Theory Guide Eq. (12-453)]
        const scalar Aref =
            constant::mathematical::pi / 4.0 * sqr(Dinj_);

        // Sp coefficient: F_D = Sp * (Uc - U)
        // [Lambert et al. (2019), Eq. (12-452)]
        const scalar SpBlob =
            nParticle_equivalent * 0.5 * rhoG * CdBlob_ * Aref * Urmag_local;

        // Effective Re with Dinj for debug log only
        const scalar ReEff_log =
            rhoG * Urmag_local * Dinj_ / (muc + VSMALL);

        dragDebugCSV
        (
            debug_, "BLOB",
            tc, userFlag, msState, tColumnBreakup, tSecStart,
            dCurr, Dinj_, Dinj_,
            Re, ReEff_log, CdBlob_, 1.0
        );

        // Su = Zero, Sp = SpBlob  →  force = SpBlob*(Uc-U) = F_D exactly
        // tMom = mass/SpBlob  →  finite and physically meaningful
        return forceSuSp(Zero, SpBlob);
    }

    // ------------------------------------------------------------------
    // Pilch-Erdman deformation / disc-drag stage.
    // FIX-DR1: viscosity-corrected weCritPE.
    // [Schmehl et al. (1998), Eq. (38); Lambert et al. (2019);
    //  Pilch & Erdman (1987), Eq. (5)]
    // ------------------------------------------------------------------
    if (tSecStart > 0.0)
    {
        const scalar WeInitPE =
            rhoG*sqr(UrelPE0)*Dparent0/(sigmaLocal + VSMALL);

        const scalar ohPE0 =
            rhoL > VSMALL
          ? p.mu()/(Foam::sqrt(rhoL*Dparent0*sigmaLocal) + VSMALL)
          : 0.0;

        const scalar weCritPE = 12.0*(1.0 + 1.077*pow(ohPE0, 1.6));

        scalar CdEff     = CdSphere(Re);
        scalar dRefEff   = dCurr;
        scalar dragScale = 1.0;
        scalar tElapsedLog = 0.0;
        scalar tDefLog     = 0.0;
        scalar tbLog       = 0.0;
        Foam::word phaseTag = "PE_POST";

        if (WeInitPE > weCritPE)
        {
            const scalar tElapsed = max(tc - tSecStart, 0.0);
            tElapsedLog = tElapsed;

            // t* [Pilch & Erdman (1987), Eq. (3); Lambert et al. (2019), Eq. (3)]
            const scalar tStar =
                (Dparent0/(UrelPE0 + VSMALL))
               *Foam::sqrt(rhoL/(rhoG + VSMALL));

            // tDef = 1.6 t* [Schmehl et al. (1998), Eq. (42)]
            const scalar tDef = 1.6*tStar;
            tDefLog = tDef;

            // FIX-DR3: (We-12) for P&E correlations.
            // [Pilch & Erdman (1987), Eqs. (8)-(12); Madabhushi (2003), Eq. (5)]
            const scalar weExcess = max(WeInitPE - 12.0, VSMALL);

            scalar tbOverTstar = 5.5;
            if      (WeInitPE < 18.0)   tbOverTstar = 6.0  *pow(weExcess, -0.25);
            else if (WeInitPE < 45.0)   tbOverTstar = 2.45 *pow(weExcess,  0.25);
            else if (WeInitPE < 351.0)  tbOverTstar = 14.1 *pow(weExcess, -0.25);
            else if (WeInitPE < 2670.0) tbOverTstar = 0.766*pow(weExcess,  0.25);

            const scalar tb = tbOverTstar*tStar;
            tbLog = tb;

            if (tElapsed < tDef)
            {
                // Deformation stage: linear Cd ramp sphere->disc.
                // [Lambert et al. (2019); Schmehl et al. (1998), after Eq. (46)]
                const scalar frac = tElapsed/(tDef + VSMALL);
                CdEff = CdSphere(Re)*(1.0 - frac) + CdDisc_*frac;

                // Dref ramp. NOTE-1: capped at We=100.
                // [Madabhushi (2003), Eq. (6); Schmehl et al. (1998), Eq. (45)]
                dRefEff = (WeInitPE < 100.0)
                  ? (1.0 + 0.19*Foam::sqrt(WeInitPE)*frac)*Dparent0
                  : (1.0 + 1.9*frac)*Dparent0;

                phaseTag = "PE_DEFORM";
            }
            else if (tElapsed < tb)
            {
                // Disc stage. [Lambert et al. (2019); Schmehl et al. (1998)]
                CdEff = CdDisc_;
                dRefEff = (WeInitPE < 100.0)
                  ? (1.0 + 0.19*Foam::sqrt(WeInitPE))*Dparent0
                  : 2.9*Dparent0;

                phaseTag = "PE_DISC";
            }
            else
            {
                // Post-breakup: sphere.
                CdEff   = CdSphere(Re);
                dRefEff = dCurr;
                phaseTag = "PE_POST";
            }
        }
        else
        {
            CdEff   = CdSphere(Re);
            dRefEff = dCurr;
            phaseTag = "PE_SUBCRIT";
        }

        // Area correction: dragScale = (Dref/Dcurr)^2
        // [Schmehl et al. (1998), Eqs. (45)-(46); Lambert et al. (2019), Eq. (9)]
        dragScale = sqr(max(dRefEff, 1.0e-12)/dCurr);

        dragDebugCSV
        (
            debug_, phaseTag,
            tc, userFlag, msState, tColumnBreakup, tSecStart,
            dCurr, dCurr, dRefEff,
            Re, Re, CdEff, dragScale,
            WeInitPE, tElapsedLog, tDefLog, tbLog
        );

        return forceSuSp
        (
            Zero,
            mass*0.75*CdEff*dragScale*Re*muc/(rhoL*sqr(dCurr))
        );
    }

    // ------------------------------------------------------------------
    // Fallback: spherical drag. Should not be reached in normal operation.
    // ------------------------------------------------------------------
    {
        const scalar CdFallback = CdSphere(Re);

        dragDebugCSV
        (
            debug_, "FALLBACK",
            tc, userFlag, msState, tColumnBreakup, tSecStart,
            dCurr, dCurr, dCurr,
            Re, Re, CdFallback, 1.0
        );

        return forceSuSp
        (
            Zero,
            mass*0.75*CdFallback*Re*muc/(rhoL*sqr(dCurr))
        );
    }
}
