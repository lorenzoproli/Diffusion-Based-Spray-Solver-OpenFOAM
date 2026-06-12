/*------------------------------------------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Lorenzo Proli, Politecnico di Torino
    SPDX-License-Identifier: GPL-3.0-or-later
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
      - All regimes use the standard DPM semi-implicit Sp formulation:
          Sp = mass * (3/4) * Cd * Re * mu / (rhoL * d^2)
        which corresponds to the particle acceleration equation from
        KIVA [Reitz (1987), Eq. (20)]:
          a = (3/4) * Cd * rhoG * |Urel| / (rhoL * d) * (Ug - Up)
        The blob regime uses CdBlob (constant) in place of CdSphere(Re).
        [Madabhushi (2003): "The drag coefficient, CD, on the droplet is
         set to 1.48"; Reitz (1987), Eq. (20); Fluent Theory Guide
         Eq. (12-447): "drag remains constant CD = 1.48"]
      - PE-active parcels use deformation / disc drag law with area
        correction dragScale = (Dref/d)^2 for the deformed cross-section.
        [Schmehl et al. (1998), Eqs. (45)-(46); Lambert et al. (2019), Eq. (9)]
      - tColumnBreakup uses local gas velocity mag(td.Uc()) and local
        gas density td.rhoc(), consistent with the breakup model

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
    CdDisc_(1.2),
    C0_(this->coeffs().getOrDefault("C0", 3.44)),
    Dinj_(this->coeffs().getOrDefault("Dinj", 0.0016)),
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
    // Column-breakup time — uses local gas velocity and density.
    // [Madabhushi (2003), Eq. (1); Lambert et al. (2019), Eq. (1)]
    // ------------------------------------------------------------------
    const scalar UgMag = mag(td.Uc());
    const scalar tColumnBreakup =
        C0_*(Dinj_/(UgMag + VSMALL))*Foam::sqrt(rhoL/(rhoG + VSMALL));

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
    // Standard DPM drag with constant CdBlob in place of CdSphere(Re).
    // This is the formulation used in KIVA [Reitz (1987), Eq. (20)] and
    // in Fluent DPM, within which Madabhushi calibrated CdBlob = 1.48.
    //
    // Sp = mass * (3/4) * CdBlob * Re * mu / (rhoL * d^2)
    //    = mass * (3/4) * CdBlob * rhoG * |Urel| / (rhoL * d)
    //
    // Acceleration: a = (3/4) * CdBlob * rhoG * |Urel| / (rhoL * d)
    // As wave stripping reduces d while nParticle stays constant,
    // a scales as 1/d — the blob accelerates as it becomes lighter.
    //
    // [Madabhushi (2003): "The drag coefficient, CD, on the droplet
    //  is set to 1.48 to simulate the motion of liquid jet in crossflow";
    //  Reitz (1987), Eq. (20); Fluent Theory Guide Eq. (12-447)]
    // ------------------------------------------------------------------
    if (userFlag < 0.5 && tSecStart <= 0.0 && tc < tColumnBreakup)
    {
        dragDebugCSV
        (
            debug_, "BLOB",
            tc, userFlag, msState, tColumnBreakup, tSecStart,
            dCurr, dCurr, dCurr,
            Re, Re, CdBlob_, 1.0
        );

        return forceSuSp
        (
            Zero,
            mass*0.75*CdBlob_*Re*muc/(rhoL*sqr(dCurr))
        );
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
