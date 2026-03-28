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
      - Blob regime uses CdBlob_ with Dinj_ as effective force diameter
      - Blob regime uses ReEff consistent with Dinj_
      - PE-active parcels use deformation / disc drag law
      - UgRef_ is read from the model coefficients dictionary and used
        consistently for tColumnBreakup evaluation

    Modifications (Thesis — Madabhushi crossflow model):
      FIX-DR1: weCritPE in the PE drag branch now uses the viscosity-corrected
               formula We_crit = 12*(1 + 1.077*Oh^1.6) consistent with
               Schmehl et al. (1998) Eq.(38), Lambert et al. (2019) onset
               criterion, and the identical formula already used in
               Madabhushi.C PATH 1 and PATH 4.
               Previously the drag used the bare value 12.0, activating
               deformed-disc drag for parcels that physically would not
               deform (viscous droplets with Oh>0), causing an over-estimate
               of the drag force on those parcels.

      FIX-DR2: Debug output restructured as CSV (dragDebug.log).
               One line per drag evaluation, separated by phase tag:
               BLOB, PE_DEFORM, PE_DISC, PE_POST, CHILD, FALLBACK.
               Readable with: awk -F, '$1=="BLOB"' dragDebug.log
               All per-tag counter limits removed: when debug=true every
               evaluation is written. Disable debug in production runs.

      FIX-DR3: weExcess for tb/t* correlations now uses (We - 12) as per
               Pilch & Erdman (1987) original calibration and Madabhushi
               (2003) Eq. (5), rather than (We - weCritPE). The viscosity-
               corrected weCritPE is retained for the admission condition
               only.
               [Pilch & Erdman (1987), Eqs. (8)-(12);
                Madabhushi (2003), Eq. (5)]

      NOTE-1:  Deformation ramp for We >= 100 uses (1 + 1.9*frac)*Dparent0,
               which equals Madabhushi (2003) Eq. (6) evaluated at We = 100
               as explicitly prescribed: "For We > 100, Eq. (6) is evaluated
               at We = 100", giving 0.19*sqrt(100)=1.9. This is also
               confirmed by the official Lambert/Fluent documentation figure.
               [Madabhushi (2003), Eq. (6); Lambert et al. (2019), Fluent doc]

      NOTE-2:  tb/t* piecewise correlations use strict-less-than (<) at
               interval boundaries rather than the less-or-equal (<=) of
               Pilch & Erdman (1987) Eqs. [8]-[12]. The physical impact is
               negligible because the correlations are continuous at the
               boundary values by construction. Documented here for
               bibliographic completeness.
               [Pilch & Erdman (1987), Eqs. [8]-[12]]

\*------------------------------------------------------------------------------------------------------------*/

#include "MadabhushiDragForce.H"
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
            // CSV header
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

        // No per-tag counter limit: all events are written when debug=true.
        // Disable debug in production runs to avoid large log files.
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
    // Returning Zero here prevents double-counting when both calcCoupled
    // and calcNonCoupled are summed by SprayParcel::calcBreakup to
    // evaluate tMom and the net aerodynamic force on the parcel.
    // [OpenFOAM ParticleForce API convention: coupled = two-way exchange
    //  with the carrier phase; non-coupled = one-way body force only]
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

    // Use the parcel's local (temperature-dependent) surface tension.
    // [FIX Bug7: replaced hardcoded sigma_ with p.sigma()]
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
    // Column-breakup time of the intact liquid column.
    // Madabhushi uses the crossflow gas velocity magnitude.
    // [Madabhushi (2003), Eq. (1); Lambert et al. (2019), Eq. (1)]
    // ------------------------------------------------------------------
    const scalar tColumnBreakup =
        C0_*(Dinj_/(UgRef_ + VSMALL))*Foam::sqrt(rhoL/(rhoRef_ + VSMALL));

    // ------------------------------------------------------------------
    // Child parcels are treated as PE-active only after an explicit latch.
    // We require:
    //   - child marker user_ > 0
    //   - breakup-state marker ms = 2
    //   - positive PE start time
    //   - positive latched PE relative velocity
    // This avoids ambiguous entry into the PE drag branch.
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
    // Secondary children not yet latched into standard PE retain
    // ordinary spherical drag.
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
    // Default drag state
    // ------------------------------------------------------------------
    scalar CdEff     = CdSphere(Re);
    scalar dRefEff   = dCurr;
    scalar dForceEff = dCurr;
    scalar dragScale = 1.0;
    scalar ReEff     = Re;

    // ------------------------------------------------------------------
    // Blob / intact liquid-column regime (core parcel only)
    //
    // CdBlob_ = 1.48: value used in Madabhushi (2003) and Lambert (2019)
    // for the intact liquid-column regime.
    // [Madabhushi (2003), jet-regime description;
    //  Lambert et al. (2019), model summary for pre-column-breakup drag]
    // ------------------------------------------------------------------
    if (userFlag < 0.5 && tSecStart <= 0.0 && tc < tColumnBreakup)
    {
        CdEff     = CdBlob_;
        dRefEff   = Dinj_;
        dForceEff = Dinj_;
        dragScale = 1.0;

        // Reynolds number rescaled consistently with the effective force
        // diameter used in the intact-column drag evaluation.
        ReEff = Re*(dForceEff/(dCurr + VSMALL));

        dragDebugCSV
        (
            debug_, "BLOB",
            tc, userFlag, msState, tColumnBreakup, tSecStart,
            dCurr, dForceEff, dRefEff,
            Re, ReEff, CdEff, dragScale
        );
    }

    // ------------------------------------------------------------------
    // Pilch-Erdman deformation / disc-drag stage
    // Applies to:
    //   - core parcel after column breakup (userFlag=0, tSecStart>0)
    //   - child parcels already latched into standard PE breakup
    //
    // FIX-DR1: weCritPE now uses the viscosity-corrected formula
    //   We_crit = 12 * (1 + 1.077 * Oh^1.6)
    // consistent with Schmehl et al. (1998) Eq.(38),
    // Lambert et al. (2019) onset criterion, and Madabhushi.C PATH 1/4.
    // Previously used the bare value 12.0, which activated deformed-disc
    // drag for viscous parcels that would not physically deform.
    // [Schmehl et al. (1998), Eq. (38);
    //  Lambert et al. (2019), PE onset criterion;
    //  Pilch & Erdman (1987), Eq. (5)]
    // ------------------------------------------------------------------
    else if (tSecStart > 0.0)
    {
        // Weber number evaluated with the latched PE-onset state.
        // [Pilch & Erdman (1987), Eq. (1);
        //  Lambert et al. (2019), Eqs. (2)-(3)]
        const scalar WeInitPE =
            rhoG*sqr(UrelPE0)*Dparent0/(sigmaLocal + VSMALL);

        // FIX-DR1: Ohnesorge number at PE onset for viscous correction.
        // Uses Dparent0 consistent with breakup model PATH 4 (Madabhushi.C).
        // [Pilch & Erdman (1987), Eq. (2); Lambert et al. (2019), Eq. (8)]
        const scalar ohPE0 =
            rhoL > VSMALL
          ? p.mu()/(Foam::sqrt(rhoL*Dparent0*sigmaLocal) + VSMALL)
          : 0.0;

        // FIX-DR1: viscosity-corrected critical Weber number.
        // [Schmehl et al. (1998), Eq. (38);
        //  Lambert et al. (2019), onset criterion;
        //  Pilch & Erdman (1987), Eq. (5)]
        const scalar weCritPE = 12.0*(1.0 + 1.077*pow(ohPE0, 1.6));

        scalar tElapsedLog = 0.0;
        scalar tDefLog     = 0.0;
        scalar tbLog       = 0.0;
        Foam::word phaseTag = "PE_POST";

        if (WeInitPE > weCritPE)
        {
            const scalar tElapsed = max(tc - tSecStart, 0.0);
            tElapsedLog = tElapsed;

            // Characteristic breakup timescale t*.
            // [Pilch & Erdman (1987), Eq. (3);
            //  Lambert et al. (2019), Eq. (3)]
            const scalar tStar =
                (Dparent0/(UrelPE0 + VSMALL))
               *Foam::sqrt(rhoL/(rhoG + VSMALL));

            // Deformation time tDef = 1.6 t*.
            // [Schmehl et al. (1998), Eq. (42);
            //  Lambert et al. (2019), Eq. (2)]
            const scalar tDef = 1.6*tStar;
            tDefLog = tDef;

            // FIX-DR3: Use (We - 12) for the breakup-time correlations,
            // consistent with the original P&E calibration and Madabhushi
            // (2003) Eq. (5). The viscosity-corrected weCritPE is used only
            // for the admission condition above.
            // [Pilch & Erdman (1987), Eqs. (8)-(12);
            //  Madabhushi (2003), Eq. (5)]
            scalar weExcess = max(WeInitPE - 12.0, VSMALL);

            // Dimensionless total breakup time tb/t* from Pilch-Erdman.
            // NOTE-2: strict-less-than boundaries vs the <= of Pilch & Erdman
            // (1987) Eqs. [8]-[12]; negligible physical impact — see header.
            // [Pilch & Erdman (1987), Eqs. (8)-(12);
            //  Schmehl et al. (1998), Eq. (43)]
            scalar tbOverTstar = 5.5;
            if      (WeInitPE < 18.0)   tbOverTstar = 6.0  *pow(weExcess, -0.25);
            else if (WeInitPE < 45.0)   tbOverTstar = 2.45 *pow(weExcess,  0.25);
            else if (WeInitPE < 351.0)  tbOverTstar = 14.1 *pow(weExcess, -0.25);
            else if (WeInitPE < 2670.0) tbOverTstar = 0.766*pow(weExcess,  0.25);

            const scalar tb = tbOverTstar*tStar;
            tbLog = tb;

            if (tElapsed < tDef)
            {
                // Deformation stage: linear CD ramp sphere -> disc.
                // [Lambert et al. (2019), deformed-droplet drag;
                //  Schmehl et al. (1998), text after Eq. (46)]
                const scalar frac = tElapsed/(tDef + VSMALL);
                CdEff = CdSphere(Re)*(1.0 - frac) + CdDisc_*frac;

                // Deformed reference diameter ramp.
                // For We < 100: D_ref(t) = (1 + 0.19*sqrt(We)*t/tDef)*D0
                // For We >= 100: D_ref(t) = (1 + 1.9*t/tDef)*D0
                //   where 1.9 = 0.19*sqrt(100), i.e. Madabhushi (2003)
                //   Eq. (6) evaluated at We=100 as prescribed: "For We>100,
                //   Eq. (6) is evaluated at We=100." Confirmed by the
                //   official Lambert/Fluent documentation figure.
                // [Madabhushi (2003), Eq. (6);
                //  Schmehl et al. (1998), Eq. (45);
                //  Lambert et al. (2019), Fluent documentation figure]
                dRefEff = (WeInitPE < 100.0)
                  ? (1.0 + 0.19*Foam::sqrt(WeInitPE)*frac)*Dparent0
                  : (1.0 + 1.9*frac)*Dparent0;

                phaseTag = "PE_DEFORM";
            }
            else if (tElapsed < tb)
            {
                // Fully deformed stage: constant disc drag.
                // [Lambert et al. (2019), disc-drag assumption;
                //  Schmehl et al. (1998), Sec. 4.3.3]
                CdEff = CdDisc_;
                dRefEff = (WeInitPE < 100.0)
                  ? (1.0 + 0.19*Foam::sqrt(WeInitPE))*Dparent0
                  : 2.9*Dparent0;

                phaseTag = "PE_DISC";
            }
            else
            {
                // Post-breakup: spherical drag again.
                CdEff   = CdSphere(Re);
                dRefEff = dCurr;
                phaseTag = "PE_POST";
            }
        }
        else
        {
            // Below viscosity-corrected PE critical Weber threshold.
            CdEff   = CdSphere(Re);
            dRefEff = dCurr;
            phaseTag = "PE_SUBCRIT";
        }

        // Force diameter stays on current parcel; area correction via dragScale.
        // [Schmehl et al. (1998), Eqs. (45)-(46);
        //  Lambert et al. (2019), Eq. (9)]
        dForceEff = dCurr;
        ReEff     = Re;
        dragScale = sqr(max(dRefEff, 1.0e-12)/dCurr);

        dragDebugCSV
        (
            debug_, phaseTag,
            tc, userFlag, msState, tColumnBreakup, tSecStart,
            dCurr, dForceEff, dRefEff,
            Re, ReEff, CdEff, dragScale,
            WeInitPE, tElapsedLog, tDefLog, tbLog
        );
    }
    else
    {
        // Fallback spherical drag
        CdEff     = CdSphere(Re);
        dRefEff   = dCurr;
        dForceEff = dCurr;
        dragScale = 1.0;
        ReEff     = Re;

        dragDebugCSV
        (
            debug_, "FALLBACK",
            tc, userFlag, msState, tColumnBreakup, tSecStart,
            dCurr, dForceEff, dRefEff,
            Re, ReEff, CdEff, dragScale
        );
    }

    return forceSuSp
    (
        Zero,
        mass*0.75*CdEff*dragScale*ReEff*muc/(rhoL*sqr(dForceEff))
    );
}
