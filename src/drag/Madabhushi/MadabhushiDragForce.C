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
\*------------------------------------------------------------------------------------------------------------*/

#include "MadabhushiDragForce.H"
#include "OFstream.H"
#include "HashTable.H"
#include <map>

namespace
{
    Foam::OFstream& dragLogFile()
    {
        static Foam::autoPtr<Foam::OFstream> logPtr;

        if (!logPtr.valid())
        {
            logPtr.reset(new Foam::OFstream("dragDebug.log"));
        }

        return logPtr();
    }


    inline void dragDebugLimited
    (
        const bool debug,
        const Foam::word& tag,
        const Foam::string& msg,
        const int maxCount = 8
    )
    {
        if (!debug)
        {
            return;
        }

        static std::map<std::string, int> counters;

        std::string key(tag.c_str());
        int count = counters[key];

        if (count < maxCount)
        {
            counters[key] = count + 1;
            dragLogFile() << "[" << tag << "] " << msg << Foam::endl;
        }
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
    return calcCoupled(p, td, dt, mass, Re, muc);
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
        C0_*(Dinj_/(UgRef_ + VSMALL))*Foam::sqrt(rhoL/(rhoG + VSMALL));

    // ------------------------------------------------------------------
    // Child parcels are treated as PE-active only after an explicit latch.
    // We require:
    //   - child marker user_ > 0
    //   - breakup-state marker ms = 2
    //   - positive PE start time
    //   - positive latched PE relative velocity
    // This avoids ambiguous entry into the PE drag branch.
    //
    // This is an implementation-state convention, not a model equation.
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
    //
    // This branch is an implementation choice consistent with the breakup
    // state machine: before explicit PE latch, children are treated as
    // ordinary droplets.
    // ------------------------------------------------------------------
    if (userFlag > 0.5 && !isStdPEChild)
    {
        const scalar CdEffChild = CdSphere(Re);

        dragDebugLimited
        (
            debug_,
            "DRAG_CHILD",
            "DRAG_CHILD "
          + Foam::string("tc=") + Foam::name(tc)
          + " user=" + Foam::name(userFlag)
          + " ms="   + Foam::name(msState)
          + " y="    + Foam::name(tSecStart)
          + " d="    + Foam::name(dCurr)
          + " Re="   + Foam::name(Re)
          + " Cd="   + Foam::name(CdEffChild),
            100
        );

        forceSuSp value
        (
            Zero,
            mass*0.75*CdEffChild*Re*muc/(rhoL*sqr(dCurr))
        );

        return value;
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
    // The liquid column is represented by a spherical computational parcel,
    // but the drag coefficient is kept fixed at CdBlob_ = 1.48 to mimic the
    // intact jet-column regime before column breakup.
    //
    // CdBlob_ = 1.48 is the value used in Madabhushi for the intact
    // liquid-column regime and is also retained in Lambert's reformulation.
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
        //
        // This is a consistent implementation rescaling because the drag
        // formula uses Dinj_ as effective aerodynamic diameter in the blob
        // regime, while the parcel diameter may already differ numerically.
        ReEff = Re*(dForceEff/(dCurr + VSMALL));

        dragDebugLimited
        (
            debug_,
            "DRAG_BLOB",
            "DRAG_BLOB "
          + Foam::string("tc=") + Foam::name(tc)
          + " tCb="       + Foam::name(tColumnBreakup)
          + " dCurr="     + Foam::name(dCurr)
          + " dForceEff=" + Foam::name(dForceEff)
          + " Re="        + Foam::name(Re)
          + " ReEff="     + Foam::name(ReEff)
          + " Cd="        + Foam::name(CdEff),
            100
        );
    }

    // ------------------------------------------------------------------
    // Pilch-Erdman deformation / disc-drag stage
    // Applies to:
    //   - core parcel after column breakup
    //   - child parcels already latched into standard PE breakup
    //
    // Staged treatment:
    //   1) deformation stage: sphere-to-disc transition
    //   2) fully deformed stage: disc drag
    //   3) after breakup completion: spherical drag again
    // ------------------------------------------------------------------
    else if (tSecStart > 0.0)
    {
        // Weber number evaluated with the latched PE-onset state
        // (UrelPE0, Dparent0), not with the current parcel state.
        // [Pilch & Erdman (1987), Eq. (1);
        //  Schmehl et al. (1998), Eq. (37);
        //  Lambert et al. (2019), Eqs. (2)-(3)]
        const scalar WeInitPE =
            rhoG*sqr(UrelPE0)*Dparent0/(sigma_ + VSMALL);

        // Baseline critical Weber threshold for PE breakup timing.
        // Viscous effects are introduced separately through Oh / Wecorr.
        // [Pilch & Erdman (1987), Eq. (5);
        //  Schmehl et al. (1998), Eq. (38);
        //  Lambert et al. (2019), onset criterion]
        const scalar weCritPE = 12.0;

        if (WeInitPE > weCritPE)
        {
            const scalar tElapsed = max(tc - tSecStart, 0.0);

            // Characteristic breakup timescale:
            // t* = (D0 / Urel,0) * sqrt(rho_l / rho_g)
            // [Pilch & Erdman (1987), Eq. (3);
            //  Schmehl et al. (1998), Eq. (41);
            //  Lambert et al. (2019), Eq. (3)]
            const scalar tStar =
                (Dparent0/(UrelPE0 + VSMALL))
               *Foam::sqrt(rhoL/(rhoG + VSMALL));

            // Deformation time:
            // tDef = 1.6 t*
            // [Schmehl et al. (1998), Eq. (42);
            //  Lambert et al. (2019), Eq. (2)]
            const scalar tDef = 1.6*tStar;

            // Numerical safeguard near threshold.
            scalar weExcess = max(WeInitPE - weCritPE, VSMALL);

            // Dimensionless total breakup time:
            // tb/t* = f(We)
            // [Pilch & Erdman (1987), Eqs. (8)-(12);
            //  Schmehl et al. (1998), Eq. (43);
            //  Lambert et al. (2019), piecewise PE timing law]
            scalar tbOverTstar = 5.5;

            if (WeInitPE < 18.0)
            {
                tbOverTstar = 6.0*pow(weExcess, -0.25);
            }
            else if (WeInitPE < 45.0)
            {
                tbOverTstar = 2.45*pow(weExcess, 0.25);
            }
            else if (WeInitPE < 351.0)
            {
                tbOverTstar = 14.1*pow(weExcess, -0.25);
            }
            else if (WeInitPE < 2670.0)
            {
                tbOverTstar = 0.766*pow(weExcess, 0.25);
            }

            const scalar tb = tbOverTstar*tStar;

            if (tElapsed < tDef)
            {
                const scalar frac = tElapsed/(tDef + VSMALL);

                // Linear drag interpolation from spherical drag to disc drag
                // during the deformation stage.
                //
                // The gradual transition itself is a closure used to model
                // the aerodynamic response of the flattening droplet.
                // [Schmehl et al. (1998), text after Eq. (46);
                //  Efficient Numerical Calculation..., Sec. 4.3.3;
                //  Lambert et al. (2019), deformed-droplet drag treatment]
                CdEff = CdSphere(Re)*(1.0 - frac) + CdDisc_*frac;

                // Deformed reference diameter during deformation.
                // For We < 100:
                //   Dmax/D0 = 1 + 0.19 sqrt(We)
                // and the code applies a linear ramp from D0 to Dmax.
                // For larger Weber numbers, deformation is capped at 2.9 D0.
                // [Schmehl et al. (1998), Eq. (45);
                //  Efficient Numerical Calculation..., Eq. (45);
                //  Lambert et al. (2019), Eq. (9)]
                if (WeInitPE < 100.0)
                {
                    dRefEff =
                        (1.0 + 0.19*Foam::sqrt(WeInitPE)*frac)*Dparent0;
                }
                else
                {
                    dRefEff =
                        (1.0 + 1.9*frac)*Dparent0;
                }
            }
            else if (tElapsed < tb)
            {
                // Fully deformed stage: droplet treated as a disc
                // with constant disc drag coefficient.
                //
                // The constant disc state after tDef is explicitly adopted
                // in the Schmehl-based drag closure to bypass the complexity
                // of the ongoing disintegration process.
                // [Schmehl et al. (1998), text in Sec. 4.3.3;
                //  Efficient Numerical Calculation..., Sec. 4.3.3;
                //  Lambert et al. (2019), deformed-disc assumption]
                CdEff = CdDisc_;

                // CdDisc_ = 1.2 corresponds to the deformed-disc drag value
                // used in this PE framework.
                // [Lambert et al. (2019), disc-drag assumption]
                dRefEff =
                    (WeInitPE < 100.0)
                  ? (1.0 + 0.19*Foam::sqrt(WeInitPE))*Dparent0
                  : 2.9*Dparent0;
            }
            else
            {
                // After the modeled breakup time is exceeded, the parcel is
                // again treated with spherical drag, consistent with the fact
                // that the tracked entity now represents a post-breakup
                // fragment-like droplet.
                CdEff   = CdSphere(Re);
                dRefEff = dCurr;
            }
        }
        else
        {
            // Below PE critical Weber threshold: no PE deformation drag.
            CdEff   = CdSphere(Re);
            dRefEff = dCurr;
        }

        // In PE-active drag we keep the force denominator based on the
        // current parcel diameter.
        //
        // This is an implementation choice consistent with the original
        // OpenFOAM drag-force structure: the projected-area correction is
        // accounted for separately through dragScale.
        dForceEff = dCurr;
        ReEff     = Re;

        // Aerodynamic area correction based on the ratio between deformed
        // reference diameter and current parcel diameter:
        // dragScale ~ (Dref / Dcurr)^2
        //
        // This directly reflects the use of the deformed projected area
        // during the flattening stage.
        // [Schmehl et al. (1998), Eq. (45) through Dmax;
        //  Efficient Numerical Calculation..., Eqs. (45)-(46);
        //  Lambert et al. (2019), Eq. (9)]
        dragScale = sqr(max(dRefEff, 1.0e-12)/dCurr);

        dragDebugLimited
        (
            debug_,
            "DRAG_PE",
            "DRAG_PE "
          + Foam::string("tc=")    + Foam::name(tc)
          + " user="               + Foam::name(userFlag)
          + " tCb="                + Foam::name(tColumnBreakup)
          + " tSecStart="          + Foam::name(tSecStart)
          + " Dparent0="           + Foam::name(Dparent0)
          + " UrelPE0="            + Foam::name(UrelPE0)
          + " WeInitPE="           + Foam::name(WeInitPE)
          + " dCurr="              + Foam::name(dCurr)
          + " dRefEff="            + Foam::name(dRefEff)
          + " dragScale="          + Foam::name(dragScale)
          + " Cd="                 + Foam::name(CdEff),
            100
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
    }

    forceSuSp value
    (
        Zero,
        mass*0.75*CdEff*dragScale*ReEff*muc/(rhoL*sqr(dForceEff))
    );

    return value;
}
