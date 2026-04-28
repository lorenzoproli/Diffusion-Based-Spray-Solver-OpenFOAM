/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Madabhushi breakup model — enhanced crossflow version

    Physical model backbone:
      - Column-breakup time from Madabhushi (2003), Eq. (1)
      - Wave-shedding stage based on Reitz KH wave model
      - Post-column breakup based on Pilch-Erdman correlations
      - Child-size correction after first ligament breakup from
        Lambert et al. (2019)

    Implementation-specific extensions:
      (A) Wave-shed children (user_=1) can latch into standard PE breakup
      (B) Post-catastrophic children (user_=2) can latch into standard PE breakup
      (C) Catastrophic breakup can spawn nChildren_ parcels
      (D) Explicit internal state flags separate breakup logic from drag logic
      (E) Child/mother velocity assignment order is controlled explicitly

    Persistent state convention:
      y       : PE start time tSecStart if > 0, otherwise no PE family active
      yDot    : before PE onset = stored initial mass; after PE onset = latched Dparent0
      KHindex : latched UrelPE0 at PE onset
      ms      : breakup-state marker / stripped-mass accumulator depending on regime

    Modifications (Thesis — Madabhushi crossflow model):
      FIX-MB1: Debug output restructured as CSV (breakupDebug.log).
               Separate log files per phase for easy filtering:
                 breakupDebug_wave.log     -> PATH 2 wave shedding events
                 breakupDebug_stage2.log   -> PATH 3 column breakup latch
                 breakupDebug_pe.log       -> PATH 4 PE pre-breakup progress
                 breakupDebug_breakup.log  -> PATH 4 actual breakup events
                 breakupDebug_stdpe.log    -> PATH 1 sentinel PE latch checks
               All per-tag counter limits removed: when debug=true every
               event is written. Disable debug in production runs to avoid
               large log files.

      FIX-MB2: (BUG-1) Condition isFirstCorePEBreakup narrowed from
               (ms < 0.0) to (ms > -1.5 && ms < -0.5).
               The previous form was semantically too broad: any future
               negative ms value would have incorrectly activated FLig.
               The only bibliographically correct trigger is ms == -1.0,
               which is the value set exclusively by PATH 3 at column-
               breakup latch.
               [Lambert et al. (2019), Eq. (11): FLig applies only to
                child droplets created immediately after column breakup]

      FIX-MB3: (BUG-2) yDot initial-mass latch now excludes sentinel
               parcels (wave-shed ms=-10, post-cat ms=-20).
               Previously the latch fired for any parcel with y<=0 and
               yDot<=0, including post-cat sentinel fragments returning
               through PATH 1. The sentinel yDot was overwritten with
               the fragment's current mass — a superfluous and potentially
               misleading state assignment, even though PATH 1 always
               returns before yDot is consumed in PATH 2.

      FIX-MB5: (BUG-4) In PATH 3, d is now updated from parcelMass and
               nParticle after the mass reabsorb, before latching yDot = d.
               Previously parcelMass was incremented by ms but d was left at
               its wave-shed value, creating an inconsistent (parcelMass, d,
               nParticle) triple. The latched Dparent0 = yDot was therefore
               underestimated, causing weInitPE and tStar to be computed with
               a diameter smaller than the physical intact-column cross-section
               at the moment of column fracture. Now d is back-calculated as
               cbrt(parcelMass / (nParticle * rho * pi/6)) so the triple is
               consistent before the PE latch.
               [Madabhushi (2003), column-breakup model description]

      FIX-MB6: Removed the previous cap min(uNormalMagRaw, UrelPE0)
               from the rim-expansion velocity. The implemented value now
               follows directly the Madabhushi/Lambert expression
               u_n = 5*Dparent0/(tb - tDef), with only VSMALL in the
               denominator for numerical protection.
               [Madabhushi (2003), text after Eq. (9);
                Lambert et al. (2019), Eq. (5)]

      FIX-MB11: Fragment #1 STD_PE relatch now uses the actual slip velocity
                of the newly generated fragment, reconstructed as
                UrelFrag0 = UFrag[0] - UgasLocal, instead of the pre-breakup
                parent Urmag. The relatch dStable formula has also been made
                identical to PATH 1: WeCrit is used only for admission, while
                dStable uses the inviscid critical Weber number 12.
                This prevents already-stable fragments from being re-latched
                due to an obsolete pre-breakup relative velocity.

      FIX-MB7: stage2 log now records the updated d (post-reabsorb) to
               reflect the actual Dparent0 used in PATH 4.

      FIX-MB9: (BUG-6) STD_PE child relatch after breakup now applies
               the same We > Wec AND d > dStable admission criteria
               used by PATH 1 for sentinel children.
               Previously, Fragment #1 of a STD_PE breakup event was
               unconditionally relatched into a new PE cycle (ms=2.0)
               regardless of whether the post-breakup diameter was
               already below the maximum stable diameter or the local
               Weber number had dropped below Wec. This caused an
               unlimited chain of PE breakup cycles on already-stable
               fragments, systematically underestimating D32 in the
               near-field region (z = 12-25 mm).
               The fix mirrors exactly the PATH 1 admission logic:
               if the post-breakup fragment satisfies We > Wec AND
               d > dStable it is relatched (ms=2.0); otherwise it
               becomes a post-cat sentinel (ms=-20.0) and is subject
               to the standard PATH 1 evaluation on subsequent steps.
               Instantaneous Urmag and d are used (consistent with
               PATH 1), not the frozen UrelPE0/Dparent0 of the current
               PE cycle.
               [Pilch & Erdman (1987), Eqs. (20) and (33):
                breakup ceases when all fragments are smaller than the
                maximum stable diameter; OpenFOAM PilchErdman.C: same
                dStable check applied before every breakup update]

      NOTE-1:  weWaveCrit = 6.0 corresponds to the bag-breakup onset
               criterion We_2 > 6.0 of Reitz (1987), Eq. (15a), which
               defines the minimum gas-side Weber number for wave-shedding
               activity in the KH wave model.
               [Reitz (1987), Eq. (15a)]

      NOTE-2:  B1PE = 0.375 and B2PE = 0.2274 are the fragment-cloud
               velocity constants from the standard OpenFOAM PilchErdman
               breakup model, themselves derived from Pilch & Erdman
               (1987) Eq. [20] for incompressible gas-liquid systems.
               [Pilch & Erdman (1987), Eq. [20];
                OpenFOAM PilchErdman.C default values]

      NOTE-3:  tb/t* piecewise correlations use strict-less-than (<) at
               interval boundaries rather than the less-or-equal (<=) of
               Pilch & Erdman (1987) Eqs. [8]-[12]. The physical impact is
               negligible because the correlations are continuous at the
               boundary values by construction. Documented here for
               bibliographic completeness.
               [Pilch & Erdman (1987), Eqs. [8]-[12]]

\*---------------------------------------------------------------------------*/

#include "Madabhushi.H"
#include "OFstream.H"
#include "mathematicalConstants.H"
#include "Pstream.H"
#include "autoPtr.H"

namespace
{
    // -------------------------------------------------------------------------
    // FIX-MB1: one log file per phase, CSV format
    // -------------------------------------------------------------------------

    Foam::OFstream& waveLogFile()
    {
        static Foam::autoPtr<Foam::OFstream> p;
        if (!p.valid())
        {
            p.reset(new Foam::OFstream("breakupDebug_wave.log"));
            p() << "tc,d,dOld,dChildTarget,massStripped,ms,weGas,lambdaKH,tauKH,"
                << "childVmag,parentUmag" << Foam::endl;
        }
        return p();
    }

    Foam::OFstream& stage2LogFile()
    {
        static Foam::autoPtr<Foam::OFstream> p;
        if (!p.valid())
        {
            p.reset(new Foam::OFstream("breakupDebug_stage2.log"));
            p() << "tc,tColumnBreakup,rhoc,UgRef,Dparent0,UrelPE0,Urmag"
                << Foam::endl;
        }
        return p();
    }

    Foam::OFstream& peLogFile()
    {
        static Foam::autoPtr<Foam::OFstream> p;
        if (!p.valid())
        {
            p.reset(new Foam::OFstream("breakupDebug_pe.log"));
            p() << "tc,tSecStart,tElapsed,tDef,tb,Dparent0,UrelPE0,"
                << "weInitPE,weCritPE,dSMD,FLigEff,isStdPEChild"
                << Foam::endl;
        }
        return p();
    }

    Foam::OFstream& breakupLogFile()
    {
        static Foam::autoPtr<Foam::OFstream> p;
        if (!p.valid())
        {
            p.reset(new Foam::OFstream("breakupDebug_breakup.log"));
            p() << "mode,tc,Dparent0,UrelPE0,tElapsed,tDef,tb,"
                << "uNormalMag,UmotherMag,dSMD,FLigEff,dFrag0,UFrag0mag"
                << Foam::endl;
        }
        return p();
    }

    Foam::OFstream& stdPeLogFile()
    {
        static Foam::autoPtr<Foam::OFstream> p;
        if (!p.valid())
        {
            p.reset(new Foam::OFstream("breakupDebug_stdpe.log"));
            p() << "event,type,tc,d,Urmag,weCurr,weCritCurr,dStable,allow"
                << Foam::endl;
        }
        return p();
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Madabhushi<CloudType>::Madabhushi
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(dict, owner, typeName),
    C0_(this->coeffDict().getOrDefault("C0", 3.44)),
    FLig_(this->coeffDict().getOrDefault("FLig", 0.4)),
    b0_(this->coeffDict().getOrDefault("b0", 0.61)),
    b1_(this->coeffDict().getOrDefault("b1", 10.0)),
    Dinj_(this->coeffDict().getOrDefault("Dinj", 0.0016)),
    nChildren_(this->coeffDict().getOrDefault("nChildren", 5)),
    minChildDiameter_(this->coeffDict().getOrDefault("minChildDiameter", 1.0e-7)),
    debug_(this->coeffDict().getOrDefault("debug", false)),
    childMsInit_(-GREAT),
    childUserInit_(0.0),
    pendingChildren_(),
    pendingChildUserFlag_(0.0),
    parentUserUpdate_(-GREAT)
{}


template<class CloudType>
Foam::Madabhushi<CloudType>::Madabhushi
(
    const Madabhushi<CloudType>& model
)
:
    BreakupModel<CloudType>(model),
    C0_(model.C0_),
    FLig_(model.FLig_),
    b0_(model.b0_),
    b1_(model.b1_),
    Dinj_(model.Dinj_),
    nChildren_(model.nChildren_),
    minChildDiameter_(model.minChildDiameter_),
    debug_(model.debug_),
    childMsInit_(-GREAT),
    childUserInit_(0.0),
    pendingChildren_(),
    pendingChildUserFlag_(0.0),
    parentUserUpdate_(-GREAT)
{}


template<class CloudType>
Foam::Madabhushi<CloudType>::~Madabhushi()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::Madabhushi<CloudType>::update
(
    const scalar dt,
    const vector& g,
    scalar& d,
    scalar& tc,
    scalar& ms,
    scalar& nParticle,
    scalar& KHindex,
    scalar& y,
    scalar& yDot,
    const scalar d0,
    const scalar rho,
    const scalar mu,
    const scalar sigma,
    const vector& U,
    const scalar rhoc,
    const scalar muc,
    const vector& Urel,
    const scalar Urmag,
    const scalar tMom,
    scalar& dChild,
    scalar& massChild
)
{
    bool addParcel = false;

    // -----------------------------------------------------------------
    // Reset mutable signaling state for possible child creation
    // -----------------------------------------------------------------
    childMsInit_ = -GREAT;
    childUserInit_ = 0.0;
    pendingChildren_.clear();
    pendingChildUserFlag_ = 0.0;
    parentUserUpdate_ = -GREAT;

    const scalar pi = constant::mathematical::pi;
    const scalar rhoPiOver6 = rho*pi/6.0;

    scalar parcelMass = nParticle*pow3(d)*rhoPiOver6;

    // -----------------------------------------------------------------
    // Before PE onset, yDot is used to store the initial mass of the
    // intact core parcel only once. After PE onset, yDot is repurposed
    // to store Dparent0.
    // This is an implementation choice, not a physical model equation.
    //
    // FIX-MB3: Sentinel parcels (wave-shed ms=-10, post-cat ms=-20) must
    // be excluded from this latch. They satisfy y<=0 and yDot<=0 when
    // they first enter update(), but they are not core parcels — writing
    // parcelMass into yDot would be a superfluous and misleading state
    // assignment. PATH 1 always returns before yDot is consumed in
    // PATH 2, so the impact is currently inert, but excluding sentinels
    // here makes the invariant explicit and safe against future changes.
    // -----------------------------------------------------------------
    const bool isWaveSentinelEarly    = (ms > -11.0 && ms < -9.0);
    const bool isPostCatSentinelEarly = (ms > -21.0 && ms < -19.0);

    if (y <= 0.0 && yDot <= 0.0
     && !isWaveSentinelEarly && !isPostCatSentinelEarly)
    {
        yDot = parcelMass;
    }

    const scalar injectorRadius = 0.5*Dinj_;

    // -----------------------------------------------------------------
    // Persistent-variable convention
    //
    // y       > 0   : tSecStart for PE-active parcels
    // y       <= 0  : no PE family active
    //
    // yDot          : Dparent0 latched at PE onset
    // KHindex       : UrelPE0 latched at PE onset
    //
    // ms encoding:
    //   Wave stage (y<=0):  accumulated stripped mass
    //   ms = -1.0:          first Madabhushi PE breakup (FLig_ applies)
    //   ms =  2.0:          standard PE child (spherical drag before latch)
    //   ms = -10.0:         sentinel: wave-shed child, needs STD_PE check
    //   ms = -20.0:         sentinel: post-cat child, needs STD_PE check
    //
    // user_ is set outside update() by calcBreakup:
    //   0.0 = core parcel
    //   1.0 = wave-shed child
    //   2.0 = post-cat child
    //
    // This state encoding is an implementation convention.
    // -----------------------------------------------------------------

    // -----------------------------------------------------------------
    // Column-breakup time for the intact liquid column.
    // Madabhushi uses the crossflow gas velocity magnitude, not the
    // instantaneous parcel-relative velocity.
    // [Madabhushi (2003), Eq. (1); Lambert et al. (2019), Eq. (1)]
    // -----------------------------------------------------------------
    // Local gas velocity magnitude: Ug = U_droplet - Urel = U_gas (lab frame)
    const scalar UgMag = mag(U - Urel);
    const scalar tColumnBreakup =
        C0_*(Dinj_/(UgMag + VSMALL))*Foam::sqrt(rho/(rhoc + VSMALL));

    tc += dt;

    // =================================================================
    // PATH 0: terminal inert droplet
    // =================================================================
    if (ms > 2.5)
    {
        nParticle = parcelMass/(pow3(d)*rhoPiOver6 + VSMALL);
        return false;
    }

    // =================================================================
    // PATH 1: sentinel children can latch into standard PE
    // =================================================================
    // ms = -10.0 : wave-shed child   (set by PATH 2 via childMsInit_)
    // ms = -20.0 : post-catastrophic child (set by PATH 4 breakup event)
    //
    // The same boolean flags computed above for the yDot guard are
    // reused here — they identify exactly the same sentinel states.
    // =================================================================

    const bool isWaveSentinel    = isWaveSentinelEarly;
    const bool isPostCatSentinel = isPostCatSentinelEarly;

    if (isWaveSentinel || isPostCatSentinel)
    {
        // Current Weber number based on the local, instantaneous child state.
        // [Pilch & Erdman (1987), Eq. (1); Lambert et al. (2019), Eqs. (8)-(10)]
        const scalar weCurr =
            rhoc*sqr(Urmag)*d/(sigma + VSMALL);

        // Ohnesorge number based on the local child state.
        // [Pilch & Erdman (1987), Eq. (2); Lambert et al. (2019), Eq. (8)]
        const scalar ohCurr =
            mu/(Foam::sqrt(rho*d*sigma + VSMALL));

        // Critical Weber number for breakup onset (viscosity corrected).
        // [Pilch & Erdman (1987), Eq. (5); Schmehl et al. (1998), Eq. (38);
        //  Lambert et al. (2019), onset criterion]
        const scalar weCritCurr =
            12.0*(1.0 + 1.077*pow(ohCurr, 1.6));

        // -----------------------------------------------------------------
        // Secondary-PE admission follows Pilch-Erdman logic, not only the
        // instantaneous Weber threshold. A child parcel can enter STD_PE
        // only if it is both supercritical and larger than the maximum stable
        // diameter predicted from the current fragment state.
        // [Pilch & Erdman (1987), Eqs. (20) and (33); OpenFOAM PilchErdman]
        // -----------------------------------------------------------------
        scalar dStable = GREAT;
        bool allowStdPE = false;

        if (weCurr > weCritCurr)
        {
            // tb/t* piecewise correlations.
            // NOTE-3: strict-less-than boundaries vs Pilch & Erdman (1987)
            // <=; negligible physical impact, documented in header NOTE-3.
            // [Pilch & Erdman (1987), Eqs. [8]-[12]]
            scalar taubBar = 5.5;

            if (weCurr < 2670.0)
            {
                if (weCurr > 351.0)
                {
                    taubBar = 0.766*pow(max(weCurr - 12.0, VSMALL), 0.25);
                }
                else if (weCurr > 45.0)
                {
                    taubBar = 14.1*pow(max(weCurr - 12.0, VSMALL), -0.25);
                }
                else if (weCurr > 18.0)
                {
                    taubBar = 2.45*pow(max(weCurr - 12.0, VSMALL), 0.25);
                }
                else if (weCurr > 12.0)
                {
                    taubBar = 6.0*pow(max(weCurr - 12.0, VSMALL), -0.25);
                }
            }

            const scalar rho12 = Foam::sqrt(max(rhoc/(rho + VSMALL), VSMALL));

            // Fragment-cloud velocity at breakup completion.
            // B1PE and B2PE are the standard OpenFOAM PilchErdman constants,
            // derived from Pilch & Erdman (1987) Eq. [20] for incompressible
            // gas-liquid systems (Cd=0.5, B=0.0758 modified correlations).
            // They are used here to keep the child-admission criterion
            // consistent with the reference OpenFOAM PilchErdman model.
            // [Pilch & Erdman (1987), Eq. [20];
            //  OpenFOAM PilchErdman.C: B1_=0.375, B2_=0.2274]
            const scalar B1PE = 0.375;
            const scalar B2PE = 0.2274;
            const scalar Vd =
                Urmag*rho12*(B1PE*taubBar + B2PE*sqr(taubBar));

            // Maximum stable diameter for the current child state.
            // P&E Eq.(33) uses We_c=12 (inviscid); viscous correction applies
            // only to the admission check above, not to d_stable.
            // [Pilch & Erdman (1987), Eq. (33)]
            // FIX-MB10: clamp to [0,1] before squaring. When Vd >= Urmag
            // (fragment fully decelerates within tb, which occurs as We->12+
            // because taubBar diverges), children emerge at near-zero relative
            // velocity -> all children stable -> dStable must be GREAT.
            // Without the clamp, sqr of a negative gives large Vd1 -> tiny
            // dStable -> spurious near-threshold breakup.
            const scalar Vd1 = max(
                sqr(max(1.0 - Vd/(Urmag + VSMALL), 0.0)), SMALL);
            dStable = 12.0*sigma/(Vd1*rhoc*sqr(Urmag) + VSMALL);

            allowStdPE = (d > dStable);
        }

        // FIX-MB1: structured CSV log for PATH 1
        if (debug_)
        {
            stdPeLogFile()
                << "CHECK,"
                << (isWaveSentinel ? "Wave" : "PostCat") << ","
                << tc << "," << d << "," << Urmag << ","
                << weCurr << "," << weCritCurr << "," << dStable << ","
                << (allowStdPE ? 1 : 0)
                << Foam::endl;
        }

        if (!allowStdPE)
        {
            // Keep the sentinel state and spherical drag. The parcel remains
            // eligible for re-evaluation at later time steps if the local
            // conditions become favourable for STD_PE onset.
            nParticle = parcelMass/(pow3(d)*rhoPiOver6 + VSMALL);
            return false;
        }

        // Latch standard PE onset using the current local state.
        y = tc;            // tSecStart
        yDot = d;          // Dparent0
        KHindex = Urmag;   // UrelPE0
        ms = 2.0;          // standard PE child

        if (debug_)
        {
            stdPeLogFile()
                << "LATCH,"
                << (isWaveSentinel ? "Wave" : "PostCat") << ","
                << tc << "," << d << "," << Urmag << ","
                << weCurr << "," << weCritCurr << "," << dStable << ",1"
                << Foam::endl;
        }

        nParticle = parcelMass/(pow3(d)*rhoPiOver6 + VSMALL);
        return false;
    }

    // =================================================================
    // PATH 2: stage 1 — intact liquid column / KH wave shedding
    // =================================================================
    if (y <= 0.0 && tc < tColumnBreakup)
    {
        // Gas-side Weber number in the wave model.
        // [Madabhushi (2003), nomenclature; Reitz (1987)]
        const scalar weGasWave =
            rhoc*sqr(Urmag)*injectorRadius/(sigma + VSMALL);

        // Liquid-side Weber number in the wave model.
        // [Madabhushi (2003), nomenclature; Reitz (1987)]
        const scalar weLiqWave =
            rho*sqr(Urmag)*injectorRadius/(sigma + VSMALL);

        // Minimum gas-side Weber number for wave-shedding activity.
        // [Reitz (1987), Eq. (15a): We_2 > 6.0 for bag-breakup onset;
        //  this threshold activates the KH wave-shedding mechanism]
        const scalar weWaveCrit = 6.0;
        const scalar strippedMassFractionThreshold = 0.05;

        const scalar reLiquid =
            rho*Urmag*injectorRadius/(mu + VSMALL);

        // Ohnesorge number Z = sqrt(We_l) / Re_l.
        // [Madabhushi (2003), nomenclature; Reitz (1987), after Eq. (5)]
        const scalar ohWave =
            Foam::sqrt(max(weLiqWave, VSMALL))/(reLiquid + VSMALL);

        // Taylor parameter T = Z * sqrt(We_g).
        // [Madabhushi (2003), nomenclature; Reitz (1987)]
        const scalar taWave =
            ohWave*Foam::sqrt(max(weGasWave, VSMALL));

        // Most unstable KH wave growth rate.
        // [Reitz (1987), Eq. (5); Madabhushi (2003), Eq. (4)]
        const scalar omegaKH =
            (0.34 + 0.38*pow(weGasWave, 1.5))
           /((1.0 + ohWave)*(1.0 + 1.4*pow(taWave, 0.6)))
           *Foam::sqrt(sigma/(rho*pow3(injectorRadius) + VSMALL));

        // Most unstable KH wavelength.
        // [Reitz (1987), Eq. (4); Madabhushi (2003), Eq. (2)]
        const scalar lambdaKH =
            9.02*injectorRadius*(1.0 + 0.45*Foam::sqrt(ohWave))
           *(1.0 + 0.4*pow(taWave, 0.7))
           /pow(1.0 + 0.87*pow(weGasWave, 1.67), 0.6);

        // Characteristic KH-wave stripping timescale.
        // [Reitz (1987), Eq. (12); Madabhushi (2003), Eq. (3)]
        const scalar tauKH =
            3.726*b1_*injectorRadius/(omegaKH*lambdaKH + VSMALL);

        // Characteristic child diameter: D_child = 2 * B0 * Lambda.
        // [Reitz (1987), Eq. (10a); Madabhushi (2003), Eq. (2)]
        const scalar dChildWaveTarget =
            max(2.0*b0_*lambdaKH, 2.0e-6);

        if
        (
            dChildWaveTarget <= Dinj_
         && dChildWaveTarget < d
         && weGasWave > weWaveCrit
        )
        {
            const scalar relaxationFraction = dt/(tauKH + VSMALL);
            const scalar dOld = d;

            // Parent-diameter relaxation toward KH child size.
            // [Reitz (1987), Eq. (11); Madabhushi (2003), Eq. (3)]
            d = (relaxationFraction*dChildWaveTarget + dOld)
              /(1.0 + relaxationFraction);

            const scalar massStrippedThisStep =
                nParticle*rhoPiOver6*(pow3(dOld) - pow3(d));

            ms += massStrippedThisStep;

            // [FIX Bug3] Sync parcelMass with the stripped mass.
            parcelMass = nParticle*pow3(d)*rhoPiOver6;

            if (parcelMass < SMALL)
            {
                parcelMass = SMALL;
            }

            const scalar initialWaveMass = max(yDot, parcelMass + ms);
            const scalar childMassTarget =
                strippedMassFractionThreshold*initialWaveMass;

            if (ms > childMassTarget && childMassTarget > SMALL)
            {
                vector streamwiseDir = U/(mag(U) + VSMALL);

                vector refAxis =
                    (mag(streamwiseDir.x()) < 0.7)
                  ? vector(1, 0, 0)
                  : vector(0, 1, 0);

                vector normal1 = (streamwiseDir ^ refAxis);
                normal1 /= (mag(normal1) + VSMALL);

                vector normal2 = (streamwiseDir ^ normal1);
                normal2 /= (mag(normal2) + VSMALL);

                const scalar phi =
                    2.0*pi*this->owner().rndGen().template sample01<scalar>();

                // Empirical angular-spread coefficient.
                // [Reitz (1987), Eq. (7); Madabhushi (2003), A1 = 0.188]
                const scalar A1 = 0.188;

                const scalar tanThetaHalf =
                    A1*lambdaKH*omegaKH/(mag(U) + VSMALL);

                const scalar v1Mag = mag(U)*tanThetaHalf*sin(phi);
                const scalar v2Mag = mag(U)*tanThetaHalf*cos(phi);

                const vector childVelocityPerturbation =
                    v1Mag*normal1 + v2Mag*normal2;

                const vector parentVelocity = U;
                vector& Umutable = const_cast<vector&>(U);

                // [FIX Bug1] Set U to child velocity for clone.
                Umutable = parentVelocity + childVelocityPerturbation;

                addParcel = true;
                dChild = dChildWaveTarget;
                massChild = childMassTarget;

                childMsInit_ = -10.0;
                childUserInit_ = 1.0;

                // [FIX Bug3] Debit child mass from stripped budget.
                ms -= massChild;

                // FIX-MB1: CSV log for wave shedding event
                if (debug_)
                {
                    waveLogFile()
                        << tc << "," << d << "," << dOld << ","
                        << dChildWaveTarget << "," << massStrippedThisStep
                        << "," << ms << "," << weGasWave << ","
                        << lambdaKH << "," << tauKH << ","
                        << mag(parentVelocity + childVelocityPerturbation) << ","
                        << mag(parentVelocity)
                        << Foam::endl;
                }

                // [FIX Bug1] Sentinel pending child to restore parent velocity.
                pendingChildren_.append
                (
                    pendingChild(-1.0, 0.0, parentVelocity)
                );
            }
        }

        nParticle = parcelMass/(pow3(d)*rhoPiOver6 + VSMALL);
        return addParcel;
    }

    // =================================================================
    // PATH 3: latch PE onset exactly once for the core parcel
    // =================================================================
    const bool isChildParcel =
        (ms > -21.0 && ms < -19.0)   // post-cat sentinel
     || (ms > -11.0 && ms < -9.0)    // wave-shed sentinel
     || (ms > 1.5   && ms < 2.5);    // standard-PE child

    if (y <= 0.0 && tc >= tColumnBreakup && !isChildParcel)
    {
        // [FIX Bug4] Reabsorb residual stripped mass.
        // The wave-shedding stage accumulates stripped mass in ms. Any
        // residual ms > 0 at column-breakup time represents mass that belongs
        // to the intact column but has not yet been emitted as a child parcel.
        // This mass is reabsorbed into the core blob before the PE latch.
        if (ms > 0.0)
        {
            parcelMass += ms;
        }

        // FIX-MB5: Update d to reflect the reabsorbed mass before latching
        // Dparent0 = yDot = d.
        //
        // Previously d was latched at its current (wave-shed) value, while
        // parcelMass had been incremented by ms. This created an inconsistency:
        // nParticle = parcelMass / ((pi/6)*rho*d^3) was inflated, so Dparent0
        // was underestimated relative to the actual blob volume. The PE stage
        // then computed weInitPE and tStar with a diameter that was too small,
        // leading to an underestimate of tb and an overestimate of the breakup
        // rate.
        //
        // The physically consistent approach is to back-calculate d from the
        // total parcelMass and the current nParticle before updating yDot:
        //   d = cbrt(parcelMass / (nParticle * rho * pi/6))
        // This ensures that (parcelMass, d, nParticle) form a consistent triple
        // at the column-breakup latch point.
        //
        // [Madabhushi (2003): "The jet is represented by spherical droplets
        //  with diameter equal to the jet" — the diameter represents the intact
        //  column cross-section, which must be consistent with the total mass.]
        d = Foam::cbrt(parcelMass/(nParticle*rhoPiOver6 + VSMALL));

        y = tc;                         // tSecStart
        yDot = d;                       // Dparent0 — now consistent with mass
        KHindex = max(Urmag, 1.0e-12);  // UrelPE0
        ms = -1.0;                      // first Madabhushi PE breakup

        // FIX-MB1: CSV log for column breakup latch — no limit (critical event)
        // FIX-MB7: log the updated d (post-reabsorb) so stage2 log reflects
        // the actual Dparent0 that will be used in PATH 4.
        if (debug_)
        {
            stage2LogFile()
                << tc << "," << tColumnBreakup << "," << rhoc << ","
                << UgMag << "," << d << "," << Urmag << "," << Urmag
                << Foam::endl;
        }
    }

    // =================================================================
    // PATH 4: PE stage
    // =================================================================
    if (y > 0.0)
    {
        const scalar tSecStart = y;
        const scalar Dparent0  = max(yDot, 1.0e-12);
        const scalar UrelPE0   = max(KHindex, 1.0e-12);
        const scalar tElapsed  = max(tc - tSecStart, 0.0);

        // Weber number frozen at PE onset.
        // [Madabhushi (2003); Lambert et al. (2019), Eqs. (2)-(3)]
        const scalar weInitPE =
            rhoc*sqr(UrelPE0)*Dparent0/(sigma + VSMALL);

        // Ohnesorge number at PE onset.
        // [Pilch & Erdman (1987), Eq. (2); Lambert et al. (2019), Eq. (8)]
        const scalar ohPE0 =
            mu/(Foam::sqrt(rho*Dparent0*sigma) + VSMALL);

        // Critical Weber number with viscous correction.
        // [Pilch & Erdman (1987), Eq. (5); Schmehl et al. (1998), Eq. (38)]
        const scalar weCritPE =
            12.0*(1.0 + 1.077*pow(ohPE0, 1.6));

        const bool isStandardPEChild = (ms > 1.5 && ms < 2.5);

        if (weInitPE > weCritPE)
        {
            // Characteristic breakup timescale t*.
            // [Pilch & Erdman (1987), Eq. (3); Lambert et al. (2019), Eq. (3)]
            const scalar tStar =
                (Dparent0/(UrelPE0 + VSMALL))
               *Foam::sqrt(rho/(rhoc + VSMALL));

            // Deformation time tDef = 1.6 t*.
            // [Schmehl et al. (1998), Eq. (42); Lambert et al. (2019), Eq. (2)]
            const scalar tDef = 1.6*tStar;

            if (tElapsed >= tDef)
            {
                // Deformed reference diameter.
                // [Schmehl et al. (1998), Eq. (45); Lambert et al. (2019), Eq. (9)]
                const scalar dRefPE =
                    (weInitPE < 100.0)
                  ? Dparent0*(1.0 + 0.19*Foam::sqrt(weInitPE))
                  : 2.9*Dparent0;

                // Ohnesorge number based on deformed reference diameter.
                // [Lambert et al. (2019), Eq. (8): Oh = mu/sqrt(rho*D_ref*sigma)]
                const scalar ohPE =
                    mu/(Foam::sqrt(rho*dRefPE*sigma) + VSMALL);

                // Corrected Weber number for viscous effects.
                // [Schmehl et al. (1998), Eq. (46); Lambert et al. (2019), Eq. (10)]
                const scalar weCorr =
                    weInitPE/(1.0 + 1.077*pow(ohPE, 1.6));

                // Target SMD after breakup.
                // [Schmehl et al. (1998), Eq. (48); Lambert et al. (2019), Eq. (7)]
                const scalar dSMD =
                    1.5
                   *(pow(max(ohPE, VSMALL), 0.2)
                   /(pow(max(weCorr, VSMALL), 0.25) + VSMALL))
                   *Dparent0;

                // FLig applies only to the first Madabhushi PE breakup of
                // the liquid column, identified by ms == -1.0.
                // This value is set exclusively by PATH 3 at the column-
                // breakup latch and is the only state for which FLig should
                // be active. The previous condition (ms < 0.0) was semantically
                // too broad: any future negative ms value would have incorrectly
                // activated FLig. The narrow range (ms > -1.5 && ms < -0.5)
                // pins the trigger to ms = -1.0 only. (FIX-MB2)
                // [Lambert et al. (2019), Eq. (11): "child droplets created
                //  immediately after column breakup"]
                const bool isFirstCorePEBreakup =
                    (ms > -1.5 && ms < -0.5) && !isStandardPEChild;

                const scalar FLigEff =
                    isFirstCorePEBreakup ? FLig_ : 1.0;

                // tb/t* from Pilch-Erdman piecewise correlations.
                // NOTE-3: strict-less-than boundaries are used here rather
                // than the less-or-equal (<=) of the original source. The
                // physical impact is negligible since the correlations are
                // continuous at the boundary values by construction.
                // [Pilch & Erdman (1987), Eqs. (8)-(12);
                //  Schmehl et al. (1998), Eq. (43);
                //  Madabhushi (2003), Eq. (5)]
                //
                // FIX-MB8: The tb/t* correlations are fitted against the
                // original Pilch & Erdman data using (We - 12) as the
                // independent variable, where 12 is the critical Weber
                // number in the original P&E formulation (without viscous
                // correction). Using (weInitPE - weCritPE) instead would
                // shift the argument of the power-law fits away from the
                // domain of their calibration for viscous liquids.
                // The viscosity-corrected weCritPE is used only for the
                // admission condition (weInitPE > weCritPE above); the
                // breakup-time correlations use the original (We - 12).
                // [Pilch & Erdman (1987), Eqs. (8)-(12);
                //  Madabhushi (2003), Eq. (5): explicitly uses (We-12)]
                scalar weExcess = max(weInitPE - 12.0, VSMALL);
                scalar tbOverTstar = 5.5;

                if      (weInitPE < 18.0)   tbOverTstar = 6.0  *pow(weExcess, -0.25);
                else if (weInitPE < 45.0)   tbOverTstar = 2.45 *pow(weExcess,  0.25);
                else if (weInitPE < 351.0)  tbOverTstar = 14.1 *pow(weExcess, -0.25);
                else if (weInitPE < 2670.0) tbOverTstar = 0.766*pow(weExcess,  0.25);

                const scalar tb = tbOverTstar*tStar;

                // FIX-MB1: CSV log for PE pre-breakup progress
                if (debug_ && !isStandardPEChild)
                {
                    peLogFile()
                        << tc << "," << tSecStart << "," << tElapsed << ","
                        << tDef << "," << tb << "," << Dparent0 << ","
                        << UrelPE0 << "," << weInitPE << "," << weCritPE << ","
                        << dSMD << "," << FLigEff << ","
                        << (isStandardPEChild ? 1 : 0)
                        << Foam::endl;
                }

                // =====================================================
                // Breakup event
                // =====================================================
                if (tElapsed >= tb)
                {
                    const vector Uparent0 = U;

                    // Reconstruct the local gas velocity at the breakup event.
                    // Urel is passed as Up - Ug, hence Ug = Up - Urel.
                    // This is later used to compute the real slip velocity of
                    // fragment #1 after its rim-expansion velocity is assigned.
                    const vector UgasLocal = Uparent0 - Urel;

                    const scalar totalMassToBreak =
                        nParticle*rhoPiOver6*pow3(d);

                    // Orthonormal basis for radial/rim expansion.
                    vector axialDir = Uparent0/(mag(Uparent0) + VSMALL);

                    vector refAxis =
                        (mag(axialDir.x()) < 0.7)
                      ? vector(1, 0, 0)
                      : vector(0, 1, 0);

                    vector normal1 = (axialDir ^ refAxis);
                    normal1 /= (mag(normal1) + VSMALL);

                    vector normal2 = (axialDir ^ normal1);
                    normal2 /= (mag(normal2) + VSMALL);

                    // Normal velocity magnitude from rim expansion.
                    // u_n = 5 * D_parent / (tb - tDef).
                    // No cap with UrelPE0 is applied: this follows the
                    // Madabhushi/Lambert formulation directly. VSMALL is kept
                    // only to avoid division by zero in degenerate timesteps.
                    // [Madabhushi (2003), text after Eq. (9);
                    //  Lambert et al. (2019), Eq. (5)]
                    const scalar uNormalMag =
                        5.0*Dparent0/(tb - tDef + VSMALL);

                    const label nFragments = max(nChildren_, label(2));

                    DynamicList<scalar> dFrag(nFragments);
                    DynamicList<scalar> massFrag(nFragments);
                    DynamicList<vector> UFrag(nFragments);

                    for (label i = 0; i < nFragments; ++i)
                    {
                        const scalar gaussianSample =
                            this->owner().rndGen().template GaussNormal<scalar>();

                        // Root-normal child-diameter distribution.
                        // [Lambert et al. (2019), Eq. (6);
                        //  Madabhushi (2003), Eq. (9)]
                        const scalar dBarChild =
                            1.2*dSMD*sqr(1.0 + 0.238*gaussianSample);

                        const scalar dTarget = FLigEff*dBarChild;

                        // Numerical guard only: avoid degenerate root-normal
                        // samples with vanishing diameter and unbounded
                        // nParticle. This is not intended as a physical floor.
                        dFrag.append(max(minChildDiameter_, min(dTarget, Dparent0)));

                        const scalar alpha =
                            2.0*pi*this->owner().rndGen().template sample01<scalar>();

                        UFrag.append
                        (
                            Uparent0
                          + uNormalMag
                           *(
                                cos(alpha)*normal1
                              + sin(alpha)*normal2
                            )
                        );
                    }

                    // -------------------------------------------------
                    // Convert sampled diameters into child parcel masses.
                    // Each child parcel represents the same number of physical
                    // droplets; parcel mass is therefore proportional to d_i^3.
                    // This conserves total mass without giving tiny fragments
                    // an artificially large nParticle.
                    // -------------------------------------------------
                    scalar sumD3 = 0.0;
                    forAll(dFrag, i)
                    {
                        sumD3 += pow3(dFrag[i]);
                    }

                    const scalar nChildPerParcel =
                        totalMassToBreak/(rhoPiOver6*sumD3 + VSMALL);

                    forAll(dFrag, i)
                    {
                        massFrag.append(nChildPerParcel*rhoPiOver6*pow3(dFrag[i]));
                    }

                    // FIX-MB1: CSV log for breakup event — no limit
                    if (debug_)
                    {
                        breakupLogFile()
                            << (isStandardPEChild ? "stdPE" : "Madabhushi") << ","
                            << tc << "," << Dparent0 << "," << UrelPE0 << ","
                            << tElapsed << "," << tDef << "," << tb << ","
                            << uNormalMag << "," << mag(Uparent0) << ","
                            << dSMD << "," << FLigEff << ","
                            << dFrag[0] << "," << mag(UFrag[0])
                            << Foam::endl;
                    }

                    // -------------------------------------------------
                    // Fragment #1 -> remains on the parent parcel
                    // -------------------------------------------------
                    d = dFrag[0];
                    parcelMass = massFrag[0];

                    const vector parentFinalVelocity = UFrag[0];

                    if (isStandardPEChild)
                    {
                        // FIX-MB11: relatch must use the actual slip velocity
                        // of fragment #1 after assigning the rim-expansion
                        // velocity, not the pre-breakup parent Urmag.
                        const vector UrelFrag0 = parentFinalVelocity - UgasLocal;
                        const scalar UrmagFrag0 = mag(UrelFrag0);

                        const scalar weCurrRelatch =
                            rhoc*sqr(UrmagFrag0)*d/(sigma + VSMALL);

                        const scalar ohCurrRelatch =
                            mu/(Foam::sqrt(rho*d*sigma) + VSMALL);

                        const scalar weCritRelatch =
                            12.0*(1.0 + 1.077*pow(ohCurrRelatch, 1.6));

                        bool allowRelatch = false;
                        scalar dStableRelatch = GREAT;

                        if (weCurrRelatch > weCritRelatch)
                        {
                            scalar taubBarRelatch = 5.5;

                            if (weCurrRelatch < 2670.0)
                            {
                                if (weCurrRelatch > 351.0)
                                {
                                    taubBarRelatch =
                                        0.766*pow
                                        (
                                            max(weCurrRelatch - 12.0, VSMALL),
                                            0.25
                                        );
                                }
                                else if (weCurrRelatch > 45.0)
                                {
                                    taubBarRelatch =
                                        14.1*pow
                                        (
                                            max(weCurrRelatch - 12.0, VSMALL),
                                            -0.25
                                        );
                                }
                                else if (weCurrRelatch > 18.0)
                                {
                                    taubBarRelatch =
                                        2.45*pow
                                        (
                                            max(weCurrRelatch - 12.0, VSMALL),
                                            0.25
                                        );
                                }
                                else if (weCurrRelatch > 12.0)
                                {
                                    taubBarRelatch =
                                        6.0*pow
                                        (
                                            max(weCurrRelatch - 12.0, VSMALL),
                                            -0.25
                                        );
                                }
                            }

                            const scalar rho12Relatch =
                                Foam::sqrt(max(rhoc/(rho + VSMALL), VSMALL));

                            // B1PE=0.375, B2PE=0.2274 [OpenFOAM PilchErdman.C]
                            const scalar VdRelatch =
                                UrmagFrag0*rho12Relatch
                               *(0.375*taubBarRelatch
                               + 0.2274*sqr(taubBarRelatch));

                            // Same stability logic used in PATH 1:
                            // clamp the velocity-history factor before squaring.
                            const scalar Vd1Relatch =
                                max
                                (
                                    sqr
                                    (
                                        max
                                        (
                                            1.0 - VdRelatch/(UrmagFrag0 + VSMALL),
                                            0.0
                                        )
                                    ),
                                    SMALL
                                );

                            // Maximum stable diameter uses We_c = 12.
                            // The viscous correction is used only in the
                            // admission criterion weCurrRelatch > weCritRelatch.
                            dStableRelatch =
                                12.0*sigma
                               /(Vd1Relatch*rhoc*sqr(UrmagFrag0) + VSMALL);

                            allowRelatch = (d > dStableRelatch);
                        }

                        if (debug_)
                        {
                            stdPeLogFile()
                                << "RELATCH_CHECK,Frag0,"
                                << tc << "," << d << "," << UrmagFrag0 << ","
                                << weCurrRelatch << "," << weCritRelatch << ","
                                << dStableRelatch << ","
                                << (allowRelatch ? 1 : 0)
                                << Foam::endl;
                        }

                        if (allowRelatch)
                        {
                            // Fragment still unstable: relatch for new PE cycle.
                            y       = tc;
                            yDot    = d;
                            KHindex = UrmagFrag0;
                            ms      = 2.0;
                        }
                        else
                        {
                            // Fragment stable or sub-critical: becomes post-cat
                            // sentinel. PATH 1 will re-evaluate on subsequent
                            // steps if local conditions change.
                            y       = 0.0;
                            yDot    = 0.0;
                            KHindex = 0.0;
                            ms      = -20.0;
                        }
                    }
                    else
                    {
                        y = 0.0;
                        yDot = 0.0;
                        KHindex = 0.0;
                        ms = -20.0;
                    }

                    // The retained parent parcel is now also a post-catastrophic
                    // fragment. Keep the user_ tag consistent for diagnostics
                    // and for any user-based branching.
                    parentUserUpdate_ = 2.0;

                    // -------------------------------------------------
                    // Fragment #2 -> returned through addParcel
                    // -------------------------------------------------
                    dChild = dFrag[1];
                    massChild = massFrag[1];

                    {
                        // [FIX Bug2] Set U to UFrag[1] for clone.
                        vector& Umutable = const_cast<vector&>(U);
                        Umutable = UFrag[1];
                    }

                    addParcel = true;
                    childMsInit_ = -20.0;
                    childUserInit_ = 2.0;

                    // [FIX Bug2] Sentinel to restore parent velocity after clone.
                    pendingChildUserFlag_ = 2.0;

                    pendingChildren_.append
                    (
                        pendingChild(-1.0, 0.0, parentFinalVelocity)
                    );

                    // Fragment #3..n -> pending children
                    for (label i = 2; i < nFragments; ++i)
                    {
                        pendingChildren_.append
                        (
                            pendingChild(dFrag[i], massFrag[i], UFrag[i])
                        );
                    }

                    nParticle = parcelMass/(pow3(d)*rhoPiOver6 + VSMALL);
                    return addParcel;
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // Default behavior: no breakup event in this time step
    // -----------------------------------------------------------------
    nParticle = parcelMass/(pow3(d)*rhoPiOver6 + VSMALL);
    return false;
}
