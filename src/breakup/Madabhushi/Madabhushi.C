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
\*---------------------------------------------------------------------------*/

#include "Madabhushi.H"
#include "OFstream.H"
#include "mathematicalConstants.H"
#include "Pstream.H"
#include "autoPtr.H"
#include <map>
#include <fstream>

namespace
{
    Foam::OFstream& breakupLogFile()
    {
        static Foam::autoPtr<Foam::OFstream> logPtr;

        if (!logPtr.valid())
        {
            logPtr.reset(new Foam::OFstream("breakupDebug.log"));
        }

        return logPtr();
    }


    inline void writeDebugLimited
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
            breakupLogFile() << "[" << tag << "] " << msg << Foam::endl;
        }
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
    sigma_(this->coeffDict().getOrDefault("sigma", 0.072)),
    nChildren_(this->coeffDict().getOrDefault("nChildren", 5)),
    UgRef_(this->coeffDict().getOrDefault("UgRef", 62.5)),
    debug_(this->coeffDict().getOrDefault("debug", false)),
    childMsInit_(-GREAT),
    childUserInit_(0.0),
    pendingChildren_(),
    pendingChildUserFlag_(0.0)
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
    sigma_(model.sigma_),
    nChildren_(model.nChildren_),
    UgRef_(model.UgRef_),
    debug_(model.debug_),
    childMsInit_(-GREAT),
    childUserInit_(0.0),
    pendingChildren_(),
    pendingChildUserFlag_(0.0)
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

    const scalar pi = constant::mathematical::pi;
    const scalar rhoPiOver6 = rho*pi/6.0;

    scalar parcelMass = nParticle*pow3(d)*rhoPiOver6;

    // -----------------------------------------------------------------
    // Before PE onset, yDot is used to store the initial mass of the
    // intact core parcel only once. After PE onset, yDot is repurposed
    // to store Dparent0.
    // This is an implementation choice, not a physical model equation.
    // -----------------------------------------------------------------
    if (y <= 0.0 && yDot <= 0.0)
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
    //   ms =  1.0:          subsequent Madabhushi PE breakup (FLig = 1)
    //   ms =  2.0:          standard PE child (spherical drag before latch)
    //   ms = -10.0:         sentinel: wave-shed child, needs STD_PE check
    //   ms = -20.0:         sentinel: post-cat child, needs STD_PE check
    //   ms =  3.0:          inert droplet, no further breakup
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
    const scalar tColumnBreakup =
        C0_*(Dinj_/(UgRef_ + VSMALL))*Foam::sqrt(rho/(rhoc + VSMALL));

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
    // ms = -10.0 : wave-shed child
    // ms = -20.0 : post-catastrophic child
    //
    // Standard PE is allowed if:
    //   We > We_crit
    //   d  > dMinStdPENumerical
    //
    // No dStable filter is used here because it was too restrictive for
    // the crossflow child population.
    // =================================================================

    const bool isWaveSentinel    = (ms > -11.0 && ms < -9.0);
    const bool isPostCatSentinel = (ms > -21.0 && ms < -19.0);

    if (isWaveSentinel || isPostCatSentinel)
    {
        // Current Weber number based on the local, instantaneous child state.
        // Same Weber-number definition as Pilch-Erdman and Lambert:
        // We = rho_g * U_rel^2 * D / sigma.
        // [Pilch & Erdman (1987), Eq. (1); Lambert et al. (2019), text and Eqs. (8)-(10)]
        const scalar weCurr =
            rhoc*sqr(Urmag)*d/(sigma_ + VSMALL);

        // Ohnesorge number based on the local child state.
        // [Pilch & Erdman (1987), Eq. (2); Lambert et al. (2019), Eq. (8)]
        const scalar ohCurr =
            mu/(Foam::sqrt(rho*d*sigma_ + VSMALL));

        // Critical Weber number for breakup onset.
        // [Pilch & Erdman (1987), Eq. (5); Schmehl et al. (1998), Eq. (38);
        //  Lambert et al. (2019), onset criterion within PE framework]
        const scalar weCritCurr =
            12.0*(1.0 + 1.077*pow(ohCurr, 1.6));

        // Numerical safeguard only, not a physical model parameter.
        const scalar dMinStdPENumerical = 1.0e-5;

        const bool allowStdPE =
        (
            weCurr > weCritCurr
         && d > dMinStdPENumerical
        );

        writeDebugLimited
        (
            debug_,
            "STD_PE_CHECK",
            "STD_PE_CHECK "
          + Foam::string("type=") + (isWaveSentinel ? "Wave" : "PostCat")
          + " tc=" + Foam::name(tc)
          + " d=" + Foam::name(d)
          + " Urmag=" + Foam::name(Urmag)
          + " WeCurr=" + Foam::name(weCurr)
          + " WeCritCurr=" + Foam::name(weCritCurr)
          + " dMinStdPE=" + Foam::name(dMinStdPENumerical)
          + " allow=" + Foam::name(allowStdPE ? 1 : 0),
            50
        );

        if (!allowStdPE)
        {
            // Child remains inert if onset conditions are not met.
            ms = 3.0;
            y = 0.0;
            yDot = 0.0;
            KHindex = 0.0;
            nParticle = parcelMass/(pow3(d)*rhoPiOver6 + VSMALL);
            return false;
        }

        // Latch standard PE onset using the current local state.
        y = tc;            // tSecStart
        yDot = d;          // Dparent0
        KHindex = Urmag;   // UrelPE0
        ms = 2.0;          // standard PE child

        writeDebugLimited
        (
            debug_,
            "STD_PE",
            "STD_PE_LATCH "
          + Foam::string("type=") + (isWaveSentinel ? "Wave" : "PostCat")
          + " tc=" + Foam::name(tc)
          + " d=" + Foam::name(d)
          + " Urmag=" + Foam::name(Urmag)
          + " WeCurr=" + Foam::name(weCurr)
          + " WeCritCurr=" + Foam::name(weCritCurr),
            50
        );

        nParticle = parcelMass/(pow3(d)*rhoPiOver6 + VSMALL);
        return false;
    }

    // =================================================================
    // PATH 2: stage 1 — intact liquid column / KH wave shedding
    // =================================================================
    if (y <= 0.0 && tc < tColumnBreakup)
    {
        // Gas-side Weber number in the wave model.
        // [Madabhushi (2003), nomenclature and wave-model section; Reitz (1987)]
        const scalar weGasWave =
            rhoc*sqr(Urmag)*injectorRadius/(sigma_ + VSMALL);

        // Liquid-side Weber number in the wave model.
        // [Madabhushi (2003), nomenclature and wave-model section; Reitz (1987)]
        const scalar weLiqWave =
            rho*sqr(Urmag)*injectorRadius/(sigma_ + VSMALL);

        // Minimum gas-side Weber number required for KH stripping onset.
        // This threshold belongs to the wave-stage model description.
        // [Madabhushi (2003), wave-shedding stage]
        const scalar weWaveCrit = 6.0;

        // Fraction of initial core mass required before emitting one
        // discrete child parcel from the continuous stripped mass budget.
        // This is an implementation choice used to discretize shedding.
        // Madabhushi states a 5% shed-mass criterion as a computational feature.
        // [Madabhushi (2003), text in wave-model section]
        const scalar strippedMassFractionThreshold = 0.05;

        // Liquid Reynolds number used in the KH wave model.
        // [Madabhushi (2003), wave-model section; Reitz (1987)]
        const scalar reLiquid =
            rho*Urmag*injectorRadius/(mu + VSMALL);

        // Ohnesorge number in Reitz/Madabhushi KH formulation:
        // Z = sqrt(We_l) / Re_l.
        // [Madabhushi (2003), nomenclature; Reitz (1987), after Eq. (5)]
        const scalar ohWave =
            Foam::sqrt(max(weLiqWave, VSMALL))/(reLiquid + VSMALL);

        // Taylor parameter in Reitz/Madabhushi KH formulation:
        // T = Z * sqrt(We_g).
        // [Madabhushi (2003), nomenclature; Reitz (1987)]
        const scalar taWave =
            ohWave*Foam::sqrt(max(weGasWave, VSMALL));

        // Most unstable KH wave growth rate.
        // [Reitz (1987), Eq. (5); Madabhushi (2003), Eq. (4)]
        const scalar omegaKH =
            (0.34 + 0.38*pow(weGasWave, 1.5))
           /((1.0 + ohWave)*(1.0 + 1.4*pow(taWave, 0.6)))
           *Foam::sqrt(sigma_/(rho*pow3(injectorRadius) + VSMALL));

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

        // Characteristic child diameter generated by the wave model:
        // D_child = 2 * B0 * Lambda, with B0 = 0.61.
        // [Reitz (1987), Eq. (10a) with B0; Madabhushi (2003), text and Eq. (2)]
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
            // [Reitz (1987), Eq. (11); Madabhushi (2003), Eq. (3)-based radius evolution]
            d = (relaxationFraction*dChildWaveTarget + dOld)
              /(1.0 + relaxationFraction);

            const scalar massStrippedThisStep =
                nParticle*rhoPiOver6*(pow3(dOld) - pow3(d));

            ms += massStrippedThisStep;

            const scalar initialWaveMass = max(yDot, parcelMass);

            // Discrete child emission criterion from accumulated stripped mass.
            // This is a computational discretization rule, not a standalone model equation.
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

                // Empirical angular-spread coefficient for initial off-axis
                // child velocity after wave shedding.
                // [Reitz (1987), Eq. (7) with A1; Madabhushi (2003), text using A1 = 0.188;
                //  Spray Trajectories of Liquid Fuel Jets in Subsonic Crossflows (1997), Eq. (10)]
                const scalar A1 = 0.188;

                // Spray half-angle relation used to assign child transverse velocity.
                // [Reitz (1987), Eq. (7); Madabhushi (2003), text after Eq. (4)]
                const scalar tanThetaHalf =
                    A1*lambdaKH*omegaKH/(mag(U) + VSMALL);

                const scalar v1Mag = mag(U)*tanThetaHalf*sin(phi);
                const scalar v2Mag = mag(U)*tanThetaHalf*cos(phi);

                const vector childVelocityPerturbation =
                    v1Mag*normal1 + v2Mag*normal2;

                const vector parentVelocity = U;
                vector& Umutable = const_cast<vector&>(U);

                // Child velocity is written before cloning so the spawned
                // parcel inherits the intended kinematics.
                Umutable = parentVelocity + childVelocityPerturbation;

                addParcel = true;
                dChild = dChildWaveTarget;
                massChild = childMassTarget;

                // Mark child as wave-shed sentinel for later STD_PE check.
                childMsInit_ = -10.0;
                childUserInit_ = 1.0;

                ms -= massChild;
                parcelMass -= massChild;

                if (parcelMass < SMALL)
                {
                    parcelMass = SMALL;
                }

                // Restore parent velocity: the core parcel remains on the
                // intact-column trajectory.
                Umutable = parentVelocity;
            }
        }

        nParticle = parcelMass/(pow3(d)*rhoPiOver6 + VSMALL);
        return addParcel;
    }

    // =================================================================
    // PATH 3: latch PE onset exactly once for the core parcel
    // =================================================================
    if (y <= 0.0 && tc >= tColumnBreakup)
    {
        y = tc;                         // tSecStart
        yDot = d;                       // Dparent0
        KHindex = max(Urmag, 1.0e-12);  // UrelPE0
        ms = -1.0;                      // first Madabhushi PE breakup

        writeDebugLimited
        (
            debug_,
            "ENTER_STAGE2",
            "ENTER_STAGE2 "
          + Foam::string("tc=") + Foam::name(tc)
          + " tColumnBreakup=" + Foam::name(tColumnBreakup)
          + " rhoc=" + Foam::name(rhoc)
          + " UgRef=" + Foam::name(UgRef_)
          + " tSecStart=" + Foam::name(y)
          + " Dparent0=" + Foam::name(yDot)
          + " UrelPE0=" + Foam::name(KHindex),
            150
        );
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
        // [Madabhushi (2003), text after column breakup; Lambert et al. (2019), Eqs. (2)-(3)]
        const scalar weInitPE =
            rhoc*sqr(UrelPE0)*Dparent0/(sigma_ + VSMALL);

        const scalar weCritPE = 12.0;

        // Standard PE children use ms = 2.0, while the Madabhushi core uses
        // negative/other values before its own first breakup.
        const bool isStandardPEChild = (ms > 1.5 && ms < 2.5);

        if (weInitPE > weCritPE)
        {
            // Characteristic breakup timescale t*.
            // [Pilch & Erdman (1987), Eq. (3); Lambert et al. (2019), Eq. (3);
            //  Madabhushi (2003), text before Eq. (5)]
            const scalar tStar =
                (Dparent0/(UrelPE0 + VSMALL))
               *Foam::sqrt(rho/(rhoc + VSMALL));

            // Deformation time tDef = 1.6 t*.
            // [Schmehl et al. (1998), Eq. (42); Lambert et al. (2019), Eq. (2);
            //  Madabhushi (2003), text before Eq. (5)]
            const scalar tDef = 1.6*tStar;

            if (tElapsed >= tDef)
            {
                // Deformed reference diameter.
                // [Schmehl et al. (1998), Eq. (45); Lambert et al. (2019), Eq. (9);
                //  Madabhushi (2003), Eq. (6)]
                const scalar dRefPE =
                    (weInitPE < 100.0)
                  ? Dparent0*(1.0 + 0.19*Foam::sqrt(weInitPE))
                  : 2.9*Dparent0;

                // Ohnesorge number based on the deformed reference diameter.
                // [Pilch & Erdman (1987), Eq. (2); Schmehl et al. (1998), Eq. (37);
                //  Lambert et al. (2019), Eq. (8)]
                const scalar ohPE =
                    mu/(Foam::sqrt(rho*dRefPE*sigma_) + VSMALL);

                // Corrected Weber number for viscous effects.
                // [Schmehl et al. (1998), Eq. (46); Lambert et al. (2019), Eq. (10);
                //  Madabhushi (2003), text after Eq. (8)]
                const scalar weCorr =
                    weInitPE/(1.0 + 1.077*pow(ohPE, 1.6));

                // Target SMD after breakup.
                // [Schmehl et al. (1998), Eq. (48); Lambert et al. (2019), Eq. (7);
                //  Madabhushi (2003), Eq. (8)]
                const scalar dSMD =
                    1.5
                   *(pow(max(ohPE, VSMALL), 0.2)
                   /(pow(max(weCorr, VSMALL), 0.25) + VSMALL))
                   *Dparent0;

                // FLig applies only to the first Madabhushi core breakup,
                // not to standard-PE children and not to DSMD itself.
                // [Lambert et al. (2019), Eq. (11)]
                const bool isFirstCorePEBreakup =
                    (ms < 0.0) && !isStandardPEChild;

                const scalar FLigEff =
                    isFirstCorePEBreakup ? FLig_ : 1.0;

                // Dimensionless total breakup time tb/t* from Pilch-Erdman
                // piecewise correlations.
                // [Pilch & Erdman (1987), Eqs. (8)-(12);
                //  Schmehl et al. (1998), Eq. (43);
                //  Madabhushi (2003), Eq. (5)]
                scalar weExcess = max(weInitPE - weCritPE, VSMALL);
                scalar tbOverTstar = 5.5;

                if (weInitPE < 18.0)
                {
                    tbOverTstar = 6.0*pow(weExcess, -0.25);
                }
                else if (weInitPE < 45.0)
                {
                    tbOverTstar = 2.45*pow(weExcess, 0.25);
                }
                else if (weInitPE < 351.0)
                {
                    tbOverTstar = 14.1*pow(weExcess, -0.25);
                }
                else if (weInitPE < 2670.0)
                {
                    tbOverTstar = 0.766*pow(weExcess, 0.25);
                }

                const scalar tb = tbOverTstar*tStar;

                if (!isStandardPEChild)
                {
                    writeDebugLimited
                    (
                        debug_,
                        "CORE_PREBREAKUP",
                        "CORE_PREBREAKUP "
                      + Foam::string("tc=") + Foam::name(tc)
                      + " tStartPE=" + Foam::name(tSecStart)
                      + " tElapsed=" + Foam::name(tElapsed)
                      + " Dparent0=" + Foam::name(Dparent0)
                      + " UrelPE0=" + Foam::name(UrelPE0)
                      + " WeInitPE=" + Foam::name(weInitPE)
                      + " tDef=" + Foam::name(tDef)
                      + " tb=" + Foam::name(tb),
                        200
                    );
                }

                // =====================================================
                // Breakup event
                // =====================================================
                if (tElapsed >= tb && d > 2.0e-6)
                {
                    const vector Uparent0 = U;
                    const scalar totalMassToBreak =
                        nParticle*rhoPiOver6*pow3(d);

                    // Build local orthonormal basis for radial/rim expansion.
                    vector axialDir = Uparent0/(mag(Uparent0) + VSMALL);

                    vector refAxis =
                        (mag(axialDir.x()) < 0.7)
                      ? vector(1, 0, 0)
                      : vector(0, 1, 0);

                    vector normal1 = (axialDir ^ refAxis);
                    normal1 /= (mag(normal1) + VSMALL);

                    vector normal2 = (axialDir ^ normal1);
                    normal2 /= (mag(normal2) + VSMALL);

                    // Characteristic normal velocity scale assigned to fragments
                    // after breakup. This corresponds to the rim-expansion-based
                    // closure used in Madabhushi/Lambert:
                    // u_n = 5 D_parent / (tb - tDef).
                    // [Lambert et al. (2019), Eq. (5); Madabhushi (2003), text after Eq. (9)]
                    const scalar uNormalMag =
                        min(5.0*Dparent0/(tb - tDef + VSMALL), UrelPE0);

                    const label nFragments = nChildren_;
                    const scalar massPerFragment =
                        totalMassToBreak/scalar(nFragments);

                    DynamicList<scalar> dFrag(nFragments);
                    DynamicList<vector> UFrag(nFragments);

                    for (label i = 0; i < nFragments; ++i)
                    {
                        const scalar gaussianSample =
                            this->owner().rndGen().template GaussNormal<scalar>();

                        // Root-normal child-diameter distribution centered on D0.5.
                        // [Schmehl et al. (1998), Eqs. (49)-(51);
                        //  Lambert et al. (2019), Eq. (6);
                        //  Madabhushi (2003), Eq. (9)]
                        const scalar dBarChild =
                            1.2*dSMD*sqr(1.0 + 0.238*gaussianSample);

                        const scalar dTarget =
                            FLigEff*dBarChild;

                        dFrag.append(max(2.0e-6, min(dTarget, Dparent0)));

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

                    const word breakupTag =
                        isStandardPEChild
                      ? "BREAKUP_DONE_STDPE"
                      : "BREAKUP_DONE_CORE";

                    writeDebugLimited
                    (
                        debug_,
                        breakupTag,
                        "BREAKUP_DONE "
                      + Foam::string("mode=")
                      + (isStandardPEChild ? "stdPE" : "Madabhushi")
                      + " nFrag=" + Foam::name(nFragments)
                      + " tc=" + Foam::name(tc)
                      + " Dparent0=" + Foam::name(Dparent0)
                      + " UrelPE0=" + Foam::name(UrelPE0)
                      + " tElapsed=" + Foam::name(tElapsed)
                      + " tDef=" + Foam::name(tDef)
                      + " tb=" + Foam::name(tb)
                      + " u_n_mag=" + Foam::name(uNormalMag)
                      + " UmotherMag=" + Foam::name(mag(Uparent0))
                      + " DSMD=" + Foam::name(dSMD)
                      + " FLigEff=" + Foam::name(FLigEff)
                      + " dFrag0=" + Foam::name(dFrag[0])
                      + " UFrag0Mag=" + Foam::name(mag(UFrag[0])),
                        200
                    );

                    // -------------------------------------------------
                    // Fragment #1 -> remains on the parent parcel
                    // -------------------------------------------------
                    d = dFrag[0];
                    parcelMass = massPerFragment;

                    {
                        vector& Umutable = const_cast<vector&>(U);
                        Umutable = UFrag[0];
                    }

                    if (isStandardPEChild)
                    {
                        // Standard-PE child remains PE-active after breakup.
                        y = tc;
                        yDot = d;
                        KHindex = Urmag;
                        ms = 2.0;
                    }
                    else
                    {
                        // Madabhushi core fragments after first catastrophic
                        // breakup become post-cat sentinels.
                        y = 0.0;
                        yDot = 0.0;
                        KHindex = 0.0;
                        ms = -20.0;
                    }

                    // -------------------------------------------------
                    // Fragment #2 -> returned through addParcel
                    // -------------------------------------------------
                    dChild = dFrag[1];
                    massChild = massPerFragment;

                    {
                        vector& Umutable = const_cast<vector&>(U);
                        Umutable = UFrag[1];
                    }

                    addParcel = true;
                    childMsInit_ = -20.0;
                    childUserInit_ = 2.0;

                    // -------------------------------------------------
                    // Fragment #3..n -> stored as pending children
                    // -------------------------------------------------
                    pendingChildUserFlag_ = 2.0;

                    for (label i = 2; i < nFragments; ++i)
                    {
                        pendingChildren_.append
                        (
                            pendingChild(dFrag[i], massPerFragment, UFrag[i])
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
