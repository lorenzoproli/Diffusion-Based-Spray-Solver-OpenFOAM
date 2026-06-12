/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2026 Lorenzo Proli, Politecnico di Torino
    SPDX-License-Identifier: GPL-3.0-or-later
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Modifications (Thesis — Madabhushi crossflow model):
    FIX-SP1: Removed erroneous modification of parent parcel user_ flag before
             cloning. The block that set this->user()=2.0 or 1.0 on the parent
             before calling new SprayParcel(*this) caused the parent core parcel
             to be reclassified as a child, disabling blob drag and preventing
             correct crossflow deflection. The child already receives the correct
             user_ via childUserInit_ (set by the breakup model) at line
             child->user() = cloud.breakup().childUserInit().
             [Implementation-state convention; see Madabhushi.C state machine]

    FIX-SP2: Added structured CSV debug output to sprayParcelDebug.log when
             the breakup model debug flag is active. One line per cloning event
             with parent and child state for post-run diagnostics.

\*---------------------------------------------------------------------------*/

#include "SprayParcel.H"
#include "BreakupModel.H"
#include "CompositionModel.H"
#include "AtomizationModel.H"
#include "OFstream.H"
#include "autoPtr.H"

// -------------------------------------------------------------------------
// FIX-SP2: structured debug log for cloning events
// -------------------------------------------------------------------------
namespace
{
    Foam::OFstream& sprayParcelLogFile()
    {
        static Foam::autoPtr<Foam::OFstream> logPtr;
        if (!logPtr.valid())
        {
            logPtr.reset(new Foam::OFstream("sprayParcelDebug.log"));
            // CSV header
            logPtr() << "event,"
                     << "parent_user,parent_ms,parent_d,parent_Umag,"
                     << "child_user,child_ms,child_d,child_Umag,"
                     << "parent_tc,parent_age"
                     << Foam::endl;
        }
        return logPtr();
    }
}

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::SprayParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::SprayParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    ParcelType::cellValueSourceCorrection(cloud, td, dt);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::SprayParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const auto& composition = cloud.composition();
    const auto& liquids = composition.liquids();

    // Check if parcel belongs to liquid core
    if (liquidCore() > 0.5)
    {
        // Liquid core parcels should not experience coupled forces
        cloud.forces().setCalcCoupled(false);
    }

    // Get old mixture composition
    scalarField X0(liquids.X(this->Y()));

    // Check if we have critical or boiling conditions
    scalar TMax = liquids.Tc(X0);
    const scalar T0 = this->T();
    const scalar pc0 = td.pc();
    if (liquids.pv(pc0, T0, X0) >= pc0*0.999)
    {
        // Set TMax to boiling temperature
        TMax = liquids.pvInvert(pc0, X0);
    }

    // Set the maximum temperature limit
    cloud.constProps().setTMax(TMax);

    // Store the parcel properties
    this->Cp() = liquids.Cp(pc0, T0, X0);
    sigma_ = liquids.sigma(pc0, T0, X0);
    const scalar rho0 = liquids.rho(pc0, T0, X0);
    this->rho() = rho0;
    const scalar mass0 = this->mass();
    mu_ = liquids.mu(pc0, T0, X0);

    ParcelType::calc(cloud, td, dt);

    if (td.keepParticle)
    {
        // Reduce the stripped parcel mass due to evaporation only when ms
        // represents a physical accumulated stripped mass, not a state flag
        if (this->ms() >= 0.0)
        {
            this->ms() -= this->ms()*(mass0 - this->mass())/(mass0 + ROOTVSMALL);
        }

        // Update Cp, sigma, density and diameter due to change in temperature
        // and/or composition
        scalar T1 = this->T();
        scalarField X1(liquids.X(this->Y()));

        this->Cp() = liquids.Cp(td.pc(), T1, X1);

        sigma_ = liquids.sigma(td.pc(), T1, X1);

        scalar rho1 = liquids.rho(td.pc(), T1, X1);
        this->rho() = rho1;

        mu_ = liquids.mu(td.pc(), T1, X1);

        scalar d1 = this->d()*cbrt(rho0/rho1);
        this->d() = d1;

        if (liquidCore() > 0.5)
        {
            calcAtomization(cloud, td, dt);

            // Preserve the total mass/volume by increasing the number of
            // particles in parcels due to breakup
            scalar d2 = this->d();
            this->nParticle() *= pow3(d1/d2);
        }
        else
        {
            calcBreakup(cloud, td, dt);
        }
    }

    // Restore coupled forces
    cloud.forces().setCalcCoupled(true);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::SprayParcel<ParcelType>::calcAtomization
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const auto& atomization = cloud.atomization();

    if (!atomization.active())
    {
        return;
    }

    const auto& composition = cloud.composition();
    const auto& liquids = composition.liquids();

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = td.rhoc()*RR*td.Tc()/td.pc();
    scalar R = RR/Wc;
    scalar Tav = atomization.Taverage(this->T(), td.Tc());

    // Calculate average gas density based on average temperature
    scalar rhoAv = td.pc()/(R*Tav);

    scalar soi = cloud.injectors().timeStart();
    scalar currentTime = cloud.db().time().value();
    const vector pos(this->position());
    const vector& injectionPos = this->position0();

    // Disregard the continuous phase when calculating the relative velocity
    // (in line with the deactivated coupled assumption)
    scalar Urel = mag(this->U());

    scalar t0 = max(0.0, currentTime - this->age() - soi);
    scalar t1 = min(t0 + dt, cloud.injectors().timeEnd() - soi);

    // This should be the vol flow rate from when the parcel was injected
    scalar volFlowRate = cloud.injectors().volumeToInject(t0, t1)/dt;

    scalar chi = 0.0;
    if (atomization.calcChi())
    {
        chi = this->chi(cloud, td, liquids.X(this->Y()));
    }

    atomization.update
    (
        dt,
        this->d(),
        this->liquidCore(),
        this->tc(),
        this->rho(),
        mu_,
        sigma_,
        volFlowRate,
        rhoAv,
        Urel,
        pos,
        injectionPos,
        cloud.pAmbient(),
        chi,
        cloud.rndGen()
    );
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::SprayParcel<ParcelType>::calcBreakup
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    auto& breakup = cloud.breakup();

    if (!breakup.active())
    {
        return;
    }

    if (breakup.solveOscillationEq())
    {
        solveTABEq(cloud, td, dt);
    }

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = td.rhoc()*RR*td.Tc()/td.pc();
    scalar R = RR/Wc;
    scalar Tav = cloud.atomization().Taverage(this->T(), td.Tc());

    // Calculate average gas density based on average temperature
    scalar rhoAv = td.pc()/(R*Tav);
    scalar muAv = td.muc();
    vector Urel = this->U() - td.Uc();
    scalar Urmag = mag(Urel);
    scalar Re = this->Re(rhoAv, this->U(), td.Uc(), this->d(), muAv);

    const typename TrackCloudType::parcelType& p =
        static_cast<const typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);
    const scalar mass = p.mass();
    const typename TrackCloudType::forceType& forces = cloud.forces();
    const forceSuSp Fcp = forces.calcCoupled(p, ttd, dt, mass, Re, muAv);
    const forceSuSp Fncp = forces.calcNonCoupled(p, ttd, dt, mass, Re, muAv);
    this->tMom() = mass/(Fcp.Sp() + Fncp.Sp() + ROOTVSMALL);

    const vector g = cloud.g().value();

    scalar parcelMassChild = 0.0;
    scalar dChild = 0.0;
    if
    (
        breakup.update
        (
            dt,
            g,
            this->d(),
            this->tc(),
            this->ms(),
            this->nParticle(),
            this->KHindex(),
            this->y(),
            this->yDot(),
            this->d0(),
            this->rho(),
            mu_,
            sigma_,
            this->U(),
            td.rhoc(),
            muAv,
            Urel,
            Urmag,
            this->tMom(),
            dChild,
            parcelMassChild
        )
    )
    {
        // Apply an optional user_ update requested by the breakup model
        // after the breakup state has been computed. This avoids changing
        // the parent user_ before the model decides the current regime, but
        // still lets the retained parent fragment be reclassified correctly.
        const scalar parentUserUpdate = cloud.breakup().parentUserUpdate();
        if (parentUserUpdate > -0.5)
        {
            this->user() = parentUserUpdate;
        }

        // ----------------------------------------------------------------
        // FIX-SP1: The block that previously modified this->user() here
        // (setting it to 1.0 or 2.0 based on ms()) has been removed.
        //
        // REASON: Modifying the parent's user_ flag BEFORE cloning it
        // caused the parent core parcel to be reclassified as a child
        // parcel. This disabled blob drag on the parent (MadabhushiDragForce
        // branches on userFlag > 0.5 to select spherical drag instead of
        // blob drag), which prevented the liquid column from acquiring
        // streamwise momentum from the crossflow and therefore suppressed
        // the correct jet deflection trajectory.
        //
        // The child already receives the correct user_ through:
        //   child->user() = cloud.breakup().childUserInit();
        // which is set by the breakup model before returning addParcel=true.
        // The parent must retain user_=0 to stay in the blob/core drag path.
        // ----------------------------------------------------------------

        // ----------------------------------------------------------------
        // FIX-SP2: structured debug log — parent state at cloning time
        // ----------------------------------------------------------------
        if (breakup.debug)
        {
            sprayParcelLogFile()
                << "CLONE_FIRST_CHILD,"
                << this->user()   << ","
                << this->ms()     << ","
                << this->d()      << ","
                << mag(this->U()) << ","
                << cloud.breakup().childUserInit() << ","
                << cloud.breakup().childMsInit()   << ","
                << dChild                          << ","
                << mag(this->U())                  << ","
                << this->tc()                      << ","
                << this->age()
                << Foam::endl;
        }

        // Add child parcel as copy of parent
        SprayParcel<ParcelType>* child = new SprayParcel<ParcelType>(*this);
        child->origId() = this->getNewParticleID();
        child->origProc() = Pstream::myProcNo();
        child->d() = dChild;
        child->d0() = dChild;
        const scalar massChild = child->mass();
        child->mass0() = massChild;
        child->nParticle() = parcelMassChild/(massChild + ROOTVSMALL);

        // Recompute relative velocity and Reynolds for the new child
        const vector UrelChild = child->U() - td.Uc();
        const scalar UrmagChild = mag(UrelChild);
        const scalar ReChild = rhoAv*UrmagChild*dChild/muAv;

        const forceSuSp Fcp =
            forces.calcCoupled(*child, ttd, dt, massChild, ReChild, muAv);
        const forceSuSp Fncp =
            forces.calcNonCoupled(*child, ttd, dt, massChild, ReChild, muAv);

        child->age() = 0.0;
        child->liquidCore() = 0.0;
        // KHindex is repurposed as UrelPE0 (latched relative velocity at
        // PE onset) by the Madabhushi breakup model. Initialising to 0.0
        // ensures the child enters PATH 1 (sentinel latch) with a clean
        // state and that PATH 4 (PE stage) is not entered prematurely
        // due to a stale non-zero value from the parent copy constructor.
        // [Implementation-state convention; see Madabhushi.C PATH 1/PATH 4]
        child->KHindex() = 0.0;
        child->y() = cloud.breakup().y0();
        child->yDot() = cloud.breakup().yDot0();
        child->tc() = 0.0;
        child->ms() = cloud.breakup().childMsInit();
        child->injector() = this->injector();
        child->tMom() = massChild/(Fcp.Sp() + Fncp.Sp() + ROOTVSMALL);
        child->user() = cloud.breakup().childUserInit();
        child->calcDispersion(cloud, td, dt);

        cloud.addParticle(child);

        // --- Multi-child support: create extra children from buffer ---
        if (cloud.breakup().hasPendingChildren())
        {
            const label nPending = cloud.breakup().nPendingChildren();

            for (label ip = 0; ip < nPending; ++ip)
            {
                scalar pd = -1.0;
                scalar pmass = 0.0;
                vector pvel = Zero;

                cloud.breakup().getPendingChild(ip, pd, pmass, pvel);

                // Sentinel d < 0: restore parent velocity
                if (pd < 0.0)
                {
                    this->U() = pvel;

                    // FIX-SP2: log parent velocity restore
                    if (breakup.debug)
                    {
                        sprayParcelLogFile()
                            << "PARENT_VEL_RESTORE,"
                            << this->user()   << ","
                            << this->ms()     << ","
                            << this->d()      << ","
                            << mag(this->U()) << ","
                            << ",,,,,"
                            << this->tc()     << ","
                            << this->age()
                            << Foam::endl;
                    }
                    continue;
                }

                // Create extra child as copy of parent
                SprayParcel<ParcelType>* extraChild =
                    new SprayParcel<ParcelType>(*this);

                extraChild->origId() = this->getNewParticleID();
                extraChild->origProc() = Pstream::myProcNo();
                extraChild->d() = pd;
                extraChild->d0() = pd;

                const scalar mExtraChild = extraChild->mass();
                extraChild->mass0() = mExtraChild;
                extraChild->nParticle() = pmass/(mExtraChild + ROOTVSMALL);

                extraChild->U() = pvel;

                extraChild->age() = 0.0;
                extraChild->liquidCore() = 0.0;
                // Same rationale as child->KHindex() = 0.0 above:
                // initialise to 0.0 so PATH 1 performs a clean latch.
                // [Implementation-state convention; see Madabhushi.C PATH 1]
                extraChild->KHindex() = 0.0;
                extraChild->y() = cloud.breakup().y0();
                extraChild->yDot() = cloud.breakup().yDot0();
                extraChild->tc() = 0.0;
                extraChild->ms() = -20.0;
                extraChild->injector() = this->injector();
                extraChild->user() = cloud.breakup().pendingChildUserFlag();

                const vector UrelExtra = extraChild->U() - td.Uc();
                const scalar UrmagExtra = mag(UrelExtra);
                const scalar ReExtra = rhoAv*UrmagExtra*pd/muAv;

                const forceSuSp FcpExtra =
                    forces.calcCoupled(*extraChild, ttd, dt, pmass, ReExtra, muAv);

                const forceSuSp FncpExtra =
                    forces.calcNonCoupled(*extraChild, ttd, dt, pmass, ReExtra, muAv);

                extraChild->tMom() = pmass/(FcpExtra.Sp() + FncpExtra.Sp() + ROOTVSMALL);
                extraChild->calcDispersion(cloud, td, dt);

                // FIX-SP2: log extra child creation
                if (breakup.debug)
                {
                    sprayParcelLogFile()
                        << "CLONE_EXTRA_CHILD,"
                        << this->user()                            << ","
                        << this->ms()                              << ","
                        << this->d()                               << ","
                        << mag(this->U())                          << ","
                        << cloud.breakup().pendingChildUserFlag()  << ","
                        << "-20.0,"
                        << pd                                      << ","
                        << mag(pvel)                               << ","
                        << this->tc()                              << ","
                        << this->age()
                        << Foam::endl;
                }

                cloud.addParticle(extraChild);
            }

            cloud.breakup().clearPendingChildren();
        }
    }
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::SprayParcel<ParcelType>::chi
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalarField& X
) const
{
    // Modifications to take account of the flash boiling on primary break-up

    const auto& composition = cloud.composition();
    const auto& liquids = composition.liquids();

    scalar chi = 0.0;
    scalar T0 = this->T();
    scalar p0 = td.pc();
    scalar pAmb = cloud.pAmbient();

    scalar pv = liquids.pv(p0, T0, X);

    forAll(liquids, i)
    {
        if (pv >= 0.999*pAmb)
        {
            // Liquid is boiling - calc boiling temperature

            const liquidProperties& liq = liquids.properties()[i];
            scalar TBoil = liq.pvInvert(p0);

            scalar hl = liq.hl(pAmb, TBoil);
            scalar iTp = liq.h(pAmb, T0) - pAmb/liq.rho(pAmb, T0);
            scalar iTb = liq.h(pAmb, TBoil) - pAmb/liq.rho(pAmb, TBoil);

            chi += X[i]*(iTp - iTb)/hl;
        }
    }

    return clamp(chi, zero_one{});
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::SprayParcel<ParcelType>::solveTABEq
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const scalar& TABCmu = cloud.breakup().TABCmu();
    const scalar& TABtwoWeCrit = cloud.breakup().TABtwoWeCrit();
    const scalar& TABComega = cloud.breakup().TABComega();

    scalar r = 0.5*this->d();
    scalar r2 = r*r;
    scalar r3 = r*r2;

    // Inverse of characteristic viscous damping time
    scalar rtd = 0.5*TABCmu*mu_/(this->rho()*r2);

    // Oscillation frequency (squared)
    scalar omega2 = TABComega*sigma_/(this->rho()*r3) - rtd*rtd;

    if (omega2 > 0)
    {
        scalar omega = sqrt(omega2);
        scalar We =
            this->We(td.rhoc(), this->U(), td.Uc(), r, sigma_)/TABtwoWeCrit;

        // Initial values for y and yDot
        scalar y0 = this->y() - We;
        scalar yDot0 = this->yDot() + y0*rtd;

        // Update distortion parameters
        scalar c = cos(omega*dt);
        scalar s = sin(omega*dt);
        scalar e = exp(-rtd*dt);

        this->y() = We + e*(y0*c + (yDot0/omega)*s);
        this->yDot() = (We - this->y())*rtd + e*(yDot0*c - omega*y0*s);
    }
    else
    {
        // Reset distortion parameters
        this->y() = 0;
        this->yDot() = 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SprayParcel<ParcelType>::SprayParcel(const SprayParcel<ParcelType>& p)
:
    ParcelType(p),
    d0_(p.d0_),
    position0_(p.position0_),
    sigma_(p.sigma_),
    mu_(p.mu_),
    liquidCore_(p.liquidCore_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_)
{}


template<class ParcelType>
Foam::SprayParcel<ParcelType>::SprayParcel
(
    const SprayParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    d0_(p.d0_),
    position0_(p.position0_),
    sigma_(p.sigma_),
    mu_(p.mu_),
    liquidCore_(p.liquidCore_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "SprayParcelIO.C"


// ************************************************************************* //
