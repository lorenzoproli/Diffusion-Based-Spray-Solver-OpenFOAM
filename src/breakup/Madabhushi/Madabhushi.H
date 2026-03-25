/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
-------------------------------------------------------------------------------
    Madabhushi breakup model with Lambert ligament extension.
    Modified to support:
      - Multi-child breakup (nChildren parcels at catastrophic event)
      - Explicit flags for Wave-shed and post-catastrophic children
      - Standard Pilch-Erdman for flagged children
\*---------------------------------------------------------------------------*/

#ifndef Madabhushi_H
#define Madabhushi_H

#include "BreakupModel.H"
#include "DynamicList.H"

namespace Foam
{

template<class CloudType>
class Madabhushi
:
    public BreakupModel<CloudType>
{
public:

    // Public data types

        //- Data for a pending child parcel (multi-child catastrophic breakup)
        struct pendingChild
        {
            scalar d;
            scalar mass;
            vector velocity;

            pendingChild()
            :
                d(0),
                mass(0),
                velocity(Zero)
            {}

            pendingChild(scalar dd, scalar mm, const vector& vv)
            :
                d(dd),
                mass(mm),
                velocity(vv)
            {}
        };


private:

    // Private data

        // Column-breakup constant
        scalar C0_;

        // Ligament weighting factor (applied ONLY at first PE breakup)
        scalar FLig_;

        // Reitz Wave constants
        scalar b0_;
        scalar b1_;

        // Injector diameter
        scalar Dinj_;

        // Surface tension
        scalar sigma_;

        //- Number of child parcels at catastrophic breakup (default 5)
        label nChildren_;

        scalar UgRef_;

        Switch debug_;

        // --- Mutable state for child signaling to calcBreakup ---

        //- ms sentinel for the first child returned by addParcel=true
        //  Replaces the hardcoded -GREAT in calcBreakup.
        //    -GREAT  = normal (OpenFOAM default)
        //    -10.0   = Wave-shed child → needs standard PE
        //    -20.0   = post-catastrophic child → needs standard PE
        //     0.0    = generic reset
        mutable scalar childMsInit_;

        //- user flag for the first child returned by addParcel=true.
        //  Replaces the hardcoded 0.0 in calcBreakup.
        //    0.0 = core parcel (Madabhushi drag)
        //    1.0 = Wave-shed child (spherical drag)
        //    2.0 = post-catastrophic child (spherical drag)
        mutable scalar childUserInit_;

        //- Buffer for extra children beyond the first (indices 2..nChildren_)
        mutable DynamicList<pendingChild> pendingChildren_;

        //- user flag for pending children
        mutable scalar pendingChildUserFlag_;


public:

    //- Runtime type information
    TypeName("Madabhushi");


    // Constructors

        Madabhushi
        (
            const dictionary& dict,
            CloudType& owner
        );

        Madabhushi(const Madabhushi<CloudType>& bum);

        virtual autoPtr<BreakupModel<CloudType>> clone() const
        {
            return autoPtr<BreakupModel<CloudType>>
            (
                new Madabhushi<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Madabhushi();


    // Member Functions

        // --- Standard BreakupModel interface ---

        virtual bool update
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
        );


        // --- Multi-child / flag interface ---
        //     (to be called by a MODIFIED calcBreakup in SprayParcel.C)

        //- Override: ms value for child parcel
        virtual scalar childMsInit() const { return childMsInit_; }

        //- Override: user value for child parcel
        virtual scalar childUserInit() const { return childUserInit_; }

        //- Override: does the model have pending children?
        virtual bool hasPendingChildren() const
        {
            return pendingChildren_.size() > 0;
        }

        //- Number of pending children in the buffer
        //  FIX-BUG3: Added virtual/override to explicitly mark this as an
        //  override of BreakupModel::nPendingChildren().
        virtual label nPendingChildren() const override { return pendingChildren_.size(); }

        //- Access pending child i (internal convenience, not virtual)
        const pendingChild& getPendingChild(label i) const
        {
            return pendingChildren_[i];
        }

        //- Override base-class 4-argument version called by SprayParcel.C
        //  FIX-BUG1: Without this override the base-class default returns
        //  velocity=Zero, causing the parent velocity restore sentinel to
        //  zero-out the parent U and all extra children to be lost.
        virtual void getPendingChild
        (
            const label i,
            scalar& d,
            scalar& mass,
            vector& velocity
        ) const override
        {
            const pendingChild& c = pendingChildren_[i];
            d        = c.d;
            mass     = c.mass;
            velocity = c.velocity;
        }

        //- user flag for pending children
        scalar pendingChildUserFlag() const { return pendingChildUserFlag_; }

        //- Clear the pending children buffer (after calcBreakup consumes them)
        //  FIX-BUG2: Removed const qualifier to match the non-const base-class
        //  virtual signature in BreakupModel.H. With const, this was a separate
        //  function (not an override) and the base-class no-op was called instead.
        virtual void clearPendingChildren() override { pendingChildren_.clear(); }
};

} // End namespace Foam

#ifdef NoRepository
    #include "Madabhushi.C"
#endif

#endif
