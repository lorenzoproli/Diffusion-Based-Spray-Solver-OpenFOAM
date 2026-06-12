/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2026 Lorenzo Proli, Politecnico di Torino
    SPDX-License-Identifier: GPL-3.0-or-later
-------------------------------------------------------------------------------
    Madabhushi custom drag force with flag-based regime selection.

    Modified: checks p.user() to distinguish core parcels from children.
      user_ = 0  →  Madabhushi deformation/disc drag (core parcels)
      user_ = 1  →  spherical drag (Wave-shed children)
      user_ = 2  →  spherical drag (post-catastrophic children)
\*---------------------------------------------------------------------------*/

#ifndef MadabhushiDragForce_H
#define MadabhushiDragForce_H

#include "ParticleForce.H"
#include "Switch.H"

namespace Foam
{

template<class CloudType>
class MadabhushiDragForce
:
    public ParticleForce<CloudType>
{
    // Private data

        scalar CdBlob_;
        scalar CdDisc_;
        scalar C0_;
        scalar Dinj_;
        const Switch debug_;


    // Private Member Functions

        //- Standard sphere drag as function of Re
        inline scalar CdSphere(const scalar Re) const
        {
            if (Re > 1000.0)
            {
                return 0.424;
            }
            return 24.0/(Re + VSMALL)
                  *(1.0 + (1.0/6.0)*pow(Re, 2.0/3.0));
        }


public:

    TypeName("MadabhushiDragForce");


    // Constructors

        MadabhushiDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        MadabhushiDragForce
        (
            const MadabhushiDragForce<CloudType>& df
        );

        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new MadabhushiDragForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MadabhushiDragForce();


    // Member Functions

        virtual forceSuSp calcCoupled
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const scalar dt,
            const scalar mass,
            const scalar Re,
            const scalar muc
        ) const;

        virtual forceSuSp calcNonCoupled
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const scalar dt,
            const scalar mass,
            const scalar Re,
            const scalar muc
        ) const;
};

} // End namespace Foam

#ifdef NoRepository
    #include "MadabhushiDragForce.C"
#endif

#endif
