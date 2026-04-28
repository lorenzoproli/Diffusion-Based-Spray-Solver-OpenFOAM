/*---------------------------------------------------------------------------*\
  ParcelPlaneSampler

  Cloud function object that records every parcel transit through one or more
  user-defined polygonal planes. For each crossing event one row is written
  to <writeTimeDir>/parcelPlaneSampler.dat with the columns:

      time   x   y   z   Ux   Uy   Uz   d   nParticle   user   ms   tc   KHindex   yState   yDotState

  Designed to reproduce the DPM-plane sampling methodology of
  Lambert et al. (ILASS 2019) — i.e. per-parcel data needed to compute
  nParticle-weighted SMD(y) and Ux(y) profiles bin-by-bin.

  Polygon/plane intersection logic is taken from Foam::ParticleCollector
  (OpenFOAM v2406, OpenCFD Ltd.).
\*---------------------------------------------------------------------------*/

#ifndef ParcelPlaneSampler_H
#define ParcelPlaneSampler_H

#include "CloudFunctionObject.H"
#include "face.H"
#include "OFstream.H"
#include "DynamicList.H"

namespace Foam
{

template<class CloudType>
class ParcelPlaneSampler
:
    public CloudFunctionObject<CloudType>
{
    typedef typename CloudType::parcelType parcelType;

    // Geometry (one polygon = one sampling plane)

        Field<point>  points_;
        List<face>    faces_;
        Field<vector> normal_;

    // Filter

        const label parcelType_;
        bool negateParcelsOppositeNormal_;

    // Per-event buffers (filled in postMove, flushed in write)

        DynamicList<scalar> evTime_;
        DynamicList<label>  evFace_;
        DynamicList<scalar> evX_;
        DynamicList<scalar> evY_;
        DynamicList<scalar> evZ_;
        DynamicList<scalar> evUx_;
        DynamicList<scalar> evUy_;
        DynamicList<scalar> evUz_;
        DynamicList<scalar> evD_;
        DynamicList<scalar> evN_;
        DynamicList<scalar> evUser_;
        DynamicList<scalar> evMs_;
        DynamicList<scalar> evTc_;
        DynamicList<scalar> evKHindex_;
        DynamicList<scalar> evYstate_;
        DynamicList<scalar> evYDotState_;

    // Output

        bool logToFile_;
        autoPtr<OFstream> outputFilePtr_;

    // Private member functions

        void makeLogFile();

        //- Initialise polygons from "polygons" list
        void initPolygons(const List<Field<point>>& polygons);

        //- Test segment p1->p2 against face facei. If it crosses, append
        //  intersection-time event to the buffers.
        void sampleSegment
        (
            const parcelType& p,
            const point& p1,
            const point& p2,
            const scalar dt
        );


protected:

    void write();


public:

    TypeName("parcelPlaneSampler");

    // Constructors

        ParcelPlaneSampler
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        ParcelPlaneSampler(const ParcelPlaneSampler<CloudType>& pps);

        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ParcelPlaneSampler<CloudType>(*this)
            );
        }

    virtual ~ParcelPlaneSampler() = default;


    // Member functions

        virtual bool postMove
        (
            parcelType& p,
            const scalar dt,
            const point& position0,
            const typename parcelType::trackingData& td
        );
};

} // End namespace Foam

#ifdef NoRepository
    #include "ParcelPlaneSampler.C"
#endif

#endif
