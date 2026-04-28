/*---------------------------------------------------------------------------*\
  ParcelPlaneSampler — implementation.
  Polygon-plane intersection logic adapted from Foam::ParticleCollector
  (OpenFOAM v2406, OpenCFD Ltd.).
\*---------------------------------------------------------------------------*/

#include "ParcelPlaneSampler.H"
#include "Pstream.H"
#include "triangle.H"
#include "cloud.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::ParcelPlaneSampler<CloudType>::makeLogFile()
{
    if (!logToFile_) return;

    if (Pstream::master())
    {
        mkDir(this->writeTimeDir());

        outputFilePtr_.reset
        (
            new OFstream(this->writeTimeDir()/(type() + ".dat"))
        );

        outputFilePtr_()
            << "# Source     : " << type() << nl
            << "# Planes     : " << faces_.size() << nl;

        outputFilePtr_() << "# Geometry   :" << nl;
        forAll(faces_, i)
        {
            outputFilePtr_()
                << "#   plane " << i
                << "  centre=" << faces_[i].centre(points_)
                << "  normal=" << normal_[i] << nl;
        }

        outputFilePtr_()
            << "#" << nl
            << "# Output (one row per parcel transit):" << nl
            << "#   time  x  y  z  Ux  Uy  Uz  d  nParticle  user  ms  tc  KHindex  yState  yDotState  faceID" << nl;
    }
}


template<class CloudType>
void Foam::ParcelPlaneSampler<CloudType>::initPolygons
(
    const List<Field<point>>& polygons
)
{
    label nPoints = 0;
    forAll(polygons, polyI)
    {
        if (polygons[polyI].size() < 3)
        {
            FatalIOErrorInFunction(this->coeffDict())
                << "polygons must have at least 3 points"
                << exit(FatalIOError);
        }
        nPoints += polygons[polyI].size();
    }

    points_.setSize(nPoints);
    faces_.setSize(polygons.size());

    label off = 0;
    forAll(faces_, fi)
    {
        const Field<point>& pp = polygons[fi];
        face f(identity(pp.size(), off));
        UIndirectList<point>(points_, f) = pp;
        faces_[fi].transfer(f);
        off += pp.size();
    }
}


template<class CloudType>
void Foam::ParcelPlaneSampler<CloudType>::sampleSegment
(
    const parcelType& p,
    const point& p1,
    const point& p2,
    const scalar dt
)
{
    forAll(faces_, fi)
    {
        const point& pf = points_[faces_[fi][0]];
        const scalar d1 = normal_[fi] & (p1 - pf);
        const scalar d2 = normal_[fi] & (p2 - pf);

        if (sign(d1) == sign(d2)) continue;  // no crossing this step

        if (negateParcelsOppositeNormal_)
        {
            // Skip parcels travelling against +normal
            const scalar Un = normal_[fi] & p.U();
            if (Un < 0) continue;
        }

        // Plane intersection point
        const scalar denom = d1 - d2;
        if (mag(denom) < SMALL) continue;
        const scalar frac = d1/denom;
        const point pHit  = p1 + frac*(p2 - p1);

        // Inside-polygon test (same as ParticleCollector)
        const face& f = faces_[fi];
        const vector areaNorm = f.areaNormal(points_);
        bool inside = true;
        for (label i = 0; i < f.size(); ++i)
        {
            const label j = f.fcIndex(i);
            const triPointRef t(pHit, points_[f[i]], points_[f[j]]);
            if ((areaNorm & t.areaNormal()) < 0)
            {
                inside = false;
                break;
            }
        }
        if (!inside) continue;

        // Time of crossing within dt (linear interp)
        const scalar tNow = this->owner().mesh().time().value();
        const scalar tHit = tNow - dt + frac*dt;

        evTime_.append(tHit);
        evFace_.append(fi);
        evX_.append(pHit.x());
        evY_.append(pHit.y());
        evZ_.append(pHit.z());
        evUx_.append(p.U().x());
        evUy_.append(p.U().y());
        evUz_.append(p.U().z());
        evD_.append(p.d());
        evN_.append(p.nParticle());
        evUser_.append(p.user());
        evMs_.append(p.ms());
        evTc_.append(p.tc());
        evKHindex_.append(p.KHindex());
        evYstate_.append(p.y());
        evYDotState_.append(p.yDot());
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::ParcelPlaneSampler<CloudType>::write()
{
    // Gather all events to master
    List<List<scalar>> aT(Pstream::nProcs());
    List<List<scalar>> aX(Pstream::nProcs());
    List<List<scalar>> aY(Pstream::nProcs());
    List<List<scalar>> aZ(Pstream::nProcs());
    List<List<scalar>> aUx(Pstream::nProcs());
    List<List<scalar>> aUy(Pstream::nProcs());
    List<List<scalar>> aUz(Pstream::nProcs());
    List<List<scalar>> aD(Pstream::nProcs());
    List<List<scalar>> aN(Pstream::nProcs());
    List<List<scalar>> aUser(Pstream::nProcs());
    List<List<scalar>> aMs(Pstream::nProcs());
    List<List<scalar>> aTc(Pstream::nProcs());
    List<List<scalar>> aKH(Pstream::nProcs());
    List<List<scalar>> aYstate(Pstream::nProcs());
    List<List<scalar>> aYDotState(Pstream::nProcs());
    List<List<label>>  aF(Pstream::nProcs());

    aT[Pstream::myProcNo()]  = evTime_;
    aX[Pstream::myProcNo()]  = evX_;
    aY[Pstream::myProcNo()]  = evY_;
    aZ[Pstream::myProcNo()]  = evZ_;
    aUx[Pstream::myProcNo()] = evUx_;
    aUy[Pstream::myProcNo()] = evUy_;
    aUz[Pstream::myProcNo()] = evUz_;
    aD[Pstream::myProcNo()]  = evD_;
    aN[Pstream::myProcNo()]  = evN_;
    aUser[Pstream::myProcNo()] = evUser_;
    aMs[Pstream::myProcNo()] = evMs_;
    aTc[Pstream::myProcNo()] = evTc_;
    aKH[Pstream::myProcNo()] = evKHindex_;
    aYstate[Pstream::myProcNo()] = evYstate_;
    aYDotState[Pstream::myProcNo()] = evYDotState_;
    aF[Pstream::myProcNo()]  = evFace_;

    Pstream::gatherList(aT);
    Pstream::gatherList(aX);
    Pstream::gatherList(aY);
    Pstream::gatherList(aZ);
    Pstream::gatherList(aUx);
    Pstream::gatherList(aUy);
    Pstream::gatherList(aUz);
    Pstream::gatherList(aD);
    Pstream::gatherList(aN);
    Pstream::gatherList(aUser);
    Pstream::gatherList(aMs);
    Pstream::gatherList(aTc);
    Pstream::gatherList(aKH);
    Pstream::gatherList(aYstate);
    Pstream::gatherList(aYDotState);
    Pstream::gatherList(aF);

    if (Pstream::master() && logToFile_)
    {
        if (!outputFilePtr_) makeLogFile();

        // Flatten and sort by time, write each row
        DynamicList<label> order;
        DynamicList<scalar> T;
        DynamicList<scalar> X, Y, Z, Ux, Uy, Uz, D, N;
        DynamicList<scalar> User, Ms, Tc, KH, Ystate, YDotState;
        DynamicList<label>  F;
        forAll(aT, proci)
        {
            forAll(aT[proci], k)
            {
                T.append(aT[proci][k]);
                X.append(aX[proci][k]);
                Y.append(aY[proci][k]);
                Z.append(aZ[proci][k]);
                Ux.append(aUx[proci][k]);
                Uy.append(aUy[proci][k]);
                Uz.append(aUz[proci][k]);
                D.append(aD[proci][k]);
                N.append(aN[proci][k]);
                User.append(aUser[proci][k]);
                Ms.append(aMs[proci][k]);
                Tc.append(aTc[proci][k]);
                KH.append(aKH[proci][k]);
                Ystate.append(aYstate[proci][k]);
                YDotState.append(aYDotState[proci][k]);
                F.append(aF[proci][k]);
            }
        }

        // Stable sort by time so the file is readable chronologically
        labelList idx(identity(T.size()));
        std::sort
        (
            idx.begin(), idx.end(),
            [&T](label a, label b){ return T[a] < T[b]; }
        );

        OFstream& os = outputFilePtr_();
        os.precision(10);
        forAll(idx, i)
        {
            const label k = idx[i];
            os  << T[k] << ' '
                << X[k] << ' ' << Y[k] << ' ' << Z[k] << ' '
                << Ux[k] << ' ' << Uy[k] << ' ' << Uz[k] << ' '
                << D[k] << ' ' << N[k] << ' '
                << User[k] << ' ' << Ms[k] << ' ' << Tc[k] << ' '
                << KH[k] << ' ' << Ystate[k] << ' ' << YDotState[k] << ' '
                << F[k] << nl;
        }
        os.flush();

        Log_<< type() << ": wrote " << T.size()
            << " parcel transit events" << endl;
    }

    // Clear local buffers (reset on every write)
    evTime_.clear(); evFace_.clear();
    evX_.clear();    evY_.clear();    evZ_.clear();
    evUx_.clear();   evUy_.clear();   evUz_.clear();
    evD_.clear();    evN_.clear();
    evUser_.clear(); evMs_.clear(); evTc_.clear();
    evKHindex_.clear(); evYstate_.clear(); evYDotState_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParcelPlaneSampler<CloudType>::ParcelPlaneSampler
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    points_(),
    faces_(),
    normal_(),
    parcelType_(this->coeffDict().template getOrDefault<label>("parcelType", -1)),
    negateParcelsOppositeNormal_
    (
        this->coeffDict().template getOrDefault<bool>
        (
            "negateParcelsOppositeNormal", false
        )
    ),
    evTime_(), evFace_(),
    evX_(), evY_(), evZ_(),
    evUx_(), evUy_(), evUz_(),
    evD_(), evN_(),
    evUser_(), evMs_(), evTc_(),
    evKHindex_(), evYstate_(), evYDotState_(),
    logToFile_(this->coeffDict().template getOrDefault<bool>("log", true)),
    outputFilePtr_()
{
    List<Field<point>> polygons(this->coeffDict().lookup("polygons"));
    initPolygons(polygons);

    vector n0(this->coeffDict().lookup("normal"));
    n0.normalise();
    normal_ = vectorField(faces_.size(), n0);

    makeLogFile();
}


template<class CloudType>
Foam::ParcelPlaneSampler<CloudType>::ParcelPlaneSampler
(
    const ParcelPlaneSampler<CloudType>& pps
)
:
    CloudFunctionObject<CloudType>(pps),
    points_(pps.points_),
    faces_(pps.faces_),
    normal_(pps.normal_),
    parcelType_(pps.parcelType_),
    negateParcelsOppositeNormal_(pps.negateParcelsOppositeNormal_),
    evTime_(), evFace_(),
    evX_(), evY_(), evZ_(),
    evUx_(), evUy_(), evUz_(),
    evD_(), evN_(),
    evUser_(), evMs_(), evTc_(),
    evKHindex_(), evYstate_(), evYDotState_(),
    logToFile_(pps.logToFile_),
    outputFilePtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ParcelPlaneSampler<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point& position0,
    const typename parcelType::trackingData& /*td*/
)
{
    if (parcelType_ != -1 && parcelType_ != p.typeId()) return true;

    sampleSegment(p, position0, p.position(), dt);
    return true;  // never remove parcels
}


// ************************************************************************* //
