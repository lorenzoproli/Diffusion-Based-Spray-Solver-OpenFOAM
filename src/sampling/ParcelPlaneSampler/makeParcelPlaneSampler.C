#include "basicSprayCloud.H"
#include "ParcelPlaneSampler.H"
#include "CloudFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    makeCloudFunctionObjectType(ParcelPlaneSampler, basicSprayCloud);
}
