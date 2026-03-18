#include "basicSprayCloud.H"
#include "MadabhushiDragForce.H"

namespace Foam
{
    makeParticleForceModelType(MadabhushiDragForce, basicSprayCloud);
}
