#include "basicSprayCloud.H" 
#include "Madabhushi.H"

// Registra semplicemente il tuo modello nella selection table esistente
namespace Foam
{
    makeBreakupModelType(Madabhushi, basicSprayCloud);
}
