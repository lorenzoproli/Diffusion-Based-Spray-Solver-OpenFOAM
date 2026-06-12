// Copyright (C) 2026 Lorenzo Proli, Politecnico di Torino
// SPDX-License-Identifier: GPL-3.0-or-later
#include "basicSprayCloud.H"
#include "Madabhushi.H"

// Registra semplicemente il tuo modello nella selection table esistente
namespace Foam
{
    makeBreakupModelType(Madabhushi, basicSprayCloud);
}
