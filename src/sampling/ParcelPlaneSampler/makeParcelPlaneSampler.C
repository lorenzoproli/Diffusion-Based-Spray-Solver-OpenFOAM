// Copyright (C) 2026 Lorenzo Proli, Politecnico di Torino
// SPDX-License-Identifier: GPL-3.0-or-later
#include "basicSprayCloud.H"
#include "ParcelPlaneSampler.H"
#include "CloudFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    makeCloudFunctionObjectType(ParcelPlaneSampler, basicSprayCloud);
}
