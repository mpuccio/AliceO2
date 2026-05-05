// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "ALICE3GlobalReconstructionWorkflow/MagneticFieldHelper.h"

#include "DetectorsBase/Propagator.h"
#include "Field/MagFieldParam.h"
#include "Field/MagneticField.h"
#include "Framework/Logger.h"

#include <TGeoGlobalMagField.h>

namespace o2::trk
{

void ensureAlice3Field(float bz)
{
  auto* globalField = TGeoGlobalMagField::Instance();
  if (globalField->GetField() == nullptr) {
    auto* field = new o2::field::MagneticField("ALICE3Mag", "ALICE 3 Magnetic Field", bz / 5.f, 0.0, o2::field::MagFieldParam::k5kGUniform);
    globalField->SetField(field);
    globalField->Lock();
    LOGP(info, "Installed ALICE3 magnetic field with Bz={} kG", bz);
  }
  auto* propagator = o2::base::Propagator::Instance(false);
  propagator->updateField();
  propagator->setNominalBz(bz);
}

} // namespace o2::trk
