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

/// \file MCVertexParam.h
/// \brief The smearing of the MC vertex list: the MC vertex reader's
///        parameter, GTRKMCVertex.

#ifndef GTRK_STUDY_MCVERTEXPARAM_H_
#define GTRK_STUDY_MCVERTEXPARAM_H_

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

#include <cstdint>

namespace o2::gtrk::study
{

/// The smearing of the MC vertex list (the MC vertex reader, which the study
/// workflow adds with `--mc-vertices`). The published variances are these values squared; its
/// time is the collision's, a box of width 0. The fixed seed makes the resulting vertex lists reproducible.
// O2ParamDef deletes the copy constructor and supplies the rest of the boilerplate.
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
struct MCVertexParam : public o2::conf::ConfigurableParamHelper<MCVertexParam> {
  double sigmaZ = 0.01;           ///< [cm]
  double sigmaTransverse = 0.002; ///< [cm], independently in x and y
  std::uint64_t randomSeed = 0x4754524b56455254ULL;

  O2ParamDef(MCVertexParam, "GTRKMCVertex");
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_MCVERTEXPARAM_H_
