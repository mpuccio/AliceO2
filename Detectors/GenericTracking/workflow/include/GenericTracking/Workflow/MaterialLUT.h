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

/// \file MaterialLUT.h
/// \brief The tracking's material LUT, a CCDB condition per detector.
///
/// | detector | CCDB path |
/// |---|---|
/// | ITS | `GLO/Param/MatLUT` (O2's) |
/// | MFT | `MFT/Param/MatLUT` (`o2-gtrk-build-matlut -d MFT`) |
/// | TRKFT3 | `TRK/Param/MatLUT` (`o2-gtrk-build-matlut -d TRKFT3`, one LUT for TRK and FT3) |
///
/// `MaterialCondition<Det>` is the condition input `Conditions<Det, true>`
/// holds (Conditions.h): the devices that track fetch it, the measurement
/// producer and the MC vertex reader (`Conditions<Det>`) hold
/// `NoMaterialCondition` in its place and fetch nothing. The LUT arrives
/// through `finaliseCCDB`, rectified as GRPGeomHelper rectifies
/// `GLO/Param/MatLUT`, and the CCDB cache owns it. Offline, from a CCDB
/// snapshot through DPL's `--condition-remap`; a program without DPL (the
/// bench) reads the same path through `o2::ccdb::BasicCCDBManager`.

#ifndef GTRK_WORKFLOWS_MATERIALLUT_H_
#define GTRK_WORKFLOWS_MATERIALLUT_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Workflow/DataSpecs.h"

#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/MatLayerCylSet.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InputRecord.h"
#include "Framework/InputSpec.h"
#include "Framework/ProcessingContext.h"
#include "Headers/DataHeader.h"

#include <fairlogger/Logger.h>

#include <concepts>
#include <string>
#include <vector>

namespace o2::gtrk::workflow
{

/// The data origin of \a Det's material LUT: the detector's, but the ITS's,
/// O2's global one, GLO.
template <DetectorLayout Det>
[[nodiscard]] o2::header::DataOrigin materialLutOrigin()
{
  return std::same_as<Det, detectors::ITS> ? o2::header::DataOrigin{"GLO"} : detectorOrigin<Det>();
}

/// The CCDB path of \a Det's material LUT.
template <DetectorLayout Det>
[[nodiscard]] std::string materialLutPath()
{
  return materialLutOrigin<Det>().template as<std::string>() + "/Param/MatLUT";
}

/// \a Det's material LUT as a condition input (see above).
template <DetectorLayout Det>
class MaterialCondition
{
 public:
  /// Declares the input: called last by the holding `Conditions`'
  /// constructor, so that the LUT's input follows the others.
  static void declare(std::vector<o2::framework::InputSpec>& inputs) { inputs.emplace_back(inputBinding<Det>(kBinding), materialLutOrigin<Det>(), "MATLUT", 0, o2::framework::Lifetime::Condition, o2::framework::ccdbParamSpec(materialLutPath<Det>())); }

  /// At the top of the device's run: delivers a new one through finaliseCCDB.
  static void checkUpdates(o2::framework::ProcessingContext& pc) { pc.inputs().get<o2::base::MatLayerCylSet*>(inputBinding<Det>(kBinding).c_str()); }

  /// Whether \a matcher is the LUT, which is then taken.
  bool finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj)
  {
    if (!(matcher == o2::framework::ConcreteDataMatcher(materialLutOrigin<Det>(), "MATLUT", 0))) {
      return false;
    }
    mLut = o2::base::MatLayerCylSet::rectifyPtrFromFile(static_cast<o2::base::MatLayerCylSet*>(obj));
    LOGP(info, "{}'s material LUT from the CCDB ({})", Det::name, materialLutPath<Det>());
    return true;
  }

  /// The CCDB cache's; null before the first timeframe.
  [[nodiscard]] const o2::base::MatLayerCylSet* lut() const { return mLut; }

 private:
  static constexpr const char* kBinding = "gtrkMatLUT";

  const o2::base::MatLayerCylSet* mLut = nullptr;
};

/// The ITS's material LUT inside the GPU reconstruction: GLO/Param/MatLUT,
/// which the reconstruction's device requests and GRPGeomHelper holds,
/// rectified. No input of its own.
struct GRPGeomMaterialCondition {
  static void declare(std::vector<o2::framework::InputSpec>& /*inputs*/) {}
  static void checkUpdates(o2::framework::ProcessingContext& /*pc*/) {}
  static bool finaliseCCDB(o2::framework::ConcreteDataMatcher& /*matcher*/, void* /*obj*/) { return false; }
  /// Null before the first timeframe.
  [[nodiscard]] static const o2::base::MatLayerCylSet* lut() { return o2::base::GRPGeomHelper::instance().getMatLUT(); }
};

/// What a `Conditions` without the material LUT holds in its place: no
/// input, nothing fetched.
struct NoMaterialCondition {
  static void declare(std::vector<o2::framework::InputSpec>& /*inputs*/) {}
  static void checkUpdates(o2::framework::ProcessingContext& /*pc*/) {}
  static bool finaliseCCDB(o2::framework::ConcreteDataMatcher& /*matcher*/, void* /*obj*/) { return false; }
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOWS_MATERIALLUT_H_
