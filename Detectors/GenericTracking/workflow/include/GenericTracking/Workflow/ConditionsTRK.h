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

/// \file ConditionsTRK.h
/// \brief ALICE 3's run conditions for a GTRK device (Conditions.h), and its
///        topology dictionaries as CCDB conditions.
///
/// ALICE 3's conditions are CCDB objects, as the ITS's: the aligned geometry
/// (`GLO/Config/GeometryAligned`, through GRPGeomHelper), the readout frames'
/// `TRKAlmiraParam` (`TRK/Config/AlmiraParam`, its singleton updated by the
/// CCDB fetch, as `ITS/Config/AlpideParam`'s) and the topology dictionaries
/// (`TRK/Calib/ClusterDictionary`, `FT3/Calib/ClusterDictionary`). The field
/// is O2's ALICE 3 field (`alice3BzKiloGauss`): O2 builds it from the macro
/// `ALICE3_MAGFIELD_MACRO` names, not from a GRP; the beam line is the nominal
/// one (the simulation's diamond). Offline, from CCDB snapshots through DPL's
/// `--condition-remap`.

#ifndef GTRK_WORKFLOWS_CONDITIONSTRK_H_
#define GTRK_WORKFLOWS_CONDITIONSTRK_H_

#include "GenericTracking/Alice/TRK.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Workflow/Conditions.h"

#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsRaw/HBFUtils.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InitContext.h"
#include "Framework/InputRecord.h"
#include "Framework/InputSpec.h"
#include "Framework/ProcessingContext.h"
#include "Headers/DataHeader.h"
#include "TRKBase/AlmiraParam.h"
#include "TRKReconstruction/TopologyDictionary.h"

#include <fairlogger/Logger.h>

#include <format>
#include <memory>
#include <stdexcept>
#include <vector>

namespace o2::gtrk::workflow
{

/// The TRK's and FT3's topology dictionaries, condition inputs of the device
/// that holds them.
class TrkDictionaryConditions
{
 public:
  explicit TrkDictionaryConditions(std::vector<o2::framework::InputSpec>& inputs)
  {
    inputs.emplace_back("trkDictionary", "TRK", "CLUSDICT", 0, o2::framework::Lifetime::Condition, o2::framework::ccdbParamSpec("TRK/Calib/ClusterDictionary"));
    inputs.emplace_back("ft3Dictionary", "FT3", "CLUSDICT", 0, o2::framework::Lifetime::Condition, o2::framework::ccdbParamSpec("FT3/Calib/ClusterDictionary"));
  }

  /// At the top of the device's run: delivers a new one through finaliseCCDB.
  static void checkUpdates(o2::framework::ProcessingContext& pc)
  {
    pc.inputs().get<alice::TrkDictionary*>("trkDictionary");
    pc.inputs().get<alice::Ft3Dictionary*>("ft3Dictionary");
  }

  /// Whether \a matcher is one of them, which is then taken. Throws
  /// std::runtime_error if its layers are not the TRKFT3 layout's barrel
  /// layers or discs.
  bool finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj)
  {
    constexpr int kBarrel = Detector<detectors::TRKFT3>::nBarrelLayers;
    if (matcher == o2::framework::ConcreteDataMatcher("TRK", "CLUSDICT", 0)) {
      mTrk = checked(static_cast<const alice::TrkDictionary*>(obj), "TRK", kBarrel);
      return true;
    }
    if (matcher == o2::framework::ConcreteDataMatcher("FT3", "CLUSDICT", 0)) {
      mFt3 = checked(static_cast<const alice::Ft3Dictionary*>(obj), "FT3", Detector<detectors::TRKFT3>::nLayers - kBarrel);
      return true;
    }
    return false;
  }

  /// The CCDB cache's; null before the first timeframe.
  [[nodiscard]] const alice::TrkDictionary* trk() const { return mTrk; }
  [[nodiscard]] const alice::Ft3Dictionary* ft3() const { return mFt3; }

 private:
  template <typename Dictionary>
  static const Dictionary* checked(const Dictionary* dictionary, const char* name, int layers)
  {
    if (dictionary->getNumberOfLayers() != layers) {
      throw std::runtime_error(std::format("{}'s topology dictionary has {} layers, the TRKFT3 layout {}", name, dictionary->getNumberOfLayers(), layers));
    }
    LOGP(info, "{}'s topology dictionary from the CCDB, its offsets {}", name, dictionary->getApplyOffsets() ? "applied" : "not applied");
    return dictionary;
  }

  const alice::TrkDictionary* mTrk = nullptr;
  const alice::Ft3Dictionary* mFt3 = nullptr;
};

/// ALICE 3's (see above).
template <bool WithMaterial>
class Conditions<detectors::TRKFT3, WithMaterial>
{
 public:
  explicit Conditions(std::vector<o2::framework::InputSpec>& inputs)
    : mRequest(std::make_shared<o2::base::GRPGeomRequest>(false, false, false, false, false, o2::base::GRPGeomRequest::Aligned, inputs)),
      mDictionaries(inputs)
  {
    inputs.emplace_back("almiraParam", origin(), "ALMIRAPARAM", 0, o2::framework::Lifetime::Condition, o2::framework::ccdbParamSpec("TRK/Config/AlmiraParam"));
    MaterialConditionOf<detectors::TRKFT3, WithMaterial>::declare(inputs);
  }

  /// The orbits per timeframe of the ALICE 3 simulation's HBFUtils: it
  /// writes no GRPECS.
  [[nodiscard]] static int orbitsPerTimeframe() { return o2::raw::HBFUtils::Instance().getNOrbitsPerTF(); }

  /// None: the beam line is the nominal one.
  static constexpr bool kMeanVertex = false;

  [[nodiscard]] static o2::header::DataOrigin origin() { return detectorOrigin<detectors::TRKFT3>(); }

  /// The field as the propagator's nominal Bz, which every single-value use
  /// reads (GridConfig.h's `nominalBzKiloGauss`).
  void init(o2::framework::InitContext& /*ic*/)
  {
    o2::base::GRPGeomHelper::instance().setRequest(mRequest);
    o2::base::Propagator::Instance(true)->setNominalBz(static_cast<float>(alice::alice3BzKiloGauss()));
    LOGP(info, "TRK: the ALICE 3 field's Bz {} kG", o2::base::Propagator::Instance()->getNominalBz());
  }

  void checkUpdates(o2::framework::ProcessingContext& pc) const
  {
    o2::base::GRPGeomHelper::instance().checkUpdates(pc);
    pc.inputs().get<o2::trk::AlmiraParam*>("almiraParam"); // the CCDB fetcher updates the singleton
    TrkDictionaryConditions::checkUpdates(pc);
    MaterialConditionOf<detectors::TRKFT3, WithMaterial>::checkUpdates(pc);
  }

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj)
  {
    if (mMaterial.finaliseCCDB(matcher, obj) || o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
      return;
    }
    if (matcher == o2::framework::ConcreteDataMatcher(origin(), "ALMIRAPARAM", 0)) {
      LOGP(info, "TRK's AlmiraParam from the CCDB");
      return;
    }
    mDictionaries.finaliseCCDB(matcher, obj);
  }

  /// The readout frames of the run's `TRKAlmiraParam` (TRK.h).
  [[nodiscard]] std::vector<FrameTiming> frameTiming() const { return alice::trkFrameTiming(); }

  /// The topology dictionaries' (TRK.h).
  [[nodiscard]] std::vector<double> measurementDepthsCM() const { return alice::measurementDepthsCM(*mDictionaries.trk(), *mDictionaries.ft3()); }

  /// The run's material LUT, the CCDB cache's; null before the first timeframe.
  [[nodiscard]] const o2::base::MatLayerCylSet* materialLut() const
    requires WithMaterial
  {
    return mMaterial.lut();
  }

 private:
  std::shared_ptr<o2::base::GRPGeomRequest> mRequest;
  TrkDictionaryConditions mDictionaries;
  [[no_unique_address]] MaterialConditionOf<detectors::TRKFT3, WithMaterial> mMaterial;
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOWS_CONDITIONSTRK_H_
