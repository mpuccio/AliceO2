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

/// \file Conditions.h
/// \brief The run's conditions for a GTRK device, per detector:
///        `Conditions<Det, WithMaterial>`, and the ITS's and MFT's
///        (ALICE 3's in ConditionsTRK.h).
///
/// Every detector's give a device, from the first timeframe on, what it
/// places its layout from (`alice::placeChips`, which reads the detector's
/// GeometryTGeo), the field as the propagator's nominal Bz (GridConfig.h's
/// `nominalBzKiloGauss`), its readout frames per layer and its layers'
/// measurement depths; they do not change within a run (a device sees one).
/// The MeanVertex does, and is the measurement producer's own input, per
/// timeframe, where the detector has one (`kMeanVertex`):
///
/// - `Conditions(inputs)` declares its condition inputs at construction;
/// - `init(ic)` in the device's init, `checkUpdates(pc)` at the top of its
///   run, `finaliseCCDB(matcher, obj)` forwarded from the device's;
/// - `frameTiming()`; `measurementDepthsCM()`, empty for 0 on every layer;
///   `origin()`, the detector's data origin; `kMeanVertex`, whether the
///   detector's beam line is the run's MeanVertex (`GLO/Calib/MeanVertex`);
/// - with `WithMaterial`, the material LUT (MaterialLUT.h) as one more
///   condition, `materialLut()`: the devices that track take it, the
///   measurement producer and the MC vertex reader (the default) fetch none.
///
/// The ITS's and MFT's are what O2's own trackers take by default: not the
/// full TGeo but the detector's lightweight `GeometryTGeo` object
/// (`ITS/Config/Geometry`, `MFT/Config/Geometry`), with its chip matrices
/// filled, which the placement reads like a loaded geometry; the field
/// (`GLO/Config/GRPMagField`) through GRPGeomHelper; and `<DET>/Config/AlpideParam`
/// for the frames. Offline, from CCDB snapshots through DPL's
/// `--condition-remap`.

#ifndef GTRK_WORKFLOWS_CONDITIONS_H_
#define GTRK_WORKFLOWS_CONDITIONS_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/MaterialLUT.h"

#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InitContext.h"
#include "Framework/InputRecord.h"
#include "Framework/InputSpec.h"
#include "Framework/ProcessingContext.h"
#include "Headers/DataHeader.h"
#include "ITSBase/GeometryTGeo.h"
#include "MFTBase/GeometryTGeo.h"

#include <fairlogger/Logger.h>

#include <concepts>
#include <format>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace o2::gtrk::workflow
{

/// The conditions a device of \a Det reads (see above); with \a WithMaterial
/// also the material LUT.
template <DetectorLayout Det, bool WithMaterial = false>
class Conditions;

/// The material LUT's condition a `Conditions` holds: the LUT's, or none.
template <DetectorLayout Det, bool WithMaterial>
using MaterialConditionOf = std::conditional_t<WithMaterial, MaterialCondition<Det>, NoMaterialCondition>;

/// The ITS's and the MFT's.
template <DetectorLayout Det, bool WithMaterial>
  requires std::same_as<Det, detectors::ITS> || std::same_as<Det, detectors::MFT>
class Conditions<Det, WithMaterial>
{
  static constexpr bool kITS = std::same_as<Det, detectors::ITS>;

 public:
  explicit Conditions(std::vector<o2::framework::InputSpec>& inputs)
    : mRequest(std::make_shared<o2::base::GRPGeomRequest>(false, true, false, true, false, o2::base::GRPGeomRequest::None, inputs))
  {
    const auto detector = origin().template as<std::string>();
    inputs.emplace_back("geomTGeo", origin(), "GEOMTGEO", 0, o2::framework::Lifetime::Condition, o2::framework::ccdbParamSpec(detector + "/Config/Geometry"));
    inputs.emplace_back("alpideParam", origin(), "ALPIDEPARAM", 0, o2::framework::Lifetime::Condition, o2::framework::ccdbParamSpec(detector + "/Config/AlpideParam"));
    MaterialConditionOf<Det, WithMaterial>::declare(inputs);
  }

  static constexpr bool kMeanVertex = true;

  [[nodiscard]] static o2::header::DataOrigin origin() { return detectorOrigin<Det>(); }

  void init(o2::framework::InitContext& /*ic*/) { o2::base::GRPGeomHelper::instance().setRequest(mRequest); }

  void checkUpdates(o2::framework::ProcessingContext& pc) const
  {
    o2::base::GRPGeomHelper::instance().checkUpdates(pc);
    pc.inputs().get<GeometryTGeoType*>("geomTGeo");
    pc.inputs().get<AlpideParam*>("alpideParam"); // the CCDB fetcher updates the singleton
    MaterialConditionOf<Det, WithMaterial>::checkUpdates(pc);
  }

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj)
  {
    // The material first: the ITS's LUT is GLO/MATLUT, GRPGeomHelper's description.
    if (mMaterial.finaliseCCDB(matcher, obj) || o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
      return;
    }
    if (matcher == o2::framework::ConcreteDataMatcher(origin(), "GEOMTGEO", 0)) {
      GeometryTGeoType::adopt(static_cast<GeometryTGeoType*>(obj));
      LOGP(info, "{} GeometryTGeo from the CCDB adopted", Det::name);
    }
  }

  /// The readout frames of the run's `DPLAlpideParam`, as O2's ITS tracker
  /// takes them (`ITSTrackingInterface::updateTimeDependentParams`): frame `k`
  /// at `k * length + delay + bias` BC, one per layer, every layer the same.
  /// Throws std::runtime_error for staggered frames: the converter reads one
  /// ROF stream.
  [[nodiscard]] std::vector<FrameTiming> frameTiming() const
  {
    const auto& param = AlpideParam::Instance();
    const FrameTiming result{.lengthBC = param.getROFLengthInBC(0), .firstStartBC = param.getROFDelayInBC(0) + param.getROFBiasInBC(0)};
    for (int layer = 1; layer < AlpideParam::getNLayers(); ++layer) {
      if (param.getROFLengthInBC(layer) != result.lengthBC || param.getROFDelayInBC(layer) + param.getROFBiasInBC(layer) != result.firstStartBC) {
        throw std::runtime_error(std::format("{}'s readout frames are staggered (layer {}): not read yet", Det::name, layer));
      }
    }
    if (result.lengthBC <= 0) {
      throw std::runtime_error(std::format("{}'s AlpideParam gives readout frames of {} BC", Det::name, result.lengthBC));
    }
    return std::vector<FrameTiming>(static_cast<std::size_t>(Detector<Det>::nLayers), result);
  }

  /// The run's orbits per timeframe, its GRPECS's, as O2's ITS tracking
  /// takes them.
  [[nodiscard]] static int orbitsPerTimeframe() { return o2::base::GRPGeomHelper::getNHBFPerTF(); }

  /// None: the ALPIDE's clusters are at the effective sensitive layer's
  /// centre, 0 on every layer.
  [[nodiscard]] std::vector<double> measurementDepthsCM() const { return {}; }

  /// The run's material LUT, the CCDB cache's; null before the first timeframe.
  [[nodiscard]] const o2::base::MatLayerCylSet* materialLut() const
    requires WithMaterial
  {
    return mMaterial.lut();
  }

 private:
  using GeometryTGeoType = std::conditional_t<kITS, o2::its::GeometryTGeo, o2::mft::GeometryTGeo>;
  using AlpideParam = o2::itsmft::DPLAlpideParam<kITS ? o2::detectors::DetID::ITS : o2::detectors::DetID::MFT>;

  std::shared_ptr<o2::base::GRPGeomRequest> mRequest;
  [[no_unique_address]] MaterialConditionOf<Det, WithMaterial> mMaterial;
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOWS_CONDITIONS_H_
