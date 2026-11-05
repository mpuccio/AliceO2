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

#include "GenericTracking/Study/MCVertexReaderSpec.h"
#include "GenericTracking/Study/MCVertexParam.h"

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Workflow/Conditions.h"
#include "GenericTracking/Workflow/ConditionsTRK.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/RecoParam.h"

#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/StringUtils.h"
#include "DetectorsRaw/HBFUtils.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Logger.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "SimulationDataFormat/DigitizationContext.h"
#include "Steer/MCKinematicsReader.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <ranges>
#include <string>
#include <vector>

using namespace o2::framework;

namespace o2::gtrk::study
{

namespace
{

/// splitmix64's finaliser: a mix, not a generator.  Cheap, stateless, and it
/// decorrelates the small consecutive integers an identity is made of.
[[nodiscard]] constexpr std::uint64_t mix64(std::uint64_t x)
{
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

/// Fold a run seed and a collision's index into one seed.
[[nodiscard]] constexpr std::uint64_t seedFrom(std::uint64_t seed, std::size_t collision) { return mix64(mix64(seed) ^ static_cast<std::uint64_t>(collision)); }

/// Reads truth only at this workflow boundary and publishes the VertexSeed
/// payload the own vertexer's list has too (GTRKReco<DET>.vertexer).
template <DetectorLayout Det>
class MCVertexReaderTask : public Task
{
 public:
  void init(InitContext& ic) final
  {
    const auto& param = MCVertexParam::Instance();
    if (!std::isfinite(param.sigmaZ) || !std::isfinite(param.sigmaTransverse) || param.sigmaZ < 0. || param.sigmaTransverse < 0.) {
      LOG(fatal) << "GTRKMCVertex: the smearing sigmas must be finite and non-negative";
      return;
    }
    const auto filename = o2::utils::Str::concat_string(o2::utils::Str::rectifyDirectory(ic.options().get<std::string>("input-dir")),
                                                        ic.options().get<std::string>("collision-context"));
    mContext.reset(o2::steer::DigitizationContext::loadFromFile(filename));
    if (!mContext) {
      LOGP(fatal, "cannot load collision context from {}", filename);
      return;
    }
    if (mContext->getEventRecords().size() != mContext->getEventParts().size()) {
      LOGP(fatal, "collision context has {} records but {} event-part lists", mContext->getEventRecords().size(), mContext->getEventParts().size());
      return;
    }
    if (!std::ranges::is_sorted(mContext->getEventRecords(), {}, [](const auto& record) { return record.getTimeNS(); })) {
      LOG(fatal) << "collision times are not sorted";
      return;
    }

    mKinematics = std::make_unique<o2::steer::MCKinematicsReader>(mContext.get());
    if (!mKinematics->isInitialized()) {
      LOGP(fatal, "cannot open MC kinematics referenced by {}", filename);
      return;
    }
  }

  void run(ProcessingContext& pc) final
  {
    const auto& param = MCVertexParam::Instance();
    const auto& timing = pc.services().get<TimingInfo>();
    const double startNS = static_cast<double>(timing.firstTForbit) * o2::constants::lhc::LHCOrbitNS;
    // The simulation's timeframes, which cut the collision context.
    const double endNS = startNS + (1.e-3 * static_cast<double>(workflow::timeframeDurationPS(o2::raw::HBFUtils::Instance().getNOrbitsPerTF())));
    const auto& records = mContext->getEventRecords();
    const auto& parts = mContext->getEventParts();
    const auto first = std::ranges::lower_bound(records, startNS, {}, [](const auto& record) { return record.getTimeNS(); });
    const auto last = std::ranges::lower_bound(first, records.end(), endNS, {}, [](const auto& record) { return record.getTimeNS(); });

    std::vector<VertexSeed> vertices;
    vertices.reserve(static_cast<std::size_t>(last - first));
    const auto varZ = static_cast<Scalar>(param.sigmaZ * param.sigmaZ);
    const auto varTransverse = static_cast<Scalar>(param.sigmaTransverse * param.sigmaTransverse);

    for (auto record = first; record != last; ++record) {
      const auto collision = static_cast<std::size_t>(record - records.begin());
      if (parts[collision].empty()) {
        LOGP(fatal, "collision {} has no MC event part", collision);
        return;
      }
      const auto& part = parts[collision].front();
      const auto& header = mKinematics->getMCEventHeader(part.sourceID, part.entryID);

      // Per-collision seeding keeps the result independent of scheduling,
      // skipped timeframes and the order in which devices are run.
      std::mt19937_64 generator(seedFrom(param.randomSeed, collision));
      std::normal_distribution<double> normal;
      const double zSmear = normal(generator);

      // The collision's time, a box of width 0.
      constexpr double kPSPerNS = 1.e3;
      const auto timePS = std::max<std::int64_t>(std::llround((record->getTimeNS() - startNS) * kPSPerNS), 0);
      const VertexSeed vertex{.timePS = static_cast<std::uint32_t>(timePS),
                              .widthPS = 0U,
                              .z = static_cast<Scalar>(header.GetZ() + param.sigmaZ * zSmear),
                              .varZ = varZ,
                              .varTransverse = varTransverse};
      vertices.push_back(vertex);
    }
    std::ranges::sort(vertices, {}, &VertexSeed::timePS);

    pc.outputs().snapshot(Output{workflow::gDataOrigin, workflow::description<Det>(workflow::kVertices), 0}, vertices);
    LOGP(info, "published {} MC vertex seeds for TF {} starting at orbit {}", vertices.size(), timing.tfCounter, timing.firstTForbit);
  }

 private:
  std::unique_ptr<o2::steer::DigitizationContext> mContext;
  std::unique_ptr<o2::steer::MCKinematicsReader> mKinematics;
};

} // namespace

template <DetectorLayout Det>
DataProcessorSpec getMCVertexReaderSpec()
{
  // The clusters only trigger each timeframe; their content is not read.
  return DataProcessorSpec{.name = workflow::deviceName<Det>("mc-vertex-reader"),
                           .inputs = Inputs{InputSpec{"clusters", workflow::gDataOrigin, workflow::description<Det>(workflow::kClusters), 0, Lifetime::Timeframe}},
                           .outputs = Outputs{OutputSpec{workflow::gDataOrigin, workflow::description<Det>(workflow::kVertices), 0, Lifetime::Timeframe}},
                           .algorithm = AlgorithmSpec{adaptFromTask<MCVertexReaderTask<Det>>()},
                           .options = Options{{"collision-context", VariantType::String, "collisioncontext.root", {"collision context and MC-event mapping"}},
                                              {"input-dir", VariantType::String, "none", {"input directory"}}}};
}

template DataProcessorSpec getMCVertexReaderSpec<detectors::ITS>();
template DataProcessorSpec getMCVertexReaderSpec<detectors::MFT>();
template DataProcessorSpec getMCVertexReaderSpec<detectors::TRKFT3>();

} // namespace o2::gtrk::study
