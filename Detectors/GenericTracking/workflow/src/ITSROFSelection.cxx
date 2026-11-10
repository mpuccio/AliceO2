// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "GenericTracking/Workflow/ITSROFSelection.h"
#include "GenericTracking/Workflow/RecoParam.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DataFormatsITSMFT/PhysTrigger.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "Framework/DataRefUtils.h"
#include "Framework/InputRecord.h"
#include "Framework/InputRecordWalker.h"
#include "Framework/ProcessingContext.h"
#include "Framework/TimingInfo.h"
#include "ITStracking/FastMultEst.h"

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>

namespace o2::gtrk::workflow
{

namespace
{
using Estimator = o2::its::FastMultEst;
using Overlap = Estimator::ROFOverlapTableN;
using Mask = Estimator::ROFMaskTableN;
using Timing = o2::itsmft::tracking::ROFTimingLayer;
using Alpide = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>;
constexpr int kLayers = Estimator::NLayers;

int triggerType()
{
  const auto type = recoParam<detectors::ITS>().selectWithTriggers;
  if (type < 0 || type > 2) {
    throw std::invalid_argument("GTRKRecoITS.selectWithTriggers must be 0 (none), 1 (ITS physics) or 2 (TRD)");
  }
  return type;
}
} // namespace

struct ITSROFSelection::State {
  struct ClusterROF {
    int layer;
    int rof;
  };
  Overlap overlap;
  Mask mask;
  std::vector<ClusterROF> clusterROFs;
};

void ITSROFSelection::addInputs(std::vector<o2::framework::InputSpec>& inputs)
{
  using namespace o2::framework;
  inputs.emplace_back("its_selection_clusters", ConcreteDataTypeMatcher{"ITS", "COMPCLUSTERS"}, Lifetime::Timeframe);
  inputs.emplace_back("its_selection_rofs", ConcreteDataTypeMatcher{"ITS", "CLUSTERSROF"}, Lifetime::Timeframe);
  if (triggerType() == 1) {
    inputs.emplace_back("its_selection_triggers", "ITS", "PHYSTRIG", 0, Lifetime::Timeframe);
  } else if (triggerType() == 2) {
    inputs.emplace_back("its_selection_triggers", "TRD", "TRKTRGRD", 0, Lifetime::Timeframe);
  }
  if (o2::its::FastMultEstConfig::Instance().isVtxMultCutRequested() && !recoParam<detectors::ITS>().vertexer) {
    throw std::invalid_argument("fastMultConfig vertex cuts require GTRKRecoITS.vertexer=true (upstream VertexSeed inputs have no contributor counts)");
  }
}

ITSROFSelection::ITSROFSelection(o2::framework::ProcessingContext& pc, std::span<const FrameTiming> timing, int orbitsPerTF)
  : mState(std::make_unique<State>())
{
  using namespace o2::framework;
  std::array<gsl::span<const o2::itsmft::CompClusterExt>, kLayers> clusters{};
  std::array<gsl::span<const o2::itsmft::ROFRecord>, kLayers> rofs{};
  std::array<bool, kLayers> seen{};
  for (const auto& ref : InputRecordWalker{pc.inputs()}) {
    const auto* header = DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
    if (header->dataOrigin != o2::header::DataOrigin{"ITS"}) {
      continue;
    }
    const auto layer = header->subSpecification;
    if (header->dataDescription == o2::header::DataDescription{"COMPCLUSTERS"}) {
      if (layer >= kLayers) {
        throw std::runtime_error("ITS selection: invalid cluster stream");
      }
      clusters[layer] = pc.inputs().get<gsl::span<o2::itsmft::CompClusterExt>>(ref);
    } else if (header->dataDescription == o2::header::DataDescription{"CLUSTERSROF"}) {
      if (layer >= kLayers) {
        throw std::runtime_error("ITS selection: invalid ROF stream");
      }
      rofs[layer] = pc.inputs().get<gsl::span<o2::itsmft::ROFRecord>>(ref);
      seen[layer] = true;
    }
  }
  const auto streams = std::ranges::count(seen, true);
  if (streams != 1 && streams != kLayers) {
    throw std::runtime_error("ITS selection requires one monolithic stream or seven layer streams");
  }
  const bool staggered = streams == kLayers;
  if (timing.size() != kLayers || !seen[0]) {
    throw std::runtime_error("ITS selection: missing layer timing or stream zero");
  }
  auto& state = *mState;
  const auto& alpide = Alpide::Instance();
  const auto& config = recoParam<detectors::ITS>();
  for (int layer = 0; layer < kLayers; ++layer) {
    const int nROFs = std::max(orbitsPerTF * (o2::constants::lhc::LHCMaxBunches / timing[layer].lengthBC), static_cast<int>(rofs[staggered ? layer : 0].size()));
    state.overlap.defineLayer(layer, {static_cast<Timing::BCType>(nROFs), static_cast<Timing::BCType>(timing[layer].lengthBC), static_cast<Timing::BCType>(alpide.getROFDelayInBC(layer)), static_cast<Timing::BCType>(alpide.getROFBiasInBC(layer)), static_cast<Timing::BCType>(config.frameMarginBC.empty() ? 0 : config.frameMarginBC[layer])});
  }
  state.overlap.init();
  state.mask = Mask{state.overlap};
  static const o2::itsmft::ChipMappingITS mapping;
  for (int stream = 0; stream < streams; ++stream) {
    const auto offset = state.clusterROFs.size();
    state.clusterROFs.resize(offset + clusters[stream].size(), {-1, -1});
    std::size_t next = 0;
    for (std::size_t r = 0; r < rofs[stream].size(); ++r) {
      const auto& record = rofs[stream][r];
      if (record.getFirstEntry() < 0 || static_cast<std::size_t>(record.getFirstEntry()) != next || record.getNEntries() < 0 || static_cast<std::size_t>(record.getNEntries()) > clusters[stream].size() - next) {
        throw std::runtime_error("ITS selection: cluster ROFs do not cover the stream in order");
      }
      for (int i = 0; i < record.getNEntries(); ++i, ++next) {
        const int layer = staggered ? stream : mapping.getLayer(clusters[stream][next].getSensorID());
        state.clusterROFs[offset + next] = {layer, static_cast<int>(r)};
      }
    }
    if (next != clusters[stream].size()) {
      throw std::runtime_error("ITS selection: cluster ROFs do not cover the stream");
    }
  }
  const auto firstOrbit = pc.services().get<TimingInfo>().firstTForbit;
  const o2::InteractionRecord firstIR{0, firstOrbit};
  std::vector<o2::itsmft::PhysTrigger> triggers;
  if (triggerType() == 1) {
    const auto input = pc.inputs().get<gsl::span<o2::itsmft::PhysTrigger>>("its_selection_triggers");
    triggers.assign(input.begin(), input.end());
  } else if (triggerType() == 2) {
    auto lastIR = firstIR;
    const auto input = pc.inputs().get<gsl::span<o2::trd::TriggerRecord>>("its_selection_triggers");
    for (const auto& trigger : input) {
      if (trigger.getBCData() >= lastIR && trigger.getNumberOfTracklets()) {
        lastIR = trigger.getBCData();
        triggers.push_back({.ir = lastIR, .data = 0});
      }
    }
  }
  // The native estimator indexes the selection array with the trigger ROF;
  // discard triggers outside this input before handing it that array.
  const auto& mult = o2::its::FastMultEstConfig::Instance();
  const int layer = mult.isMultCutRequested() ? std::clamp(mult.cutMultClusLayer, 0, kLayers - 1) : state.overlap.getView().getClock();
  const auto count = rofs[staggered ? layer : 0].size();
  std::erase_if(triggers, [&](const auto& trigger) {
    const auto bc = trigger.ir.differenceInBC(firstIR);
    return bc < 0 || bc > std::numeric_limits<Timing::BCType>::max() || static_cast<std::size_t>(state.overlap.getLayer(layer).getROF(static_cast<Timing::BCType>(bc))) >= count;
  });
  Estimator estimator;
  const auto selected = estimator.selectROFs(rofs, clusters, triggers, firstOrbit, staggered, state.overlap.getView(), state.mask);
  LOGP(info, "ITS ROF selection: {} of {} ROFs on selection layer {}", selected, count, layer);
}

ITSROFSelection::~ITSROFSelection() = default;

bool ITSROFSelection::selectWithVertices(std::span<const VertexSeed> seeds, std::span<const std::uint32_t> contributors)
{
  if (!o2::its::FastMultEstConfig::Instance().isVtxMultCutRequested()) {
    return false;
  }
  if (seeds.size() != contributors.size()) {
    throw std::runtime_error("ITS vertex selection: contributor counts do not match the vertices");
  }
  std::vector<o2::its::Vertex> vertices;
  vertices.reserve(seeds.size());
  for (std::size_t v = 0; v < seeds.size(); ++v) {
    auto& vertex = vertices.emplace_back();
    vertex.setNContributors(static_cast<ushort>(std::min<std::uint32_t>(contributors[v], std::numeric_limits<ushort>::max())));
    const auto start = bcOfPS(seeds[v].timePS);
    const auto end = bcOfPS(static_cast<std::int64_t>(seeds[v].timePS) + seeds[v].widthPS);
    vertex.setTimeStamp({static_cast<std::uint32_t>(start), static_cast<std::uint32_t>(end - start)});
  }
  Estimator{}.selectROFsWithVertices(vertices, mState->overlap.getView(), mState->mask);
  return true;
}

std::vector<std::uint8_t> ITSROFSelection::clockMask() const
{
  const auto clock = mState->overlap.getView().getClock();
  std::vector<std::uint8_t> result(mState->overlap.getLayer(clock).mNROFsTF);
  const auto mask = mState->mask.getView();
  for (std::size_t r = 0; r < result.size(); ++r) {
    result[r] = mask.isROFEnabled(clock, static_cast<int>(r));
  }
  return result;
}

std::optional<ITSROFSelection::Measurements> ITSROFSelection::filter(MeasurementView source, std::span<const std::uint32_t> offsets, std::span<const std::uint32_t> originalRows) const
{
  const auto mask = mState->mask.getView();
  const auto accepted = [&](std::uint32_t row) {
    const auto cluster = source.clusterIndex[row];
    if (cluster >= mState->clusterROFs.size()) {
      throw std::runtime_error("ITS selection: measurement cluster index out of range");
    }
    const auto [layer, rof] = mState->clusterROFs[cluster];
    return mask.isROFEnabled(layer, rof);
  };
  std::uint32_t firstRejected = 0U;
  while (firstRejected < offsets.back() && accepted(firstRejected)) {
    ++firstRejected;
  }
  if (firstRejected == offsets.back()) {
    return std::nullopt;
  }
  std::vector<std::uint32_t> rows;
  for (std::uint32_t r = 0; r < firstRejected; ++r) {
    rows.push_back(r);
  }
  for (std::uint32_t r = firstRejected + 1U; r < offsets.back(); ++r) {
    if (accepted(r)) {
      rows.push_back(r);
    }
  }
  Measurements selected;
  resizeMeasurementColumns(selected.columns, rows.size());
  selected.originalRows.reserve(rows.size());
  for (std::size_t i = 0; i < rows.size(); ++i) {
    Measurement measurement;
    measurement.apply(source[rows[i]], [](auto& to, const auto& from) -> auto& { to = from; return to; });
    writeMeasurement(selected.columns, i, measurement);
    selected.originalRows.push_back(originalRows.empty() ? rows[i] : originalRows[rows[i]]);
  }
  selected.offsets.reserve(offsets.size());
  auto at = rows.begin();
  for (const auto offset : offsets) {
    at = std::lower_bound(at, rows.end(), offset);
    selected.offsets.push_back(static_cast<std::uint32_t>(at - rows.begin()));
  }
  return selected;
}

} // namespace o2::gtrk::workflow
