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

/// \file testTrackletStage.cxx
/// \brief Tracklet generation against a whole-layer brute-force reference.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TrackletStage
#include <boost/test/unit_test.hpp>
#include <numbers>

#include "GenericTracking/Base/LayerLinks.h"
#include "ToyBarrel3.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackletStage.h"

#include "GPUCommonMath.h"
#include "GridTestSupport.h"
#include "StageOutputComparison.h"
#include "SerialReferences.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace o2::gtrk
{

namespace
{

/// Three barrel layers with an extra 0 -> 2 link. The table order makes the
/// outgoing links of layer 0 non-contiguous and therefore observable. Some
/// trajectory on 0 -> 2 passes beyond layer 1's ends, so it skips nothing by
/// the table: whether a pair on it leaves a hole is the pair's.
///
/// Layer 2 is five modules long, |z| 5 cm, so the fixture's targets lie on
/// it: the hole test's per-source check reads the target's extent.
struct StageToy {
  static constexpr const char* name = "StageToy";
  static constexpr std::array<BarrelLayerSpec, 3> barrelSpecs{{
    detectors::ToyBarrel3::barrelSpecs[0],
    detectors::ToyBarrel3::barrelSpecs[1],
    {.r = 8., .nStaves = 4, .nModulesPerStave = 5, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = detectors::ToyBarrel3::sensor},
  }};
  static constexpr auto discSpecs = detectors::ToyBarrel3::discSpecs;
  static constexpr auto barrelLayers = placedLayers(barrelSpecs);
  static constexpr auto discLayers = placedLayers(discSpecs);
  static constexpr std::array<LayerLink, 3> layerLinks{{
    {.from = 0, .to = 2, .sourceMin = -10., .sourceMax = 10., .skipsLayer = false},
    {.from = 1, .to = 2, .sourceMin = -10., .sourceMax = 10., .skipsLayer = false},
    {.from = 0, .to = 1, .sourceMin = -2., .sourceMax = 2., .skipsLayer = false},
  }};
  static constexpr int nSensors = countSensors(barrelLayers, discLayers);
};

[[maybe_unused]] const bool kStageToyPublished = publishToy<StageToy>();

using Frame = TimeFrame<StageToy>;
using Sources = LayerRowMap<StageToy>;
using Output = StageOutput<TrackletT>;
using Cpu = CpuStage<TrackletT>;
/// The tracklet stage alone: chains of two clusters, so the iteration takes
/// every link (`takenLinks`).
constexpr int kChainClusters = 2;
constexpr std::size_t kReferenceMeasurements = 9;
constexpr std::size_t kExtraSources = 64;
constexpr std::size_t kMeasurements = kReferenceMeasurements + kExtraSources;
/// Every row at 110 ns, in its toy tick's box [92, 133) ns.
constexpr double kRowTimeNS = 110.;
const TimeEstPS kRowBox = test::toyTimeBox(kRowTimeNS);

struct InputPoint {
  LayerID layer = kInvalidLayerID;
  Scalar r = 0.F;
  Scalar phi = 0.F;
  Scalar z = 0.F;
};

SearchWindowParams windowParameters()
{
  SearchWindowParams result;
  result.minPtGeV = 1.;
  result.bzKiloGauss = 20.;
  result.materialX2X0 = 0.;
  result.positionSigma = 0.01;
  result.nSigma = 3.;
  return result;
}

GridParams gridParameters()
{
  GridParams result;
  result.blockWidth = 20'000;
  return test::withToyFrames(result);
}

class Fixture
{
 public:
  static constexpr VertexConstraint kBeamConstraint{.z = 0., .sigmaZ = 0.5, .sigmaTransverse = 0.001};

  Fixture()
    : geometry(test::cappedGrid<StageToy>(20., gridParameters(), windowParameters(), kBeamConstraint, 8, 4)),
      offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U)
  {
    constexpr std::array<InputPoint, kReferenceMeasurements> reference{{
      {.layer = 0, .r = 2.F, .phi = 0.F, .z = 1.F},
      {.layer = 0, .r = 2.F, .phi = 0.4F, .z = -0.5F},
      {.layer = 0, .r = 2.F, .phi = 1.F, .z = 20.F},
      {.layer = 1, .r = 4.F, .phi = 0.003F, .z = 2.F},
      {.layer = 1, .r = 4.F, .phi = (2 * std::numbers::pi_v<Scalar>)-0.003F, .z = 2.F},
      {.layer = 1, .r = 4.F, .phi = 0.F, .z = 0.2F},
      {.layer = 2, .r = 8.F, .phi = 0.005F, .z = 4.F},
      {.layer = 2, .r = 8.F, .phi = (2 * std::numbers::pi_v<Scalar>)-0.005F, .z = 4.F},
      {.layer = 2, .r = 8.F, .phi = 0.F, .z = 0.F},
    }};
    std::array<InputPoint, kMeasurements> input{};
    std::ranges::copy(reference, input.begin());
    for (std::size_t i = 0; i < kExtraSources; ++i) {
      input[kReferenceMeasurements + i] = {
        .layer = 0, .r = 2.F, .phi = (2 * std::numbers::pi_v<Scalar>)*static_cast<Scalar>(i) / static_cast<Scalar>(kExtraSources), .z = 20.F};
    }

    frame = test::frameOf(geometry, input.size(), [&](std::size_t i) { return test::FramePoint{.layer = input[i].layer, .x = input[i].r * std::cos(input[i].phi), .y = input[i].r * std::sin(input[i].phi), .z = input[i].z, .r = input[i].r, .phi = input[i].phi}; }, [&](std::size_t i) { return kRowBox; }, [&](std::size_t i) { return test::PointVariances{1e-4F, 1e-4F, 1e-4F}; }, columns, offsets);
    used.assign(frame->size(), kFree);
    std::array<MeasurementRange, Sources::kLayerCount> rows{};
    for (std::size_t layer = 0; layer < rows.size(); ++layer) {
      rows[layer] = frame->grid().layerRows(static_cast<LayerID>(layer));
    }
    sources = Sources{rows, {true, true, false}};
  }

  GridGeometry<StageToy> geometry;
  std::vector<std::uint32_t> offsets;
  MeasurementStore columns;
  std::unique_ptr<Frame> frame;
  std::vector<std::uint8_t> used; ///< the used mask: every row free
  Sources sources;
  VertexSeed beamLine = vertexSeed(kBeamConstraint);
  std::array<VertexSeed, 4> vertices{{
    {.timePS = 0U, .z = 0.F},
    {.timePS = 109'500U, .z = 2.F},
    {.timePS = 110'500U, .z = 0.F},
    {.timePS = 200'000U, .z = 0.F},
  }};
  std::array<std::uint32_t, Sources::kLayerCount> maxCandidateWidthPS{kRowBox.width, kRowBox.width, kRowBox.width};
};

[[nodiscard]] bool vertexUsed(const IterationParameters<StageToy>& parameters, const auto& source, const VertexSeed& vertex)
{
  return vertexUsable(parameters, timeBox(source), vertex);
}

template <typename Stage>
std::unique_ptr<Output> runStage(const Fixture& fixture, const Stage& stage)
{
  auto output = std::make_unique<Output>();
  output->setMemoryPool(std::make_shared<BoundedMemoryResource>());
  output->prepare(fixture.sources.size());
  BOOST_REQUIRE(runSerial(*output, stage) == StageStatus::kComplete);
  return output;
}

std::vector<std::uint32_t> bruteForceTargets(const Fixture& fixture, const IterationParameters<StageToy>& parameters, std::span<const VertexSeed> vertices,
                                             Sources::Index sourceElement, std::size_t linkIndex)
{
  const auto sourceAt = fixture.sources.row(sourceElement);
  const auto source = fixture.frame->measurement(sourceAt.row);
  const auto& link = StageToy::layerLinks[linkIndex];
  std::vector<std::uint32_t> result;
  const Scalar sourceSecond = secondCoordinate<StageToy>(sourceAt.layer, source.r, source.z);
  if (link.from != sourceAt.layer || sourceSecond < parameters.sourceMin[linkIndex] || sourceSecond > parameters.sourceMax[linkIndex]) {
    return result;
  }
  const auto sourcePoint = pairPoint(source);
  const auto targets = fixture.frame->grid().layerRows(link.to);
  for (auto targetRow = targets.begin; targetRow < targets.end; ++targetRow) {
    const auto target = fixture.frame->measurement(targetRow);
    if (!timeBox(source).isCompatible(timeBox(target))) {
      continue;
    }
    // Without skips, no pair whose segment crosses layer 1 (LinkHoles.h).
    if (segmentLeavesHole(parameters.holeLayers[linkIndex], parameters.barrelLayers, parameters.layerSpans, source.r, source.z, target.r, target.z)) {
      continue;
    }
    for (const auto& vertex : vertices) {
      if (!vertexUsed(parameters, source, vertex)) {
        continue;
      }
      const Scalar scattering = sourceScattering<StageToy>(parameters.scattering, parameters.pairCut, vertex, sourceAt.layer, sourcePoint);
      if (pairCompatible(parameters.pairCut, vertex, sourcePoint, scattering, pairPoint(target))) {
        result.push_back(targetRow);
        break;
      }
    }
  }
  return result;
}

void verify(VertexSource vertexSource, std::span<const VertexSeed> suppliedVertices = {}, bool checkSpecialCases = true, bool skipLinks = true)
{
  Fixture fixture;
  const IterationSettings settings{.windows = windowParameters(), .vertexSource = vertexSource, .minClusters = kChainClusters, .skipLinks = skipLinks};
  const std::array<double, 3> material{};
  const auto parameters = iterationParameters<StageToy>(settings, Fixture::kBeamConstraint, scatteringTable<StageToy>(material));
  const std::span<const VertexSeed> vertexList = suppliedVertices.empty() ? std::span<const VertexSeed>{fixture.vertices} : suppliedVertices;
  const std::span<const VertexSeed> selected = vertexSource == VertexSource::kBeamLine ? std::span<const VertexSeed>{&fixture.beamLine, 1U} : vertexList;
  const auto vertexWindows = vertexWindowParameters<StageToy>(parameters, settings.windows, vertexList);

  const TrackletStage<StageToy> stage{
    parameters, fixture.sources, *fixture.frame, fixture.geometry.allBlocks(), fixture.beamLine, vertexList, fixture.maxCandidateWidthPS,
    vertexWindows, fixture.used.data()};
  Output output;
  output.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  output.prepare(fixture.sources.size());
  BOOST_REQUIRE((runSerial(output, stage) == StageStatus::kComplete));

  const auto view = std::as_const(output).view();
  std::vector<Tracklet> expectedOrder;
  for (Sources::Index sourceElement = 0; sourceElement < fixture.sources.size(); ++sourceElement) {
    const auto sourceAt = fixture.sources.row(sourceElement);
    const auto source = fixture.frame->measurement(sourceAt.row);
    std::size_t outputAt = output.offsets()[sourceElement];

    for (auto at = outLinks<StageToy>.begin(sourceAt.layer); at < outLinks<StageToy>.end(sourceAt.layer); ++at) {
      const std::size_t linkIndex = outLinks<StageToy>.links[at];
      const auto brute = bruteForceTargets(fixture, parameters, selected, sourceElement, linkIndex);
      std::vector<std::uint32_t> emitted;
      while (outputAt < output.offsets()[sourceElement + 1U] && view[outputAt].link == linkIndex) {
        emitted.push_back(view[outputAt].target);
        ++outputAt;
      }
      auto sortedEmitted = emitted;
      std::ranges::sort(sortedEmitted);
      auto sortedBrute = brute;
      std::ranges::sort(sortedBrute);
      BOOST_REQUIRE(sortedEmitted.size() == sortedBrute.size());
      for (std::size_t i = 0; i < sortedBrute.size(); ++i) {
        BOOST_TEST(sortedEmitted[i] == sortedBrute[i]);
      }

      const Scalar sourceSecond = secondCoordinate<StageToy>(sourceAt.layer, source.r, source.z);
      if (sourceSecond < parameters.sourceMin[linkIndex] || sourceSecond > parameters.sourceMax[linkIndex]) {
        continue;
      }
      const auto& link = StageToy::layerLinks[linkIndex];
      const auto window = candidateWindow(parameters.windows[linkIndex], timeBox(source), source.phi, source.r, source.z, fixture.maxCandidateWidthPS[link.to]);
      fixture.frame->grid().visitCandidates(link.to, window, fixture.geometry.allBlocks(), [&](std::uint32_t targetRow) {
        if (std::ranges::find(brute, targetRow) != brute.end()) {
          const auto target = fixture.frame->measurement(targetRow);
          const Scalar dx = target.x - source.x;
          const Scalar dy = target.y - source.y;
          const Scalar dz = target.z - source.z;
          expectedOrder.push_back({.source = sourceAt.row,
                                   .target = targetRow,
                                   .link = static_cast<std::uint16_t>(linkIndex),
                                   .lambda = o2::gpu::GPUCommonMath::ATan2(dz, o2::gpu::GPUCommonMath::Sqrt((dx * dx) + (dy * dy))),
                                   .phi = o2::gpu::GPUCommonMath::ATan2(dy, dx)});
        }
      });
    }
    BOOST_TEST(outputAt == output.offsets()[sourceElement + 1U]);
  }

  BOOST_REQUIRE(output.size() == expectedOrder.size());
  for (std::size_t i = 0; i < expectedOrder.size(); ++i) {
    BOOST_TEST(view[i].source == expectedOrder[i].source);
    BOOST_TEST(view[i].target == expectedOrder[i].target);
    BOOST_TEST(view[i].link == expectedOrder[i].link);
    BOOST_TEST(view[i].lambda == expectedOrder[i].lambda);
    BOOST_TEST(view[i].phi == expectedOrder[i].phi);
  }

  // This source's window wraps around phi = 0. The high-phi cell is visited
  // before the low-phi cell, even though its global rows come later.
  if (checkSpecialCases) {
    const auto wrappedSource = static_cast<std::uint32_t>(std::ranges::find(fixture.columns.clusterIndex, 0U) - fixture.columns.clusterIndex.begin());
    const auto lowTarget = static_cast<std::uint32_t>(std::ranges::find(fixture.columns.clusterIndex, 3U) - fixture.columns.clusterIndex.begin());
    const auto highTarget = static_cast<std::uint32_t>(std::ranges::find(fixture.columns.clusterIndex, 4U) - fixture.columns.clusterIndex.begin());
    const auto wrappedElement = fixture.sources.element(0, wrappedSource);
    std::vector<std::uint32_t> wrappedOrder;
    for (auto i = output.offsets()[wrappedElement]; i < output.offsets()[wrappedElement + 1U]; ++i) {
      if (view[i].link == 2U) {
        wrappedOrder.push_back(view[i].target);
      }
    }
    BOOST_REQUIRE(wrappedOrder.size() == 2U);
    BOOST_TEST(wrappedOrder[0] == highTarget);
    BOOST_TEST(wrappedOrder[1] == lowTarget);

    constexpr std::size_t kFirstPassingVertex = 2U;
    if (vertexSource == VertexSource::kVertexList && selected.size() > kFirstPassingVertex) {
      const auto sourcePoint = pairPoint(fixture.frame->measurement(wrappedSource));
      const auto targetPoint = pairPoint(fixture.frame->measurement(lowTarget));
      for (std::size_t i = 0; i < kFirstPassingVertex; ++i) {
        if (!vertexUsed(parameters, fixture.frame->measurement(wrappedSource), selected[i])) {
          continue;
        }
        const Scalar scattering = sourceScattering<StageToy>(parameters.scattering, parameters.pairCut, selected[i], LayerID{0}, sourcePoint);
        BOOST_TEST(!pairCompatible(parameters.pairCut, selected[i], sourcePoint, scattering, targetPoint));
      }
      BOOST_REQUIRE(vertexUsed(parameters, fixture.frame->measurement(wrappedSource), selected[kFirstPassingVertex]));
      const Scalar passingScattering =
        sourceScattering<StageToy>(parameters.scattering, parameters.pairCut, selected[kFirstPassingVertex], LayerID{0}, sourcePoint);
      BOOST_TEST(pairCompatible(parameters.pairCut, selected[kFirstPassingVertex], sourcePoint, passingScattering, targetPoint));
      bool found = false;
      for (std::size_t i = 0; i < output.size(); ++i) {
        found = found || (view[i].source == wrappedSource && view[i].target == lowTarget && view[i].link == 2U);
      }
      BOOST_TEST(found);
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(BeamLineMatchesBruteForceAndSerialOrder) { verify(VertexSource::kBeamLine); }

BOOST_AUTO_TEST_CASE(VertexListMatchesBruteForceAndAcceptsAtLeastOneVertex) { verify(VertexSource::kVertexList); }

/// Without skips, 0 -> 2 makes only the pairs past layer 1's ends: from
/// (2, 1) to (8, 4) the segment is at z 2 at r 4, beyond layer 1's |z| 1;
/// from (2, -0.5) to (8, 0) it crosses layer 1 at z -1/3.
BOOST_AUTO_TEST_CASE(WithoutSkipsAPairCrossingALayerIsNotATracklet)
{
  verify(VertexSource::kBeamLine, {}, true, false);
  verify(VertexSource::kVertexList, {}, true, false);
  const std::array<double, 3> material{};
  const auto table = scatteringTable<StageToy>(material);
  const auto withSkips = iterationParameters<StageToy>({.windows = windowParameters(), .minClusters = kChainClusters}, Fixture::kBeamConstraint, table);
  const auto without =
    iterationParameters<StageToy>({.windows = windowParameters(), .minClusters = kChainClusters, .skipLinks = false}, Fixture::kBeamConstraint, table);
  BOOST_TEST(withSkips.holeLayers[0] == 0U);
  BOOST_TEST(without.holeLayers[0] == 0b010U);
  BOOST_TEST(without.holeLayers[1] == 0U);
  BOOST_TEST(without.holeLayers[2] == 0U);
  const auto& spans = without.layerSpans;
  BOOST_TEST(!segmentLeavesHole(without.holeLayers[0], without.barrelLayers, spans, 2.F, 1.F, 8.F, 4.F));
  BOOST_TEST(segmentLeavesHole(without.holeLayers[0], without.barrelLayers, spans, 2.F, -0.5F, 8.F, 0.F));
}

BOOST_AUTO_TEST_CASE(VertexListWithOverlappingIntervalsMatchesBruteForce)
{
  const std::array<VertexSeed, 2> vertices{{
    {.timePS = 109'500U, .z = 0.F, .varZ = 1.e-4F},
    {.timePS = 110'500U, .z = 0.01F, .varZ = 1.e-4F},
  }};
  verify(VertexSource::kVertexList, vertices, false);
}

BOOST_AUTO_TEST_CASE(VertexListWithDisjointIntervalsMatchesBruteForce)
{
  const std::array<VertexSeed, 2> vertices{{
    {.timePS = 109'500U, .z = 0.F, .varZ = 1.e-4F},
    {.timePS = 110'500U, .z = 1.F, .varZ = 1.e-4F},
  }};
  verify(VertexSource::kVertexList, vertices, false);
}

BOOST_AUTO_TEST_CASE(VertexListBeyondCacheMatchesBruteForceAndSerialOrder)
{
  const std::array<VertexSeed, 4> vertices{{
    {.timePS = 109'400U, .z = 1.2F},
    {.timePS = 109'700U, .z = 1.2F},
    {.timePS = 110'000U, .z = 0.F},
    {.timePS = 110'300U, .z = -1.F},
  }};
  verify(VertexSource::kVertexList, vertices);
}

/// A halo source works only within the gate's window of box starts; a core
/// source always does. Every fixture row's box starts at 92 ns.
BOOST_AUTO_TEST_CASE(TheSourceGateSkipsHaloSourcesOutsideItsWindow)
{
  Fixture fixture;
  const IterationSettings settings{.windows = windowParameters(), .vertexSource = VertexSource::kBeamLine, .minClusters = kChainClusters};
  const std::array<double, 3> material{};
  const auto parameters = iterationParameters<StageToy>(settings, Fixture::kBeamConstraint, scatteringTable<StageToy>(material));
  const auto vertexWindows = vertexWindowParameters<StageToy>(parameters, settings.windows, {});
  std::array<MeasurementRange, Sources::kLayerCount> everything{};
  std::array<MeasurementRange, Sources::kLayerCount> nothing{};
  for (std::size_t layer = 0; layer < everything.size(); ++layer) {
    everything[layer] = fixture.frame->grid().layerRows(static_cast<LayerID>(layer));
  }
  const auto run = [&](const SourceGate& gate) {
    const TrackletStage<StageToy> stage{parameters,
                                        fixture.sources,
                                        *fixture.frame,
                                        fixture.geometry.allBlocks(),
                                        fixture.beamLine,
                                        {},
                                        fixture.maxCandidateWidthPS,
                                        vertexWindows,
                                        fixture.used.data(),
                                        gate};
    return runStage(fixture, stage);
  };

  const auto ungated = run({});
  BOOST_REQUIRE(ungated->size() > 0U);
  BOOST_TEST(sameStageOutput(*run({.core = everything.data(), .earliestPS = 200'000, .latestPS = 300'000}), *ungated));
  BOOST_TEST(sameStageOutput(*run({.core = nothing.data(), .earliestPS = 0, .latestPS = 120'000}), *ungated));
  BOOST_TEST(run({.core = nothing.data(), .earliestPS = 200'000, .latestPS = 300'000})->size() == 0U);
}

/// An earlier iteration's exclusive rows are neither sources nor targets; its
/// shared ones are both.
BOOST_AUTO_TEST_CASE(TheUsedMaskRemovesExclusiveRows)
{
  Fixture fixture;
  const IterationSettings settings{.windows = windowParameters(), .vertexSource = VertexSource::kBeamLine, .minClusters = kChainClusters};
  const std::array<double, 3> material{};
  const auto parameters = iterationParameters<StageToy>(settings, Fixture::kBeamConstraint, scatteringTable<StageToy>(material));
  const auto vertexWindows = vertexWindowParameters<StageToy>(parameters, settings.windows, {});
  const auto run = [&](const std::vector<std::uint8_t>& used) {
    const TrackletStage<StageToy> stage{
      parameters, fixture.sources, *fixture.frame, fixture.geometry.allBlocks(), fixture.beamLine, {}, fixture.maxCandidateWidthPS, vertexWindows, used.data()};
    return runStage(fixture, stage);
  };
  const auto free = run(fixture.used);
  BOOST_REQUIRE(free->size() > 0U);
  const auto view = free->view();
  const auto source = view.source[0];
  const auto target = view.target[0];
  const auto count = [&](const Output& output, auto keep) {
    std::size_t result = 0U;
    for (std::uint32_t i = 0; i < output.size(); ++i) {
      result += keep(output.view().source[i], output.view().target[i]) ? 1U : 0U;
    }
    return result;
  };
  const auto all = [](std::uint32_t, std::uint32_t) { return true; };

  // Shared: nothing changes.
  auto used = fixture.used;
  used[source] = kShared;
  used[target] = kShared;
  BOOST_TEST(sameStageOutput(*run(used), *free));

  // An exclusive source: its tracklets go, and only they.
  used = fixture.used;
  used[source] = kExclusive;
  const auto withoutSource = run(used);
  const auto fromSource = count(*free, [&](std::uint32_t from, std::uint32_t) { return from == source; });
  BOOST_TEST(fromSource > 0U);
  BOOST_TEST(count(*withoutSource, [&](std::uint32_t from, std::uint32_t to) { return from == source || to == source; }) == 0U);
  BOOST_TEST(count(*withoutSource, all) == count(*free, [&](std::uint32_t from, std::uint32_t to) { return from != source && to != source; }));

  // An exclusive target: the tracklets to it go, and only they.
  used = fixture.used;
  used[target] = kExclusive;
  const auto withoutTarget = run(used);
  BOOST_TEST(count(*withoutTarget, [&](std::uint32_t from, std::uint32_t to) { return from == target || to == target; }) == 0U);
  BOOST_TEST(count(*withoutTarget, all) == count(*free, [&](std::uint32_t from, std::uint32_t to) { return from != target && to != target; }));
}

BOOST_AUTO_TEST_CASE(ForcedVertexCapacityOneKeepsTheSameTracklets)
{
  Fixture fixture;
  const IterationSettings settings{.windows = windowParameters(), .vertexSource = VertexSource::kVertexList, .minClusters = kChainClusters};
  const std::array<double, 3> material{};
  const auto parameters = iterationParameters<StageToy>(settings, Fixture::kBeamConstraint, scatteringTable<StageToy>(material));
  // Beyond the beam spot's 1.5 cm reach and compatible in time: not used, so
  // no window has to hold it and the capacity cannot change the result.
  std::array<VertexSeed, 5> vertices{};
  std::ranges::copy(fixture.vertices, vertices.begin());
  vertices.back() = {.timePS = 110'200U, .z = 5.F};
  std::ranges::sort(vertices, {}, &VertexSeed::timePS);
  const auto vertexWindows = vertexWindowParameters<StageToy>(parameters, settings.windows, vertices);
  const TrackletStage<StageToy> normal{
    parameters, fixture.sources, *fixture.frame, fixture.geometry.allBlocks(), fixture.beamLine, vertices, fixture.maxCandidateWidthPS,
    vertexWindows, fixture.used.data()};
  const TrackletStage<StageToy, 1U> forced{
    parameters, fixture.sources, *fixture.frame, fixture.geometry.allBlocks(), fixture.beamLine, vertices, fixture.maxCandidateWidthPS,
    vertexWindows, fixture.used.data()};
  BOOST_TEST(!vertexUsable(parameters, kRowBox, vertices[2]));

  const auto normalOutput = runStage(fixture, normal);
  const auto forcedOutput = runStage(fixture, forced);
  BOOST_TEST(sameStageOutput(*normalOutput, *forcedOutput));
}

/// The beam-line stage makes a beam-line iteration's tracklets, and serves no
/// other iteration.
BOOST_AUTO_TEST_CASE(TheBeamLineStageKeepsTheSameTracklets)
{
  Fixture fixture;
  const std::array<double, 3> material{};
  const IterationSettings settings{.windows = windowParameters(), .vertexSource = VertexSource::kBeamLine, .minClusters = kChainClusters};
  const auto parameters = iterationParameters<StageToy>(settings, Fixture::kBeamConstraint, scatteringTable<StageToy>(material));
  const auto vertexWindows = vertexWindowParameters<StageToy>(parameters, settings.windows, {});
  const auto run = [&](const auto& stage) {
    return runStage(fixture, stage);
  };
  const TrackletStage<StageToy> general{
    parameters, fixture.sources, *fixture.frame, fixture.geometry.allBlocks(), fixture.beamLine, {}, fixture.maxCandidateWidthPS, vertexWindows, fixture.used.data()};
  const TrackletStageFor<StageToy, true> beamLine{
    parameters, fixture.sources, *fixture.frame, fixture.geometry.allBlocks(), fixture.beamLine, {}, fixture.maxCandidateWidthPS, vertexWindows, fixture.used.data()};
  const auto expected = run(general);
  BOOST_REQUIRE(expected->size() > 0U);
  BOOST_TEST(sameStageOutput(*run(beamLine), *expected));

  const IterationSettings listSettings{.windows = windowParameters(), .vertexSource = VertexSource::kVertexList, .minClusters = kChainClusters};
  const auto listParameters = iterationParameters<StageToy>(listSettings, Fixture::kBeamConstraint, scatteringTable<StageToy>(material));
  const auto listWindows = vertexWindowParameters<StageToy>(listParameters, listSettings.windows, fixture.vertices);
  BOOST_CHECK_THROW((TrackletStageFor<StageToy, true>{listParameters, fixture.sources, *fixture.frame, fixture.geometry.allBlocks(), fixture.beamLine,
                                                      fixture.vertices, fixture.maxCandidateWidthPS, listWindows, fixture.used.data()}),
                    std::invalid_argument);
}

/// A displaced iteration: its stage makes
/// exactly the pairs its test admits -- against every row of the target, as
/// the windows must hold them all -- and more than the prompt stage, whose
/// beam point admits fewer; neither stage serves the other's iteration, and
/// a displaced iteration needs the beam line and a z bound.
BOOST_AUTO_TEST_CASE(TheDisplacedStageMakesThePairsItsTestAdmits)
{
  Fixture fixture;
  const std::array<double, 3> material{};
  const auto table = scatteringTable<StageToy>(material);
  auto windows = windowParameters();
  windows.impactMax = 1.;
  windows.impactZMax = 10.;
  const IterationSettings settings{.windows = windows, .vertexSource = VertexSource::kBeamLine, .minClusters = kChainClusters};
  const auto parameters = iterationParameters<StageToy>(settings, Fixture::kBeamConstraint, table);
  const IterationSettings promptSettings{.windows = windowParameters(), .vertexSource = VertexSource::kBeamLine, .minClusters = kChainClusters};
  const auto prompt = iterationParameters<StageToy>(promptSettings, Fixture::kBeamConstraint, table);
  const auto vertexWindows = vertexWindowParameters<StageToy>(parameters, settings.windows, {});
  for (std::size_t link = 0; link < StageToy::layerLinks.size(); ++link) {
    BOOST_TEST(parameters.windows[link].phiHalfWidth >= prompt.windows[link].phiHalfWidth);
    BOOST_TEST(std::isinf(parameters.sourceMax[link]));
  }
  const auto run = [&](const auto& stage) {
    return runStage(fixture, stage);
  };
  const TrackletStageFor<StageToy, true, true> stage{
    parameters, fixture.sources, *fixture.frame, fixture.geometry.allBlocks(), fixture.beamLine, {}, fixture.maxCandidateWidthPS, vertexWindows, fixture.used.data()};
  const auto output = run(stage);
  const auto view = std::as_const(*output).view();
  std::size_t admitted = 0U;
  for (Sources::Index sourceElement = 0; sourceElement < fixture.sources.size(); ++sourceElement) {
    const auto sourceAt = fixture.sources.row(sourceElement);
    const auto source = fixture.frame->measurement(sourceAt.row);
    const auto sourcePoint = pairPoint(source);
    std::vector<std::pair<std::size_t, std::uint32_t>> expected;
    for (auto at = outLinks<StageToy>.begin(sourceAt.layer); at < outLinks<StageToy>.end(sourceAt.layer); ++at) {
      const std::size_t link = outLinks<StageToy>.links[at];
      const auto targets = fixture.frame->grid().layerRows(StageToy::layerLinks[link].to);
      const Scalar scattering = sourceScattering<StageToy>(parameters.scattering, parameters.pairCut, fixture.beamLine, sourceAt.layer, sourcePoint);
      for (auto targetRow = targets.begin; targetRow < targets.end; ++targetRow) {
        const auto target = fixture.frame->measurement(targetRow);
        if (!timeBox(source).isCompatible(timeBox(target)) ||
            segmentLeavesHole(parameters.holeLayers[link], parameters.barrelLayers, parameters.layerSpans, source.r, source.z, target.r, target.z)) {
          continue;
        }
        const auto result = displacedPairOf(parameters, fixture.beamLine, sourcePoint, scattering, pairPoint(target), 0.F);
        if (result.transverse && result.longitudinal) {
          expected.emplace_back(link, targetRow);
        }
      }
    }
    std::vector<std::pair<std::size_t, std::uint32_t>> emitted;
    for (auto at = output->offsets()[sourceElement]; at < output->offsets()[sourceElement + 1U]; ++at) {
      emitted.emplace_back(view[at].link, view[at].target);
    }
    std::ranges::sort(expected);
    std::ranges::sort(emitted);
    BOOST_TEST((emitted == expected));
    admitted += expected.size();
  }
  BOOST_TEST(admitted == output->size());
  const TrackletStageFor<StageToy, true> promptStage{
    prompt, fixture.sources, *fixture.frame, fixture.geometry.allBlocks(), fixture.beamLine, {}, fixture.maxCandidateWidthPS, vertexWindows, fixture.used.data()};
  BOOST_TEST(output->size() > run(promptStage)->size());

  BOOST_CHECK_THROW((TrackletStageFor<StageToy, true>{parameters,
                                                      fixture.sources,
                                                      *fixture.frame,
                                                      fixture.geometry.allBlocks(),
                                                      fixture.beamLine,
                                                      {},
                                                      fixture.maxCandidateWidthPS,
                                                      vertexWindows,
                                                      fixture.used.data()}),
                    std::invalid_argument);
  BOOST_CHECK_THROW((TrackletStageFor<StageToy, true, true>{prompt,
                                                            fixture.sources,
                                                            *fixture.frame,
                                                            fixture.geometry.allBlocks(),
                                                            fixture.beamLine,
                                                            {},
                                                            fixture.maxCandidateWidthPS,
                                                            vertexWindows,
                                                            fixture.used.data()}),
                    std::invalid_argument);
  auto noZ = settings;
  noZ.windows.impactZMax = 0.;
  BOOST_CHECK_THROW(static_cast<void>(iterationParameters<StageToy>(noZ, Fixture::kBeamConstraint, table)), std::invalid_argument);
  auto list = settings;
  list.vertexSource = VertexSource::kVertexList;
  BOOST_CHECK_THROW(static_cast<void>(iterationParameters<StageToy>(list, Fixture::kBeamConstraint, table)), std::invalid_argument);
}

/// A stage bound to copies of its inputs -- what a device run binds to --
/// gives the same tracklets, reads the copies and nothing else, and leaves
/// the stage it was bound from untouched.
BOOST_AUTO_TEST_CASE(ABoundStageReadsItsBindingsOnly)
{
  Fixture fixture;
  for (const auto vertexSource : {VertexSource::kBeamLine, VertexSource::kVertexList}) {
    const IterationSettings settings{.windows = windowParameters(), .vertexSource = vertexSource, .minClusters = kChainClusters};
    const std::array<double, 3> material{};
    const auto parameters = iterationParameters<StageToy>(settings, Fixture::kBeamConstraint, scatteringTable<StageToy>(material));
    const auto vertexWindows = vertexWindowParameters<StageToy>(parameters, settings.windows, fixture.vertices);
    const std::array<MeasurementRange, Sources::kLayerCount> nothing{};
    const SourceGate gate{.core = nothing.data(), .earliestPS = 0, .latestPS = 120'000};
    const TrackletStage<StageToy> stage{parameters,
                                        fixture.sources,
                                        *fixture.frame,
                                        fixture.geometry.allBlocks(),
                                        fixture.beamLine,
                                        fixture.vertices,
                                        fixture.maxCandidateWidthPS,
                                        vertexWindows,
                                        fixture.used.data(),
                                        gate};
    const auto run = [&](const auto& generator) {
      auto output = std::make_unique<Output>();
      output->setMemoryPool(std::make_shared<BoundedMemoryResource>());
      output->prepare(fixture.sources.size());
      BOOST_REQUIRE(runSerial(*output, generator) == StageStatus::kComplete);
      return output;
    };
    const auto reference = run(stage);
    BOOST_REQUIRE(reference->size() > 0U);

    // Copies of everything the stage reads, at other addresses.
    const auto parametersCopy = std::make_unique<IterationParameters<StageToy>>(parameters);
    const auto sourcesCopy = std::make_unique<Sources>(fixture.sources);
    const auto frameCopy = std::make_unique<Frame>(*fixture.frame);
    const auto windowsCopy = std::make_unique<VertexWindowParameters<StageToy>>(vertexWindows);
    // The vertices the stage reads, as its iteration chose: the beam line or the list.
    const std::vector<VertexSeed> verticesCopy = vertexSource == VertexSource::kBeamLine
                                                   ? std::vector<VertexSeed>{fixture.beamLine}
                                                   : std::vector<VertexSeed>(fixture.vertices.begin(), fixture.vertices.end());
    const auto coreCopy = nothing;
    auto usedCopy = fixture.used;
    const TrackletBindings<StageToy> bindings{.parameters = parametersCopy.get(),
                                              .sources = sourcesCopy.get(),
                                              .frame = frameCopy.get(),
                                              .used = usedCopy.data(),
                                              .vertexWindows = windowsCopy.get(),
                                              .vertices = verticesCopy.data(),
                                              .gateCore = coreCopy.data()};
    const auto bound = stage.bind(bindings);
    BOOST_TEST(sameStageOutput(*run(bound), *reference));

    // The bound stage reads the copied mask, the original the fixture's.
    std::ranges::fill(usedCopy, kExclusive);
    BOOST_TEST(run(bound)->size() == 0U);
    BOOST_TEST(sameStageOutput(*run(stage), *reference));

    // A binding without an input the stage reads is refused.
    auto missing = bindings;
    missing.gateCore = nullptr;
    BOOST_CHECK_THROW(static_cast<void>(stage.bind(missing)), std::invalid_argument);
    missing = bindings;
    missing.vertices = nullptr;
    BOOST_CHECK_THROW(static_cast<void>(stage.bind(missing)), std::invalid_argument);
  }
}

BOOST_AUTO_TEST_CASE(CpuStageMatchesTheSerialReference)
{
  Fixture fixture;
  BOOST_REQUIRE(fixture.sources.size() > 64U);
  for (const auto vertexSource : {VertexSource::kBeamLine, VertexSource::kVertexList}) {
    const IterationSettings settings{.windows = windowParameters(), .vertexSource = vertexSource, .minClusters = kChainClusters};
    const std::array<double, 3> material{};
    const auto parameters = iterationParameters<StageToy>(settings, Fixture::kBeamConstraint, scatteringTable<StageToy>(material));
    const auto vertexWindows = vertexWindowParameters<StageToy>(parameters, settings.windows, fixture.vertices);
    const TrackletStage<StageToy> stage{
      parameters, fixture.sources, *fixture.frame, fixture.geometry.allBlocks(), fixture.beamLine, fixture.vertices, fixture.maxCandidateWidthPS,
      vertexWindows, fixture.used.data()};

    Output serial;
    serial.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    serial.prepare(fixture.sources.size());
    BOOST_REQUIRE((runSerial(serial, stage) == StageStatus::kComplete));
    BOOST_REQUIRE(serial.size() > 0U);

    for (const int threads : {1, 2, 4}) {
      for (const Output::Index chunk : {1U, 3U, 8U, 17U, 64U}) {
        Output parallel;
        const auto parallelPool = std::make_shared<BoundedMemoryResource>();
        parallel.setMemoryPool(parallelPool);
        parallel.prepare(fixture.sources.size());
        detail::TaskArena cpuArena{threads};
        Cpu cpu{cpuArena};
        cpu.setMemoryPool(parallelPool);
        const auto status = cpu.run(parallel, chunk, stage);
        BOOST_REQUIRE_MESSAGE(status == StageStatus::kComplete,
                              "vertex source " << static_cast<int>(vertexSource) << ", threads " << threads << ", chunk " << chunk);
        BOOST_TEST(sameStageOutput(parallel, serial), "vertex source " << static_cast<int>(vertexSource) << ", threads " << threads << ", chunk " << chunk);
      }
    }
  }
}

} // namespace o2::gtrk
