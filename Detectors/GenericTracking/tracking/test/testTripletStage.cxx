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

/// \file testTripletStage.cxx
/// \brief Triplet generation against a brute-force tracklet-pair reference.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TripletStage
#include <boost/test/unit_test.hpp>
#include <numbers>

#include "ToyBarrel4Disc8.h"
#include "ToyBarrel3.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Tracking/TrackletPaths.h"
#include "GenericTracking/Tracking/TripletPermutation.h"
#include "GenericTracking/Tracking/TripletStage.h"

#include "GridTestSupport.h"
#include "Oracles.h"
#include "StageOutputComparison.h"
#include "SerialReferences.h"
#include "UniformMaterial.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace o2::gtrk
{

namespace
{

using Det = detectors::ToyBarrel3;
using Frame = TimeFrame<Det>;
using Scratch = TimeFrameScratch<Det>;
using Sources = LayerRowMap<Det>;
using Tracklets = StageOutput<TrackletT>;
using Triplets = StageOutput<TripletT>;
using Permutation = TripletPermutation<Det>;
using Controller = Tracker<Det, UniformMaterial>;

/// The largest variances per layer over \a slice's halo rows.
template <DetectorLayout Layout>
std::array<TripletTargetBounds, static_cast<std::size_t>(Detector<Layout>::nLayers)> haloTargetBounds(const TimeFrame<Layout>& frame,
                                                                                                      const typename TimeFrame<Layout>::TimeSlice& slice)
{
  return boundsOver(frame, std::span<const MeasurementRange, HaloBounds<Layout>::kLayerCount>{slice.halo}).target;
}

constexpr std::size_t kTracks = 4U;
constexpr std::size_t kNoisePerLayer = 3U;
constexpr std::size_t kMeasurements = (3U * kTracks) + (3U * kNoisePerLayer);

struct InputPoint {
  LayerID layer = kInvalidLayerID;
  Scalar x = 0.F;
  Scalar y = 0.F;
  Scalar z = 0.F;
  Scalar r = 0.F;
  Scalar phi = 0.F;
};

SearchWindowParams windowParameters()
{
  SearchWindowParams result;
  result.minPtGeV = 1.;
  result.bzKiloGauss = 20.;
  result.materialX2X0 = 0.;
  result.positionSigma = 0.003;
  result.nSigma = 3.;
  return result;
}

GridParams gridParameters()
{
  GridParams result;
  result.blockWidth = 20'000;
  return test::withToyFrames(result);
}

/// The toy's three layers make triplets and no road: its chains are 3
/// clusters, or no link would be taken (`takenLinks`).
IterationSettings iterationSettings() { return {.windows = windowParameters(), .vertexSource = VertexSource::kBeamLine, .minClusters = 3}; }

InputPoint helixPoint(std::size_t track, LayerID layer)
{
  constexpr std::array<double, kTracks> kPhi{{0.10, 1.35, 3.05, 4.75}};
  constexpr std::array<double, kTracks> kSlope{{0.04, -0.06, 0.08, -0.03}};
  constexpr std::array<double, kTracks> kPt{{1.2, 1.8, 3.0, 1.4}};
  const double radius = layerDistance<Det>(layer);
  const double curvature = 1. / curvatureRadius(kPt[track], 20.);
  const double arc = 2. * std::asin(0.5 * curvature * radius) / curvature;
  const double direction = kPhi[track] + (curvature * arc);
  const double x = (std::sin(direction) - std::sin(kPhi[track])) / curvature;
  const double y = (std::cos(kPhi[track]) - std::cos(direction)) / curvature;
  double phi = std::atan2(y, x);
  if (phi < 0.) {
    phi += 2 * std::numbers::pi_v<double>;
  }
  return {.layer = layer,
          .x = static_cast<Scalar>(x),
          .y = static_cast<Scalar>(y),
          .z = static_cast<Scalar>(kSlope[track] * arc),
          .r = static_cast<Scalar>(radius),
          .phi = static_cast<Scalar>(phi)};
}

class Fixture
{
 public:
  static constexpr VertexConstraint kBeamConstraint{.z = 0., .sigmaZ = 0.02, .sigmaTransverse = 0.001};

  Fixture()
    : geometry(test::cappedGrid<Det>(20., gridParameters(), windowParameters(), kBeamConstraint, 16, 8)),
      offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U)
  {
    std::array<InputPoint, kMeasurements> input{};
    for (std::size_t track = 0; track < kTracks; ++track) {
      for (std::size_t layer = 0; layer < 3U; ++layer) {
        input[logicalIndex(track, layer)] = helixPoint(track, static_cast<LayerID>(layer));
      }
    }
    for (std::size_t layer = 0; layer < 3U; ++layer) {
      const double radius = layerDistance<Det>(static_cast<LayerID>(layer));
      for (std::size_t noise = 0; noise < kNoisePerLayer; ++noise) {
        const double phi = 0.55 + (1.85 * static_cast<double>(noise)) + (0.23 * static_cast<double>(layer));
        const std::size_t index = (3U * kTracks) + (layer * kNoisePerLayer) + noise;
        input[index] = {.layer = static_cast<LayerID>(layer),
                        .x = static_cast<Scalar>(radius * std::cos(phi)),
                        .y = static_cast<Scalar>(radius * std::sin(phi)),
                        .z = static_cast<Scalar>((noise % 2U == 0U ? 0.65 : -0.55) + (0.1 * static_cast<double>(layer))),
                        .r = static_cast<Scalar>(radius),
                        .phi = static_cast<Scalar>(phi)};
      }
    }

    frame = test::frameOf(geometry, input.size(), [&](std::size_t i) { return input[i]; }, [&](std::size_t i) { return test::toyTimeBox(10.); }, [&](std::size_t i) { const auto scale = static_cast<Scalar>(1U + (i % 3U)); return test::PointVariances{scale * 1.e-6F, scale * 2.e-6F, scale * 3.e-6F}; }, columns, offsets);
  }

  [[nodiscard]] static constexpr std::size_t logicalIndex(std::size_t track, std::size_t layer) { return (3U * track) + layer; }

  [[nodiscard]] std::uint32_t row(std::size_t logical) const
  {
    return static_cast<std::uint32_t>(std::ranges::find(columns.clusterIndex, static_cast<std::uint32_t>(logical)) - columns.clusterIndex.begin());
  }

  GridGeometry<Det> geometry;
  std::vector<std::uint32_t> offsets;
  MeasurementStore columns;
  std::unique_ptr<Frame> frame;
};

using Mixed = detectors::ToyBarrel4Disc8;
using MixedFrame = TimeFrame<Mixed>;
using MixedSources = LayerRowMap<Mixed>;
using MixedScratch = TimeFrameScratch<Mixed>;
using MixedController = Tracker<Mixed, UniformMaterial>;
constexpr std::size_t kMixedMeasurements = 10U;

struct MixedInputPoint {
  LayerID layer = kInvalidLayerID;
  double timeNS = 0.;
  Scalar x = 0.F;
  Scalar y = 0.F;
  Scalar z = 0.F;
  Scalar r = 0.F;
  Scalar phi = 0.F;
};

IterationSettings mixedIterationSettings() { return {.windows = windowParameters(), .vertexSource = VertexSource::kBeamLine}; }

MixedInputPoint mixedHelixPoint(LayerID layer, double phiAtOrigin, double slope, double ptGeV, double timeNS)
{
  const double curvature = 1. / curvatureRadius(ptGeV, 20.);
  double arc = 0.;
  if (Detector<Mixed>::isBarrelLayer(layer)) {
    const double radius = layerDistance<Mixed>(layer);
    arc = 2. * std::asin(0.5 * curvature * radius) / curvature;
  } else {
    arc = discZ<Mixed>(layer) / slope;
  }
  const double direction = phiAtOrigin + (curvature * arc);
  const double x = (std::sin(direction) - std::sin(phiAtOrigin)) / curvature;
  const double y = (std::cos(phiAtOrigin) - std::cos(direction)) / curvature;
  double phi = std::atan2(y, x);
  if (phi < 0.) {
    phi += 2 * std::numbers::pi_v<double>;
  }
  return {.layer = layer,
          .timeNS = timeNS,
          .x = static_cast<Scalar>(x),
          .y = static_cast<Scalar>(y),
          .z = static_cast<Scalar>(slope * arc),
          .r = static_cast<Scalar>(std::hypot(x, y)),
          .phi = static_cast<Scalar>(phi)};
}

std::uint16_t mixedLink(LayerID from, LayerID to)
{
  const auto link = linkOf<Mixed>(from, to);
  if (link == kNoLink) {
    BOOST_FAIL("no link " << int(from) << " -> " << int(to));
  }
  return static_cast<std::uint16_t>(link);
}

class MixedFixture
{
 public:
  static constexpr VertexConstraint kBeamConstraint{.z = 0., .sigmaZ = 0.02, .sigmaTransverse = 0.001};
  static constexpr std::size_t kCentral1 = 0U;
  static constexpr std::size_t kCentral2 = 1U;
  static constexpr std::size_t kCentral3 = 2U;
  static constexpr std::size_t kForward1 = 3U;
  static constexpr std::size_t kForward2 = 4U;
  static constexpr std::size_t kForward3 = 5U;
  static constexpr std::size_t kForwardDisc = 6U;
  static constexpr std::size_t kLate1 = 7U;
  static constexpr std::size_t kLate2 = 8U;
  static constexpr std::size_t kLate3 = 9U;

  MixedFixture()
    : geometry(test::cappedGrid<Mixed>(80., gridParameters(), windowParameters(), kBeamConstraint, 32, 16)),
      offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U)
  {
    // Flight times between the layers: a nanosecond or less, one tick for all.
    // From layer 2 the forward track continues to layer 3 and to the +22 cm
    // disc, and from the disc to layer 3: its points lie on its helix; that
    // it would cross the disc first does not matter to the stages.
    constexpr double time12 = 0.1;
    constexpr double time23 = 0.1;
    constexpr double time28 = 0.2;
    const std::array<MixedInputPoint, kMixedMeasurements> input{{
      mixedHelixPoint(1, 0.2, 0.15, 2., 10.),
      mixedHelixPoint(2, 0.2, 0.15, 2., 10. + time12),
      mixedHelixPoint(3, 0.2, 0.15, 2., 10. + time12 + time23),
      mixedHelixPoint(1, 2.2, 1.3, 2.5, 10.),
      mixedHelixPoint(2, 2.2, 1.3, 2.5, 10. + time12),
      mixedHelixPoint(3, 2.2, 1.3, 2.5, 10. + time12 + time23),
      mixedHelixPoint(8, 2.2, 1.3, 2.5, 10. + time12 + time28),
      mixedHelixPoint(1, 4.5, 0.25, 2., 10.),
      mixedHelixPoint(2, 4.5, 0.25, 2., 35. + time12),
      mixedHelixPoint(3, 4.5, 0.25, 2., 60. + time12 + time23),
    }};

    frame = test::frameOf(geometry, input.size(), [&](std::size_t i) { return input[i]; }, [&](std::size_t i) { return test::toyTimeBox(input[i].timeNS); }, [&](std::size_t i) { return test::PointVariances{9.e-6F, 9.e-6F, 9.e-6F}; }, columns, offsets);
  }

  [[nodiscard]] std::uint32_t row(std::size_t logical) const
  {
    return static_cast<std::uint32_t>(std::ranges::find(columns.clusterIndex, static_cast<std::uint32_t>(logical)) - columns.clusterIndex.begin());
  }

  GridGeometry<Mixed> geometry;
  std::vector<std::uint32_t> offsets;
  MeasurementStore columns;
  std::unique_ptr<MixedFrame> frame;
};

constexpr std::size_t kSliceTracks = 4U;
constexpr std::size_t kSliceMeasurements = 3U * kSliceTracks;

class SliceFixture
{
 public:
  static constexpr VertexConstraint kBeamConstraint = Fixture::kBeamConstraint;

  SliceFixture()
    : geometry(test::cappedGrid<Det>(140., gridParameters(), windowParameters(), kBeamConstraint, 16, 8)),
      offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U)
  {
    constexpr double flight01 = 0.2; // [ns], a few centimetres at c
    constexpr double flight12 = 0.2;
    // Track 1 crosses the 40 ns block edge: its first hit in the 25-50 ns
    // tick, whose box starts at 17 ns, its others in the next, starting at
    // 42 ns (the toys' boxes, GridTestSupport.h).
    const std::array<double, kSliceTracks> start{{10., 50. - (0.5 * flight01), 70., 110.}};

    std::array<InputPoint, kSliceMeasurements> input{};
    std::array<double, kSliceMeasurements> times{};
    for (std::size_t track = 0; track < kSliceTracks; ++track) {
      for (std::size_t layer = 0; layer < 3U; ++layer) {
        const std::size_t logical = logicalIndex(track, layer);
        input[logical] = helixPoint(track, static_cast<LayerID>(layer));
        times[logical] = start[track] + (layer == 0U ? 0. : flight01) + (layer == 2U ? flight12 : 0.);
      }
    }

    frame = test::frameOf(geometry, input.size(), [&](std::size_t i) { return input[i]; }, [&](std::size_t i) { return test::toyTimeBox(times[i]); }, [&](std::size_t i) { return test::PointVariances{1.e-6F, 2.e-6F, 3.e-6F}; }, columns, offsets);
  }

  [[nodiscard]] static constexpr std::size_t logicalIndex(std::size_t track, std::size_t layer) { return (3U * track) + layer; }

  [[nodiscard]] std::uint32_t row(std::size_t logical) const
  {
    return static_cast<std::uint32_t>(std::ranges::find(columns.clusterIndex, static_cast<std::uint32_t>(logical)) - columns.clusterIndex.begin());
  }

  GridGeometry<Det> geometry;
  std::vector<std::uint32_t> offsets;
  MeasurementStore columns;
  std::unique_ptr<Frame> frame;
};

/// Compare a tracker result with the triplet-stage output that preceded its
/// level sweep. All stage structure and semantic triplet fields must match;
/// only the level column is expected to differ.
bool sameTripletOutputIgnoringLevelSweep(const Triplets& lhs, const Triplets& rhs)
{
  if (!sameStageShape(lhs, rhs)) {
    return false;
  }
  const auto a = lhs.view();
  const auto b = rhs.view();
  const std::size_t size = lhs.size();
  return std::memcmp(a.first, b.first, size * sizeof(*a.first)) == 0 && std::memcmp(a.second, b.second, size * sizeof(*a.second)) == 0 &&
         std::memcmp(a.curvature, b.curvature, size * sizeof(*a.curvature)) == 0 &&
         std::memcmp(a.varCurvature, b.varCurvature, size * sizeof(*a.varCurvature)) == 0;
}

/// The path bounds at 3 clusters and any road end: no path is too short, so
/// the permutation drops only the first tracklets whose middle has no
/// tracklets and the stage applies its own cuts alone.
template <DetectorLayout Layout>
TrackletPathBounds everyPath(TrackletPaths<Layout>& paths, IterationParameters<Layout> parameters, const LayerRowMap<Layout>& sources,
                             const Tracklets& tracklets)
{
  parameters.minClusters = 3U;
  parameters.roadEndLayers = 0U;
  paths.setMemoryPool(std::make_shared<BoundedMemoryResource>());
  paths.build(parameters, sources, tracklets, SerialForEach{});
  return paths.bounds();
}

/// The permutation is every tracklet whose middle is a source with tracklets,
/// in arena order; returns how many tracklets it excludes because their
/// middle is not a source.
template <DetectorLayout Layout>
std::size_t checkPermutation(const IterationParameters<Layout>& parameters, const LayerRowMap<Layout>& sources, const Tracklets& tracklets,
                             const TripletPermutation<Layout>& permutation)
{
  const auto view = tracklets.view();
  const auto offsets = tracklets.offsets();
  std::vector<std::uint32_t> expected;
  std::size_t notSources = 0U;
  for (std::uint32_t first = 0; first < tracklets.size(); ++first) {
    const LayerID middleLayer = parameters.linkTargets[view.link[first]];
    if (!sources.contains(middleLayer, view.target[first])) {
      ++notSources;
      continue;
    }
    const auto element = sources.element(middleLayer, view.target[first]);
    if (offsets[element] != offsets[element + 1U]) {
      expected.push_back(first);
    }
  }
  BOOST_REQUIRE(permutation.size() == expected.size());
  for (std::size_t index = 0; index < expected.size(); ++index) {
    BOOST_TEST(permutation.indices()[index] == expected[index]);
  }
  return notSources;
}

template <DetectorLayout Layout>
std::uint32_t permutationElement(const TripletPermutation<Layout>& permutation, std::uint32_t first)
{
  const auto found = std::ranges::find(permutation.indices(), first);
  return static_cast<std::uint32_t>(found - permutation.indices().begin());
}

/// Each first tracklet's triplets in canonical order: by the second
/// tracklet's link in the middle layer's out-link order, then by lambda.
template <DetectorLayout Layout>
void checkOutputOrder(const IterationParameters<Layout>& parameters, const Tracklets& tracklets, const TripletPermutation<Layout>& permutation,
                      const Triplets& triplets)
{
  const auto trackletView = tracklets.view();
  const auto tripletView = triplets.view();
  for (std::uint32_t element = 0; element < permutation.size(); ++element) {
    const std::uint32_t first = permutation.indices()[element];
    const LayerID middleLayer = parameters.linkTargets[trackletView.link[first]];
    const auto begin = triplets.offsets()[element];
    const auto end = triplets.offsets()[element + 1U];
    std::uint32_t previousLinkOrder = 0U;
    Scalar previousLambda = 0.F;
    for (auto index = begin; index < end; ++index) {
      BOOST_TEST(tripletView.first[index] == first);
      const std::size_t secondLink = trackletView.link[tripletView.second[index]];
      const auto linkBegin = parameters.outLinkOffsets[middleLayer];
      const auto linkEnd = parameters.outLinkOffsets[static_cast<std::size_t>(middleLayer) + 1U];
      std::uint32_t linkOrder = linkEnd - linkBegin;
      for (auto out = linkBegin; out < linkEnd; ++out) {
        if (parameters.outLinkIndices[out] == secondLink) {
          linkOrder = out - linkBegin;
          break;
        }
      }
      BOOST_REQUIRE(linkOrder < linkEnd - linkBegin);
      const Scalar lambda = trackletView.lambda[tripletView.second[index]];
      if (index != begin) {
        BOOST_TEST(previousLinkOrder <= linkOrder);
        if (previousLinkOrder == linkOrder) {
          BOOST_TEST(previousLambda <= lambda);
        }
      }
      previousLinkOrder = linkOrder;
      previousLambda = lambda;
    }
  }
}

template <DetectorLayout Layout>
struct TripletStageChecks {
  using Parameters = IterationParameters<Layout>;
  using SourcesT = LayerRowMap<Layout>;
  using FrameT = TimeFrame<Layout>;
  using PermutationT = TripletPermutation<Layout>;

  TripletStageChecks(const Parameters& parameters, const SourcesT& sources, const FrameT& frame, const Tracklets& tracklets,
                     const TripletTargetBounds* targetBounds)
    : mParameters(&parameters), mSources(&sources), mFrame(&frame), mTracklets(&tracklets), mTargetBounds(targetBounds)
  {
  }

  TrackletPaths<Layout> paths;
  TrackletPathBounds bounds;
  PermutationT permutation;
  Triplets serial;
  std::size_t excludedFirstTracklets = 0U;

  void run()
  {
    bounds = everyPath(paths, *mParameters, *mSources, *mTracklets);
    permutation.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    BOOST_REQUIRE(permutation.build(*mParameters, *mSources, *mTracklets, bounds, SerialForEach{}, 256U) == StageStatus::kComplete);
    excludedFirstTracklets = checkPermutation(*mParameters, *mSources, *mTracklets, permutation);

    serial.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    serial.prepare(permutation.size());
    BOOST_REQUIRE(runSerial(serial, stage(permutation)) == StageStatus::kComplete);
    checkOutputOrder(*mParameters, *mTracklets, permutation, serial);
  }

  void checkBackends() const
  {
    const auto theStage = stage(permutation);
    for (const int threads : {1, 2, 4}) {
      for (const Triplets::Index chunk : {1U, 2U, 3U, 7U, 16U}) {
        Triplets parallel;
        const auto parallelPool = std::make_shared<BoundedMemoryResource>();
        parallel.setMemoryPool(parallelPool);
        parallel.prepare(permutation.size());
        detail::TaskArena cpuArena{threads};
        CpuStage<TripletT> cpu{cpuArena};
        cpu.setMemoryPool(parallelPool);
        BOOST_REQUIRE(cpu.run(parallel, chunk, theStage) == StageStatus::kComplete);
        BOOST_TEST(sameStageOutput(parallel, serial), "threads " << threads << ", chunk " << chunk);
      }
    }
  }

 private:
  [[nodiscard]] TripletStage<Layout> stage(const PermutationT& order) const
  {
    return {*mParameters, *mSources, *mFrame, *mTracklets, mTargetBounds, order.data(), bounds};
  }

  const Parameters* mParameters = nullptr;
  const SourcesT* mSources = nullptr;
  const FrameT* mFrame = nullptr;
  const Tracklets* mTracklets = nullptr;
  const TripletTargetBounds* mTargetBounds = nullptr;
};

/// Compare source-ordered stage output with an exhaustive continuation scan.
/// The return value counts candidates that pass the lambda window but fail a
/// later triplet cut.
template <DetectorLayout Layout>
std::size_t checkAgainstBruteForce(const IterationParameters<Layout>& parameters, const LayerRowMap<Layout>& sources, const TimeFrame<Layout>& frame,
                                   const Tracklets& tracklets, const TripletTargetBounds* targetBounds, const TripletPermutation<Layout>& sourceOrder,
                                   const Triplets& serial)
{
  const auto trackletView = tracklets.view();
  const auto tripletView = serial.view();
  std::size_t postLambdaRejected = 0U;
  for (std::uint32_t firstIndex = 0; firstIndex < tracklets.size(); ++firstIndex) {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> brute;
    const std::size_t firstLink = trackletView.link[firstIndex];
    const LayerID middleLayer = parameters.linkTargets[firstLink];
    const std::uint32_t middleRow = trackletView.target[firstIndex];
    if (sources.contains(middleLayer, middleRow)) {
      const auto first = timedPointAt(frame.measurements(), trackletView.source[firstIndex]);
      const auto middle = timedPointAt(frame.measurements(), middleRow);
      const auto middleElement = sources.element(middleLayer, middleRow);
      for (auto candidateIndex = tracklets.offsets()[middleElement]; candidateIndex < tracklets.offsets()[middleElement + 1U]; ++candidateIndex) {
        const std::size_t secondLink = trackletView.link[candidateIndex];
        const LayerID targetLayer = parameters.linkTargets[secondLink];
        const Scalar halfWidth = tripletLambdaHalfWidth(parameters, secondLink, middleLayer, first.point, middle.point, targetBounds[targetLayer]);
        if (trackletView.lambda[candidateIndex] < trackletView.lambda[firstIndex] - halfWidth ||
            trackletView.lambda[candidateIndex] > trackletView.lambda[firstIndex] + halfWidth) {
          continue;
        }
        const auto third = timedPointAt(frame.measurements(), trackletView.target[candidateIndex]);
        TripletFitResult fit;
        if (tripletCompatible(parameters, middleLayer, first, middle, third, fit)) {
          brute.emplace_back(firstIndex, candidateIndex);
        } else {
          ++postLambdaRejected;
        }
      }
    }

    std::vector<std::pair<std::uint32_t, std::uint32_t>> emitted;
    const auto sourceElement = permutationElement(sourceOrder, firstIndex);
    if (sourceElement < sourceOrder.size()) {
      for (auto index = serial.offsets()[sourceElement]; index < serial.offsets()[sourceElement + 1U]; ++index) {
        emitted.emplace_back(tripletView.first[index], tripletView.second[index]);
      }
    }
    std::ranges::sort(brute);
    std::ranges::sort(emitted);
    BOOST_REQUIRE(emitted.size() == brute.size());
    for (std::size_t index = 0; index < brute.size(); ++index) {
      BOOST_TEST(emitted[index].first == brute[index].first);
      BOOST_TEST(emitted[index].second == brute[index].second);
    }
  }
  return postLambdaRejected;
}

std::uint32_t findTracklet(const Tracklets& tracklets, std::uint32_t source, std::uint32_t target, std::uint16_t link)
{
  const auto view = tracklets.view();
  for (std::uint32_t index = 0; index < tracklets.size(); ++index) {
    if (view.source[index] == source && view.target[index] == target && view.link[index] == link) {
      return index;
    }
  }
  return tracklets.size();
}

using SliceTriplet = std::array<std::uint32_t, 5U>;

struct SliceRun {
  std::vector<SliceTriplet> triplets;
  std::uint64_t tracklets = 0U;
  std::uint64_t haloSourceTracklets = 0U;
};

SliceRun runTripletSlices(const SliceFixture& fixture, std::int32_t widthBlocks, std::int32_t firstBoundary)
{
  const std::array<double, Sources::kLayerCount> material{};
  const auto scattering = scatteringTable<Det>(material);
  const auto parameters = iterationParameters<Det>(iterationSettings(), SliceFixture::kBeamConstraint, scattering);
  const std::int64_t haloPS = sliceReach(frameBounds(*fixture.frame), fixture.geometry).haloPS();

  Controller tracker{iterationSettings(), SliceFixture::kBeamConstraint, scattering, std::make_shared<BoundedMemoryResource>(), 2U, 2};
  const auto blockCount = fixture.geometry.blockCount();
  const auto width = std::max(1, widthBlocks);
  const auto shift = firstBoundary > 0 ? std::min(firstBoundary, width - 1) : 0;
  SliceRun result;
  for (std::int32_t begin = 0; begin < blockCount;) {
    const auto firstWidth = begin == 0 && shift > 0 ? shift : width;
    const auto end = std::min(blockCount, begin + firstWidth);
    const Frame::TimeSlice slice{*fixture.frame, {.begin = begin, .end = end}, haloPS};
    BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
    const auto trackletResult = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
    BOOST_REQUIRE(trackletResult.status == StageStatus::kComplete);
    BOOST_REQUIRE(trackletResult.tracklets != nullptr);
    const auto& tracklets = *trackletResult.tracklets;
    const auto trackletView = tracklets.view();
    result.tracklets += tracklets.size();
    for (std::uint32_t index = 0; index < tracklets.size(); ++index) {
      const LayerID sourceLayer = parameters.linkSources[trackletView.link[index]];
      const auto core = slice.core[sourceLayer];
      result.haloSourceTracklets += trackletView.source[index] >= core.begin && trackletView.source[index] < core.end ? 0U : 1U;
    }

    const Sources sources = Scratch::makeSources(slice.halo);
    TrackletPaths<Det> paths;
    const auto bounds = everyPath(paths, parameters, sources, tracklets);
    Permutation permutation;
    permutation.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    BOOST_REQUIRE(permutation.build(parameters, sources, tracklets, bounds, SerialForEach{}, 256U) == StageStatus::kComplete);
    const auto targetBounds = haloTargetBounds(*fixture.frame, slice);
    const TripletStage<Det> stage{parameters, sources, *fixture.frame, tracklets, targetBounds.data(), permutation.data(), bounds};
    Triplets triplets;
    triplets.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    triplets.prepare(permutation.size());
    BOOST_REQUIRE(runSerial(triplets, stage) == StageStatus::kComplete);
    // The core-owned triplets, first measurement in the core, partition the
    // triplets over the slices; halo-owned ones are built for the chains.
    const auto tripletView = triplets.view();
    for (std::uint32_t index = 0; index < triplets.size(); ++index) {
      const std::uint32_t first = tripletView.first[index];
      const std::uint32_t second = tripletView.second[index];
      BOOST_REQUIRE(trackletView.target[first] == trackletView.source[second]);
      const auto firstCore = slice.core[parameters.linkSources[trackletView.link[first]]];
      if (trackletView.source[first] < firstCore.begin || trackletView.source[first] >= firstCore.end) {
        continue;
      }
      result.triplets.push_back({trackletView.source[first], trackletView.target[first], trackletView.target[second],
                                 std::bit_cast<std::uint32_t>(tripletView.curvature[index]), std::bit_cast<std::uint32_t>(tripletView.varCurvature[index])});
    }
    begin = end;
  }
  std::ranges::sort(result.triplets);
  return result;
}

void checkSliceTriplets(const std::vector<SliceTriplet>& actual, const std::vector<SliceTriplet>& expected)
{
  BOOST_REQUIRE(actual.size() == expected.size());
  for (std::size_t triplet = 0; triplet < expected.size(); ++triplet) {
    for (std::size_t field = 0; field < expected[triplet].size(); ++field) {
      BOOST_TEST(actual[triplet][field] == expected[triplet][field]);
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(TripletsMatchBruteForceAndAllTrueHelicesAreFound)
{
  const Fixture fixture;
  const Frame::TimeSlice slice{*fixture.frame, fixture.geometry.allBlocks(), 0};
  const std::array<double, Sources::kLayerCount> material{};
  const auto scattering = scatteringTable<Det>(material);
  Controller tracker{iterationSettings(), Fixture::kBeamConstraint, scattering, std::make_shared<BoundedMemoryResource>(), 2U, 2};
  BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
  const auto trackletResult = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
  BOOST_REQUIRE(trackletResult.status == StageStatus::kComplete);
  BOOST_REQUIRE(trackletResult.tracklets != nullptr);
  const auto& tracklets = *trackletResult.tracklets;
  BOOST_REQUIRE(tracklets.size() > 0U);

  const Sources sources = Scratch::makeSources(slice.halo);
  const auto parameters = iterationParameters<Det>(iterationSettings(), Fixture::kBeamConstraint, scattering);
  const auto targetBounds = haloTargetBounds(*fixture.frame, slice);
  BOOST_TEST(targetBounds[2].varZ == 3.e-6F);
  BOOST_TEST(targetBounds[2].varR == 6.e-6F);
  BOOST_TEST(targetBounds[2].varRPhi == 9.e-6F);
  TripletStageChecks<Det> checks{parameters, sources, *fixture.frame, tracklets, targetBounds.data()};
  checks.run();
  BOOST_TEST(checks.excludedFirstTracklets > 0U);
  const auto& sourceOrder = checks.permutation;
  const auto& serial = checks.serial;

  const auto tripletView = serial.view();
  static_cast<void>(checkAgainstBruteForce(parameters, sources, *fixture.frame, tracklets, targetBounds.data(), sourceOrder, serial));

  std::size_t trueTriplets = 0U;
  for (std::size_t track = 0; track < kTracks; ++track) {
    const auto first = findTracklet(tracklets, fixture.row(Fixture::logicalIndex(track, 0U)), fixture.row(Fixture::logicalIndex(track, 1U)), 0U);
    const auto second = findTracklet(tracklets, fixture.row(Fixture::logicalIndex(track, 1U)), fixture.row(Fixture::logicalIndex(track, 2U)), 1U);
    BOOST_REQUIRE(first < tracklets.size());
    BOOST_REQUIRE(second < tracklets.size());
    bool found = false;
    for (std::uint32_t index = 0; index < serial.size(); ++index) {
      found = found || (tripletView.first[index] == first && tripletView.second[index] == second);
    }
    BOOST_TEST(found);
    trueTriplets += found ? 1U : 0U;
  }
  BOOST_TEST(trueTriplets == kTracks);
  BOOST_TEST_MESSAGE("fixture triplets " << serial.size() << ", true triplets " << trueTriplets);

  checks.checkBackends();
}

BOOST_AUTO_TEST_CASE(MultipleContinuationsDiscPathAndPostLambdaRejection)
{
  const MixedFixture fixture;
  const MixedFrame::TimeSlice slice{*fixture.frame, fixture.geometry.allBlocks(), 0};
  const std::array<double, MixedSources::kLayerCount> material{};
  const auto scattering = scatteringTable<Mixed>(material);
  MixedController tracker{mixedIterationSettings(), MixedFixture::kBeamConstraint, scattering, std::make_shared<BoundedMemoryResource>(), 2U, 2};
  BOOST_REQUIRE(tracker.beginTimeframe(*fixture.frame) == StageStatus::kComplete);
  const auto trackletResult = test::runSlice(tracker, 0U, *fixture.frame, slice, {});
  BOOST_REQUIRE(trackletResult.status == StageStatus::kComplete);
  BOOST_REQUIRE(trackletResult.tracklets != nullptr);
  const auto& tracklets = *trackletResult.tracklets;

  const auto link12 = mixedLink(1, 2);
  const auto link23 = mixedLink(2, 3);
  const auto link28 = mixedLink(2, 8);
  const auto link83 = mixedLink(8, 3);
  const auto centralFirst = findTracklet(tracklets, fixture.row(MixedFixture::kCentral1), fixture.row(MixedFixture::kCentral2), link12);
  const auto centralSecond = findTracklet(tracklets, fixture.row(MixedFixture::kCentral2), fixture.row(MixedFixture::kCentral3), link23);
  const auto forwardFirst = findTracklet(tracklets, fixture.row(MixedFixture::kForward1), fixture.row(MixedFixture::kForward2), link12);
  const auto forwardThird = findTracklet(tracklets, fixture.row(MixedFixture::kForward2), fixture.row(MixedFixture::kForward3), link23);
  const auto forwardDisc = findTracklet(tracklets, fixture.row(MixedFixture::kForward2), fixture.row(MixedFixture::kForwardDisc), link28);
  const auto discNext = findTracklet(tracklets, fixture.row(MixedFixture::kForwardDisc), fixture.row(MixedFixture::kForward3), link83);
  const auto lateFirst = findTracklet(tracklets, fixture.row(MixedFixture::kLate1), fixture.row(MixedFixture::kLate2), link12);
  const auto lateSecond = findTracklet(tracklets, fixture.row(MixedFixture::kLate2), fixture.row(MixedFixture::kLate3), link23);
  const std::array<std::pair<const char*, std::uint32_t>, 8> requiredTracklets{{
    {"central first", centralFirst},
    {"central second", centralSecond},
    {"forward first", forwardFirst},
    {"forward third", forwardThird},
    {"forward disc", forwardDisc},
    {"disc next", discNext},
    {"late first", lateFirst},
    {"late second", lateSecond},
  }};
  for (const auto& [name, index] : requiredTracklets) {
    BOOST_REQUIRE_MESSAGE(index < tracklets.size(), "missing " << name);
  }

  const MixedSources sources = MixedScratch::makeSources(slice.halo);
  const auto parameters = iterationParameters<Mixed>(mixedIterationSettings(), MixedFixture::kBeamConstraint, scattering);
  const auto targetBounds = haloTargetBounds(*fixture.frame, slice);
  TripletStageChecks<Mixed> checks{parameters, sources, *fixture.frame, tracklets, targetBounds.data()};
  checks.run();
  // No link leaves layer 3: the tracklets ending on it start no triplet.
  std::size_t endingOnLayer3 = 0U;
  for (std::uint32_t index = 0; index < tracklets.size(); ++index) {
    endingOnLayer3 += parameters.linkTargets[tracklets.view().link[index]] == 3U ? 1U : 0U;
  }
  BOOST_TEST(endingOnLayer3 > 0U);
  BOOST_TEST(checks.excludedFirstTracklets == endingOnLayer3);
  const auto& sourceOrder = checks.permutation;
  const auto& serial = checks.serial;
  // The tracker runs the same stage itself, with the path bounds at the
  // iteration's own length.
  BOOST_REQUIRE(trackletResult.triplets != nullptr);
  {
    TrackletPaths<Mixed> paths;
    paths.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    paths.build(parameters, sources, tracklets, SerialForEach{});
    TripletPermutation<Mixed> order;
    order.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    BOOST_REQUIRE(order.build(parameters, sources, tracklets, paths.bounds(), SerialForEach{}, 256U) == StageStatus::kComplete);
    const TripletStage<Mixed> stage{parameters, sources, *fixture.frame, tracklets, targetBounds.data(), order.data(), paths.bounds()};
    Triplets own;
    own.setMemoryPool(std::make_shared<BoundedMemoryResource>());
    own.prepare(order.size());
    BOOST_REQUIRE(runSerial(own, stage) == StageStatus::kComplete);
    BOOST_TEST(sameTripletOutputIgnoringLevelSweep(*trackletResult.triplets, own));
  }
  const auto tripletView = serial.view();
  const auto trackletView = tracklets.view();
  const std::size_t postLambdaRejected = checkAgainstBruteForce(parameters, sources, *fixture.frame, tracklets, targetBounds.data(), sourceOrder, serial);
  BOOST_TEST(postLambdaRejected > 0U);

  const auto hasTriplet = [&](std::uint32_t first, std::uint32_t second) {
    for (std::uint32_t index = 0; index < serial.size(); ++index) {
      if (tripletView.first[index] == first && tripletView.second[index] == second) {
        return true;
      }
    }
    return false;
  };

  BOOST_TEST(hasTriplet(centralFirst, centralSecond));
  BOOST_TEST(hasTriplet(forwardFirst, forwardThird));
  BOOST_TEST(hasTriplet(forwardFirst, forwardDisc));
  BOOST_TEST(hasTriplet(forwardDisc, discNext));

  std::uint32_t forwardThirdAt = serial.size();
  std::uint32_t forwardDiscAt = serial.size();
  const auto forwardElement = permutationElement(sourceOrder, forwardFirst);
  BOOST_REQUIRE(forwardElement < sourceOrder.size());
  for (auto index = serial.offsets()[forwardElement]; index < serial.offsets()[forwardElement + 1U]; ++index) {
    if (tripletView.second[index] == forwardThird) {
      forwardThirdAt = index;
    }
    if (tripletView.second[index] == forwardDisc) {
      forwardDiscAt = index;
    }
  }
  BOOST_TEST(forwardThirdAt < forwardDiscAt);

  const auto first = timedPointAt(fixture.frame->measurements(), trackletView.source[lateFirst]);
  const auto middle = timedPointAt(fixture.frame->measurements(), trackletView.target[lateFirst]);
  const auto third = timedPointAt(fixture.frame->measurements(), trackletView.target[lateSecond]);
  const Scalar halfWidth = tripletLambdaHalfWidth(parameters, link23, LayerID{2}, first.point, middle.point, targetBounds[LayerID{3}]);
  BOOST_TEST(trackletView.lambda[lateSecond] >= trackletView.lambda[lateFirst] - halfWidth);
  BOOST_TEST(trackletView.lambda[lateSecond] <= trackletView.lambda[lateFirst] + halfWidth);
  BOOST_TEST(!tripletTimeCompatible(first, third));
  BOOST_TEST(!hasTriplet(lateFirst, lateSecond));

  checks.checkBackends();
}

BOOST_AUTO_TEST_CASE(TripletsAreInvariantUnderSliceWidthAndOrigin)
{
  const SliceFixture fixture;
  const auto whole = runTripletSlices(fixture, fixture.geometry.blockCount(), 0);
  const auto oneBlock = runTripletSlices(fixture, 1, 0);
  const auto twoBlocks = runTripletSlices(fixture, 2, 0);
  const auto shifted = runTripletSlices(fixture, 2, 1);

  BOOST_REQUIRE(whole.triplets.size() == kSliceTracks);
  checkSliceTriplets(oneBlock.triplets, whole.triplets);
  checkSliceTriplets(twoBlocks.triplets, whole.triplets);
  checkSliceTriplets(shifted.triplets, whole.triplets);

  const bool duplicateRows = std::ranges::adjacent_find(whole.triplets, [](const SliceTriplet& lhs, const SliceTriplet& rhs) {
                               return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
                             }) != whole.triplets.end();
  BOOST_TEST(!duplicateRows);

  const std::uint32_t crossingFirst = fixture.row(SliceFixture::logicalIndex(1U, 0U));
  const std::uint32_t crossingMiddle = fixture.row(SliceFixture::logicalIndex(1U, 1U));
  const std::uint32_t crossingThird = fixture.row(SliceFixture::logicalIndex(1U, 2U));
  const bool crossingFound = std::ranges::any_of(
    whole.triplets, [&](const SliceTriplet& triplet) { return triplet[0] == crossingFirst && triplet[1] == crossingMiddle && triplet[2] == crossingThird; });
  BOOST_TEST(crossingFound);

  const std::int64_t haloPS = sliceReach(frameBounds(*fixture.frame), fixture.geometry).haloPS();
  const Frame::TimeSlice left{*fixture.frame, {.begin = 0, .end = 2}, haloPS};
  const Frame::TimeSlice right{*fixture.frame, {.begin = 2, .end = 4}, haloPS};
  BOOST_TEST(crossingFirst >= left.core[0].begin);
  BOOST_TEST(crossingFirst < left.core[0].end);
  BOOST_TEST(crossingMiddle >= right.core[1].begin);
  BOOST_TEST(crossingMiddle < right.core[1].end);
  BOOST_TEST(crossingMiddle >= left.halo[1].begin);
  BOOST_TEST(crossingMiddle < left.halo[1].end);

  for (const auto* run : {&oneBlock, &twoBlocks, &shifted}) {
    BOOST_REQUIRE(run->tracklets > 0U);
    BOOST_TEST(run->haloSourceTracklets > 0U);
  }
  const auto fraction = [](const SliceRun& run) { return static_cast<double>(run.haloSourceTracklets) / static_cast<double>(run.tracklets); };
  BOOST_TEST_MESSAGE("halo-source tracklet fractions: one-block " << fraction(oneBlock) << ", two-block " << fraction(twoBlocks) << ", shifted "
                                                                  << fraction(shifted));
}

} // namespace o2::gtrk
