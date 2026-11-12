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

/// \file Iteration.h
/// \brief What one tracking iteration is configured with, and everything its
///        kernels read that follows from it.
///
/// Iterations run one after another over the same timeframe, each with looser
/// cuts over what the previous ones left unclaimed. Three lifetimes are kept
/// apart (settings, parameters, timeframe context):
///
///  - `IterationSettings`, what a person configures, per run;
///  - `IterationParameters<Det>`, what follows from the settings, the geometry
///    and the material, built once at initialisation for every iteration;
///  - the timeframe context -- measurements, vertex seeds, used mask -- which
///    is not here.
///
/// Nothing in `IterationParameters` depends on the timeframe, so a device
/// backend uploads the whole array once. It holds C arrays rather than
/// pointers for the same reason.

#ifndef GTRK_TRACKING_ITERATION_H_
#define GTRK_TRACKING_ITERATION_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/TrackFitter.h"
#include "GenericTracking/Tracking/LinkHoles.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/SearchWindow.h"

#include "GPUCommonMath.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>

namespace o2::gtrk
{

/// What a person configures for one iteration. Host-side.
struct IterationSettings {
  /// The window physics. The minimum pT, nSigma and a displaced iteration's
  /// bounds are what an iteration sets; the field, the mass, the uniform
  /// material that sizes the grid (`materialX2X0`) and the resolution are the
  /// run's.
  SearchWindowParams windows;
  VertexSource vertexSource = VertexSource::kBeamLine;
  /// The shortest road this iteration fits, in clusters: at least 4 for a
  /// road, two triplets (`minRoadStartLevel`);
  /// 3 for a pass that stops at triplets, the vertexing pass. It also names
  /// the links the iteration takes: only those on a path of this many layers
  /// (`takenLinks`).
  int minClusters = 4;
  /// The kinks a fit pass may take; 0 for none.
  int maxKinks = 1;
  FitReference fitReference = FitReference::kFilteredState;
  /// The triplet and neighbour cuts' kink allowance: the factor on theta0 at pT_min a triplet's middle
  /// point may have scattered by, beyond the nominal cuts, to be kept as a
  /// kinked triplet; 1 for none.
  double kinkSigma = 1.;
  /// The graph's own kink allowance, the triplet and neighbour cuts' alone:
  /// above 0 it replaces `kinkSigma`
  /// there, which then sets only the road fit's and the follow's; 0, the
  /// default, `kinkSigma` everywhere.
  double graphKinkSigma = 0.;
  /// Whether this iteration makes tracklets that leave a hole: whose
  /// segment crosses another layer's envelope (LinkHoles.h). Off, it takes
  /// no link that skips a layer (`LayerLink::skipsLayer`) and on the others no pair
  /// that crosses one, so a track beyond a shorter layer's end is still made.
  bool skipLinks = true;
  /// The layers this iteration's links join, a bit per layer; 0 for every
  /// layer. A link is taken only when both its layers are in the mask: the
  /// vertexing pass takes the three innermost.
  std::uint64_t layerMask = 0U;
  /// The layers a road may end on, a bit per layer; 0 for any. With a mask a
  /// road starts at every triplet ending on one of them, so a road's prefix
  /// is a road of its own when it ends there (O2's ITS `StartLayerMask`, whose
  /// roads grow inward from their outermost cell); without one only at the
  /// triplets nothing extends.
  std::uint64_t roadEndLayers = 0U;
  /// The layers a road may start on, its innermost measurement's, a bit per
  /// layer; 0 for any. The roads counted are the
  /// same chains; one whose innermost measurement lies on another layer is
  /// none, as one starting outside the slice's core.
  std::uint64_t roadStartLayers = 0U;
  /// The lowest pT a track of each cluster count may have [GeV], by the
  /// count; 0 for none. O2's ITS `TrackingParameters::MinPt`, by track
  /// length. The final track's: where the iteration
  /// follows, after the follow, on the track's count then; else on the road.
  static constexpr std::size_t kClusterSlots = 16U;
  std::array<double, kClusterSlots> minPtGeVByClusters{};
  /// How many of a track's innermost clusters this iteration's selection and
  /// follow let it share; 0, strict: the
  /// follow does not go inward otherwise.
  std::uint32_t sharedInnermost = 0U;
  /// Whether this iteration follows its tracks after their refit, outward
  /// and then inward, as a CKF does (`FollowStage`):
  /// the layers every link of the layout reaches, in the order the track
  /// crosses them, the first with a row within the gate taken; a layer
  /// without one is a hole where the track crosses one of its sensors.
  bool follow = true;
  /// The holes a followed track may cross in each direction; past them the
  /// direction ends.
  int followHoles = 0;
  /// Whether the road's outermost row is decided again by the follow:
  /// the road fitted without it, the follow's
  /// first row takes its place where it is another; where the follow takes
  /// none, the end stays. Without `follow`, the end alone.
  bool followEnd = false;
  /// The innermost clusters a road the road fit refuses may be fitted again
  /// without, one more each time, down to `minClusters`; 0 for none.
  int innerTrim = 0;
  /// Which of a barrel plane's two crossings the follow's refit takes
  /// (`CrossingRule`).
  CrossingRule crossingRule = CrossingRule::kFirstAhead;
  /// Whether the follow's refit keeps the longest accepted prefix of the
  /// rows it would refuse (`follow::Refit`'s Prefix).
  bool followPrefix = false;
  /// Whether the follow takes the smallest chi2 where more than one row
  /// passes a step, instead of the lookahead's confirmed row (`FollowStage`'s Greedy).
  bool followGreedy = false;
  /// The unconfirmed drop's chi2 below which a direction's unconfirmed last row stays:
  /// the unconfirmed drop removes it only at or
  /// above; 0, the default, none stays.
  double followKeepChi2 = 0.;
  /// The end decisions on the last row's own update (`FollowStage`'s
  /// EndDensity), each 0, the
  /// default, off: inward, the sqrt|V| [cm^2] above which a nominal last row
  /// at the layout's end -- no inward link -- is refused; outward, the
  /// Gaussian density g [cm^-2] above which the unconfirmed drop keeps a nominal last row.
  double followEndMaxSqrtDet = 0.;
  double followKeepDensity = 0.;
  /// Whether the inward walk's last row, where a nominal step took it, is
  /// refused above the gate against that step's prediction updated with the
  /// beam line (`FollowStage`'s BeamEnd); requires the follow and the beam constraint. False, the
  /// default, off.
  bool followBeamEnd = false;
  /// The Delta chi2 below which the follow trims a track's last outward row
  /// where it is one of two rows or more passing the gate: the other's
  /// chi2 less its own (`follow::TrimAmbiguous`); 0, the default, off.
  double followTrimMargin = 0.;
  /// The clusters a track of this iteration must hold after its follow; 0
  /// for none. A shorter one is
  /// dropped and its rows given back.
  int minClustersAfterFollow = 0;
  /// Whether each road is updated with the beam line on its innermost surface
  /// and refused above the update's gate;
  /// the constrained state starts the follow inward, and the refits leave
  /// the beam line out.
  bool beamConstraint = false;
  /// Whether every earlier iteration's track is followed again with this
  /// iteration's settings, before this iteration makes its seeds, the rows
  /// it takes its own;
  /// requires the follow. A length after the follow applies to this
  /// iteration's own tracks only.
  bool followEarlier = false;
  /// Whether this iteration takes the previous one's graph instead of
  /// building its own (one graph for both rounds):
  /// the previous iteration's slices build the graph once with this
  /// iteration's kink allowance; the roads of nominal triplets the nominal
  /// fit accepts are the previous iteration's, every other road is fitted
  /// with this iteration's settings and is this iteration's. Requires the
  /// previous iteration's graph settings (the windows, the vertex source,
  /// minClusters, skipLinks, layerMask, roadEndLayers, roadStartLayers and the pT floors by
  /// cluster count), and no kink allowance there and one here.
  bool sharedGraph = false;
  /// Whether the iteration tests its pairs and triplets against their
  /// measurements' cluster shapes: each measurement's slope interval
  /// (`MeasurementT::slopeLow`, `slopeHigh`) must hold the candidate's local
  /// slope there (ShapeCut.h).
  /// False, the default, off.
  bool shapeTest = false;

  [[nodiscard]] bool operator==(const IterationSettings&) const = default;
};

/// One link's geometry as the triplet stage bounds it:
/// the shortest 3D distance between measurements on its two layers, and the
/// longest transverse chord a tracklet on it can have.
struct TripletLinkBounds {
  Scalar lengthMin = 0.F; ///< [cm], rounded down
  Scalar chordMax = 0.F;  ///< [cm], rounded up
};

/// The bounds of the link \a from -> \a to, whose tracklets lie within
/// \a phiHalfWidth of their source's azimuth.
///
/// Two points at radii `r_1`, `r_2` and heights `z_1`, `z_2` are at least
/// `hypot(r_2 - r_1, z_2 - z_1)` apart, so the gaps between the layers' radial
/// and z extents bound the length from below. Their transverse chord is
/// `d^2 = (r_2 - r_1)^2 + 2 r_1 r_2 (1 - cos(dphi))`, bounded above by the
/// largest radial difference and the window's azimuth.
template <DetectorLayout Det>
[[nodiscard]] inline TripletLinkBounds tripletLinkBounds(LayerID from, LayerID to, double phiHalfWidth)
{
  const double fromLow = radiusMin<Det>(from);
  const double fromHigh = radiusMax<Det>(from);
  const double toLow = radiusMin<Det>(to);
  const double toHigh = radiusMax<Det>(to);
  const double gapR = std::max({toLow - fromHigh, fromLow - toHigh, 0.});
  const double gapZ = std::max({zMin<Det>(to) - zMax<Det>(from), zMin<Det>(from) - zMax<Det>(to), 0.});
  const double lengthMin = std::hypot(gapR, gapZ);
  auto narrowedLength = static_cast<Scalar>(lengthMin);
  if (static_cast<double>(narrowedLength) > lengthMin) {
    narrowedLength = std::nextafter(narrowedLength, 0.F);
  }
  const double radialMax = std::max(std::abs(toHigh - fromLow), std::abs(fromHigh - toLow));
  // Half a turn at most: a prompt window stays within a quarter, a displaced
  // one may not.
  const double azimuth = std::min(phiHalfWidth, std::numbers::pi_v<double>);
  const double chordMax = std::sqrt((radialMax * radialMax) + (2. * fromHigh * toHigh * (1. - std::cos(azimuth))));
  return {.lengthMin = narrowedLength, .chordMax = conservativeScalarBound(chordMax)};
}

/// Everything one iteration's kernels read.
template <DetectorLayout Det>
// NOLINTNEXTLINE(clang-analyzer-optin.performance.Padding): the fields follow the configuration; their order is the designated initialisers'
struct IterationParameters {
  static constexpr std::size_t kLinkCount = nLinks<Det>;
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);
  static_assert(kLinkCount <= std::size_t{std::numeric_limits<std::uint16_t>::max()} + 1U);

  PairCut pairCut;
  VertexSource vertexSource = VertexSource::kBeamLine;
  ScatteringTable<Det> scattering;
  /// One per entry of `Det::layerLinks`, in its order. They choose the grid
  /// cells a source's query visits and are derived about the beam spot; the
  /// exact pair cut is then evaluated against each vertex seed.
  LayerPairWindow windows[kLinkCount]{};
  /// Float copies of the host geometry's source intervals, narrowed once when
  /// the parameters cross into tracking code.
  Scalar sourceMin[kLinkCount]{};
  Scalar sourceMax[kLinkCount]{};
  TripletLinkBounds tripletBounds[kLinkCount]{};
  /// Whether a passive scatterer can lie between the link's layers: the pair
  /// test then takes `betweenScattering`.
  bool passiveBetween[kLinkCount]{};
  /// The largest `betweenScattering` a pair on the link can take: the link's
  /// `betweenKinkBound` squared, rounded up; 0 without a passive between. The
  /// tracklet stage tests a pair at it first and computes the pair's own term
  /// only for the pairs it keeps: the test only widens with the term, so the
  /// same pairs pass.
  Scalar betweenBound[kLinkCount]{};
  /// Device-facing copy of the detector's link topology. Runtime indexing of
  /// the detector's constexpr host tables is not valid in device code.
  std::uint32_t outLinkOffsets[kLayerCount + 1U]{};
  std::uint16_t outLinkIndices[kLinkCount]{};
  std::uint32_t inLinkOffsets[kLayerCount + 1U]{};
  std::uint16_t inLinkIndices[kLinkCount]{};
  LayerID linkSources[kLinkCount]{};
  LayerID linkTargets[kLinkCount]{};
  /// Per link, the layers a pair on it can cross (`linkHoleLayers`): in an
  /// iteration that takes no skips a tracklet crossing one is not made;
  /// 0 in one that takes them.
  std::uint64_t holeLayers[kLinkCount]{};
  /// The passives a pair on the link can cross, a bit each (`passivesBetween`):
  /// the ones `betweenScattering` tries.
  std::uint64_t betweenPassives[kLinkCount]{};
  /// The barrel layers, a bit per layer, and every layer's span, for that test.
  std::uint64_t barrelLayers = 0U;
  LayerSpan layerSpans[kLayerCount]{};
  /// The beam spot's z and its n-sigma reach [cm]: the range the tracking is
  /// to be efficient in, and the only one a vertex-list iteration uses.
  Scalar vertexZCentre = 0.F;
  Scalar vertexZReach = 0.F;
  /// A displaced iteration's impact bound and its closest approach's z bound
  /// about `vertexZCentre` [cm] (`SearchWindowParams::impactMax`); 0 for a prompt one.
  Scalar impactMax = 0.F;
  Scalar impactZMax = 0.F;
  /// The adjacent-triplet fit's chi2 bound: three degrees of freedom at the
  /// confidence of `nSigma`, from `chiSquareQuantile`.
  Scalar neighbourChi2Max = 0.F;
  /// A proposed cluster's update chi2 bound, two coordinates.
  Scalar clusterChi2Max = 0.F;
  /// The road fit's chi2 bound by degrees of freedom: `2 N - 5` for `N`
  /// clusters, two fewer per kink, at the confidence of `nSigma`
  /// (`roadChi2Cut`); a box gives the time none.
  Scalar roadChi2MaxByNdf[(2 * maxChainLayers<Det>)-4]{};
  /// The kinks a fit pass may take.
  std::uint8_t maxKinks = 0U;
  FitReference fitReference = FitReference::kFilteredState;
  /// `IterationSettings::followEnd`, `innerTrim`
  /// and `crossingRule`.
  bool followEnd = false;
  std::uint8_t innerTrim = 0U;
  CrossingRule crossingRule = CrossingRule::kFirstAhead;
  /// The update chi2 beyond which a step is a kink: two degrees of freedom at
  /// `nSigma`.
  Scalar kinkChi2Min = 0.F;
  /// The kink allowance's square, the factor on the middle layer's scattering
  /// variance a kinked triplet is tested with; 1 for none (`IterationSettings`).
  Scalar kinkScale2 = 1.F;
  /// The triplet and neighbour cuts' kink allowance's square
  /// (`IterationSettings::graphKinkSigma`); `kinkScale2` without one of its own.
  Scalar graphKinkScale2 = 1.F;
  /// The lowest level a road can start at: a road from level `L` has `L + 2`
  /// clusters, so `minClusters - 2`, and at least 2.
  std::uint32_t minRoadStartLevel = 2U;
  /// `IterationSettings::minClusters`: the clusters a chain must reach to be
  /// used, a road's or the vertexing pass's triplet's. The stages drop the
  /// tracklets and triplets on no path that long (TrackletPaths.h).
  std::uint32_t minClusters = 4U;
  /// The layers a road may end on; 0 for any (`IterationSettings`).
  std::uint64_t roadEndLayers = 0U;
  /// The layers a road may start on; 0 for any (`IterationSettings`).
  std::uint64_t roadStartLayers = 0U;
  /// The lowest pT of a track by its cluster count [GeV]; 0 for none
  /// (`IterationSettings::minPtGeVByClusters`).
  Scalar minPtByClusters[maxChainLayers<Det> + 1U]{};
  /// `IterationSettings`' shareable innermost clusters and follow switch.
  std::uint32_t sharedInnermost = 0U;
  bool follow = true;
  /// `IterationSettings::sharedGraph`.
  bool sharedGraph = false;
  /// `IterationSettings::beamConstraint`,
  /// with its terms below.
  bool beamConstraint = false;
  /// `IterationSettings::followEarlier`.
  bool followEarlier = false;
  /// `IterationSettings::followPrefix` and `followGreedy`.
  bool followPrefix = false;
  bool followGreedy = false;
  /// `IterationSettings::followKeepChi2`.
  Scalar followKeepChi2 = 0.F;
  /// The resolution term of the follow's window, `positionSigma` squared,
  /// as the pair windows take it [cm^2].
  Scalar followPositionVariance = 0.F;
  /// `IterationSettings::minClustersAfterFollow`; 0 for none.
  std::uint32_t minClustersAfterFollow = 0U;
  /// The beam-line update's terms (`beamConstraint`): the beam line's
  /// position in the fit's frame [cm] and its transverse width squared
  /// [cm^2], and the update's gate, `chiSquareQuantile(1, nSigma)`.
  Scalar beamX = 0.F;
  Scalar beamY = 0.F;
  Scalar beamVarTransverse = 0.F;
  Scalar beamChi2Max = 0.F;
  /// `IterationSettings::followHoles`, and the layout's every link out of
  /// each layer and into it, which the follow walks outward and inward.
  std::uint32_t followHoles = 0U;
  std::uint32_t allOutLinkOffsets[kLayerCount + 1U]{};
  std::uint16_t allOutLinkIndices[kLinkCount]{};
  std::uint32_t allInLinkOffsets[kLayerCount + 1U]{};
  std::uint16_t allInLinkIndices[kLinkCount]{};
  /// `IterationSettings::followEndMaxSqrtDet` and `followKeepDensity`.
  Scalar followEndMaxSqrtDet = 0.F;
  Scalar followKeepDensity = 0.F;
  /// `IterationSettings::followTrimMargin`.
  Scalar followTrimMargin = 0.F;
  /// `IterationSettings::followBeamEnd`, set only with the follow and the beam
  /// constraint; read on the host only to choose the follow's instantiation
  /// (`FollowStage`'s BeamEnd).
  bool followBeamEnd = false;
  /// `IterationSettings::shapeTest`, read on the host only to choose the
  /// tracklet and triplet stages' instantiations.
  bool shapeTest = false;
};

/// Whether an iteration uses \a vertex for a source with time box \a source:
/// compatible in time, and -- from a vertex list -- lying
/// with its own n sigma inside the beam spot's reach. Every window the query
/// can take is derived over that reach, so the tracklets do not depend on
/// which window a source takes; a vertex beyond it is outside the range the
/// tracking is to be efficient in.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool vertexUsable(const IterationParameters<Det>& parameters, const TimeEstPS& source, const VertexSeed& vertex)
{
  if (parameters.vertexSource == VertexSource::kVertexList && o2::gpu::CAMath::Abs(vertex.z - parameters.vertexZCentre) + (parameters.pairCut.nSigma * o2::gpu::CAMath::Sqrt(vertex.varZ)) > parameters.vertexZReach) {
    return false;
  }
  return vertexTimeCompatible(source, TimeEstPS{.start = vertex.timePS, .width = vertex.widthPS});
}

/// Per-link host bounds for the intervals centred on individual vertices.
/// They depend on one timeframe's vertex list and are kept separate from the
/// iteration parameters uploaded once per run.
template <DetectorLayout Det>
struct VertexWindowParameters {
  static constexpr std::size_t kLinkCount = nLinks<Det>;
  VertexPairWindow windows[kLinkCount]{};
};

/// The bounds for one timeframe's vertex list, over the reach `vertexUsable`
/// admits and with the pair cut's scattering table. \a settings are the ones
/// \a parameters were derived from.
template <DetectorLayout Det>
[[nodiscard]] inline VertexWindowParameters<Det> vertexWindowParameters(const IterationParameters<Det>& parameters, const SearchWindowParams& settings, std::span<const VertexSeed> vertices)
{
  double maxSigmaZ = 0.;
  double maxSigmaTransverse = 0.;
  for (const auto& vertex : vertices) {
    maxSigmaZ = std::max(maxSigmaZ, std::sqrt(static_cast<double>(vertex.varZ)));
    maxSigmaTransverse = std::max(maxSigmaTransverse, std::sqrt(static_cast<double>(vertex.varTransverse)));
  }
  const auto centre = static_cast<double>(parameters.vertexZCentre);
  const auto reach = static_cast<double>(parameters.vertexZReach);
  const Interval vertexZ{.low = centre - reach, .high = centre + reach};
  const VertexConstraint envelope{.z = centre, .sigmaZ = settings.nSigma > 0. ? reach / settings.nSigma : 0., .sigmaTransverse = maxSigmaTransverse};

  VertexWindowParameters<Det> result;
  for (std::size_t link = 0; link < VertexWindowParameters<Det>::kLinkCount; ++link) {
    const auto& pair = Det::layerLinks[link];
    const double kink = sourceKinkBound<Det>(parameters.scattering, settings, pair.from, envelope);
    const double between = betweenKinkBound<Det>(parameters.scattering, settings, pair.from, pair.to, envelope);
    result.windows[link] = vertexPairWindow<Det>(pair.from, pair.to, settings, vertexZ, maxSigmaZ, kink, between);
  }
  return result;
}

/// The links an iteration with \a settings takes: those its skip rule and
/// its mask allow, and of them only the ones on a path of `minClusters`
/// layers or more through the allowed links. No road through any other link
/// reaches `minClusters`, so no tracklet, triplet or neighbour is made on it.
/// The one predicate for the stages' link lists and for the grid's
/// sizing (GridSizing.h), so that the two cannot disagree.
template <DetectorLayout Det>
[[nodiscard]] inline std::array<bool, nLinks<Det>> takenLinks(const IterationSettings& settings)
{
  const auto inMask = [&](LayerID layer) { return settings.layerMask == 0U || ((settings.layerMask >> layer) & 1U) != 0U; };
  std::array<bool, nLinks<Det>> allowed{};
  for (std::size_t link = 0; link < nLinks<Det>; ++link) {
    const auto& pair = Det::layerLinks[link];
    allowed[link] = (settings.skipLinks || !pair.skipsLayer) && inMask(pair.from) && inMask(pair.to);
  }
  // The most layers a path of allowed links visits ending at each layer, and
  // leaving it, the layer itself counted in both.
  constexpr auto kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);
  std::array<int, kLayers> ending{};
  std::array<int, kLayers> leaving{};
  ending.fill(1);
  leaving.fill(1);
  for (const auto layer : layerOrder<Det>) {
    for (auto i = outLinks<Det>.begin(layer); i < outLinks<Det>.end(layer); ++i) {
      const auto link = outLinks<Det>.links[i];
      if (allowed[link]) {
        auto& to = ending[Det::layerLinks[link].to];
        to = std::max(to, ending[layer] + 1);
      }
    }
  }
  for (auto order = kLayers; order-- > 0;) {
    const auto layer = layerOrder<Det>[order];
    for (auto i = outLinks<Det>.begin(layer); i < outLinks<Det>.end(layer); ++i) {
      const auto link = outLinks<Det>.links[i];
      if (allowed[link]) {
        leaving[layer] = std::max(leaving[layer], leaving[Det::layerLinks[link].to] + 1);
      }
    }
  }
  std::array<bool, nLinks<Det>> taken{};
  for (std::size_t link = 0; link < nLinks<Det>; ++link) {
    const auto& pair = Det::layerLinks[link];
    taken[link] = allowed[link] && ending[pair.from] + leaving[pair.to] >= settings.minClusters;
  }
  return taken;
}

/// Derive one iteration's parameters. \a scattering is built once per run
/// from the per-layer material and is the same for every iteration. \a alike,
/// where given, is an earlier iteration's parameters from the same beam spot
/// and table with equal `settings.windows`: the links' windows and bounds
/// depend on those alone and are taken from it instead of derived again.
template <DetectorLayout Det>
[[nodiscard]] inline IterationParameters<Det> iterationParameters(const IterationSettings& settings, const VertexConstraint& beamSpot, const ScatteringTable<Det>& scattering, const IterationParameters<Det>* alike = nullptr)
{
  if (settings.windows.impactMax > 0. && (settings.vertexSource != VertexSource::kBeamLine || !(settings.windows.impactZMax > 0.))) {
    throw std::invalid_argument("a displaced iteration (impactMax > 0) needs the beam line as its vertex and a z bound (impactZMax > 0)");
  }
  IterationParameters<Det> parameters;
  parameters.pairCut = pairCut(settings.windows);
  parameters.vertexSource = settings.vertexSource;
  parameters.scattering = scattering;
  // The link lists the stages walk, without the links this iteration does
  // not take: no tracklet is made on them, and no triplet crosses them.
  const auto keep = [&](const std::array<bool, nLinks<Det>>& taken, const LinksByLayer<Det>& all, std::uint32_t* offsets, std::uint16_t* indices) {
    std::uint32_t kept = 0U;
    for (std::size_t layer = 0; layer < IterationParameters<Det>::kLayerCount; ++layer) {
      offsets[layer] = kept;
      for (auto i = all.begin(static_cast<LayerID>(layer)); i < all.end(static_cast<LayerID>(layer)); ++i) {
        if (taken[all.links[i]]) {
          indices[kept++] = static_cast<std::uint16_t>(all.links[i]);
        }
      }
    }
    offsets[IterationParameters<Det>::kLayerCount] = kept;
  };
  const auto taken = takenLinks<Det>(settings);
  keep(taken, outLinks<Det>, parameters.outLinkOffsets, parameters.outLinkIndices);
  keep(taken, inLinks<Det>, parameters.inLinkOffsets, parameters.inLinkIndices);
  // The follow walks every link of the layout, with the window of its own
  // prediction.
  parameters.follow = settings.follow;
  parameters.followEarlier = settings.follow && settings.followEarlier;
  parameters.followPrefix = settings.follow && settings.followPrefix;
  parameters.followGreedy = settings.follow && settings.followGreedy;
  parameters.followKeepChi2 = static_cast<Scalar>(std::max(settings.followKeepChi2, 0.));
  parameters.followEndMaxSqrtDet = static_cast<Scalar>(std::max(settings.followEndMaxSqrtDet, 0.));
  parameters.followKeepDensity = static_cast<Scalar>(std::max(settings.followKeepDensity, 0.));
  parameters.followTrimMargin = settings.follow ? static_cast<Scalar>(std::max(settings.followTrimMargin, 0.)) : 0.F;
  parameters.followBeamEnd = settings.follow && settings.beamConstraint && settings.followBeamEnd;
  parameters.shapeTest = settings.shapeTest;
  parameters.followHoles = static_cast<std::uint32_t>(std::max(settings.followHoles, 0));
  parameters.followEnd = settings.followEnd;
  if (settings.follow || settings.followEnd) {
    std::array<bool, nLinks<Det>> every{};
    every.fill(true);
    keep(every, outLinks<Det>, parameters.allOutLinkOffsets, parameters.allOutLinkIndices);
    keep(every, inLinks<Det>, parameters.allInLinkOffsets, parameters.allInLinkIndices);
    parameters.followPositionVariance = static_cast<Scalar>(settings.windows.positionSigma * settings.windows.positionSigma);
  }
  parameters.minClustersAfterFollow = static_cast<std::uint32_t>(std::max(settings.minClustersAfterFollow, 0));
  for (std::size_t layer = 0; layer < IterationParameters<Det>::kLayerCount; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    parameters.layerSpans[layer] = layerSpan<Det>(id);
  }
  parameters.barrelLayers = barrelLayerBits<Det>();
  // A displaced iteration's tracks need not come from the beam spot the
  // links' intervals were traced from: the whole layer.
  const bool displaced = settings.windows.impactMax > 0.;
  for (std::size_t link = 0; link < IterationParameters<Det>::kLinkCount; ++link) {
    const auto& pair = Det::layerLinks[link];
    parameters.linkSources[link] = pair.from;
    parameters.linkTargets[link] = pair.to;
    if (alike != nullptr) {
      parameters.windows[link] = alike->windows[link];
      parameters.passiveBetween[link] = alike->passiveBetween[link];
      parameters.betweenBound[link] = alike->betweenBound[link];
      parameters.tripletBounds[link] = alike->tripletBounds[link];
    } else {
      parameters.windows[link] = layerPairWindow<Det>(pair.from, pair.to, settings.windows, beamSpot, scattering);
      const double between = betweenKinkBound<Det>(scattering, settings.windows, pair.from, pair.to, beamSpot);
      parameters.passiveBetween[link] = between > 0.;
      parameters.betweenBound[link] = conservativeScalarBound(between * between);
      parameters.tripletBounds[link] = tripletLinkBounds<Det>(pair.from, pair.to, static_cast<double>(parameters.windows[link].phiHalfWidth));
    }
    parameters.betweenPassives[link] = passivesBetween<Det>(scattering, pair.from, pair.to);
    parameters.holeLayers[link] = settings.skipLinks ? 0U : linkHoleLayers<Det>(pair.from, pair.to);
    parameters.sourceMin[link] = displaced ? -std::numeric_limits<Scalar>::infinity() : static_cast<Scalar>(pair.sourceMin);
    parameters.sourceMax[link] = displaced ? std::numeric_limits<Scalar>::infinity() : static_cast<Scalar>(pair.sourceMax);
  }
  parameters.vertexZCentre = static_cast<Scalar>(beamSpot.z);
  parameters.vertexZReach = static_cast<Scalar>(settings.windows.nSigma * beamSpot.sigmaZ);
  parameters.impactMax = static_cast<Scalar>(settings.windows.impactMax);
  parameters.impactZMax = static_cast<Scalar>(settings.windows.impactZMax);
  parameters.neighbourChi2Max = static_cast<Scalar>(chiSquareQuantile(3, settings.windows.nSigma));
  parameters.clusterChi2Max = static_cast<Scalar>(chiSquareQuantile(2, settings.windows.nSigma));
  parameters.beamConstraint = settings.beamConstraint;
  parameters.sharedGraph = settings.sharedGraph;
  parameters.beamX = static_cast<Scalar>(beamSpot.x);
  parameters.beamY = static_cast<Scalar>(beamSpot.y);
  parameters.beamVarTransverse = static_cast<Scalar>(beamSpot.sigmaTransverse * beamSpot.sigmaTransverse);
  parameters.beamChi2Max = static_cast<Scalar>(chiSquareQuantile(1, settings.windows.nSigma));
  for (int ndf = 1; ndf <= (2 * static_cast<int>(maxChainLayers<Det>)) - kNBoundParams; ++ndf) {
    parameters.roadChi2MaxByNdf[ndf] = static_cast<Scalar>(chiSquareQuantile(ndf, settings.windows.nSigma));
  }
  parameters.maxKinks = static_cast<std::uint8_t>(std::max(settings.maxKinks, 0));
  parameters.fitReference = settings.fitReference;
  parameters.innerTrim = static_cast<std::uint8_t>(std::clamp(settings.innerTrim, 0, static_cast<int>(maxChainLayers<Det>)));
  parameters.crossingRule = settings.crossingRule;
  parameters.kinkChi2Min = static_cast<Scalar>(chiSquareQuantile(2, settings.windows.nSigma));
  {
    const double sigma = std::max(settings.kinkSigma, 1.);
    parameters.kinkScale2 = static_cast<Scalar>(sigma * sigma);
    const double graphSigma = settings.graphKinkSigma > 0. ? std::max(settings.graphKinkSigma, 1.) : sigma;
    parameters.graphKinkScale2 = static_cast<Scalar>(graphSigma * graphSigma);
  }
  parameters.minRoadStartLevel = static_cast<std::uint32_t>(std::max(settings.minClusters, 4) - 2);
  parameters.minClusters = static_cast<std::uint32_t>(std::max(settings.minClusters, 0));
  parameters.roadEndLayers = settings.roadEndLayers;
  parameters.roadStartLayers = settings.roadStartLayers;
  constexpr std::size_t kCounts = std::min(maxChainLayers<Det> + 1U, IterationSettings::kClusterSlots);
  for (std::size_t count = 0; count < kCounts; ++count) {
    parameters.minPtByClusters[count] = static_cast<Scalar>(settings.minPtGeVByClusters[count]);
  }
  parameters.sharedInnermost = settings.sharedInnermost;
  return parameters;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_ITERATION_H_
