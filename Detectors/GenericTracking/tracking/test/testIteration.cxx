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

/// \file testIteration.cxx
/// \brief One iteration's parameters, derived from its settings.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Iteration
#include <boost/test/unit_test.hpp>
#include "CommonConstants/LHCConstants.h"
#include "GridTestSupport.h"

#include "GenericTracking/Study/Validation.h"
#include "GenericTracking/Base/LayerLinks.h"
#include "ToyBarrel4Disc8.h"
#include "ToyBarrel8.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/SlicePlan.h"

#include "UniformMaterial.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <vector>

namespace o2::gtrk
{

namespace
{

using detectors::ToyBarrel4Disc8;

/// A pair's transverse and longitudinal terms together, as the tracklet
/// stage takes them in two halves (`pairCompatible`).
PairTerms pairTerms(const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target)
{
  const auto pair = pairChord(source, target);
  auto terms = transversePairTerms(cut, vertex, source, sourceScatter, target, pair, 0.F);
  addLongitudinalPairTerms(terms, cut, vertex, source, sourceScatter, target, pair, 0.F);
  return terms;
}

/// An iteration at \a minPtGeV taking every link: chains of two clusters.
IterationParameters<ToyBarrel4Disc8> parametersAt(double minPtGeV)
{
  IterationSettings settings;
  settings.minClusters = 2;
  settings.windows.minPtGeV = minPtGeV;
  settings.vertexSource = VertexSource::kVertexList;
  return iterationParameters<ToyBarrel4Disc8>(settings, VertexConstraint{}, scatteringTable<ToyBarrel4Disc8>(uniformLayerMaterial<ToyBarrel4Disc8>(0.01)));
}

PairPoint helixPoint(double curvature, double arc, double cotTheta, double vertexZ)
{
  const double turn = curvature * arc;
  const double x = std::sin(turn) / curvature;
  const double y = (1. - std::cos(turn)) / curvature;
  return {
    .x = static_cast<Scalar>(x),
    .y = static_cast<Scalar>(y),
    .z = static_cast<Scalar>(vertexZ + (cotTheta * arc)),
    .r = static_cast<Scalar>(std::hypot(x, y)),
    .varZ = 1.e-4F,
    .varR = 1.e-6F,
    .varRPhi = 1.e-6F,
  };
}

/// The vertex range the per-vertex intervals are checked over [cm].
constexpr double kVertexRange = 5.;

/// \a steep puts a barrel source at the end of its stave farthest from the
/// vertex, where the path factor of the scattering is largest.
void verifyVertexInterval(double minPtGeV, Scalar vertexZ, LayerID from, LayerID to, bool steep)
{
  SearchWindowParams params;
  params.minPtGeV = minPtGeV;
  const auto table = scatteringTable<ToyBarrel4Disc8>(uniformLayerMaterial<ToyBarrel4Disc8>(params.materialX2X0));
  const auto cut = pairCut(params);
  const VertexSeed vertex{.z = vertexZ, .varZ = 1.e-4F, .varTransverse = 4.e-6F};
  const auto pair = layerPairWindow<ToyBarrel4Disc8>(from, to, params, VertexConstraint{}, table);
  const auto bounds = vertexPairWindow<ToyBarrel4Disc8>(from, to, params, {.low = -kVertexRange, .high = kVertexRange}, 0.01, sourceKinkBound<ToyBarrel4Disc8>(table, params, from, VertexConstraint{}));

  const bool sourceBarrel = Detector<ToyBarrel4Disc8>::isBarrelLayer(from);
  const bool targetBarrel = Detector<ToyBarrel4Disc8>::isBarrelLayer(to);
  const auto sourceR = static_cast<Scalar>(sourceBarrel ? layerDistance<ToyBarrel4Disc8>(from) : 0.5 * (secondMin<ToyBarrel4Disc8>(from) + secondMax<ToyBarrel4Disc8>(from)));
  auto sourceZ = static_cast<Scalar>(sourceBarrel ? static_cast<double>(vertexZ) + (0.5 * static_cast<double>(sourceR)) : discZ<ToyBarrel4Disc8>(from));
  if (sourceBarrel && steep) {
    sourceZ = static_cast<Scalar>(vertexZ < 0.F ? zMax<ToyBarrel4Disc8>(from) : zMin<ToyBarrel4Disc8>(from));
  }
  if (sourceBarrel && !targetBarrel) {
    const double targetR = 0.5 * (secondMin<ToyBarrel4Disc8>(to) + secondMax<ToyBarrel4Disc8>(to));
    const auto vertexValue = static_cast<double>(vertexZ);
    sourceZ = static_cast<Scalar>(vertexValue + (static_cast<double>(sourceR) * (discZ<ToyBarrel4Disc8>(to) - vertexValue) / targetR));
  }
  const auto variance = static_cast<Scalar>(params.positionSigma * params.positionSigma);
  const PairPoint source{.x = sourceR, .z = sourceZ, .r = sourceR, .varZ = variance, .varR = variance, .varRPhi = variance};
  const Scalar scattering = sourceScattering<ToyBarrel4Disc8>(table, cut, vertex, from, source);
  const auto interval = vertexSecondRange(pair, bounds, source.r, source.z, vertex.z, vertex.varZ, cut.nSigma);

  std::size_t accepted = 0U;
  constexpr int kSecondSteps = 400;
  constexpr int kPhiSteps = 40;
  const auto sourceRadius = static_cast<double>(sourceR);
  const auto sourcePositionZ = static_cast<double>(sourceZ);
  const auto vertexPositionZ = static_cast<double>(vertexZ);
  const double lineSecond = targetBarrel ? vertexPositionZ + ((sourcePositionZ - vertexPositionZ) * layerDistance<ToyBarrel4Disc8>(to) / sourceRadius) : sourceRadius * (discZ<ToyBarrel4Disc8>(to) - vertexPositionZ) / (sourcePositionZ - vertexPositionZ);
  for (int secondStep = 0; secondStep <= kSecondSteps + 1; ++secondStep) {
    const double second = secondStep == kSecondSteps + 1 ? lineSecond : secondMin<ToyBarrel4Disc8>(to) + ((secondMax<ToyBarrel4Disc8>(to) - secondMin<ToyBarrel4Disc8>(to)) * static_cast<double>(secondStep) / kSecondSteps);
    for (int phiStep = -kPhiSteps; phiStep <= kPhiSteps; ++phiStep) {
      const auto phi = static_cast<double>(pair.phiHalfWidth) * phiStep / kPhiSteps;
      const auto targetR = static_cast<Scalar>(targetBarrel ? layerDistance<ToyBarrel4Disc8>(to) : second);
      const auto targetZ = static_cast<Scalar>(targetBarrel ? second : discZ<ToyBarrel4Disc8>(to));
      const PairPoint target{
        .x = targetR * static_cast<Scalar>(std::cos(phi)),
        .y = targetR * static_cast<Scalar>(std::sin(phi)),
        .z = targetZ,
        .r = targetR,
        .varZ = variance,
        .varR = variance,
        .varRPhi = variance,
      };
      if (!pairCompatible(cut, vertex, source, scattering, target)) {
        continue;
      }
      ++accepted;
      BOOST_TEST_CONTEXT("pT " << minPtGeV << " vertex " << vertexZ << " link " << int(from) << " -> " << int(to) << " steep " << steep) { BOOST_TEST(secondRangeContains(interval, static_cast<Scalar>(second))); }
    }
  }
  BOOST_TEST(accepted > 0U, "pT " << minPtGeV << " vertex " << vertexZ << " link " << int(from) << " -> " << int(to));
}

} // namespace

BOOST_AUTO_TEST_CASE(WithoutSkipLinksAnIterationTakesOnlyTheLinksThatSkipNoLayer)
{
  IterationSettings settings;
  settings.minClusters = 2; // every link
  const auto all = iterationParameters<ToyBarrel4Disc8>(settings, VertexConstraint{}, scatteringTable<ToyBarrel4Disc8>(uniformLayerMaterial<ToyBarrel4Disc8>(0.01)));
  settings.skipLinks = false;
  const auto direct = iterationParameters<ToyBarrel4Disc8>(settings, VertexConstraint{}, scatteringTable<ToyBarrel4Disc8>(uniformLayerMaterial<ToyBarrel4Disc8>(0.01)));
  std::size_t skipping = 0U;
  for (std::size_t link = 0; link < nLinks<ToyBarrel4Disc8>; ++link) {
    BOOST_TEST_CONTEXT("link " << link)
    {
      BOOST_TEST(study::linkTaken(all, link));
      BOOST_TEST(study::linkTaken(direct, link) == !ToyBarrel4Disc8::layerLinks[link].skipsLayer);
      std::size_t in = 0U;
      const auto target = static_cast<std::size_t>(ToyBarrel4Disc8::layerLinks[link].to);
      for (auto at = direct.inLinkOffsets[target]; at < direct.inLinkOffsets[target + 1U]; ++at) {
        in += direct.inLinkIndices[at] == link ? 1U : 0U;
      }
      BOOST_TEST(in == (ToyBarrel4Disc8::layerLinks[link].skipsLayer ? 0U : 1U));
    }
    skipping += ToyBarrel4Disc8::layerLinks[link].skipsLayer ? 1U : 0U;
  }
  BOOST_TEST(skipping > 0U);
  BOOST_TEST(direct.outLinkOffsets[IterationParameters<ToyBarrel4Disc8>::kLayerCount] == nLinks<ToyBarrel4Disc8> - skipping);
}

BOOST_AUTO_TEST_CASE(PerVertexIntervalsContainEveryAcceptedPair)
{
  struct Case {
    LayerID from;
    LayerID to;
  };
  constexpr std::array<Case, 3> cases{{{.from = 2, .to = 3}, {.from = 8, .to = 9}, {.from = 0, .to = 8}}};
  for (const double pt : {0.1, 1.}) {
    for (const auto vertex : {static_cast<Scalar>(-kVertexRange), static_cast<Scalar>(kVertexRange)}) {
      for (const auto& pair : cases) {
        for (const bool steep : {false, true}) {
          verifyVertexInterval(pt, vertex, pair.from, pair.to, steep);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(PerVertexIntervalUsesThatVertexSigma)
{
  const SearchWindowParams params;
  constexpr LayerID kSource = 2;
  constexpr LayerID kTarget = 3;
  const auto pair = layerPairWindow<ToyBarrel4Disc8>(kSource, kTarget, params, VertexConstraint{});
  const auto bounds = vertexPairWindow<ToyBarrel4Disc8>(kSource, kTarget, params, {.low = -kVertexRange, .high = kVertexRange}, 0.1, uniformSourceKink(params));
  const auto sourceR = static_cast<Scalar>(layerDistance<ToyBarrel4Disc8>(kSource));
  const auto narrow = vertexSecondRange(pair, bounds, sourceR, 10.F, 0.F, 1.e-4F, static_cast<Scalar>(params.nSigma));
  const auto wide = vertexSecondRange(pair, bounds, sourceR, 10.F, 0.F, 1.e-2F, static_cast<Scalar>(params.nSigma));
  BOOST_TEST(wide.halfWidth > narrow.halfWidth);
}

BOOST_AUTO_TEST_CASE(ALooserIterationOpensEveryWindow)
{
  const auto tight = parametersAt(1.);
  const auto loose = parametersAt(0.1);
  BOOST_TEST(loose.pairCut.twoRMin < tight.pairCut.twoRMin);
  BOOST_TEST(loose.pairCut.invBetaMomentum2 > tight.pairCut.invBetaMomentum2);
  for (std::size_t link = 0; link < nLinks<ToyBarrel4Disc8>; ++link) {
    BOOST_TEST(loose.windows[link].phiHalfWidth >= tight.windows[link].phiHalfWidth, "link " << link);
    BOOST_TEST(loose.windows[link].secondHalfWidth >= tight.windows[link].secondHalfWidth, "link " << link);
  }
}

BOOST_AUTO_TEST_CASE(TheBeamLineIsCompatibleWithEveryTime)
{
  const auto seed = vertexSeed(VertexConstraint{.z = 1.5, .sigmaZ = 6., .sigmaTransverse = 0.002});
  const TimeEstPS box{.start = seed.timePS, .width = seed.widthPS};
  BOOST_TEST(vertexTimeCompatible(TimeEstPS{.start = 0U, .width = 1U}, box));
  BOOST_TEST(vertexTimeCompatible(TimeEstPS{.start = std::numeric_limits<std::uint32_t>::max() - 1U, .width = 1U}, box));
}

namespace
{

/// A grid of ToyBarrel8's eight layers in a chain, with \a frames per layer.
GridGeometry<detectors::ToyBarrel8> chainGrid(const std::array<FrameTiming, 8>& frames)
{
  GridParams params;
  std::ranges::copy(frames, std::begin(params.frames));
  const std::vector<GridLayerWindow> windows(8U);
  return {psOfBC(o2::constants::lhc::LHCMaxBunches), params, windows};
}

/// Bounds with every layer's widest box \a widthPS[layer] and the rows of
/// every layer exactly its frames.
HaloBounds<detectors::ToyBarrel8> framedBounds(const std::array<std::uint32_t, 8>& widthPS)
{
  HaloBounds<detectors::ToyBarrel8> bounds;
  bounds.widthPS = widthPS;
  bounds.rowsAreFrames.fill(true);
  return bounds;
}

} // namespace

/// Rows that are not their frames, the same widths on every layer: each link
/// reaches that width both ways, `R+ = R-` one width per link of the longest
/// chain, the gate `2R`, and the halo one link more for the targets of the
/// tracklets starting at the gate's edge.
BOOST_AUTO_TEST_CASE(RowsNotTheirFramesReachTheWiderBox)
{
  constexpr std::uint32_t kWidthPS = test::kToyBoxWidthPS;
  std::array<std::uint32_t, 8> widths{};
  widths.fill(kWidthPS);
  auto bounds = framedBounds(widths);
  bounds.rowsAreFrames.fill(false);
  std::array<FrameTiming, 8> frames{};
  frames.fill(test::kToyFrames);
  const auto reach = sliceReach(bounds, chainGrid(frames));
  const std::int64_t chain = static_cast<std::int64_t>(maxChainLayers<detectors::ToyBarrel8> - 1) * kWidthPS;
  BOOST_TEST(reach.pairForwardPS == kWidthPS);
  BOOST_TEST(reach.pairBackwardPS == kWidthPS);
  BOOST_TEST(reach.forwardPS == chain);
  BOOST_TEST(reach.backwardPS == chain);
  BOOST_TEST(reach.gatePS() == 2 * chain);
  BOOST_TEST(reach.haloPS() == (2 * chain) + kWidthPS);
}

/// Layers whose rows are their frames, and whose frames are the same, hold a
/// pair in one frame: only the link between the 22-BC and the 81-BC layers
/// reaches, 22 BC forward and 81 BC
/// backward.
BOOST_AUTO_TEST_CASE(LayersSharingFramesDoNotReach)
{
  constexpr FrameTiming kShort{.lengthBC = 22, .firstStartBC = 0};
  constexpr FrameTiming kLong{.lengthBC = 81, .firstStartBC = 0};
  const auto shortPS = static_cast<std::uint32_t>(psOfBC(22));
  const auto longPS = static_cast<std::uint32_t>(psOfBC(81));
  const auto grid = chainGrid({kShort, kShort, kShort, kLong, kLong, kLong, kLong, kLong});
  const std::array<std::uint32_t, 8> widths{shortPS, shortPS, shortPS, longPS, longPS, longPS, longPS, longPS};
  const auto reach = sliceReach(framedBounds(widths), grid);
  BOOST_TEST(reach.pairForwardPS == shortPS);
  BOOST_TEST(reach.pairBackwardPS == longPS);
  BOOST_TEST(reach.forwardPS == shortPS);
  BOOST_TEST(reach.backwardPS == longPS);
  BOOST_TEST(reach.gatePS() == std::int64_t{shortPS} + longPS);
  BOOST_TEST(reach.haloPS() == std::int64_t{shortPS} + (2 * std::int64_t{longPS}));

  // The same frames one BC later on layer 5: its two links reach.
  auto shifted = chainGrid({kShort, kShort, kShort, kLong, kLong, {.lengthBC = 81, .firstStartBC = 1}, kLong, kLong});
  const auto shiftedReach = sliceReach(framedBounds(widths), shifted);
  BOOST_TEST(shiftedReach.forwardPS == std::int64_t{shortPS} + (2 * std::int64_t{longPS}));
  BOOST_TEST(shiftedReach.backwardPS == 3 * std::int64_t{longPS});
}

/// A layer whose rows are not exactly its frames -- widened by a margin the
/// producer applied, whatever the configuration says now -- reaches on both
/// of its links.
BOOST_AUTO_TEST_CASE(RowsThatAreNotTheirFramesReach)
{
  constexpr FrameTiming kLong{.lengthBC = 81, .firstStartBC = 0};
  const auto widePS = static_cast<std::uint32_t>(psOfBC(81 + 60));
  std::array<FrameTiming, 8> frames{};
  frames.fill(kLong);
  std::array<std::uint32_t, 8> widths{};
  widths.fill(widePS);
  auto bounds = framedBounds(widths);
  BOOST_TEST(sliceReach(bounds, chainGrid(frames)).gatePS() == 0);
  bounds.rowsAreFrames[5] = false;
  const auto reach = sliceReach(bounds, chainGrid(frames));
  BOOST_TEST(reach.pairForwardPS == widePS);
  BOOST_TEST(reach.pairBackwardPS == widePS);
  BOOST_TEST(reach.forwardPS == 2 * std::int64_t{widePS});
  BOOST_TEST(reach.backwardPS == 2 * std::int64_t{widePS});
}

BOOST_AUTO_TEST_CASE(TheLongitudinalCutFollowsTheHelixArc)
{
  constexpr double kVertexZ = 0.7;
  constexpr double kCotTheta = 0.8;
  constexpr double kCurvature = 0.006;
  constexpr Scalar kNSigma = 3.F;
  const VertexSeed vertex{.z = static_cast<Scalar>(kVertexZ)};
  const PairCut cut{.twoRMin = 300.F, .nSigma = kNSigma, .maxCurvature = 2.F / 300.F};

  for (const double charge : {-1., 1.}) {
    const auto source = helixPoint(charge * kCurvature, 50., kCotTheta, kVertexZ);
    auto target = helixPoint(charge * kCurvature, 105., kCotTheta, kVertexZ);
    const auto terms = pairTerms(cut, vertex, source, 0.F, target);
    BOOST_TEST(transverseCompatible(cut, terms));
    BOOST_TEST(longitudinalCompatible(cut, terms));
    BOOST_TEST(pairCompatible(cut, vertex, source, 0.F, target));

    const Scalar dzSource = source.z - vertex.z;
    const Scalar dzTarget = target.z - vertex.z;
    const Scalar radialResidual = (source.r * dzTarget) - (target.r * dzSource);
    const Scalar radialVariance = (dzTarget * dzTarget * source.varR) + (target.r * target.r * source.varZ) + (dzSource * dzSource * target.varR) + (source.r * source.r * target.varZ);
    BOOST_TEST(radialResidual * radialResidual > kNSigma * kNSigma * radialVariance);

    target.z += 0.2F;
    BOOST_TEST(!pairCompatible(cut, vertex, source, 0.F, target));
  }
}

BOOST_AUTO_TEST_CASE(TheSecondCoordinateWindowCoversHelixPitch)
{
  SearchWindowParams params;
  params.minPtGeV = 1.;
  const VertexConstraint vertex{.sigmaZ = 0.01};
  constexpr LayerID kSourceLayer = 2;
  constexpr LayerID kTargetLayer = 3;
  const double sourceRadius = layerDistance<ToyBarrel4Disc8>(kSourceLayer);
  const double targetRadius = layerDistance<ToyBarrel4Disc8>(kTargetLayer);
  const double curvature = 1. / curvatureRadius(params.minPtGeV, params.bzKiloGauss);
  const double sourceArc = transverseArcAtRadius(curvature, sourceRadius);
  const double targetArc = transverseArcAtRadius(curvature, targetRadius);
  constexpr double kSourceZ = 12.;
  const double targetZ = kSourceZ * targetArc / sourceArc;
  const double sourcePhi = std::asin(0.5 * curvature * sourceRadius);
  const double targetPhi = std::asin(0.5 * curvature * targetRadius);

  const auto pair = layerPairWindow<ToyBarrel4Disc8>(kSourceLayer, kTargetLayer, params, vertex);
  const auto window = candidateWindow(pair, TimeEstPS{.start = 0U, .width = 1U}, static_cast<Scalar>(sourcePhi), static_cast<Scalar>(sourceRadius), static_cast<Scalar>(kSourceZ), 1U);
  BOOST_TEST(helixSecondSpread<ToyBarrel4Disc8>(kSourceLayer, kTargetLayer, params, vertex) > 0.);
  BOOST_TEST(spatiallyCompatible(window, static_cast<Scalar>(targetPhi), static_cast<Scalar>(targetZ)));
}

/// Two disc points at one radius, 30 cm apart in z: no track from the
/// vertex does that above 1 GeV, but with a radial error of centimetres the
/// chord is unresolved and the circle through them reaches curvature 2/r. The
/// cut holds the curvature to pT_min's, where the pitch test still has a
/// finite error and rejects the pair; a helix between the same discs passes.
BOOST_AUTO_TEST_CASE(AnUnresolvedDiscChordDoesNotPassTheLongitudinalCut)
{
  SearchWindowParams params;
  params.minPtGeV = 1.;
  const auto cut = pairCut(params);
  const VertexSeed vertex{.varZ = 1.e-4F};
  constexpr Scalar kRadius = 87.F;
  constexpr Scalar kVarR = 1.F;
  const PairPoint source{.x = kRadius, .z = 130.F, .r = kRadius, .varZ = 1.e-4F, .varR = kVarR, .varRPhi = 1.e-4F};
  const PairPoint sameRadius{.x = kRadius, .y = 0.01F, .z = 160.F, .r = kRadius, .varZ = 1.e-4F, .varR = kVarR, .varRPhi = 1.e-4F};

  const auto terms = pairTerms(cut, vertex, source, 0.F, sameRadius);
  BOOST_TEST(transverseCompatible(cut, terms));
  BOOST_TEST(std::isfinite(terms.varCollinearity));
  BOOST_TEST(!longitudinalCompatible(cut, terms));

  const double curvature = 0.5 / curvatureRadius(params.minPtGeV, params.bzKiloGauss);
  const double sourceArc = 2. * std::asin(0.5 * curvature * static_cast<double>(kRadius)) / curvature;
  const double cotTheta = 130. / sourceArc;
  const auto helixSource = helixPoint(curvature, sourceArc, cotTheta, 0.);
  const auto helixTarget = helixPoint(curvature, 160. / cotTheta, cotTheta, 0.);
  BOOST_TEST(pairCompatible(cut, vertex, helixSource, 0.F, helixTarget));
}

/// The phi window holds the widest azimuth the transverse test accepts, for a
/// steep source at the end of its layer and the vertex at the far end of its
/// range, where the scattering sum's path factor is largest.
BOOST_AUTO_TEST_CASE(ThePhiWindowHoldsTheTransverseTest)
{
  const auto table = scatteringTable<ToyBarrel4Disc8>(uniformLayerMaterial<ToyBarrel4Disc8>(0.02));
  const VertexConstraint beamSpot{};
  struct Case {
    LayerID from;
    LayerID to;
  };
  for (const double pt : {0.1, 1.}) {
    SearchWindowParams params;
    params.minPtGeV = pt;
    const auto cut = pairCut(params);
    for (const auto [from, to] : {Case{.from = 0, .to = 1}, Case{.from = 1, .to = 2}, Case{.from = 2, .to = 3}, Case{.from = 6, .to = 5}, Case{.from = 0, .to = 8}}) {
      const auto window = layerPairWindow<ToyBarrel4Disc8>(from, to, params, beamSpot, table);
      const double reach = params.nSigma * beamSpot.sigmaZ;
      const VertexSeed vertex{.z = static_cast<Scalar>(-reach), .varTransverse = 4.e-6F};
      const bool barrelSource = Detector<ToyBarrel4Disc8>::isBarrelLayer(from);
      const auto sourceR = static_cast<Scalar>(radiusMin<ToyBarrel4Disc8>(from));
      const auto sourceZ = static_cast<Scalar>(barrelSource ? zMax<ToyBarrel4Disc8>(from) : zMin<ToyBarrel4Disc8>(from));
      // Between two discs the target radius follows the seed's through the z ratio.
      const double rho = (zMin<ToyBarrel4Disc8>(to) + reach) / (static_cast<double>(sourceZ) + reach);
      const auto targetR = static_cast<Scalar>(Detector<ToyBarrel4Disc8>::isBarrelLayer(from) || Detector<ToyBarrel4Disc8>::isBarrelLayer(to) ? radiusMax<ToyBarrel4Disc8>(to) : rho * static_cast<double>(sourceR));
      const auto variance = static_cast<Scalar>(params.positionSigma * params.positionSigma);
      const PairPoint source{.x = sourceR, .z = sourceZ, .r = sourceR, .varRPhi = variance};
      const Scalar scatter = sourceScattering<ToyBarrel4Disc8>(table, cut, vertex, from, source);
      double widest = 0.;
      constexpr int kSteps = 15000;
      for (int step = 0; step < kSteps; ++step) {
        const double phi = 1.e-4 * step;
        const PairPoint target{
          .x = targetR * static_cast<Scalar>(std::cos(phi)),
          .y = targetR * static_cast<Scalar>(std::sin(phi)),
          .r = targetR,
          .varRPhi = variance,
        };
        if (transverseCompatible(cut, transversePairTerms(cut, vertex, source, scatter, target))) {
          widest = phi;
        }
      }
      BOOST_TEST_CONTEXT("pT " << pt << " link " << int(from) << " -> " << int(to))
      {
        BOOST_TEST(widest > 0.);
        BOOST_TEST(widest <= static_cast<double>(window.phiHalfWidth));
      }
    }
  }
}

/// A looper's two legs are found separately: a pair with the circle's turning
/// point between its measurements fails, though its curvature is allowed.
BOOST_AUTO_TEST_CASE(APairDoesNotSpanATurningPoint)
{
  SearchWindowParams params;
  params.minPtGeV = 0.1;
  const auto cut = pairCut(params);
  const VertexSeed vertex{};
  constexpr double kCurvature = 0.05; // above pT_min, R = 20 cm
  const auto source = helixPoint(kCurvature, 0.5 / kCurvature, 0., 0.);
  const auto outgoing = helixPoint(kCurvature, 2.5 / kCurvature, 0., 0.);
  const auto returning = helixPoint(kCurvature, 4. / kCurvature, 0., 0.);
  BOOST_TEST(returning.r > source.r);
  BOOST_TEST(transverseCompatible(cut, transversePairTerms(cut, vertex, source, 0.F, outgoing)));
  BOOST_TEST(!transverseCompatible(cut, transversePairTerms(cut, vertex, source, 0.F, returning)));
}

/// A vertex's box meets the measurement's, as in
/// O2's `isVertexCompatible`, its end inclusive -- a point vertex at a box's
/// start belongs to the box, one at its end to the next.
BOOST_AUTO_TEST_CASE(AVertexBelongsToTheBoxItMeets)
{
  const TimeEstPS frame{.start = 1000U, .width = 1000U};
  BOOST_TEST(vertexTimeCompatible(frame, TimeEstPS{.start = 1000U, .width = 0U}));
  BOOST_TEST(vertexTimeCompatible(frame, TimeEstPS{.start = 1999U, .width = 0U}));
  BOOST_TEST(!vertexTimeCompatible(frame, TimeEstPS{.start = 2000U, .width = 0U}));
  BOOST_TEST(!vertexTimeCompatible(frame, TimeEstPS{.start = 999U, .width = 0U}));
  // A found vertex's box spanning the edge belongs to both frames.
  BOOST_TEST(vertexTimeCompatible(frame, TimeEstPS{.start = 1900U, .width = 200U}));
  BOOST_TEST(vertexTimeCompatible(TimeEstPS{.start = 2000U, .width = 1000U}, TimeEstPS{.start = 1900U, .width = 200U}));
  BOOST_TEST(vertexTimeCompatible(frame, TimeEstPS{.start = 500U, .width = 500U}));
}

/// A point time -- a vertex -- in the picosecond two touching frames may
/// leave between them lies in neither: containment is quantised to 1 ps.
BOOST_AUTO_TEST_CASE(APointInTheGapBetweenFramesLiesInNeither)
{
  constexpr std::int64_t kGapFrameBC = 198;
  std::int64_t gaps = 0;
  for (std::int64_t k = 0; k < 600; ++k) {
    const auto box = frameBox(k * kGapFrameBC, (k + 1) * kGapFrameBC, 0);
    const auto next = frameBox((k + 1) * kGapFrameBC, (k + 2) * kGapFrameBC, 0);
    if (next.lower() == box.upper()) {
      continue;
    }
    ++gaps;
    const TimeEstPS point{.start = box.upper(), .width = 0U};
    BOOST_TEST(!vertexTimeCompatible(box, point));
    BOOST_TEST(!vertexTimeCompatible(next, point));
  }
  BOOST_TEST(gaps > 0);
}

} // namespace o2::gtrk
