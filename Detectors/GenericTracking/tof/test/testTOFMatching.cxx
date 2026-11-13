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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE GTrkTOFMatching
#include <boost/test/unit_test.hpp>

#include "GenericTracking/TOF/Transport.h"
#include "GenericTracking/TOF/Clusters.h"
#include "DataFormatsIOTOF/Cluster.h"
#include "DetectorsBase/GeometryManager.h"
#include "IOTOFBase/GeometryTGeo.h"
#include "IOTOFBase/IOTOFBaseParam.h"
#include "GenericTracking/TOF/Workflow.h"
#include "UniformMaterial.h"
#include "CommonConstants/LHCConstants.h"
#include <TMemFile.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <TGeoMatrix.h>
#include <numbers>
#include <vector>

using namespace o2::gtrk;
using namespace o2::gtrk::tof;
namespace
{
constexpr std::size_t pion = 2;
constexpr std::uint64_t largeBC = 8000000000ULL;
Track coarseTrack(SensorID id)
{
  Track track;
  track.state.surface = id;
  track.boxStart = {largeBC, 0.F};
  track.boxWidthNS = 100.;
  track.valid = true;
  return track;
}
Hit timedHit(int index, Layer layer, double time)
{
  Hit hit;
  hit.index = index;
  hit.layer = layer;
  hit.time = shifted({largeBC, 0.F}, time);
  return hit;
}
TrackingSurface plane(float x)
{
  TrackingSurface surface;
  surface.origin[0] = x;
  surface.axis0[1] = 1.F;
  surface.axis1[2] = 1.F;
  surface.normal[0] = 1.F;
  surface.bounds = {1000.F, 1000.F};
  return surface;
}
Track physicalTrack(float field)
{
  Track track = coarseTrack(0);
  track.surface = plane(1.F);
  track.state.par[kTheta] = std::numbers::pi_v<float> / 2.F;
  track.state.par[kQOverP] = 0.5F;
  for (int i = 0; i < kNBoundParams; ++i) {
    track.state.variance(i) = i < 2 ? 1e-4F : 1e-6F;
  }
  track.outwardSurface = plane(80.F);
  track.outwardState = track.state;
  BOOST_REQUIRE(propagate(Propagator{.bzKiloGauss = field}, track.outwardState, track.surface, track.outwardSurface, PropagationDirection::kForward).status == PropagationStatus::kOk);
  return track;
}
} // namespace

BOOST_AUTO_TEST_CASE(SplitTimingRetainsPicosecondsAcrossSecondsAndBCBoundaries)
{
  Timestamp a{largeBC, 24.9F};
  const auto b = shifted(a, 0.1);
  BOOST_TEST(b.bc == a.bc + 1);
  BOOST_TEST(std::abs(differenceNS(b, a) - 0.1) < 0.000002);
  const auto c = shifted(b, -0.101);
  BOOST_TEST(std::abs(differenceNS(c, a) + 0.001) < 0.000002);
  BOOST_CHECK_THROW(shifted({0, 0.F}, -1.), std::out_of_range);
  BOOST_TEST(std::abs(differenceNS(a, b) + differenceNS(b, a)) < 1e-12);
}

BOOST_AUTO_TEST_CASE(OuterThenGuidedInnerReservesHitsBeforeRemainingTracks)
{
  std::vector<Track> tracks{coarseTrack(0), coarseTrack(1)};
  std::vector<Hit> hits{timedHit(7, Layer::Outer, 8.), timedHit(11, Layer::Inner, 6.), timedHit(19, Layer::Inner, 6.3)};
  const auto predictor = [](const Track& track, const Hit& hit) {
    Prediction p;
    p.spatialChi2 = track.state.surface == 0 ? 2.F : 0.F;
    if (hit.layer == Layer::Outer && track.state.surface != 0) {
      return Prediction{};
    }
    p.flightNS[pion] = hit.layer == Layer::Outer ? 3. : 1.;
    p.valid[pion] = true;
    return p;
  };
  auto output = match(tracks, hits, predictor);
  BOOST_REQUIRE_EQUAL(output.size(), 2);
  BOOST_TEST(output[0].outer == 7);
  BOOST_TEST(output[0].inner == 11);
  BOOST_TEST(output[1].outer == -1);
  BOOST_TEST(output[1].inner == 19);
  BOOST_TEST(output[0].innerCandidates == 1);
  BOOST_TEST(output[1].innerCandidates == 1);
  BOOST_CHECK(output[0].hypotheses[pion].valid);
  BOOST_TEST(std::abs(differenceNS(output[0].hypotheses[pion].time, {largeBC, 5.F})) < 0.000002);
  BOOST_TEST(std::abs(output[0].hypotheses[pion].sigmaNS - 0.02 / std::sqrt(2.)) < 1e-7);
  BOOST_TEST(std::abs(output[1].hypotheses[pion].sigmaNS - 0.02) < 1e-7);
}

BOOST_AUTO_TEST_CASE(OuterOnlySurvivesAnIncompatibleInnerTime)
{
  std::vector<Track> tracks{coarseTrack(0)};
  std::vector<Hit> hits{timedHit(1, Layer::Outer, 8.), timedHit(2, Layer::Inner, 7.)};
  const auto output = match(tracks, hits, [](const Track&, const Hit& hit) {
    Prediction p;
    p.spatialChi2 = 0.F;
    p.valid[pion] = true;
    p.flightNS[pion] = hit.layer == Layer::Outer ? 3. : 1.;
    return p;
  });
  BOOST_TEST(output[0].outer == 1);
  BOOST_TEST(output[0].inner == -1);
  BOOST_TEST(output[0].hypotheses[pion].valid);
  BOOST_TEST(std::abs(differenceNS(output[0].hypotheses[pion].time, {largeBC, 5.F})) < 0.000002);
}

BOOST_AUTO_TEST_CASE(CommonTrackErrorCancelsInTheGuidedGate)
{
  std::vector<Track> tracks{coarseTrack(0)};
  tracks[0].state.variance(kQOverP) = 1.F;
  std::vector<Hit> hits{timedHit(0, Layer::Outer, 8.), timedHit(1, Layer::Inner, 6.3)};
  const auto predictor = [](const Track&, const Hit& hit) {
    Prediction p;
    p.spatialChi2 = 0.F;
    p.valid[pion] = true;
    p.flightNS[pion] = hit.layer == Layer::Outer ? 3. : 1.;
    p.derivative[pion][kQOverP] = 10.; // Large but exactly shared uncertainty.
    return p;
  };
  auto rejected = match(tracks, hits, predictor);
  BOOST_TEST(rejected[0].inner == -1);
  hits[1].time = shifted({largeBC, 0.F}, 6.);
  auto accepted = match(tracks, hits, predictor);
  BOOST_TEST(accepted[0].inner == 1);
  BOOST_TEST(std::abs(accepted[0].hypotheses[pion].sigmaNS - std::sqrt(100.0002)) < 1e-5);
}

BOOST_AUTO_TEST_CASE(DeterministicArbitrationAndCoarseTimeAcceptance)
{
  std::vector<Track> tracks{coarseTrack(0), coarseTrack(1), coarseTrack(2)};
  tracks[2].valid = false;
  std::vector<Hit> hits{timedHit(9, Layer::Inner, 6.), timedHit(3, Layer::Inner, 6.), timedHit(20, Layer::Inner, 106.)};
  auto predictor = [](const Track&, const Hit&) {
    Prediction p;
    p.spatialChi2 = 0.F;
    p.valid[pion] = true;
    p.flightNS[pion] = 1.;
    return p;
  };
  auto first = match(tracks, hits, predictor);
  std::reverse(hits.begin(), hits.end());
  auto second = match(tracks, hits, predictor);
  BOOST_TEST(first[0].inner == 3);
  BOOST_TEST(first[1].inner == 9);
  BOOST_TEST(first[2].inner == -1);
  BOOST_TEST(first[0].innerCandidates == 2);
  BOOST_TEST(first[0].inner == second[0].inner);
  BOOST_TEST(first[1].inner == second[1].inner);
  hits.push_back(hits.front());
  BOOST_CHECK_THROW(match(tracks, hits, predictor), std::invalid_argument);
  BOOST_TEST(match(tracks, {}, predictor).size() == tracks.size());
}

BOOST_AUTO_TEST_CASE(FlightTimeMatchesAnalyticHelixForSignedAndZeroFields)
{
  const UniformMaterial vacuum;
  for (const float field : {0.F, -20.F, 20.F}) {
    const auto track = physicalTrack(field);
    Hit hit = timedHit(1, Layer::Outer, 8.);
    hit.surface = plane(92.F);
    auto truth = track.state;
    BOOST_REQUIRE(propagate(Propagator{.bzKiloGauss = field}, truth, track.surface, hit.surface, PropagationDirection::kForward).status == PropagationStatus::kOk);
    hit.loc0 = truth.loc0();
    hit.loc1 = truth.loc1();
    hit.varLoc0 = hit.varLoc1 = 1e-4F;
    const auto prediction = predict(track, hit, Propagator{.bzKiloGauss = field}, vacuum);
    BOOST_TEST(prediction.spatialChi2 < 1e-5F);
    const double curvature = 0.5 * field * o2::constants::math::B2C;
    const double path = field == 0.F ? 91. : std::asin(curvature * 91.) / curvature;
    for (std::size_t mass = 0; mass < kHypotheses; ++mass) {
      BOOST_REQUIRE(prediction.valid[mass]);
      const double expected = path * std::sqrt(1. + kMassGeV[mass] * kMassGeV[mass] / 4.) / kLightSpeedCMNS;
      BOOST_TEST(std::abs(prediction.flightNS[mass] - expected) < 0.00001); // 0.01 ps
    }
  }
}

BOOST_AUTO_TEST_CASE(MaterialChangesFlightTimeAndTransportFailureIsExplicit)
{
  const auto track = physicalTrack(20.F);
  const Propagator prop{.bzKiloGauss = 20.F};
  const auto target = plane(92.F);
  const auto vacuum = transport(track.state, track.surface, target, prop, UniformMaterial{}, PropagationDirection::kForward, Parameters{});
  const auto material = transport(track.state, track.surface, target, prop, UniformMaterial{.x2X0PerCm = 0.001F, .rho = 0.05F}, PropagationDirection::kForward, Parameters{}, true);
  BOOST_REQUIRE(vacuum.ok);
  BOOST_REQUIRE(material.ok);
  BOOST_TEST(material.flightNS > vacuum.flightNS);
  BOOST_TEST(std::abs(material.state.qOverP()) > std::abs(vacuum.state.qOverP()));
  auto stopped = transport(track.state, track.surface, target, prop, UniformMaterial{}, PropagationDirection::kForward, Parameters{.maxPathCM = 10.});
  BOOST_CHECK(!stopped.ok);
  auto curling = track;
  curling.state.par[kQOverP] = 100.F;
  BOOST_CHECK(!transport(curling.state, curling.surface, target, prop, UniformMaterial{}, PropagationDirection::kForward, Parameters{}).ok);
  Hit missed = timedHit(0, Layer::Outer, 8.);
  missed.surface = target;
  missed.surface.bounds = {0.001F, 0.001F};
  BOOST_CHECK(!predict(track, missed, prop, UniformMaterial{}).valid[pion]);
}

BOOST_AUTO_TEST_CASE(MassHypothesesAndMatchRecordsSurviveROOTIO)
{
  std::vector<Match> output(1);
  output[0].outer = 17;
  auto& h = output[0].hypotheses[pion];
  h.time = {largeBC, 12.123F};
  h.valid = true;
  h.sigmaNS = 0.02F;
  TMemFile file("tof-matches.root", "RECREATE");
  file.WriteObject(&output, "matches");
  const auto* restored = file.Get<std::vector<Match>>("matches");
  BOOST_REQUIRE(restored != nullptr);
  BOOST_REQUIRE_EQUAL(restored->size(), 1);
  BOOST_TEST(restored->front().outer == 17);
  BOOST_TEST(restored->front().hypotheses[pion].time.bc == largeBC);
  BOOST_TEST(restored->front().hypotheses[pion].time.fineNS == h.time.fineNS);
  BOOST_CHECK(!restored->front().hypotheses[0].valid);
}

BOOST_AUTO_TEST_CASE(WorkflowDeclaresTheMatchingContract)
{
  const auto spec = getMatcherSpec();
  BOOST_TEST(spec.name == "gtrk-tof-matcher");
  BOOST_TEST(spec.outputs.size() == 1);
  BOOST_TEST(getMatchWriterSpec().name == "gtrk-tof-match-writer");
}

BOOST_AUTO_TEST_CASE(ClusterAdapterUsesRealSegmentedSensorPlacements)
{
  const char* geometryFile = std::getenv("GTRK_TOF_GEOMETRY");
  if (!geometryFile) {
    BOOST_TEST_MESSAGE("Set GTRK_TOF_GEOMETRY to validate the adapter against an ALICE3 geometry");
    return;
  }
  o2::base::GeometryManager::loadGeometry(geometryFile, false, true);
  auto* geometry = o2::iotof::GeometryTGeo::Instance();
  geometry->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  const auto& specs = o2::iotof::ChipSpecificsParam::Instance();
  std::vector<o2::iotof::Cluster> clusters;
  for (int layer = 0; layer < 2; ++layer) {
    const auto chip = geometry->getIOTOFFirstChipIndex(layer);
    clusters.emplace_back(100, 50, 2, 2, 9, 0, chip, 24.123F, largeBC);
  }
  auto invalid = clusters.front();
  invalid.setPatternID(o2::iotof::Cluster::InvalidPatternID);
  clusters.insert(clusters.begin() + 1, invalid);
  const auto hits = convertClusters(clusters, *geometry, 0.020F);
  BOOST_REQUIRE_EQUAL(hits.size(), 2);
  BOOST_TEST(hits[0].index == 0);
  BOOST_TEST(hits[1].index == 2);
  BOOST_CHECK(hits[0].layer == Layer::Inner);
  BOOST_CHECK(hits[1].layer == Layer::Outer);
  const double local[3]{0.5 * (specs.ActiveMatrixSizeRows() - specs.PassiveEdgeTop + specs.PassiveEdgeReadOut - specs.PitchRow) - 100.5 * specs.PitchRow,
                        0., 50.5 * specs.PitchCol + 0.5 * (specs.PitchCol - specs.ActiveMatrixSizeCols())};
  for (const auto& hit : hits) {
    double expected[3];
    const auto* matrix = geometry->extractMatrixSensor(hit.surface.id);
    BOOST_REQUIRE(matrix != nullptr);
    matrix->LocalToMaster(local, expected);
    Scalar actual[3];
    surfaceToGlobal(hit.surface, hit.loc0, hit.loc1, actual);
    for (int i = 0; i < 3; ++i) {
      BOOST_TEST(std::abs(actual[i] - expected[i]) < 0.00005);
    }
    BOOST_TEST(hit.time.bc == largeBC);
    BOOST_TEST(hit.time.fineNS == clusters.front().getTime());
  }
}

BOOST_AUTO_TEST_CASE(PhysicalTwoLayerMatchingRecoversTheReferenceTime)
{
  const auto track = physicalTrack(-20.F);
  const Propagator prop{.bzKiloGauss = -20.F};
  const UniformMaterial vacuum;
  std::vector<Hit> hits;
  for (const auto layer : {Layer::Outer, Layer::Inner}) {
    auto hit = timedHit(layer == Layer::Outer ? 9 : 2, layer, 0.);
    hit.surface = plane(layer == Layer::Outer ? 92.F : 21.F);
    const auto crossing = transport(track.state, track.surface, hit.surface, prop, vacuum, PropagationDirection::kForward, Parameters{});
    BOOST_REQUIRE(crossing.ok);
    hit.loc0 = crossing.state.loc0();
    hit.loc1 = crossing.state.loc1();
    hit.varLoc0 = hit.varLoc1 = 1e-4F;
    hit.time = shifted({largeBC, 5.F}, crossing.flightNS);
    hits.push_back(hit);
  }
  const std::vector<Track> tracks{track};
  const auto output = match(tracks, hits, [&](const Track& t, const Hit& h) { return predict(t, h, prop, vacuum); });
  BOOST_TEST(output[0].outer == 9);
  BOOST_TEST(output[0].inner == 2);
  BOOST_REQUIRE(output[0].hypotheses[pion].valid);
  BOOST_TEST(std::abs(differenceNS(output[0].hypotheses[pion].time, {largeBC, 5.F})) < 0.000002);
}

BOOST_AUTO_TEST_CASE(InnerTOFCanPrecedeTheTrackReferencePoint)
{
  const auto track = physicalTrack(0.F);
  auto hit = timedHit(4, Layer::Inner, 4.9);
  hit.surface = plane(0.5F);
  hit.varLoc0 = hit.varLoc1 = 1e-4F;
  const auto p = predict(track, hit, Propagator{.bzKiloGauss = 0.F}, UniformMaterial{});
  BOOST_REQUIRE(p.valid[pion]);
  BOOST_TEST(p.flightNS[pion] < 0.);
  BOOST_TEST(std::abs(p.flightNS[pion] + 0.5 * std::sqrt(1. + kMassGeV[pion] * kMassGeV[pion] / 4.) / kLightSpeedCMNS) < 1e-7);
}

BOOST_AUTO_TEST_CASE(OuterTimingDoesNotRepeatTheSpatialBoundsGate)
{
  auto track = physicalTrack(0.F);
  // The outward fit correctly predicts the outer sensor. The independently
  // fitted inward state has a different direction; its long flight-time
  // projection can fall outside that chip after crossing the tracker.
  track.state.par[kPhi] = 0.01F;
  auto hit = timedHit(0, Layer::Outer, 8.);
  hit.surface = plane(92.F);
  hit.surface.bounds = {0.1F, 0.1F};
  hit.varLoc0 = hit.varLoc1 = 1e-4F;
  const Propagator prop{.bzKiloGauss = 0.F};
  const auto flight = transport(track.state, track.surface, hit.surface, prop, UniformMaterial{}, PropagationDirection::kForward, Parameters{});
  BOOST_REQUIRE(flight.ok);
  BOOST_CHECK(!flight.inside);
  const auto prediction = predict(track, hit, prop, UniformMaterial{});
  BOOST_TEST(prediction.spatialChi2 < 1e-5F);
  BOOST_REQUIRE(prediction.valid[pion]);
  BOOST_TEST(std::abs(prediction.flightNS[pion] - flight.flightNS) < 1e-7);
}
