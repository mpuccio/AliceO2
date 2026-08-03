// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// M5a native-vs-legacy ITS refit A/B characterization harness
// (doc/design/0001-descriptor-driven-operation-boundary.md Sec 6).
//
// Reuses the existing, already-committed native refit primitives
// (driveRefitLeg<CylinderCylinder>/nativeRefitTrackCylinderCylinder,
// NativeCylinderCylinderRefitDriver.h) and the frozen legacy ITS refit path
// (o2::its::track::fitTrack/refitTrack/refitTrackSeed, ITStracking/
// TrackHelpers.h) exactly as testNativeCylinderCylinderRefitDriver.cxx
// (Gate 3 Slice B) already does, over an identical set of independently
// constructed (never one derived from the other) seed/measurement fixtures.
//
// This file does not replace or duplicate that earlier test's own coverage
// or assertions -- it stays in the repository unchanged, as the accepted
// per-scenario oracle. What this harness adds:
//   1. a single data-driven scenario table (easy to extend for future
//      physics QA) instead of one hand-written test case per scenario;
//   2. a machine-readable JSON characterization report
//      (itsNativeVsLegacyRefitCharacterization.json, written to the test's
//      own working directory -- never the source tree) recording, per
//      scenario: both engines' success/failure classification and failure
//      reason, the structural invariants (cluster pattern/hit mask,
//      timestamp) that must match exactly when both succeed, and the
//      numeric parameter/chi2 deltas as characterization data -- never as a
//      pass/fail tolerance.
//
// Exactly like testNativeCylinderCylinderRefitDriver.cxx, this harness
// asserts structural parity (success/failure classification, cluster
// pattern, timestamp) exactly, and reports numeric differences without
// inventing a bound: the two engines' per-hit Kalman update arithmetic are
// independently implemented and are already known (that earlier test's own
// SingleHitZeroMaterialCovarianceDivergenceCharacterization) to disagree by
// a non-trivial amount even at exactly zero material. This harness does not
// re-litigate that finding; it extends its evidence base.
//
// Rollback-behavior asymmetry (recorded, not asserted as parity): the native
// driver is fully transactional on failure (every output parameter is
// written only once, at the very end of a successful call -- see
// NativeCylinderCylinderRefitDriver.h), which this harness verifies with a
// poison-sentinel check. The frozen legacy fitTrack is documented (Gate 3
// Slice B evidence, see testNativeCylinderCylinderRefitDriver.cxx's header
// comment) to mutate its working TrackParCov in place as it walks each
// layer, so a leg that fails partway through still leaves that partial
// progress in the caller's TrackITSInternal -- this harness does not assert
// legacy is transactional, because it verifiably is not.
//
// Does not wire native refit into production, does not change any
// cell/tracklet/refit formula, does not reintroduce TransitionPolicyTag
// publicly, does not replace LegacyTrackerScratch, and does not touch any
// workflow/default/output contract.

#define BOOST_TEST_MODULE ITSMFTNativeVsLegacyRefitCharacterizationHarness
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "DetectorsBase/Propagator.h"
#include "ITSMFTTracking/NativeCylinderCylinderRefitDriver.h"
#include "ITStracking/Cell.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/TrackHelpers.h"
#include "ITStracking/TrackITSInternal.h"

using namespace o2::itsmft::tracking;

namespace
{

// See testNativeCylinderCylinderRefitDriver.cxx's identical fixture for why
// this is required: every propagateToX call below uses the explicit-bz
// overload, so the singleton's own (absent) field state is never read.
struct PropagatorUninitializedFixture {
  PropagatorUninitializedFixture() { o2::base::Propagator::Instance(true); }
};

} // namespace

BOOST_TEST_GLOBAL_FIXTURE(PropagatorUninitializedFixture);

namespace
{

constexpr int NLayers = ITSNLayers;
constexpr float Bz = 5.f;
constexpr float LocalX = 2.5f;
constexpr float Alpha = 0.3f;

std::string failureReasonName(OperationFailureReason reason)
{
  switch (reason) {
    case OperationFailureReason::SourceFamilyMismatch:
      return "SourceFamilyMismatch";
    case OperationFailureReason::PredictedChi2Failure:
      return "PredictedChi2Failure";
    case OperationFailureReason::LegAcceptanceFailure:
      return "LegAcceptanceFailure";
    case OperationFailureReason::MinPtFailure:
      return "MinPtFailure";
    case OperationFailureReason::ReseedNotSupported:
      return "ReseedNotSupported";
    default:
      return "Other(" + std::to_string(static_cast<int>(reason)) + ")";
  }
}

// One shared numeric description of "NLayers barrel measurements", built
// independently into both a legacy (Cluster/TrackingFrameInfo) and a native
// (SurfaceMeasurement/SurfaceDescriptor) representation -- same convention
// as testNativeCylinderCylinderRefitDriver.cxx's SharedFixture, so a bug in
// legacy::exportBarrelTrackParCov itself would still show up as a mismatch
// rather than being built into both sides identically.
struct SharedFixture {
  struct LegacyLayer {
    o2::its::Cluster cluster;
    o2::its::TrackingFrameInfo tfInfo;
  };
  std::array<LegacyLayer, NLayers> legacyLayers{};
  std::array<std::vector<SurfaceMeasurement>, NLayers> nativeStorage{};
  LayerMeasurementSpans<NLayers> layerMeasurements{};
  std::vector<SurfaceDescriptor> catalogSurfaces;
  SurfaceCatalogView catalog{};
  const o2::its::TrackingFrameInfo* tfInfoPtrs[NLayers];
  const o2::its::Cluster* clusterPtrs[NLayers];

  SharedFixture(float xOverX0, float arealDensityGPerCm2, int nHitLayers)
  {
    catalogSurfaces.resize(NLayers);
    for (int layer = 0; layer < NLayers; ++layer) {
      const float u = 0.8f + 0.01f * static_cast<float>(layer);
      const float v = -0.45f - 0.02f * static_cast<float>(layer);
      constexpr float covUU = 0.04f;
      constexpr float covUV = 0.f;
      constexpr float covVV = 0.09f;

      catalogSurfaces[layer].id = SurfaceId{static_cast<uint16_t>(layer)};
      catalogSurfaces[layer].detectorSurfaceIndex = static_cast<uint16_t>(layer);
      catalogSurfaces[layer].kind = SurfaceKind::Cylinder;
      catalogSurfaces[layer].material = NominalSurfaceMaterial{xOverX0, arealDensityGPerCm2};

      if (layer < nHitLayers) {
        legacyLayers[layer].cluster = o2::its::Cluster(0.f, 0.f, 0.f, layer);
        legacyLayers[layer].tfInfo = o2::its::TrackingFrameInfo(0.f, 0.f, 0.f, LocalX, Alpha, {u, v}, {covUU, covUV, covVV});

        SurfaceMeasurement m{};
        m.frame.q = LocalX;
        m.frame.frameAngle = Alpha;
        m.frame.u = u;
        m.frame.v = v;
        m.covariance = {covUU, covUV, covVV};
        m.cluster = ClusterRef{ClusterSourceId{0}, 0};
        m.surface = SurfaceId{static_cast<uint16_t>(layer)};
        nativeStorage[layer] = {m};
      } else {
        nativeStorage[layer] = {}; // hole
      }
      layerMeasurements[layer] = gsl::span<const SurfaceMeasurement>(nativeStorage[layer]);
      tfInfoPtrs[layer] = &legacyLayers[layer].tfInfo;
      clusterPtrs[layer] = &legacyLayers[layer].cluster;
    }
    catalog = SurfaceCatalogView{catalogSurfaces.data(), static_cast<uint32_t>(catalogSurfaces.size())};
  }
};

SurfaceKinematicState makeInitialState()
{
  SurfaceKinematicState state{};
  state.parameters[0] = 1.25f;
  state.parameters[1] = -0.75f;
  state.parameters[2] = 0.2f;
  state.parameters[3] = -0.35f;
  state.parameters[4] = 0.8f;
  state.referenceCoordinate = LocalX;
  state.alpha = Alpha;
  state.family = StateFamily::Barrel;
  state.absCharge = 1;
  state.pid = o2::track::PID::Pion;
  return state;
}

o2::track::TrackParCovF makeInitialTrackParCovF()
{
  const auto state = makeInitialState();
  o2::track::TrackParCovF::params_t parameters{};
  o2::track::TrackParCovF::covMat_t covariance{};
  for (uint8_t i = 0; i < 5; ++i) {
    parameters[i] = state.parameters[i];
  }
  covariance.fill(0.f); // irrelevant: both refit paths reset covariance before leg A
  return o2::track::TrackParCovF{state.referenceCoordinate, state.alpha, parameters, covariance, state.absCharge, state.pid};
}

TrackSeedN<NLayers> makeNativeSeed(int nHitLayers)
{
  TrackSeedN<NLayers> seed{};
  seed.state() = makeInitialState();
  uint16_t mask = 0;
  for (int layer = 0; layer < nHitLayers; ++layer) {
    seed.getClusters()[layer] = 0;
    mask |= static_cast<uint16_t>(uint16_t(1) << layer);
  }
  seed.setHitLayerMask(LayerMask{mask});
  return seed;
}

o2::its::TrackSeed<NLayers> makeLegacySeed(int nHitLayers)
{
  o2::its::TrackSeed<NLayers> seed{};
  static_cast<o2::track::TrackParCovF&>(seed) = makeInitialTrackParCovF();
  o2::its::LayerMask mask{};
  for (int layer = 0; layer < nHitLayers; ++layer) {
    seed.getClusters()[layer] = 0;
    mask.set(layer);
  }
  seed.setHitLayerMask(mask);
  return seed;
}

// One scenario: a named point in the coverage matrix the M5a task requires
// (single-hit, multi-hit, holes, zero/nonzero material, shift-reference,
// repeat-refit-out, and the two failure/rollback modes below).
struct Scenario {
  std::string name;
  int nHitLayers{NLayers};
  float xOverX0{0.f};
  float arealDensityGPerCm2{0.f};
  float maxChi2ClusterAttachment{60.f};
  float maxChi2NDF{30.f};
  bool shiftRefToCluster{false};
  bool repeatRefitOut{false};
  float minPtThreshold{0.f}; // 0 disables the MinPt check on both sides
  bool expectBothSucceed{true};
};

struct ScenarioResult {
  std::string scenario;
  bool legacyOk{false};
  bool nativeOk{false};
  std::string nativeFailureReason;
  bool patternMatches{false};
  bool timestampMatches{false};
  bool numericComparable{false};
  double legacyQ2Pt{0.}, nativeQ2Pt{0.}, deltaQ2Pt{0.};
  double legacyChi2{0.}, nativeChi2{0.}, deltaChi2{0.};
};

std::vector<ScenarioResult>& allResults()
{
  static std::vector<ScenarioResult> results;
  return results;
}

ScenarioResult runScenario(const Scenario& scenario)
{
  SharedFixture fixture(scenario.xOverX0, scenario.arealDensityGPerCm2, scenario.nHitLayers);
  std::vector<float> minPt(NLayers + 1, 0.f);
  if (scenario.minPtThreshold > 0.f) {
    minPt[NLayers - scenario.nHitLayers] = scenario.minPtThreshold;
  }
  std::vector<float> layerxX0(NLayers, scenario.xOverX0);
  std::vector<float> layerRadii(NLayers, LocalX);

  ScenarioResult result{};
  result.scenario = scenario.name;

  // --- Legacy oracle: o2::its::track::refitTrackSeed (frozen chain). ---
  const auto legacySeed = makeLegacySeed(scenario.nHitLayers);
  const o2::its::track::TrackFitContext<NLayers> fitCtx{
    fixture.tfInfoPtrs, layerxX0.data(), NLayers, Bz,
    scenario.maxChi2ClusterAttachment, scenario.maxChi2NDF,
    o2::base::Propagator::Instance(), o2::base::PropagatorF::MatCorrType::USEMatCorrNONE,
    scenario.shiftRefToCluster, scenario.repeatRefitOut};
  o2::its::TrackITSInternal<NLayers> internalTrack;
  result.legacyOk = o2::its::track::refitTrackSeed<NLayers>(legacySeed, internalTrack, fitCtx, fixture.clusterPtrs,
                                                            layerRadii.data(), minPt.data(), /*reseedIfShorter=*/0);
  o2::its::TrackITSExt legacyTrack{};
  if (result.legacyOk) {
    legacyTrack = o2::its::makeTrackITSExt(internalTrack);
  }

  // --- Native driver: driveRefitLeg<CylinderCylinder>/
  // nativeRefitTrackCylinderCylinder (unwired). ---
  const auto nativeSeed = makeNativeSeed(scenario.nHitLayers);
  SurfaceKinematicState paramIn{};
  paramIn.parameters[0] = -777.f; // poison sentinel: must stay untouched on failure
  SurfaceKinematicState paramOut{};
  paramOut.parameters[0] = -888.f;
  float chi2 = -999.f;
  OperationFailureReason reason{};
  result.nativeOk = nativeRefitTrackCylinderCylinder<NLayers>(
    nativeSeed, fixture.layerMeasurements, fixture.catalog, Bz, scenario.shiftRefToCluster,
    scenario.maxChi2ClusterAttachment, scenario.maxChi2NDF, scenario.repeatRefitOut,
    gsl::span<const float>(minPt), /*reseedIfShorter=*/0, paramIn, paramOut, chi2, reason);
  result.nativeFailureReason = result.nativeOk ? "" : failureReasonName(reason);

  o2::its::TrackITSExt nativeTrack{};
  if (result.nativeOk) {
    BOOST_REQUIRE(exportNativeRefitToTrackITSExt<NLayers>(nativeSeed, paramIn, paramOut, chi2, nativeTrack));
  } else {
    // Native's own transactionality contract: untouched on any failure.
    BOOST_CHECK_EQUAL(paramIn.parameters[0], -777.f);
    BOOST_CHECK_EQUAL(paramOut.parameters[0], -888.f);
    BOOST_CHECK_EQUAL(chi2, -999.f);
  }

  // --- Structural invariants: asserted exactly. ---
  BOOST_CHECK_MESSAGE(result.legacyOk == result.nativeOk,
                      scenario.name << ": success/failure classification differs (legacy=" << result.legacyOk
                                    << " native=" << result.nativeOk << ")");
  BOOST_CHECK_EQUAL(result.legacyOk, scenario.expectBothSucceed);

  if (result.legacyOk && result.nativeOk) {
    result.patternMatches = (legacyTrack.getPattern() == nativeTrack.getPattern());
    result.timestampMatches = (legacyTrack.getTimeStamp().getTimeStamp() == nativeTrack.getTimeStamp().getTimeStamp());
    BOOST_CHECK_MESSAGE(result.patternMatches, scenario.name << ": cluster/hit-mask pattern differs");
    BOOST_CHECK_MESSAGE(result.timestampMatches, scenario.name << ": timestamp differs");

    result.numericComparable = true;
    result.legacyQ2Pt = legacyTrack.getParamIn().getQ2Pt();
    result.nativeQ2Pt = nativeTrack.getParamIn().getQ2Pt();
    result.deltaQ2Pt = std::abs(result.legacyQ2Pt - result.nativeQ2Pt);
    result.legacyChi2 = legacyTrack.getChi2();
    result.nativeChi2 = nativeTrack.getChi2();
    result.deltaChi2 = std::abs(result.legacyChi2 - result.nativeChi2);
    BOOST_TEST_MESSAGE(scenario.name << ": legacy Q2Pt=" << result.legacyQ2Pt << " native Q2Pt=" << result.nativeQ2Pt
                                     << " |delta|=" << result.deltaQ2Pt << "; legacy chi2=" << result.legacyChi2
                                     << " native chi2=" << result.nativeChi2 << " |delta|=" << result.deltaChi2
                                     << " (characterization only, not a tolerance assertion)");
  } else {
    // Both engines rejected this scenario: no numeric comparison is
    // meaningful. Native's own rollback is checked above; legacy's fitTrack
    // is documented non-transactional (see this file's header comment) and
    // is not asserted untouched here.
    BOOST_CHECK_MESSAGE(!result.legacyOk && !result.nativeOk,
                        scenario.name << ": expected both engines to reject this scenario");
  }

  allResults().push_back(result);
  return result;
}

std::vector<Scenario> buildScenarioTable()
{
  std::vector<Scenario> scenarios;
  scenarios.push_back({"SingleHitZeroMaterial", 1, 0.f, 0.f, 60.f, 30.f, false, false, 0.f, false});
  // A lone hit always fails fitTrack's/driveRefitLeg's own trailing
  // nCl>=3-shaped NDF check (2*1-5 < 0) on both sides -- this is a
  // structural, not a numeric, rejection, so expectBothSucceed is false.
  scenarios.push_back({"MultiHitAllLayersZeroMaterial", NLayers, 0.f, 0.f, 60.f, 30.f, false, false, 0.f, true});
  scenarios.push_back({"HolesTwoOfSevenMissing", NLayers - 2, 0.f, 0.f, 60.f, 30.f, false, false, 0.f, true});
  scenarios.push_back({"NonzeroMaterial", NLayers, 0.01f, 0.001f, 60.f, 30.f, false, false, 0.f, true}); // O(ITS Si layer) xOverX0
  scenarios.push_back({"ShiftReferenceToClusterOn", NLayers, 0.f, 0.f, 60.f, 30.f, true, false, 0.f, true});
  scenarios.push_back({"RepeatRefitOutOn", NLayers, 0.f, 0.f, 60.f, 30.f, false, true, 0.f, true});
  scenarios.push_back({"FailureChi2GateRejectsFourthHit", NLayers, 0.f, 0.f, 1.e-8f, 30.f, false, false, 0.f, false});
  scenarios.push_back({"FailureMinPtRejectsAfterLegB", NLayers, 0.f, 0.f, 60.f, 30.f, false, false, 1.e6f, false});
  return scenarios;
}

void writeJsonReport()
{
  std::ofstream out("itsNativeVsLegacyRefitCharacterization.json");
  out << "{\n";
  out << "  \"harness\": \"ITSMFTNativeVsLegacyRefitCharacterizationHarness\",\n";
  out << "  \"designNote\": \"Detectors/ITSMFT/common/tracking/doc/design/0001-descriptor-driven-operation-boundary.md\",\n";
  out << "  \"note\": \"Structural fields (patternMatches/timestampMatches) are asserted exactly by this test; numeric fields are characterization only, not a tolerance-bounded pass/fail claim.\",\n";
  out << "  \"scenarios\": [\n";
  for (size_t i = 0; i < allResults().size(); ++i) {
    const auto& r = allResults()[i];
    out << "    {\n";
    out << "      \"scenario\": \"" << r.scenario << "\",\n";
    out << "      \"legacyOk\": " << (r.legacyOk ? "true" : "false") << ",\n";
    out << "      \"nativeOk\": " << (r.nativeOk ? "true" : "false") << ",\n";
    out << "      \"nativeFailureReason\": \"" << r.nativeFailureReason << "\",\n";
    out << "      \"numericComparable\": " << (r.numericComparable ? "true" : "false") << ",\n";
    if (r.numericComparable) {
      out << "      \"patternMatches\": " << (r.patternMatches ? "true" : "false") << ",\n";
      out << "      \"timestampMatches\": " << (r.timestampMatches ? "true" : "false") << ",\n";
      out << "      \"legacyQ2Pt\": " << r.legacyQ2Pt << ",\n";
      out << "      \"nativeQ2Pt\": " << r.nativeQ2Pt << ",\n";
      out << "      \"deltaQ2Pt\": " << r.deltaQ2Pt << ",\n";
      out << "      \"legacyChi2\": " << r.legacyChi2 << ",\n";
      out << "      \"nativeChi2\": " << r.nativeChi2 << ",\n";
      out << "      \"deltaChi2\": " << r.deltaChi2 << "\n";
    } else {
      out << "      \"bothRejected\": " << (!r.legacyOk && !r.nativeOk ? "true" : "false") << "\n";
    }
    out << "    }" << (i + 1 < allResults().size() ? ",\n" : "\n");
  }
  out << "  ]\n}\n";
  BOOST_TEST_MESSAGE("wrote itsNativeVsLegacyRefitCharacterization.json (" << allResults().size() << " scenarios)");
}

} // namespace

BOOST_AUTO_TEST_CASE(RunScenarioMatrixAndWriteCharacterizationReport)
{
  for (const auto& scenario : buildScenarioTable()) {
    runScenario(scenario);
  }
  writeJsonReport();
}
