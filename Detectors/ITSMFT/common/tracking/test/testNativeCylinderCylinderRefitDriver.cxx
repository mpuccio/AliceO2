// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Gate 3 Slice B: independent legacy-oracle coverage for
// nativeRefitTrackCylinderCylinder/exportNativeRefitToTrackITSExt
// (NativeCylinderCylinderRefitDriver.h), comparing the new unwired native
// driver directly against the frozen o2::its::track::refitTrackSeed/
// refitTrack/fitTrack chain (ITStracking/TrackHelpers.h) it is designed to
// replace. Both fixtures are built independently (not one derived from the
// other via the same export adapter under test), from the same numeric
// values, so a bug in legacy::exportBarrelTrackParCov itself would still be
// caught by a parity mismatch.
//
// Every parity fixture uses reseedIfShorter = 0 (never triggers legacy's
// seedTrackForRefit mid-track geometric reseed) since the native driver
// does not reproduce that reseed formula in this slice (see
// NativeCylinderCylinderRefitDriver.h's own documented scope limitation).
//
// IMPORTANT, empirically established finding (SingleHitZeroMaterial
// CovarianceDivergenceCharacterization below): even at exactly zero
// material, the already-existing, already-committed native
// refitHit<CylinderCylinder>/barrel::update() primitive and the frozen
// legacy TrackParametrizationWithError::update() are independently
// implemented, mathematically-equivalent-in-theory Kalman update
// formulations that do NOT agree bit-for-bit, or even to simple floating-
// point-rounding precision -- a single hit already shows a ~1.6% difference
// in the fitted Y-Y covariance element, which compounds across a multi-hit
// leg into a visible (~0.2%) difference in final chi2/parameters. This is
// not a Slice B defect (verified: a single driveRefitLeg<CylinderCylinder>
// leg run in isolation reproduces the exact same divergence a full
// nativeRefitTrackCylinderCylinder call does, so Slice B's own leg-
// sequencing/orchestration code adds no additional divergence beyond what
// the pre-existing per-hit primitive already exhibits) and it is not this
// slice's place to alter refitHit/barrel::update, which are separately
// owned, already-tested, already-accepted primitives. Consequently, the
// tests below assert *structural* parity only (success/failure outcome,
// cluster pattern, timestamp -- all of which are unaffected by the update
// formula and are expected, and shown, to match exactly) and report
// resulting parameter/covariance/chi2 differences as characterization, not
// as a bit-identical or tolerance-bounded pass/fail claim -- exactly the
// same treatment already applied to the nonzero-material case below, and
// for the same reason: asserting a tolerance here would misrepresent an
// unevaluated (and, on this evidence, non-trivial) numerical difference
// between two independent implementations as validated parity.

#define BOOST_TEST_MODULE ITSMFTNativeCylinderCylinderRefitDriver
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <cstring>
#include <limits>
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

// The frozen legacy refit chain reads o2::base::Propagator::Instance().
// Constructing it with uninitialized=true (before any other call) skips its
// TGeoGlobalMagField/FairRunAna field lookup, which has neither in a plain
// unit-test process; every propagateToX call in this file uses the
// explicit-bz overload, so the singleton's own (absent) field state is
// never read.
struct PropagatorUninitializedFixture {
  PropagatorUninitializedFixture() { o2::base::Propagator::Instance(true); }
};

} // namespace

BOOST_TEST_GLOBAL_FIXTURE(PropagatorUninitializedFixture);

namespace
{

constexpr int NLayers = ITSNLayers;
constexpr float Bz = 5.f;
// Same local x ("radius") and alpha for every layer -- a degenerate but
// valid barrel fixture, matching the established convention already used by
// testDriveRefitLeg.cxx's own LegTraits<CylinderCylinder> (frame.q/
// frameAngle identical for every SurfaceId; only cluster/surface identity
// and the measured u/v coordinates vary per layer).
constexpr float LocalX = 2.5f;
constexpr float Alpha = 0.3f;

struct LegacyFixtureLayer {
  o2::its::Cluster cluster;
  o2::its::TrackingFrameInfo tfInfo;
};

// One shared numeric description of "NLayers barrel measurements", built
// independently into both a legacy (Cluster/TrackingFrameInfo) and a native
// (SurfaceMeasurement/SurfaceDescriptor) representation.
struct SharedFixture {
  std::array<LegacyFixtureLayer, NLayers> legacyLayers{};
  std::array<std::vector<SurfaceMeasurement>, NLayers> nativeStorage{};
  LayerMeasurementSpans<NLayers> layerMeasurements{};
  std::vector<SurfaceDescriptor> catalogSurfaces;
  SurfaceCatalogView catalog{};
  float xOverX0{0.f};
  float arealDensityGPerCm2{0.f};

  explicit SharedFixture(float materialXOverX0 = 0.f, float materialArealDensity = 0.f, int nHitLayers = NLayers)
    : xOverX0(materialXOverX0), arealDensityGPerCm2(materialArealDensity)
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
        // Legacy: Cluster (global xyz -- unused by the linRef fitTrack
        // branch this fixture exercises, filled for completeness/API
        // validity only) + TrackingFrameInfo.
        legacyLayers[layer].cluster = o2::its::Cluster(0.f, 0.f, 0.f, layer);
        legacyLayers[layer].tfInfo = o2::its::TrackingFrameInfo(0.f, 0.f, 0.f, LocalX, Alpha, {u, v}, {covUU, covUV, covVV});

        // Native: matching SurfaceMeasurement.
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
        nativeStorage[layer] = {}; // hole: empty span, seed leaves this layer unattached
      }
      layerMeasurements[layer] = gsl::span<const SurfaceMeasurement>(nativeStorage[layer]);
    }
    catalog = SurfaceCatalogView{catalogSurfaces.data(), static_cast<uint32_t>(catalogSurfaces.size())};
  }

  const o2::its::TrackingFrameInfo* tfInfoPtrs[NLayers];
  const o2::its::Cluster* clusterPtrs[NLayers];

  void fillLegacyPointerArrays()
  {
    for (int layer = 0; layer < NLayers; ++layer) {
      tfInfoPtrs[layer] = &legacyLayers[layer].tfInfo;
      clusterPtrs[layer] = &legacyLayers[layer].cluster;
    }
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

TrackSeedN<NLayers> makeNativeSeed(int nHitLayers = NLayers)
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

o2::its::TrackSeed<NLayers> makeLegacySeed(int nHitLayers = NLayers)
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

struct LegacyOracleResult {
  bool ok{false};
  o2::its::TrackITSExt track{};
};

LegacyOracleResult runLegacyOracle(SharedFixture& fixture, int nHitLayers, float maxChi2ClusterAttachment,
                                   float maxChi2NDF, bool shiftRefToCluster, bool repeatRefitOut,
                                   const std::vector<float>& minPt, o2::base::PropagatorF::MatCorrType matCorrType)
{
  fixture.fillLegacyPointerArrays();
  const auto legacySeed = makeLegacySeed(nHitLayers);
  std::vector<float> layerxX0(NLayers, fixture.xOverX0);
  std::vector<float> layerRadii(NLayers, LocalX);

  const o2::its::track::TrackFitContext<NLayers> fitCtx{
    fixture.tfInfoPtrs, layerxX0.data(), NLayers, Bz,
    maxChi2ClusterAttachment, maxChi2NDF,
    o2::base::Propagator::Instance(), matCorrType, shiftRefToCluster, repeatRefitOut};

  o2::its::TrackITSInternal<NLayers> internalTrack;
  LegacyOracleResult result{};
  result.ok = o2::its::track::refitTrackSeed<NLayers>(legacySeed, internalTrack, fitCtx, fixture.clusterPtrs,
                                                      layerRadii.data(), minPt.data(), /*reseedIfShorter=*/0);
  if (result.ok) {
    result.track = o2::its::makeTrackITSExt(internalTrack);
  }
  return result;
}

struct NativeDriverResult {
  bool ok{false};
  OperationFailureReason reason{};
  o2::its::TrackITSExt track{};
};

NativeDriverResult runNativeDriver(SharedFixture& fixture, int nHitLayers, float maxChi2ClusterAttachment,
                                   float maxChi2NDF, bool shiftRefToCluster, bool repeatRefitOut,
                                   const std::vector<float>& minPt)
{
  const auto nativeSeed = makeNativeSeed(nHitLayers);
  NativeDriverResult result{};
  SurfaceKinematicState paramIn{}, paramOut{};
  float chi2{0.f};
  result.ok = nativeRefitTrackCylinderCylinder<NLayers>(
    nativeSeed, fixture.layerMeasurements, fixture.catalog, Bz, shiftRefToCluster, maxChi2ClusterAttachment,
    maxChi2NDF, repeatRefitOut, gsl::span<const float>(minPt), paramIn, paramOut, chi2, result.reason);
  if (result.ok) {
    BOOST_REQUIRE(exportNativeRefitToTrackITSExt<NLayers>(nativeSeed, paramIn, paramOut, chi2, result.track));
  }
  return result;
}

// Structural parity: outcome, cluster pattern, and timestamp are all
// refit-engine-independent (unaffected by which Kalman update formula ran)
// and are asserted exactly. Parameter/covariance/chi2 values are reported
// via BOOST_TEST_MESSAGE as characterization, not asserted -- see this
// file's header comment for why (SingleHitZeroMaterialCovarianceDivergence
// CharacterizationSameEngineFinding below is the evidence).
void expectStructuralParityAndReportNumericDifferences(SharedFixture& fixture, int nHitLayers,
                                                       float maxChi2ClusterAttachment, float maxChi2NDF,
                                                       bool shiftRefToCluster, bool repeatRefitOut,
                                                       const std::vector<float>& minPt)
{
  const auto legacy = runLegacyOracle(fixture, nHitLayers, maxChi2ClusterAttachment, maxChi2NDF, shiftRefToCluster,
                                      repeatRefitOut, minPt, o2::base::PropagatorF::MatCorrType::USEMatCorrNONE);
  const auto native = runNativeDriver(fixture, nHitLayers, maxChi2ClusterAttachment, maxChi2NDF, shiftRefToCluster,
                                      repeatRefitOut, minPt);

  BOOST_REQUIRE(legacy.ok);
  BOOST_REQUIRE(native.ok);
  BOOST_CHECK_EQUAL(legacy.track.getPattern(), native.track.getPattern());
  BOOST_CHECK_EQUAL(legacy.track.getTimeStamp().getTimeStamp(), native.track.getTimeStamp().getTimeStamp());

  const float legacyQ2Pt = legacy.track.getParamIn().getQ2Pt();
  const float nativeQ2Pt = native.track.getParamIn().getQ2Pt();
  const float legacyChi2 = legacy.track.getChi2();
  const float nativeChi2 = native.track.getChi2();
  BOOST_TEST_MESSAGE("Zero-material numeric characterization: legacy paramIn.Q2Pt=" << legacyQ2Pt
                                                                                    << " native paramIn.Q2Pt=" << nativeQ2Pt << " |delta|=" << std::abs(legacyQ2Pt - nativeQ2Pt)
                                                                                    << "; legacy chi2=" << legacyChi2 << " native chi2=" << nativeChi2
                                                                                    << " |delta|=" << std::abs(legacyChi2 - nativeChi2));
}

} // namespace

// Empirically established (see this file's header comment): a single
// driveRefitLeg<CylinderCylinder> hit at exactly zero material already
// diverges from the frozen legacy fitTrack's single-hit result by ~1.6% in
// the fitted Y-Y covariance element, even though both engines' material
// step is independently confirmed a no-op at {xOverX0=0,
// arealDensityGPerCm2=0} (FamilyMaterialOperations.cxx's isNoopMaterial
// early-return; TrackParametrizationWithError::correctForMaterial's
// `xrho!=0`/`x2x0!=0` guards). The divergence is therefore in the Kalman
// update/rotate/propagate arithmetic itself, not in material handling --
// two independently-implemented, mathematically-equivalent-in-theory
// formulations that do not agree bit-for-bit. This test pins the observed
// values as a citable record; it intentionally does not assert a bound on
// the difference (see header comment).
BOOST_AUTO_TEST_CASE(SingleHitZeroMaterialCovarianceDivergenceCharacterization)
{
  SharedFixture fixture(0.f, 0.f, 1); // only layer 0 hit
  fixture.fillLegacyPointerArrays();

  auto legacySeed = makeLegacySeed(1);
  o2::its::TrackITSInternal<NLayers> legacyInternal;
  legacyInternal.paramIn = static_cast<const o2::track::TrackParCov&>(legacySeed);
  for (int l = 0; l < NLayers; ++l) {
    legacyInternal.setClusterIndex(l, legacySeed.getCluster(l));
  }
  o2::track::TrackPar legacyLinRef{legacyInternal.paramIn};
  o2::its::track::resetTrackCovariance(legacyInternal.paramIn);
  legacyInternal.setChi2(0.f);
  std::vector<float> layerxX0(NLayers, 0.f);
  const o2::its::track::TrackFitContext<NLayers> fitCtx{
    fixture.tfInfoPtrs, layerxX0.data(), NLayers, Bz, 60.f, 30.f,
    o2::base::Propagator::Instance(), o2::base::PropagatorF::MatCorrType::USEMatCorrNONE, false, false};
  // A lone hit always fails fitTrack's own trailing nCl>=3-shaped NDF check
  // (2*1-5 < 0), matching driveRefitLeg's own contract of never asserting
  // that check itself -- the per-hit rotate/propagate/chi2-gate/update
  // arithmetic already ran and mutated legacyInternal by this point
  // regardless (legacy fitTrack is not transactional), which is exactly
  // what this test reads.
  o2::its::track::fitTrack(legacyInternal, legacyInternal.paramIn, 0, 1, 1, o2::constants::math::VeryBig, 0, fitCtx,
                           &legacyLinRef);

  auto nativeSeed = makeNativeSeed(1);
  SurfaceKinematicState stateA = nativeSeed.state();
  SurfaceLinearizationReference linRefA{};
  BOOST_REQUIRE(makeLinearizationReference(stateA, linRefA));
  resetCylinderCylinderCovarianceForRefit(stateA);
  float chi2A = 0.f;
  uint32_t acceptedA = 0;
  std::array<SurfaceMeasurement, NLayers> slotsBufferA{};
  const auto slotsA = assembleRefitLegSlots<NLayers>(nativeSeed, fixture.layerMeasurements, 0, 1, 1, slotsBufferA);
  OperationFailureReason reason{};
  BOOST_REQUIRE(driveRefitLeg<TransitionPolicyTag::CylinderCylinder>(
    stateA, linRefA, chi2A, acceptedA, slotsA, fixture.catalog, Bz,
    material::MaterialTraversalDirection::AlongMomentum, false, 60.f, reason));

  BOOST_TEST_MESSAGE("legacy 1-hit chi2=" << legacyInternal.getChi2() << " y=" << legacyInternal.paramIn.getY()
                                          << " cov[YY]=" << legacyInternal.paramIn.getCov()[0]);
  BOOST_TEST_MESSAGE("native 1-hit chi2=" << chi2A << " y=" << stateA.parameters[0]
                                          << " cov[YY]=" << stateA.covariance[0]);
  // Both must at least be finite and physically sane (no crash, no NaN) --
  // that much genuinely is asserted, distinct from a numeric-equality claim.
  BOOST_CHECK(std::isfinite(chi2A));
  BOOST_CHECK(std::isfinite(stateA.parameters[0]));
  BOOST_CHECK(std::isfinite(legacyInternal.getChi2()));
  BOOST_CHECK(std::isfinite(legacyInternal.paramIn.getY()));
}

// --- Zero-material structural parity (see header comment for what "parity"
// means here), ShiftRefToCluster off, RepeatRefitOut off ---

BOOST_AUTO_TEST_CASE(ZeroMaterialStructuralParityAllHitsNoShiftNoRepeat)
{
  SharedFixture fixture(0.f, 0.f, NLayers);
  expectStructuralParityAndReportNumericDifferences(fixture, NLayers, /*maxChi2ClusterAttachment=*/60.f,
                                                    /*maxChi2NDF=*/30.f, /*shiftRefToCluster=*/false,
                                                    /*repeatRefitOut=*/false, std::vector<float>(NLayers + 1, 0.f));
}

BOOST_AUTO_TEST_CASE(ZeroMaterialStructuralParityShiftReferenceOn)
{
  SharedFixture fixture(0.f, 0.f, NLayers);
  expectStructuralParityAndReportNumericDifferences(fixture, NLayers, 60.f, 30.f, /*shiftRefToCluster=*/true, false,
                                                    std::vector<float>(NLayers + 1, 0.f));
}

BOOST_AUTO_TEST_CASE(ZeroMaterialStructuralParityRepeatRefitOutSucceeds)
{
  SharedFixture fixture(0.f, 0.f, NLayers);
  expectStructuralParityAndReportNumericDifferences(fixture, NLayers, 60.f, 30.f, false, /*repeatRefitOut=*/true,
                                                    std::vector<float>(NLayers + 1, 0.f));
}

BOOST_AUTO_TEST_CASE(ZeroMaterialStructuralParityWithHoles)
{
  // Only 5 of 7 layers hit (layers 5, 6 are holes on both sides).
  SharedFixture fixture(0.f, 0.f, 5);
  expectStructuralParityAndReportNumericDifferences(fixture, 5, 60.f, 30.f, false, false,
                                                    std::vector<float>(NLayers + 1, 0.f));
}

// --- Failure/rollback parity ------------------------------------------------

BOOST_AUTO_TEST_CASE(BothFailIdenticallyWhenChi2GateRejectsFourthHit)
{
  // A tiny maxChi2ClusterAttachment rejects hit 4 (nCl>=3 gate) on both
  // sides identically -- leg A fails on both, so neither ever reaches leg B.
  SharedFixture fixture(0.f, 0.f, NLayers);
  const auto legacy = runLegacyOracle(fixture, NLayers, /*maxChi2ClusterAttachment=*/1.e-8f, 30.f, false, false,
                                      std::vector<float>(NLayers + 1, 0.f), o2::base::PropagatorF::MatCorrType::USEMatCorrNONE);
  const auto native = runNativeDriver(fixture, NLayers, 1.e-8f, 30.f, false, false, std::vector<float>(NLayers + 1, 0.f));
  BOOST_CHECK(!legacy.ok);
  BOOST_CHECK(!native.ok);
  BOOST_CHECK(native.reason == OperationFailureReason::PredictedChi2Failure);
}

BOOST_AUTO_TEST_CASE(MinPtFailureRejectsAfterLegB)
{
  // A MinPt threshold far above the fixture's actual fitted Pt must reject
  // after leg B, matching refitTrackSeed's own trailing MinPt check.
  SharedFixture fixture(0.f, 0.f, NLayers);
  std::vector<float> minPt(NLayers + 1, 0.f);
  minPt[NLayers - NLayers] = 1.e6f; // params.MinPt[NLayers - nClAttached], nClAttached == NLayers here

  const auto legacy = runLegacyOracle(fixture, NLayers, 60.f, 30.f, false, false, minPt,
                                      o2::base::PropagatorF::MatCorrType::USEMatCorrNONE);
  const auto native = runNativeDriver(fixture, NLayers, 60.f, 30.f, false, false, minPt);
  BOOST_CHECK(!legacy.ok);
  BOOST_CHECK(!native.ok);
  BOOST_CHECK(native.reason == OperationFailureReason::MinPtFailure);
}

// --- Export-failure handling (native driver's own contract, no legacy side) ---

BOOST_AUTO_TEST_CASE(ExportFailsClosedForWrongFamilyStateAndLeavesOutputUntouched)
{
  TrackSeedN<NLayers> seed{};
  SurfaceKinematicState wrongFamily{};
  wrongFamily.family = StateFamily::Forward; // never produced by this driver; precondition violation
  o2::its::TrackITSExt sentinel{};
  sentinel.setChi2(-999.f); // recognizable poison value
  auto outTrack = sentinel;

  const bool ok = exportNativeRefitToTrackITSExt<NLayers>(seed, wrongFamily, wrongFamily, 0.f, outTrack);
  BOOST_CHECK(!ok);
  BOOST_CHECK_EQUAL(outTrack.getChi2(), -999.f); // untouched on failure
}

// --- Nonzero-material characterization (NOT a parity claim) -----------------

BOOST_AUTO_TEST_CASE(NonzeroMaterialCharacterizationRecorded)
{
  // The native refitHit<CylinderCylinder> kernel applies a PID/absCharge-
  // aware, direction-signed material correction; the frozen legacy fitTrack
  // applies a fixed-positive-sign, non-PID-aware formula
  // (correctForMaterial(*linRef, xOverX0, xOverX0*Radl*Rho, true)) --
  // TransitionPolicyOperations.h documents this as a deliberate, already-
  // accepted divergence, not a bug. This test records the resulting
  // numerical difference for a realistic ITS-scale material budget as
  // characterization only: it asserts both refits succeed, but intentionally
  // does NOT assert any bound on the difference between them -- there is no
  // agreed tolerance, and picking one here would misrepresent an
  // unevaluated physics difference as validated parity.
  SharedFixture fixture(0.01f, 0.001f, NLayers); // O(ITS Si layer) xOverX0
  const auto legacy = runLegacyOracle(fixture, NLayers, 60.f, 30.f, false, false, std::vector<float>(NLayers + 1, 0.f),
                                      o2::base::PropagatorF::MatCorrType::USEMatCorrNONE);
  const auto native = runNativeDriver(fixture, NLayers, 60.f, 30.f, false, false, std::vector<float>(NLayers + 1, 0.f));

  BOOST_REQUIRE(legacy.ok);
  BOOST_REQUIRE(native.ok);

  const float legacyQ2Pt = legacy.track.getParamIn().getQ2Pt();
  const float nativeQ2Pt = native.track.getParamIn().getQ2Pt();
  const float legacyChi2 = legacy.track.getChi2();
  const float nativeChi2 = native.track.getChi2();
  BOOST_TEST_MESSAGE("Nonzero-material characterization (xOverX0=0.01, arealDensityGPerCm2=0.001): "
                     << "legacy paramIn.Q2Pt=" << legacyQ2Pt << " native paramIn.Q2Pt=" << nativeQ2Pt
                     << " |delta|=" << std::abs(legacyQ2Pt - nativeQ2Pt) << "; legacy chi2=" << legacyChi2
                     << " native chi2=" << nativeChi2 << " |delta|=" << std::abs(legacyChi2 - nativeChi2));
  // No BOOST_CHECK on the deltas above: characterization only, per task
  // instruction not to select an arbitrary tolerance or call it parity.
}
