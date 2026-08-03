// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// M2b (GenericTrackingEngineMigration.md; ADR 0007) focused tests for
// MultiSourceTimeFrameLoader's participant-count-generic atomic loading
// transaction (loadEvent()) and its fixed-position ITS/MFT compatibility
// wrapper (loadITSAndMFT()). Covers:
//  - a synthetic three-participant transaction, proving loadEvent() has no
//    hidden two-source limit and no source-0/1 branch (it is exercised here
//    with three fake LoadTargets that know nothing about ITS/MFT);
//  - a failure in the *last* participant's staged backfill leaves the
//    normalized frame and every earlier participant's target completely
//    uncommitted -- the loader-level foundation the composition/participant
//    level tests (testCombinedTrackingComposition.cxx) already build their
//    own sidecar/scratch/publication-state proofs on top of;
//  - loadITSAndMFT() still reproduces source-qualified per-surface storage
//    and correct compact per-layer backfill on both real scratches;
//  - a source-level guard: MultiSourceTimeFrameLoader's own generic parts
//    (LoadTarget, LoadTargetImpl, AtomicLoadBinding, loadEvent()) contain no
//    ITS/MFT/fixed-source-0-or-1 code, only loadITSAndMFT()/
//    resetITSAndMFTEvent() (explicitly exempted) do.

#define BOOST_TEST_MODULE ITSMFT MultiSourceTimeFrameLoader
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <fstream>
#include <string>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/DecodedCluster.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/LegacyTrackerScratch.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/SurfaceMeasurementAdapters.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

// Same construction as testMultiSourceLoading.cxx's FakeClusterDecoder,
// narrowed to what this file's fixtures need: sensorID == detector-local
// layer, one measurement per cluster, no corruption modes.
class FakeClusterDecoder final : public ClusterDecoder
{
 public:
  explicit FakeClusterDecoder(o2::detectors::DetID::ID detector, SurfaceKind kind = SurfaceKind::Cylinder)
    : mDetector(detector), mKind(kind)
  {
  }

  o2::itsmft::ioutils::SurfaceMeasurementDecodeResult decode(
    const CompClusterExt& cluster,
    BoundedPatternCursor& patterns,
    const TopologyDictionary* dict,
    gsl::span<const SurfaceId> layerToSurface,
    ClusterSourceId source,
    uint32_t externalIndex,
    uint32_t sourceROF,
    bool) const override
  {
    const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dict);
    if (!clusterData.ok()) {
      o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
      result.error = clusterData.error;
      return result;
    }

    o2::itsmft::ioutils::SurfaceMeasurementDecodeResult result;
    const auto sensorID = cluster.getSensorID();
    const int layer = sensorID;
    result.layer = layer;
    if (layer < 0 || static_cast<size_t>(layer) >= layerToSurface.size()) {
      return result;
    }
    result.layerMapped = true;
    result.kind = mKind;

    DecodedCluster decoded{};
    decoded.global = {static_cast<float>(sensorID), static_cast<float>(cluster.getRow()), static_cast<float>(cluster.getCol())};
    decoded.cylinderFrame = {10.f + sensorID, 1.f, 2.f, 0.1f};
    decoded.rowColumnCovariance = {clusterData.sig2Row, 0.f, clusterData.sig2Col};
    decoded.shape = clusterData.shape;
    decoded.sensor = static_cast<uint32_t>(sensorID);
    decoded.layer = layer;

    const auto surface = layerToSurface[layer];
    const DetectorSensorId sensor{static_cast<uint32_t>(mDetector), decoded.sensor};
    const ClusterRef clusterRef{source, externalIndex};
    result.measurement = mKind == SurfaceKind::Disk
                           ? makeDiskSurfaceMeasurement(decoded, sensor, surface, clusterRef, sourceROF)
                           : makeCylinderSurfaceMeasurement(decoded, sensor, surface, clusterRef, sourceROF);
    return result;
  }

 private:
  o2::detectors::DetID::ID mDetector;
  SurfaceKind mKind;
};

struct BuiltLayout {
  DetectorLayout layout;
  std::vector<SurfaceDescriptor> surfaces;

  DetectorLayoutView getView() const noexcept
  {
    const auto masks = computeSurfaceKindMasks(surfaces);
    return layout.getView(surfaces, masks.first, masks.second);
  }
};

// One surface per synthetic source, all tagged the same detector (the
// FakeLoadTarget-based tests below never inspect detector identity, so a
// single tag keeps every source's own decoder/catalog consistent).
BuiltLayout makeThreeSurfaceLayout()
{
  SparseTrackingTopology topology{3};
  topology.finalize();
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.push_back(SurfaceDescriptor{SurfaceId{0}, 0, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{1}, 0, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  surfaces.push_back(SurfaceDescriptor{SurfaceId{2}, 0, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder});
  return BuiltLayout{DetectorLayout{surfaces, std::move(topology)}, std::move(surfaces)};
}

constexpr std::array<unsigned char, 3> onePixelPattern{1, 1, 0x80};

std::vector<unsigned char> makePatternBytes(size_t nClusters)
{
  std::vector<unsigned char> bytes;
  bytes.reserve(nClusters * onePixelPattern.size());
  for (size_t i = 0; i < nClusters; ++i) {
    bytes.insert(bytes.end(), onePixelPattern.begin(), onePixelPattern.end());
  }
  return bytes;
}

const TopologyDictionary& dict()
{
  static const TopologyDictionary d;
  return d;
}

// One single-cluster, single-ROF, single-surface source, `surfaceIndex`
// selecting which of the layout's own SurfaceIds it targets. `detector`/
// `kind` must agree with that SurfaceId's own catalog tag.
struct OneClusterSource {
  std::vector<CompClusterExt> clusters{{0, 1, CompCluster::InvalidPatternID, 0}};
  std::vector<unsigned char> patterns{makePatternBytes(1)};
  std::vector<ROFRecord> rofs{ROFRecord{{100, 5}, 0, 0, 1}};
  std::array<SurfaceId, 1> layerToSurface;
  o2::detectors::DetID::ID detector;
  FakeClusterDecoder decoder;

  explicit OneClusterSource(uint16_t surfaceIndex, o2::detectors::DetID::ID det = o2::detectors::DetID::ITS,
                            SurfaceKind kind = SurfaceKind::Cylinder)
    : layerToSurface{SurfaceId{surfaceIndex}}, detector(det), decoder(det, kind)
  {
  }

  ClusterSourceInput asInput(ClusterSourceId id) const
  {
    ClusterSourceInput input{};
    input.id = id;
    input.detector = detector;
    input.clusters = clusters;
    input.patterns = patterns;
    input.rofs = rofs;
    input.dictionary = &dict();
    input.layerToSurface = layerToSurface;
    input.timing = ROFTimingConfig{40, 0, 0, 0};
    input.decoder = &decoder;
    return input;
  }
};

// A LoadTarget test double that knows nothing about ITS/MFT or any legacy
// scratch: stage()/commit() are pure call/outcome recorders. Proves
// loadEvent()'s own generic machinery (staging order, all-succeed-then-
// commit-all, stop-and-leave-everything-untouched-on-first-failure) without
// needing a real LegacyTrackerScratch<NLayers> behind every binding.
class FakeLoadTarget final : public MultiSourceTimeFrameLoader::LoadTarget
{
 public:
  explicit FakeLoadTarget(bool stageSucceeds) : mStageSucceeds(stageSucceeds) {}

  LoadSourcesResult stage(const ClusterSourceInput& source, SurfaceCatalogView, const o2::InteractionRecord&) override
  {
    ++mStageCalls;
    if (!mStageSucceeds) {
      return {MultiSourceLoadError::OtherMalformedInput, source.id};
    }
    return {};
  }

  void commit() noexcept override { ++mCommitCalls; }

  int stageCalls() const noexcept { return mStageCalls; }
  int commitCalls() const noexcept { return mCommitCalls; }

 private:
  bool mStageSucceeds;
  int mStageCalls{0};
  int mCommitCalls{0};
};

} // namespace

BOOST_AUTO_TEST_CASE(ThreeParticipantSyntheticTransactionProvesNoHiddenTwoSourceLimit)
{
  // Three sources, three independent fake targets -- if loadEvent() had a
  // residual two-source assumption (e.g. only ever staging/committing
  // bindings[0]/bindings[1]), the third target's call counts below would
  // stay at zero even though the whole transaction reports success.
  const auto layout = makeThreeSurfaceLayout();
  const OneClusterSource sourceA{0};
  const OneClusterSource sourceB{1};
  const OneClusterSource sourceC{2};

  FakeLoadTarget targetA{true};
  FakeLoadTarget targetB{true};
  FakeLoadTarget targetC{true};
  const std::array<MultiSourceTimeFrameLoader::AtomicLoadBinding, 3> bindings{
    MultiSourceTimeFrameLoader::AtomicLoadBinding{sourceA.asInput(ClusterSourceId{0}), targetA},
    MultiSourceTimeFrameLoader::AtomicLoadBinding{sourceB.asInput(ClusterSourceId{1}), targetB},
    MultiSourceTimeFrameLoader::AtomicLoadBinding{sourceC.asInput(ClusterSourceId{2}), targetC}};

  TimeFrame frame;
  const auto result = MultiSourceTimeFrameLoader::loadEvent(
    frame, gsl::span<const MultiSourceTimeFrameLoader::AtomicLoadBinding>{bindings}, layout.getView().getSurfaceCatalogView(), {50, 5});
  BOOST_REQUIRE(result.ok());

  BOOST_CHECK_EQUAL(targetA.stageCalls(), 1);
  BOOST_CHECK_EQUAL(targetB.stageCalls(), 1);
  BOOST_CHECK_EQUAL(targetC.stageCalls(), 1);
  BOOST_CHECK_EQUAL(targetA.commitCalls(), 1);
  BOOST_CHECK_EQUAL(targetB.commitCalls(), 1);
  BOOST_CHECK_EQUAL(targetC.commitCalls(), 1);

  // The normalized frame itself carries all three sources' measurements --
  // loadSources() was already source-count-agnostic before this milestone,
  // and this is loadEvent()'s own proof that it drives it over the whole
  // binding list, not just a fixed leading pair.
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{0}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{1}).size(), 1u);
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{2}).size(), 1u);
}

BOOST_AUTO_TEST_CASE(FailureInLastParticipantsStagedBackfillLeavesEverythingAtBaseline)
{
  // First two targets' stage() would succeed; the third's fails. Since
  // loadEvent() only commits after *every* binding's stage() has already
  // returned ok(), nothing must be committed at all: not the normalized
  // frame, and not targetA/targetB even though their own stage() calls
  // individually succeeded before the transaction as a whole failed. This
  // is the loader-level foundation the composition-level tests
  // (testCombinedTrackingComposition.cxx's LoadFailureResetsWholeCombined
  // TFExactlyOnceAndInvalidatesPublication/CompatibilitySidecarGettersReflec
  // tSealAndReset) build their own real-sidecar/real-scratch/publication-
  // state proofs on top of: a participant's scratch/sidecar can only change
  // via its own commit()/track(), and neither runs here.
  const auto layout = makeThreeSurfaceLayout();
  const OneClusterSource sourceA{0};
  const OneClusterSource sourceB{1};
  const OneClusterSource sourceC{2};

  FakeLoadTarget targetA{true};
  FakeLoadTarget targetB{true};
  FakeLoadTarget targetC{false};
  const std::array<MultiSourceTimeFrameLoader::AtomicLoadBinding, 3> bindings{
    MultiSourceTimeFrameLoader::AtomicLoadBinding{sourceA.asInput(ClusterSourceId{0}), targetA},
    MultiSourceTimeFrameLoader::AtomicLoadBinding{sourceB.asInput(ClusterSourceId{1}), targetB},
    MultiSourceTimeFrameLoader::AtomicLoadBinding{sourceC.asInput(ClusterSourceId{2}), targetC}};

  TimeFrame frame;
  // Deliberately populate the frame's shared common-result storage before
  // the failing call, as a baseline sentinel: loadEvent() only ever touches
  // the normalized frame (commitNormalizedFrame()), never CommonTrack/
  // TrackClusterIndices storage, so this must survive untouched too.
  frame.getTrackClusterIndices().push_back(TrackClusterReference{SurfaceId{0}, SurfaceMeasurementIndex{0}});

  const auto result = MultiSourceTimeFrameLoader::loadEvent(
    frame, gsl::span<const MultiSourceTimeFrameLoader::AtomicLoadBinding>{bindings}, layout.getView().getSurfaceCatalogView(), {50, 5});
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.source == ClusterSourceId{2});

  // Every target's stage() ran (staging order), but commit() never did --
  // on any target, including the two that individually staged successfully.
  BOOST_CHECK_EQUAL(targetA.stageCalls(), 1);
  BOOST_CHECK_EQUAL(targetB.stageCalls(), 1);
  BOOST_CHECK_EQUAL(targetC.stageCalls(), 1);
  BOOST_CHECK_EQUAL(targetA.commitCalls(), 0);
  BOOST_CHECK_EQUAL(targetB.commitCalls(), 0);
  BOOST_CHECK_EQUAL(targetC.commitCalls(), 0);

  // The normalized frame was never committed: still completely empty, not
  // partially populated from sourceA/sourceB's already-decoded measurements.
  BOOST_CHECK_EQUAL(frame.getNormalizedFrame().getTotalMeasurements(), 0u);
  // The pre-existing baseline sentinel is untouched.
  BOOST_REQUIRE_EQUAL(frame.getTrackClusterIndices().size(), 1u);
  BOOST_CHECK(frame.getTrackClusterIndices()[0].surface == SurfaceId{0});
}

namespace
{

std::vector<SurfaceId> orderedRange(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

// One cluster on layer 0 only; `layerToSurface` still spans every one of
// this detector's NLayers surfaces (layers 1..NLayers-1 stay empty) --
// LegacyTrackerScratch<NLayers>::loadNormalizedSource() sizes its own
// legacy per-layer containers from layerToSurface.size(), so it must match
// NLayers exactly, unlike the single-surface synthetic sources above (which
// never touch a real LegacyTrackerScratch<NLayers> at all).
ClusterSourceInput makeSingleClusterInput(ClusterSourceId id, o2::detectors::DetID::ID det, SurfaceKind kind,
                                          const std::vector<SurfaceId>& layerToSurface,
                                          std::vector<CompClusterExt>& clustersOut, std::vector<unsigned char>& patternsOut,
                                          std::vector<ROFRecord>& rofsOut, FakeClusterDecoder& decoder)
{
  clustersOut = {{0, 1, CompCluster::InvalidPatternID, 0}};
  patternsOut = makePatternBytes(1);
  rofsOut = {ROFRecord{{100, 5}, 0, 0, 1}};

  ClusterSourceInput input{};
  input.id = id;
  input.detector = det;
  input.clusters = clustersOut;
  input.patterns = patternsOut;
  input.rofs = rofsOut;
  input.dictionary = &dict();
  input.layerToSurface = layerToSurface;
  input.timing = ROFTimingConfig{40, 0, 0, 0};
  input.decoder = &decoder;
  return input;
}

} // namespace

BOOST_AUTO_TEST_CASE(TwoParticipantITSMFTWrapperReproducesSourceQualificationAndCompactBackfills)
{
  // loadITSAndMFT() itself, directly -- not through any application-layer
  // composition -- still backfills both real scratches (the actual production
  // combined 17-surface ITS+MFT static catalog, ITS global ids 0..6, MFT
  // 7..16) and preserves source-qualified per-surface storage, now that it
  // is a thin wrapper over loadEvent().
  const auto catalogView = SurfaceCatalogView{kITSMFTCombinedStaticSurfaceCatalog.data(),
                                              static_cast<uint32_t>(kITSMFTCombinedStaticSurfaceCatalog.size())};
  const auto itsLayerToSurface = orderedRange(0, ITSNLayers);
  const auto mftLayerToSurface = orderedRange(ITSNLayers, MFTNLayers);

  TimeFrame frame;
  LegacyTrackerScratchITS itsScratch;
  LegacyTrackerScratchMFT mftScratch;

  FakeClusterDecoder itsDecoder{o2::detectors::DetID::ITS, SurfaceKind::Cylinder};
  FakeClusterDecoder mftDecoder{o2::detectors::DetID::MFT, SurfaceKind::Disk};
  std::vector<CompClusterExt> itsClusters, mftClusters;
  std::vector<unsigned char> itsPatterns, mftPatterns;
  std::vector<ROFRecord> itsRofs, mftRofs;
  const auto itsInput = makeSingleClusterInput(ClusterSourceId{0}, o2::detectors::DetID::ITS, SurfaceKind::Cylinder,
                                               itsLayerToSurface, itsClusters, itsPatterns, itsRofs, itsDecoder);
  const auto mftInput = makeSingleClusterInput(ClusterSourceId{1}, o2::detectors::DetID::MFT, SurfaceKind::Disk,
                                               mftLayerToSurface, mftClusters, mftPatterns, mftRofs, mftDecoder);

  const auto result = MultiSourceTimeFrameLoader::loadITSAndMFT(frame, itsScratch, mftScratch, itsInput, mftInput, catalogView, {50, 5});
  BOOST_REQUIRE(result.ok());

  // Compact per-layer backfill landed independently on each scratch.
  BOOST_CHECK_EQUAL(itsScratch.getTotalClusters(), 1);
  BOOST_CHECK_EQUAL(mftScratch.getTotalClusters(), 1);
  BOOST_CHECK_EQUAL(itsScratch.getNrof(0), 1);
  BOOST_CHECK_EQUAL(mftScratch.getNrof(0), 1);

  // Source qualification: surface 0 (ITS layer 0)'s one measurement is
  // ITS/source 0, surface 7 (MFT layer 0)'s is MFT/source 1.
  const auto onIts = frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{0});
  const auto onMft = frame.getNormalizedFrame().getSurfaceMeasurements(SurfaceId{ITSNLayers});
  BOOST_REQUIRE_EQUAL(onIts.size(), 1u);
  BOOST_REQUIRE_EQUAL(onMft.size(), 1u);
  BOOST_CHECK(onIts[0].sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::ITS));
  BOOST_CHECK(onIts[0].cluster.source == ClusterSourceId{0});
  BOOST_CHECK(onMft[0].sensor.detector == static_cast<uint32_t>(o2::detectors::DetID::MFT));
  BOOST_CHECK(onMft[0].cluster.source == ClusterSourceId{1});

  // The fixed-position guard itself: a mismatched id/detector pairing is
  // rejected before any staging happens, exactly as before this milestone.
  TimeFrame rejectedFrame;
  LegacyTrackerScratchITS rejectedItsScratch;
  LegacyTrackerScratchMFT rejectedMftScratch;
  auto wrongIdInput = itsInput;
  wrongIdInput.id = ClusterSourceId{5};
  const auto rejected = MultiSourceTimeFrameLoader::loadITSAndMFT(rejectedFrame, rejectedItsScratch, rejectedMftScratch, wrongIdInput, mftInput,
                                                                  catalogView, {50, 5});
  BOOST_CHECK(rejected.error == MultiSourceLoadError::UnsupportedDetector);
  BOOST_CHECK(rejected.source == ClusterSourceId{5});
  BOOST_CHECK_EQUAL(rejectedItsScratch.getTotalClusters(), 0);
  BOOST_CHECK_EQUAL(rejectedFrame.getNormalizedFrame().getTotalMeasurements(), 0u);
}

namespace
{

std::vector<std::string> readLines(const std::string& path)
{
  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path);
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(input, line)) {
    lines.push_back(line);
  }
  return lines;
}

bool isFullLineComment(const std::string& line)
{
  const auto pos = line.find_first_not_of(" \t");
  return pos != std::string::npos && line.compare(pos, 2, "//") == 0;
}

// Marks every line from the first occurrence of `marker` through the line
// where a running `openChar`/`closeChar` depth (counted from that same
// line onward) returns to zero as exempt -- robust to the exact multi-line
// formatting of the one ITS/MFT-specific declaration/definition each call
// below is used to skip, without needing to hand-maintain a line range.
void markExempt(const std::vector<std::string>& lines, const std::string& marker,
                char openChar, char closeChar, std::vector<bool>& exempt)
{
  for (size_t i = 0; i < lines.size(); ++i) {
    if (lines[i].find(marker) == std::string::npos) {
      continue;
    }
    int depth = 0;
    bool opened = false;
    for (size_t j = i; j < lines.size(); ++j) {
      exempt[j] = true;
      for (char c : lines[j]) {
        if (c == openChar) {
          ++depth;
          opened = true;
        } else if (c == closeChar) {
          --depth;
        }
      }
      if (opened && depth <= 0) {
        break;
      }
    }
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(GenericLoaderHeaderAndSourceHaveNoITSMFTOrFixedSourceSpecialization)
{
  // Source-level guard: MultiSourceTimeFrameLoader's generic parts --
  // LoadTarget, LoadTargetImpl, AtomicLoadBinding, loadEvent() -- must
  // contain no ITS/MFT-specific code and no ClusterSourceId{0}/{1} literal.
  // The one ITS/MFT-specific compatibility wrapper (loadITSAndMFT()/
  // resetITSAndMFTEvent(), plus the two explicit LoadTargetImpl<ITSNLayers>/
  // <MFTNLayers> instantiations they need) is explicitly exempted below --
  // it is expected, documented legacy-position knowledge, not a leak. Full-
  // line comments (this file's own doc comments above, and every doc
  // comment inside the production files, freely explain *why* ITS/MFT are
  // avoided) are skipped entirely: this is a code scan, not a prose scan.
  const std::string testFile = __FILE__;
  const auto testDirectory = testFile.substr(0, testFile.find_last_of('/'));
  const std::string headerFile = testDirectory + "/../include/ITSMFTTracking/MultiSourceTimeFrameLoader.h";
  const std::string sourceFile = testDirectory + "/../src/MultiSourceTimeFrameLoader.cxx";

  // Precise identifiers for real detector/layer-count/fixed-source
  // specialization -- not the bare substrings "ITS"/"MFT", which would
  // also match this library's own always-present namespace/header-guard/
  // include-path boilerplate ("ITSMFTTracking", "o2::itsmft::tracking",
  // "ALICEO2_ITSMFT_...") on every line of every file in this library.
  // "Source 0/1 specialization" in this codebase's own established
  // vocabulary (loadITSAndMFT()'s doc) means "ITS is always assigned
  // ClusterSourceId 0, MFT is always assigned ClusterSourceId 1" -- a
  // detector-position mapping, fully captured by scanning for the DetID
  // tokens below (loadITSAndMFT()'s own guard combines a source-id check
  // with a DetID check). A bare ClusterSourceId{0}/{1} alone is not
  // inherently detector-position specialization: LoadTargetImpl<NLayers>::
  // stage() legitimately numbers its own internal, single-element,
  // disposable-TimeFrame source as 0 (any value would do -- there is only
  // ever one source in that throwaway container), unrelated to this
  // transaction's real, externally supplied source ids.
  const std::vector<std::string> forbidden = {"DetID::ITS", "DetID::MFT", "LegacyTrackerScratchITS", "LegacyTrackerScratchMFT",
                                              "ITSNLayers", "MFTNLayers"};

  {
    auto lines = readLines(headerFile);
    std::vector<bool> exempt(lines.size(), false);
    markExempt(lines, "loadITSAndMFT(", '(', ')', exempt);
    markExempt(lines, "resetITSAndMFTEvent(", '(', ')', exempt);
    for (size_t i = 0; i < lines.size(); ++i) {
      if (exempt[i] || isFullLineComment(lines[i])) {
        continue;
      }
      if (lines[i].find("extern template class") != std::string::npos) {
        continue; // LoadTargetImpl<ITSNLayers>/<MFTNLayers> explicit declarations
      }
      for (const auto& token : forbidden) {
        BOOST_CHECK_MESSAGE(lines[i].find(token) == std::string::npos,
                            headerFile << ":" << (i + 1) << " contains forbidden token '" << token << "': " << lines[i]);
      }
    }
  }
  {
    auto lines = readLines(sourceFile);
    std::vector<bool> exempt(lines.size(), false);
    markExempt(lines, "MultiSourceTimeFrameLoader::loadITSAndMFT(", '{', '}', exempt);
    markExempt(lines, "MultiSourceTimeFrameLoader::resetITSAndMFTEvent(", '{', '}', exempt);
    for (size_t i = 0; i < lines.size(); ++i) {
      if (exempt[i] || isFullLineComment(lines[i])) {
        continue;
      }
      if (lines[i].find("template class MultiSourceTimeFrameLoader::LoadTargetImpl") != std::string::npos) {
        continue; // LoadTargetImpl<ITSNLayers>/<MFTNLayers> explicit instantiations
      }
      for (const auto& token : forbidden) {
        BOOST_CHECK_MESSAGE(lines[i].find(token) == std::string::npos,
                            sourceFile << ":" << (i + 1) << " contains forbidden token '" << token << "': " << lines[i]);
      }
    }
  }
}
