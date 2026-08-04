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
///
/// \file TrackingInterface.h
/// \brief Shared DPL-facing CA tracking interface for ITS and MFT
///

#ifndef ALICEO2_ITSMFT_TRACKING_TRACKINGINTERFACE_H_
#define ALICEO2_ITSMFT_TRACKING_TRACKINGINTERFACE_H_

#include <limits>
#include <memory>
#include <optional>
#include <vector>

#include <gsl/span>

#include "CommonDataFormat/IRFrame.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/DetectorTraits.h"
#ifndef GPUCA_GPUCODE
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/ClockTimingPublicationView.h"
#include "ITSMFTTracking/DetectorLayoutSet.h"
#include "ITSMFTTracking/TimeFrameLoadFailure.h"
#endif
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/LegacyTrackerScratch.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/ROFLookupTables.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

namespace o2::itsmft::tracking
{

// Host-only, immutable workflow boundary. Raw ROFRecord input remains owned
// by the workflow; this export never exposes scratch or overlap-table state.
struct CommonTrackPublicationExport {
  o2::detectors::DetID::ID detector{};
  ClusterSourceId source{};
  ClockTimingPublicationView clock;
  gsl::span<const SurfaceId> orderedSurfaces;
};

template <int NLayers>
struct MFTPublicationCompatibilityOwner {
};

template <>
struct MFTPublicationCompatibilityOwner<o2::mft::constants::mft::LayersNumber> {
  MFTPublicationCompatibility sidecar;
};

template <int NLayers>
struct ITSSharedClusterCompatibilityOwner {
};

template <>
struct ITSSharedClusterCompatibilityOwner<ITSNLayers> {
  ITSSharedClusterCompatibility sidecar;
};

// M6e1 (doc/design/0002-m6-generic-workspace-migration.md): the same
// compile-time-defaulted scratch/binding seam M6d added to
// Tracker<NLayers,...>/TrackerTraits<NLayers,...>/LegacyCATrackingParticipant
// <NLayers,...>, applied here so this class's own standalone-MFT
// instantiation can own SurfaceTrackingScratch/SurfacePlanBinding while its
// ITS instantiation (default template arguments) stays byte-for-byte on
// LegacyTrackerScratch<7>/DetectorTraversalBinding, unchanged.
template <int NLayers, typename ScratchT = LegacyTrackerScratch<NLayers>, typename BindingT = DetectorTraversalBinding>
class ITSMFTTrackingInterface : private MFTPublicationCompatibilityOwner<NLayers>, private ITSSharedClusterCompatibilityOwner<NLayers>
{
 public:
  static_assert(NLayers == ITSNLayers || NLayers == o2::mft::constants::mft::LayersNumber,
                "ITSMFTTrackingInterface supports ITS (7) and MFT (10) layer counts only");
  static constexpr o2::detectors::DetID::ID DetId = detIdFromNLayers<NLayers>();

  using ScratchN = ScratchT;
  using TrackerN = Tracker<NLayers, ScratchT, BindingT>;
  using TrackerTraitsN = TrackerTraits<NLayers, ScratchT, BindingT>;
  // Derived via ScratchN, not named directly: both LegacyTrackerScratch<NLayers>
  // and SurfaceTrackingScratch already expose these three nested aliases
  // (LegacyTrackerScratch<NLayers>'s own o2::its::ROFOverlapTable<NLayers>/
  // ROFVertexLookupTable<NLayers>/ROFMaskTable<NLayers>; SurfaceTrackingScratch's
  // own hardcoded-to-MFTNLayers equivalents) -- so this line needs no
  // NLayers-vs-ScratchT branch and is unchanged in effect for ITS.
  using ROFOverlapTableN = typename ScratchN::ROFOverlapTableN;
  using ROFVertexLookupTableN = typename ScratchN::ROFVertexLookupTableN;
  using ROFMaskTableN = typename ScratchN::ROFMaskTableN;
  using BoundedMemoryResourceN = BoundedMemoryResource;

  ITSMFTTrackingInterface(bool useMC, o2::itsmft::TrackingMode::Type mode, bool overrideBeamEst);
#ifndef GPUCA_GPUCODE
  // clusterDecoder == nullptr (the default) selects the production adapter
  // GeometryClusterDecoder<DetId>. Tests inject a fake ClusterDecoder here
  // instead of adding a virtual factory hook to the class: explicit,
  // constructor-owned lifetime, no dispatch added to the processing path.
  ITSMFTTrackingInterface(bool useMC, o2::itsmft::TrackingMode::Type mode, bool overrideBeamEst,
                          std::unique_ptr<ClusterDecoder> clusterDecoder);
#endif

  void setTrackingMode(o2::itsmft::TrackingMode::Type mode) { mTrackingMode = mode; }
  void setClusterDictionary(const o2::itsmft::TopologyDictionary* dict) { mDict = dict; }
  void setMeanVertex(const o2::dataformats::MeanVertexObject* v) { mMeanVertex = v; }

  void initialise();

  /// Phase 2 (TimeFrame load) + Phase 3 (CA tracking).
  ///
  /// - Success (load and tracking both completed): returns elapsed tracking
  ///   ms (>= 0). The loaded/tracked event is retained in the TimeFrame;
  ///   call clearTimeFrame() once outputs have been extracted.
  /// - Recoverable failure (malformed per-TF loading input, or a per-TF
  ///   resource-exhaustion exception) with DropTFUponFailure=true: the
  ///   TimeFrame is fully wiped first, then this returns exactly
  ///   kDroppedTimeFrameResult (see CATracker.h); nothing was published for
  ///   this TF. onTrackingFinished()/onTimeFrameLoaded() are not both
  ///   called in this case -- see loadTimeFrame()/runTracking() for which
  ///   ran.
  /// - Every other failure -- structural/configuration (bad or stale
  ///   catalog, non-uniform per-layer ROF timing, a missing dictionary, a
  ///   structural LoadSourcesResult), unclassified, or a recoverable
  ///   failure with DropTFUponFailure=false: the TimeFrame is fully wiped,
  ///   then the failure is rethrown (by its original type for an
  ///   unclassified exception; see TimeFrameLoadFailure.h for the loading
  ///   boundary's typed exceptions and CATracker.h for the CA-tracking
  ///   ones). No output is created on this path.
  float processTimeFrame(gsl::span<const o2::itsmft::ROFRecord> rofs,
                         gsl::span<const o2::itsmft::CompClusterExt> clusters,
                         gsl::span<const unsigned char> patterns,
                         const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels,
                         gsl::span<const o2::dataformats::IRFrame> irFrames = {});

  // Owner-level reset clears detector-local compatibility state with the
  // shared frame. Scratch-only reset deliberately does neither.
  void resetEvent()
  {
    mPublicationClock.reset();
    if constexpr (DetId == o2::detectors::DetID::ITS) {
      static_cast<ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
    }
    if constexpr (DetId == o2::detectors::DetID::MFT) {
      static_cast<MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
    }
    resetTimeFrameEvent(mFrame, mScratch);
  }

  TimeFrame& getTimeFrame() { return mFrame; }
  const TimeFrame& getTimeFrame() const { return mFrame; }
  ScratchN& getScratch() { return mScratch; }
  const ScratchN& getScratch() const { return mScratch; }
  const MFTPublicationCompatibility* getMFTPublicationCompatibility() const noexcept
  {
    if constexpr (DetId == o2::detectors::DetID::MFT) {
      return &static_cast<const MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar;
    }
    return nullptr;
  }
  const ITSSharedClusterCompatibility* getITSSharedClusterCompatibility() const noexcept
  {
    if constexpr (DetId == o2::detectors::DetID::ITS) {
      return &static_cast<const ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar;
    }
    return nullptr;
  }
  const std::vector<o2::itsmft::TrackingParameters>& getTrackingParameters() const { return mTrackParams; }
  std::optional<CommonTrackPublicationExport> getCommonTrackPublicationExport() const
  {
    if (!mPublicationClock || !mPlan || !isActive()) {
      return std::nullopt;
    }
    return CommonTrackPublicationExport{DetId, ClusterSourceId{0}, *mPublicationClock,
                                        gsl::span<const SurfaceId>{mPlan->getConfigurationKey().orderedSurfaces}};
  }
  bool isActive() const { return !mTrackParams.empty(); }
  // Actual tbb::task_arena concurrency the tracker was constructed with (0
  // if initialiseTracker() has not run yet, e.g. mTrackParams is empty).
  // Exposed for testing initialiseTracker()'s DetId-dependent nThreads
  // source (ITSCommonCATrackerParam for ITS, TrackerParamConfig<MFT> for
  // MFT unchanged) -- see TrackingInterface.cxx.
  int getTrackerNThreads() const { return mTrackerTraits ? mTrackerTraits->getNThreads() : 0; }

 protected:
  virtual void onTimeFrameLoaded() {}
  virtual void onTrackingFinished(float elapsedMs) {}

 private:
  void resolveTrackingParameters();
  void initialiseMemoryPool();
  void initialiseTracker();
  void loadTimeFrame(gsl::span<const o2::itsmft::ROFRecord> rofs,
                     gsl::span<const o2::itsmft::CompClusterExt> clusters,
                     gsl::span<const unsigned char> patterns,
                     const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels,
                     gsl::span<const o2::dataformats::IRFrame> irFrames);
  float runTracking();
  // Returns the single source-level ROFTimingConfig derived from per-layer
  // DPLAlpideParam values; throws TimeFrameLoadException if those values are
  // not uniform across layers (see ROFTimingUniformity.h). Also configures
  // mTimeFrame's ROF overlap/vertex-lookup tables as a side effect, as before.
  ROFTimingConfig configureROFLookupTables();
  void configureBeamPosition();
  void configureTrackingTopology();
  void configureROFMask(gsl::span<const o2::itsmft::ROFRecord> rofs,
                        gsl::span<const o2::dataformats::IRFrame> irFrames);
  void validateROFInput(gsl::span<const o2::itsmft::ROFRecord> rofs) const;

  bool mUseMC = false;
  bool mOverrideBeamEstimation = false;
  o2::itsmft::TrackingMode::Type mTrackingMode = o2::itsmft::TrackingMode::Unset;
  std::vector<o2::itsmft::TrackingParameters> mTrackParams;
  std::shared_ptr<BoundedMemoryResourceN> mMemoryPool;
  std::unique_ptr<TrackerTraitsN> mTrackerTraits;
  std::unique_ptr<TrackerN> mTracker;
#ifndef GPUCA_GPUCODE
  std::unique_ptr<ClusterDecoder> mClusterDecoder;
  // This interface's one immutable plan, built once in initialiseTracker()
  // from the compile-time-selected static per-detector catalog
  // (StaticDetectorCatalogs.h). Never rebuilt, invalidated, or forwarded
  // through TimeFrame -- runTracking()/loadTimeFrame() pass it explicitly to
  // Tracker::adoptDetectorLayoutSet()/TimeFrame::loadNormalizedSource().
  std::optional<DetectorLayoutSet> mPlan;
  std::optional<ClockTimingPublicationView> mPublicationClock;
  // M6e1: bind-once, built in initialiseTracker() alongside mPlan, adopted
  // onto mTracker immediately after. Only ever populated for the MFT
  // instantiation (if constexpr (DetId == MFT), see TrackingInterface.cxx) --
  // stays null for ITS, exactly today's existing behavior of never adopting
  // a binding at all (TrackerTraits<NLayers,...>::adoptDetectorTraversalBinding()
  // is optional; null means the identity-mapping fallback its own doc
  // comment already documents, unaffected by this slice for ITS).
  std::unique_ptr<BindingT> mBinding;
#endif
  const o2::itsmft::TopologyDictionary* mDict = nullptr;
  const o2::dataformats::MeanVertexObject* mMeanVertex = nullptr;
  // Gate 4 B3.1 lifetime contract: mFrame declared before mScratch, so C++'s
  // reverse-declaration-order destruction tears the scratch down first --
  // the permanent, non-templated TimeFrame always outlives the temporary
  // per-detector legacy scratch. Neither owns or stores a reference to the
  // other (see LegacyTrackerScratch.h); this owner is what binds both
  // together, via adoptScratch()/adoptFrame() on mTracker/mTrackerTraits and
  // explicit parameters everywhere else.
  TimeFrame mFrame;
  ScratchN mScratch;
  int mMFTROFrameLengthInBC = 0;
  bool mMFTTriggered = false;
};

using ITSMFTTrackingInterfaceITS = ITSMFTTrackingInterface<ITSNLayers>;
// M6e1: the standalone MFT common-CA workflow's own interface (o2-mft-ca-
// tracker-workflow, CATrackerSpec.cxx) now owns SurfaceTrackingScratch/
// SurfacePlanBinding instead of LegacyTrackerScratch<10>/
// DetectorTraversalBinding, mirroring M6d's combined-participant migration.
using ITSMFTTrackingInterfaceMFT = ITSMFTTrackingInterface<o2::mft::constants::mft::LayersNumber, SurfaceTrackingScratch, SurfacePlanBinding>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKINGINTERFACE_H_ */
