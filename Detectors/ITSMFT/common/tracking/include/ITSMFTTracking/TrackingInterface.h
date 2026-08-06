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
#ifndef GPUCA_GPUCODE
#include "ITSMFTTracking/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/ClockTimingPublicationView.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/TimeFrameLoadFailure.h"
#include "ITSMFTTracking/TrackingOperationAdapter.h"
#endif
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITStracking/BoundedAllocator.h"
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

// Host-facing standalone workflow interface over the same
// SurfaceTrackingScratch/SurfacePlanBinding model used by the combined
// participant. NLayers remains an algorithm/storage parameter until M7d.
template <int NLayers>
class ITSMFTTrackingInterface
#ifndef GPUCA_GPUCODE
  : private TrackingOperationAdapter
#endif
{
 public:
  static_assert(NLayers == ITSNLayers || NLayers == o2::mft::constants::mft::LayersNumber,
                "ITSMFTTrackingInterface supports ITS (7) and MFT (10) layer counts only");
  static constexpr o2::detectors::DetID::ID DetId = detIdFromNLayers<NLayers>();

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

#ifndef GPUCA_GPUCODE
  // The workflow owns the detector timing tables and compatibility sidecar.
  // The interface borrows this adapter for the lifetime of the workflow task.
  void bindPublicationAdapter(DetectorPublicationAdapter<NLayers>& adapter) noexcept { mDetectorPublicationAdapter = &adapter; }
  // Optional borrowed view source for callers that reuse one workflow context
  // across event resets. No table or timing storage is retained here.
  void bindROFViews(const RuntimeROFViews& views) noexcept { mBoundROFViews = &views; }
#endif

  void initialise();

  /// Phase 2 (TimeFrame load) + Phase 3 (CA tracking).
  ///
  /// - Success (load and tracking both completed): returns elapsed tracking
  ///   ms (>= 0). The loaded/tracked event is retained in the TimeFrame;
  ///   call clearTimeFrame() once outputs have been extracted.
  /// - Recoverable failure (malformed per-TF loading input, or a per-TF
  ///   resource-exhaustion exception) with DropTFUponFailure=true: the
  ///   TimeFrame is fully wiped first, then this returns exactly
  ///   kDroppedTimeFrameResult (see Tracker.h); nothing was published for
  ///   this TF. onTrackingFinished()/onTimeFrameLoaded() are not both
  ///   called in this case -- see loadTimeFrame()/runTracking() for which
  ///   ran.
  /// - Every other failure -- structural/configuration (bad or stale
  ///   catalog, non-uniform per-layer ROF timing, a missing dictionary, a
  ///   structural LoadSourcesResult), unclassified, or a recoverable
  ///   failure with DropTFUponFailure=false: the TimeFrame is fully wiped,
  ///   then the failure is rethrown (by its original type for an
  ///   unclassified exception; see TimeFrameLoadFailure.h for the loading
  ///   boundary's typed exceptions and Tracker.h for the CA-tracking
  ///   ones). No output is created on this path.
  float processTimeFrame(gsl::span<const o2::itsmft::ROFRecord> rofs,
                         gsl::span<const o2::itsmft::CompClusterExt> clusters,
                         gsl::span<const unsigned char> patterns,
                         const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels,
                         gsl::span<const o2::dataformats::IRFrame> irFrames = {},
                         RuntimeROFViews rofViews = {});

  // Owner-level reset clears adapter state and invokes the frame-owned generic
  // event reset exactly once.
  void resetEvent()
  {
    if (mDetectorPublicationAdapter != nullptr) {
      mDetectorPublicationAdapter->reset();
    }
    mFrame.resetEvent();
  }

  TimeFrame& getTimeFrame() { return mFrame; }
  const TimeFrame& getTimeFrame() const { return mFrame; }
  SurfaceTrackingScratch& getScratch() { return mFrame.getWorkspace(ClusterSourceId{0}); }
  const SurfaceTrackingScratch& getScratch() const { return mFrame.getWorkspace(ClusterSourceId{0}); }
  const MFTPublicationCompatibility* getMFTPublicationCompatibility() const noexcept
  {
    return mDetectorPublicationAdapter == nullptr ? nullptr : mDetectorPublicationAdapter->getMFTPublicationCompatibility();
  }
  const ITSSharedClusterCompatibility* getITSSharedClusterCompatibility() const noexcept
  {
    return mDetectorPublicationAdapter == nullptr ? nullptr : mDetectorPublicationAdapter->getITSSharedClusterCompatibility();
  }
  const std::vector<o2::itsmft::TrackingParameters>& getTrackingParameters() const { return mFrame.getTrackingParameters(); }
  std::optional<CommonTrackPublicationExport> getCommonTrackPublicationExport(const ClockTimingPublicationView& clock) const
  {
    if (!mFrame.isConfigured() || !isActive()) {
      return std::nullopt;
    }
    return CommonTrackPublicationExport{DetId, ClusterSourceId{0}, clock,
                                        gsl::span<const SurfaceId>{mFrame.getGraph(0).getOrderedSurfaces()}};
  }
  // Test and adapter convenience: derive a value export from the currently
  // borrowed runtime view. The interface still owns neither the tables nor
  // the clock state.
  std::optional<CommonTrackPublicationExport> getCommonTrackPublicationExport() const
  {
    if (getScratch().getTotalClusters() == 0) {
      return std::nullopt;
    }
    const auto views = getScratch().getROFViews();
    if (views.overlap.mLayerCount <= 0) {
      return std::nullopt;
    }
    return getCommonTrackPublicationExport(ClockTimingPublicationView{views.overlap.getClockLayer()});
  }
  bool isActive() const { return mFrame.isConfigured() && !mFrame.getTrackingParameters().empty(); }
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
  std::vector<o2::itsmft::TrackingParameters> resolveTrackingParameters();
  void initialiseTracker(const std::vector<o2::itsmft::TrackingParameters>& parameters);
  void loadTimeFrame(gsl::span<const o2::itsmft::ROFRecord> rofs,
                     gsl::span<const o2::itsmft::CompClusterExt> clusters,
                     gsl::span<const unsigned char> patterns,
                     const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels,
                     gsl::span<const o2::dataformats::IRFrame> irFrames,
                     RuntimeROFViews rofViews);
  float runTracking();
  void configureBeamPosition();
  void validateROFInput(gsl::span<const o2::itsmft::ROFRecord> rofs) const;

#ifndef GPUCA_GPUCODE
  bool refitSeed(const TrackSeed& seed,
                 const TrackingParameters& params,
                 float bz,
                 SurfaceTrackingScratch& scratch,
                 gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                 SurfaceCatalogView surfaceCatalog,
                 ClusterSourceId expectedSource,
                 TrackingCandidate& candidate) override;
  bool completeAccepted(gsl::span<const TrackingCandidate> candidates,
                        const TrackingParameters& params,
                        const SurfaceTrackingScratch& scratch,
                        bool final) override;
  void resetAdapterState() noexcept override;
#endif

  bool mUseMC = false;
  bool mOverrideBeamEstimation = false;
  o2::itsmft::TrackingMode::Type mTrackingMode = o2::itsmft::TrackingMode::Unset;
  std::unique_ptr<TrackerTraits> mTrackerTraits;
  std::unique_ptr<Tracker> mTracker;
#ifndef GPUCA_GPUCODE
  DetectorPublicationAdapter<NLayers>* mDetectorPublicationAdapter = nullptr;
  const RuntimeROFViews* mBoundROFViews = nullptr;
  std::unique_ptr<ClusterDecoder> mClusterDecoder;
#endif
  const o2::itsmft::TopologyDictionary* mDict = nullptr;
  const o2::dataformats::MeanVertexObject* mMeanVertex = nullptr;
  TimeFrame mFrame;
};

// M6e2: the standalone ITS common-CA workflow's own interface (o2-its-ca-
// tracker-workflow, CATrackerSpec.cxx) now owns SurfaceTrackingScratch/
// SurfacePlanBinding too, completing the same migration M6e1 already did for
// standalone MFT.
using ITSMFTTrackingInterfaceITS = ITSMFTTrackingInterface<ITSNLayers>;
// M6e1: the standalone MFT common-CA workflow's own interface (o2-mft-ca-
// tracker-workflow, CATrackerSpec.cxx) now owns SurfaceTrackingScratch/
// SurfacePlanBinding, mirroring the combined participant.
using ITSMFTTrackingInterfaceMFT = ITSMFTTrackingInterface<o2::mft::constants::mft::LayersNumber>;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TRACKINGINTERFACE_H_ */
