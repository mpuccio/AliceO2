// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_TEST_COMBINEDTRACKINGTESTSUPPORT_H_
#define ALICEO2_ITSMFT_TRACKING_TEST_COMBINEDTRACKINGTESTSUPPORT_H_

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/detail/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/detail/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/detail/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/DetectorTrackingOperationAdapterSupport.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITStracking/ROFLookupTables.h"

namespace o2::itsmft::tracking::test
{

// Test-only application-plan fixture. It composes the production Tracker,
// TrackerTraits, refit functions, plans, and bindings directly.
inline SurfaceCatalogView combinedCatalogView()
{
  return {kITSMFTCombinedStaticSurfaceCatalog.data(), static_cast<uint32_t>(kITSMFTCombinedStaticSurfaceCatalog.size())};
}

inline std::vector<SurfaceId> orderedSurfaceRange(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

inline SurfaceMask surfaceRangeMask(uint16_t first, uint16_t count)
{
  SurfaceMask result;
  for (uint16_t i = 0; i < count; ++i) {
    result.set(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

inline TrackerInitialization makeCombinedConfiguration(const TrackingParameters& itsParams,
                                                       const TrackingParameters& mftParams)
{
  TrackerInitialization configuration;
  configuration.catalog = combinedCatalogView();
  configuration.memoryPool = std::make_shared<BoundedMemoryResource>();
  TrackerIterationConfiguration iteration;
  const auto itsSurfaces = orderedSurfaceRange(0, ITSNLayers);
  const auto mftSurfaces = orderedSurfaceRange(ITSNLayers, MFTNLayers);
  SurfaceGraphSubgraph itsSubgraph;
  itsSubgraph.orderedSurfaces.assign(itsSurfaces.begin(), itsSurfaces.end());
  itsSubgraph.maxHoles = itsParams.MaxHoles;
  itsSubgraph.holeSurfaces = positionalSurfaceMask(itsParams.HoleLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));
  itsSubgraph.seedingSurfaces = positionalSurfaceMask(itsParams.StartLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));

  SurfaceGraphSubgraph mftSubgraph;
  mftSubgraph.orderedSurfaces.assign(mftSurfaces.begin(), mftSurfaces.end());
  mftSubgraph.maxHoles = mftParams.MaxHoles;
  mftSubgraph.holeSurfaces = positionalSurfaceMask(mftParams.HoleLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));
  mftSubgraph.seedingSurfaces = positionalSurfaceMask(mftParams.StartLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));

  iteration.graphSubgraphs = {std::move(itsSubgraph), std::move(mftSubgraph)};
  iteration.parameters = {itsParams, mftParams};
  iteration.bindings.push_back(SurfacePlanBinding::Declaration{ClusterSourceId{0},
                                                               surfaceRangeMask(0, ITSNLayers), itsSurfaces,
                                                               SurfaceKind::Cylinder});
  iteration.bindings.push_back(SurfacePlanBinding::Declaration{ClusterSourceId{1},
                                                               surfaceRangeMask(ITSNLayers, MFTNLayers), mftSurfaces,
                                                               SurfaceKind::Disk});
  configuration.iterations.push_back(std::move(iteration));
  return configuration;
}

class CombinedTrackingPlan
{
 public:
  CombinedTrackingPlan(std::vector<TrackingParameters> itsParams, std::vector<TrackingParameters> mftParams)
  {
    if (itsParams.size() != 1 || mftParams.size() != 1) {
      throw std::invalid_argument{"combined test application plan requires one iteration per detector"};
    }

    mConfiguration = makeCombinedConfiguration(itsParams[0], mftParams[0]);
    mITSPublicationAdapter.adoptITSSharedClusterCompatibility(&mITSCompatibility);
    mMFTPublicationAdapter.adoptMFTPublicationCompatibility(&mMFTCompatibility);
    mITSTracker = std::make_unique<Tracker>(&detail::refitDetectorSeed<o2::detectors::DetID::ITS>, ClusterSourceId{0});
    mMFTTracker = std::make_unique<Tracker>(&detail::refitDetectorSeed<o2::detectors::DetID::MFT>, ClusterSourceId{1});
    mITSTraits = std::make_unique<TrackerTraits>();
    mMFTTraits = std::make_unique<TrackerTraits>();
  }

  CombinedTrackingPlan(const CombinedTrackingPlan&) = delete;
  CombinedTrackingPlan& operator=(const CombinedTrackingPlan&) = delete;

  void adoptFrame(TimeFrame& frame)
  {
    mFrame = &frame;
    const auto result = mITSTracker->initialize(frame, mConfiguration);
    if (!result.ok()) {
      throw std::runtime_error{"combined test application plan failed to configure the TimeFrame"};
    }
    mITSTraits->setMemoryPool(frame.getMemoryPool());
    mMFTTraits->setMemoryPool(frame.getMemoryPool());
  }
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mITSTraits->setMemoryPool(pool);
    mMFTTraits->setMemoryPool(pool);
  }
  void setBz(float bz)
  {
    mFrame->setBz(bz);
  }
  void setNThreads(int n)
  {
    mITSTraits->setNThreads(n, mITSArena);
    mMFTTraits->setNThreads(n, mMFTArena);
  }

  Tracker& itsTracker() noexcept { return *mITSTracker; }
  Tracker& mftTracker() noexcept { return *mMFTTracker; }
  TrackingResult runITS()
  {
    auto result = mITSTracker->run(*mFrame, *mITSTraits);
    if (result.outcome == TrackingOutcome::Success) {
      const auto& params = mFrame->getTrackingParameters(ClusterSourceId{0});
      const auto& scratch = mFrame->getWorkspace(ClusterSourceId{0});
      for (std::size_t i = 0; i < params.size(); ++i) {
        const auto& candidates = mITSTraits->acceptedTracksForSharedStatus();
        if (i >= result.acceptedTrackCounts.size() || result.acceptedTrackCounts[i] > candidates.size() ||
            !mITSPublicationAdapter.completeAccepted(
              gsl::span<const TrackingCandidate>{candidates.data(), result.acceptedTrackCounts[i]}, params[i], scratch, i + 1 == params.size())) {
          throw std::runtime_error{"failed to seal ITS tracking compatibility"};
        }
      }
    } else {
      mITSPublicationAdapter.reset();
    }
    return result;
  }
  TrackingResult runMFT()
  {
    auto result = mMFTTracker->run(*mFrame, *mMFTTraits);
    if (result.outcome == TrackingOutcome::Success) {
      const auto& params = mFrame->getTrackingParameters(ClusterSourceId{1});
      const auto& scratch = mFrame->getWorkspace(ClusterSourceId{1});
      for (std::size_t i = 0; i < params.size(); ++i) {
        const auto& candidates = mMFTTraits->acceptedTracksForSharedStatus();
        if (i >= result.acceptedTrackCounts.size() || result.acceptedTrackCounts[i] > candidates.size() ||
            !mMFTPublicationAdapter.completeAccepted(
              gsl::span<const TrackingCandidate>{candidates.data(), result.acceptedTrackCounts[i]}, params[i], scratch, i + 1 == params.size())) {
          throw std::runtime_error{"failed to seal MFT tracking compatibility"};
        }
      }
    } else {
      mMFTPublicationAdapter.reset();
    }
    return result;
  }
  RuntimeROFViews getITSROFViews() const noexcept { return getITSScratch().getROFViews(); }
  RuntimeROFViews getMFTROFViews() const noexcept { return getMFTScratch().getROFViews(); }
  void clearPublicationSidecars() noexcept
  {
    mITSPublicationAdapter.reset();
    mMFTPublicationAdapter.reset();
  }

  std::optional<LoadSourcesResult> validateSources(const ClusterSourceInput& itsSource,
                                                   const ClusterSourceInput& mftSource) const noexcept
  {
    if (itsSource.id != ClusterSourceId{0} || itsSource.detector != o2::detectors::DetID::ITS) {
      return LoadSourcesResult{MultiSourceLoadError::UnsupportedDetector, itsSource.id};
    }
    if (mftSource.id != ClusterSourceId{1} || mftSource.detector != o2::detectors::DetID::MFT) {
      return LoadSourcesResult{MultiSourceLoadError::UnsupportedDetector, mftSource.id};
    }
    return std::nullopt;
  }

  SurfaceCatalogView catalogView() const noexcept { return combinedCatalogView(); }
  std::optional<bool> dropTFUponFailureFor(ClusterSourceId source) const noexcept
  {
    if (source == ClusterSourceId{0}) {
      return mFrame != nullptr && !mFrame->getTrackingParameters(ClusterSourceId{0}).empty()
               ? std::optional<bool>{mFrame->getTrackingParameters(ClusterSourceId{0})[0].DropTFUponFailure}
               : std::nullopt;
    }
    if (source == ClusterSourceId{1}) {
      return mFrame != nullptr && !mFrame->getTrackingParameters(ClusterSourceId{1}).empty()
               ? std::optional<bool>{mFrame->getTrackingParameters(ClusterSourceId{1})[0].DropTFUponFailure}
               : std::nullopt;
    }
    return std::nullopt;
  }
  void configureRofTables(const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource)
  {
    auto configure = [](auto& overlap, auto& vertex, auto& mask, const auto& timing, uint32_t nROFs, int layers) {
      o2::its::LayerTiming layerTiming{};
      layerTiming.mNROFsTF = nROFs;
      layerTiming.mROFLength = timing.rofLength;
      layerTiming.mROFDelay = timing.rofDelay;
      layerTiming.mROFBias = timing.rofBias;
      layerTiming.mROFAddTimeErr = timing.rofAddTimeErr;
      for (int layer = 0; layer < layers; ++layer) {
        overlap.defineLayer(layer, layerTiming);
        vertex.defineLayer(layer, layerTiming);
      }
      overlap.init();
      vertex.init();
      mask = std::remove_cvref_t<decltype(mask)>{overlap};
      mask.resetMask();
      for (int layer = 0; layer < layers; ++layer) {
        mask.setROFsEnabled(layer, 0, static_cast<int>(nROFs), 1);
      }
    };
    configure(mITSROFOverlapTable, mITSROFVertexLookupTable, mITSMultiplicityMask, itsSource.timing, static_cast<uint32_t>(itsSource.rofs.size()), ITSNLayers);
    configure(mMFTROFOverlapTable, mMFTROFVertexLookupTable, mMFTMultiplicityMask, mftSource.timing, static_cast<uint32_t>(mftSource.rofs.size()), MFTNLayers);
    mFrame->getWorkspace(ClusterSourceId{0}).setROFViews({mITSROFOverlapTable.getView(), mITSROFVertexLookupTable.getView(), mITSMultiplicityMask.getView(), mITSUPCMask.getView()});
    mFrame->getWorkspace(ClusterSourceId{1}).setROFViews({mMFTROFOverlapTable.getView(), mMFTROFVertexLookupTable.getView(), mMFTMultiplicityMask.getView(), mMFTUPCMask.getView()});
  }

  const SurfaceTrackingScratch& getITSScratch() const noexcept { return mFrame->getWorkspace(ClusterSourceId{0}); }
  const SurfaceTrackingScratch& getMFTScratch() const noexcept { return mFrame->getWorkspace(ClusterSourceId{1}); }
  gsl::span<const SurfaceId> getITSOrderedSurfaces() const noexcept { return mFrame->getBinding(0, ClusterSourceId{0})->getOrderedSurfaces(); }
  gsl::span<const SurfaceId> getMFTOrderedSurfaces() const noexcept { return mFrame->getBinding(0, ClusterSourceId{1})->getOrderedSurfaces(); }
  const ITSSharedClusterCompatibility& getITSSharedClusterCompatibility() const noexcept
  {
    return mITSCompatibility;
  }
  const MFTPublicationCompatibility& getMFTPublicationCompatibility() const noexcept
  {
    return mMFTCompatibility;
  }
  SurfaceGraphView getITSLayoutView() const noexcept { return mFrame != nullptr && mFrame->isConfigured() ? mFrame->getGraph(0).getView() : SurfaceGraphView{}; }
  SurfaceGraphView getMFTLayoutView() const noexcept { return getITSLayoutView(); }

 private:
  TrackerInitialization mConfiguration;
  TimeFrame* mFrame = nullptr;
  std::unique_ptr<Tracker> mITSTracker;
  std::unique_ptr<Tracker> mMFTTracker;
  std::unique_ptr<TrackerTraits> mITSTraits;
  std::unique_ptr<TrackerTraits> mMFTTraits;
  DetectorPublicationAdapter<ITSNLayers> mITSPublicationAdapter;
  DetectorPublicationAdapter<MFTNLayers> mMFTPublicationAdapter;
  ITSSharedClusterCompatibility mITSCompatibility;
  MFTPublicationCompatibility mMFTCompatibility;
  o2::its::ROFOverlapTable<ITSNLayers> mITSROFOverlapTable;
  o2::its::ROFVertexLookupTable<ITSNLayers> mITSROFVertexLookupTable;
  o2::its::ROFMaskTable<ITSNLayers> mITSMultiplicityMask;
  o2::its::ROFMaskTable<ITSNLayers> mITSUPCMask;
  o2::its::ROFOverlapTable<MFTNLayers> mMFTROFOverlapTable;
  o2::its::ROFVertexLookupTable<MFTNLayers> mMFTROFVertexLookupTable;
  o2::its::ROFMaskTable<MFTNLayers> mMFTMultiplicityMask;
  o2::its::ROFMaskTable<MFTNLayers> mMFTUPCMask;
  std::shared_ptr<tbb::task_arena> mITSArena;
  std::shared_ptr<tbb::task_arena> mMFTArena;
};

} // namespace o2::itsmft::tracking::test

#endif // ALICEO2_ITSMFT_TRACKING_TEST_COMBINEDTRACKINGTESTSUPPORT_H_
