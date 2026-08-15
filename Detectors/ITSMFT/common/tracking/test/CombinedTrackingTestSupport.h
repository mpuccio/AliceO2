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

#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/detail/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/detail/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/detail/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/detail/DetectorRefitSupport.h"
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
  auto itsDefinition = makeSurfaceLayoutChain(
    itsSurfaces, itsParams.MaxHoles,
    positionalSurfaceMask(itsParams.HoleLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size())),
    positionalSurfaceMask(itsParams.StartLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size())));
  auto mftDefinition = makeSurfaceLayoutChain(
    mftSurfaces, mftParams.MaxHoles,
    positionalSurfaceMask(mftParams.HoleLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size())),
    positionalSurfaceMask(mftParams.StartLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size())));
  SurfaceLayoutDefinition definition;
  definition.orderedSurfaces = std::move(itsDefinition.orderedSurfaces);
  const auto offset = static_cast<uint16_t>(definition.orderedSurfaces.size());
  definition.orderedSurfaces.insert(definition.orderedSurfaces.end(), mftDefinition.orderedSurfaces.begin(), mftDefinition.orderedSurfaces.end());
  definition.componentOffsets = {0, offset};
  definition.maxHoles = itsDefinition.maxHoles;
  definition.holeSurfaces = itsDefinition.holeSurfaces | mftDefinition.holeSurfaces;
  definition.seedingSurfaces = surfaceRangeMask(0, ITSNLayers + MFTNLayers);
  iteration.layout = std::move(definition);
  const auto combine = [&] {
    auto parameters = itsParams;
    parameters.NLayers = ITSNLayers + MFTNLayers;
    const auto concatenate = [](auto& output, const auto& prefix, const auto& suffix) {
      output = prefix;
      output.insert(output.end(), suffix.begin(), suffix.end());
    };
    concatenate(parameters.AddTimeError, itsParams.AddTimeError, mftParams.AddTimeError);
    concatenate(parameters.LayerZ, itsParams.LayerZ, mftParams.LayerZ);
    concatenate(parameters.LayerRadii, itsParams.LayerRadii, mftParams.LayerRadii);
    concatenate(parameters.LayerxX0, itsParams.LayerxX0, mftParams.LayerxX0);
    concatenate(parameters.LayerResolution, itsParams.LayerResolution, mftParams.LayerResolution);
    concatenate(parameters.SystError2Row, itsParams.SystError2Row, mftParams.SystError2Row);
    concatenate(parameters.SystError2Col, itsParams.SystError2Col, mftParams.SystError2Col);
    parameters.LayerColHalfExtent = itsParams.LayerColHalfExtent.empty() ? itsParams.LayerZ : itsParams.LayerColHalfExtent;
    const auto& mftColExtent = mftParams.LayerColHalfExtent.empty() ? mftParams.LayerZ : mftParams.LayerColHalfExtent;
    parameters.LayerColHalfExtent.insert(parameters.LayerColHalfExtent.end(), mftColExtent.begin(), mftColExtent.end());
    parameters.StartLayerMask = LayerMask{(uint32_t{1} << (ITSNLayers + MFTNLayers)) - 1u};
    parameters.HoleLayerMask = LayerMask{itsParams.HoleLayerMask.value() |
                                         (mftParams.HoleLayerMask.value() << ITSNLayers)};
    return parameters;
  };
  iteration.parameters = combine();
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
    mTracker = std::make_unique<Tracker>(&detail::refitSurfaceSeed);
    mTraits = std::make_unique<TrackerTraits>();
  }

  CombinedTrackingPlan(const CombinedTrackingPlan&) = delete;
  CombinedTrackingPlan& operator=(const CombinedTrackingPlan&) = delete;

  void adoptFrame(TimeFrame& frame)
  {
    mFrame = &frame;
    const auto result = mTracker->initialize(frame, mConfiguration);
    if (!result.ok()) {
      throw std::runtime_error{"combined test application plan failed to configure the TimeFrame"};
    }
  }
  void setBz(float bz)
  {
    mFrame->setBz(bz);
  }
  void setNThreads(int n)
  {
    mTraits->setNThreads(n, mArena);
  }

  Tracker& itsTracker() noexcept { return *mTracker; }
  Tracker& mftTracker() noexcept { return *mTracker; }
  TrackingResult runITS()
  {
    auto result = mTracker->run(*mFrame, *mTraits);
    mLastResult = result;
    if (result.outcome == TrackingOutcome::Success) {
      const auto& params = mFrame->getTrackingParameters();
      const auto& scratch = mFrame->getWorkspace();
      for (std::size_t i = 0; i < params.size(); ++i) {
        const auto& candidates = scratch.getTraversalWorkspace(i).acceptedTracks;
        if (i >= result.acceptedTrackCounts.size() || result.acceptedTrackCounts[i] != candidates.size()) {
          throw std::runtime_error{"failed to seal ITS tracking compatibility"};
        }
        std::vector<TrackingCandidate> selected;
        for (std::size_t index = 0; index < result.acceptedTrackCounts[i]; ++index) {
          if (candidates[index].track.innerState.kind == SurfaceKind::Cylinder) {
            selected.push_back(candidates[index]);
          }
        }
        if (!mITSPublicationAdapter.completeAccepted(
              selected, params[i], scratch, i + 1 == params.size())) {
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
    if (!mLastResult) {
      runITS();
    }
    auto result = *mLastResult;
    if (result.outcome == TrackingOutcome::Success) {
      const auto& params = mFrame->getTrackingParameters();
      const auto& scratch = mFrame->getWorkspace();
      for (std::size_t i = 0; i < params.size(); ++i) {
        const auto& candidates = scratch.getTraversalWorkspace(i).acceptedTracks;
        if (i >= result.acceptedTrackCounts.size() || result.acceptedTrackCounts[i] != candidates.size()) {
          throw std::runtime_error{"failed to seal MFT tracking compatibility"};
        }
        std::vector<TrackingCandidate> selected;
        for (std::size_t index = 0; index < result.acceptedTrackCounts[i]; ++index) {
          if (candidates[index].track.innerState.kind == SurfaceKind::Disk) {
            selected.push_back(candidates[index]);
          }
        }
        if (!mMFTPublicationAdapter.completeAccepted(
              selected, params[i], scratch, i + 1 == params.size())) {
          throw std::runtime_error{"failed to seal MFT tracking compatibility"};
        }
      }
    } else {
      mMFTPublicationAdapter.reset();
    }
    return result;
  }
  RuntimeROFViews getITSROFViews() const noexcept { return {mITSROFOverlapTable.getView(), mITSROFVertexLookupTable.getView(), mITSMultiplicityMask.getView(), mITSUPCMask.getView()}; }
  RuntimeROFViews getMFTROFViews() const noexcept { return {mMFTROFOverlapTable.getView(), mMFTROFVertexLookupTable.getView(), mMFTMultiplicityMask.getView(), mMFTUPCMask.getView()}; }
  void clearPublicationSidecars() noexcept
  {
    mITSPublicationAdapter.reset();
    mMFTPublicationAdapter.reset();
    mLastResult.reset();
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
      return mFrame != nullptr && !mFrame->getTrackingParameters().empty()
               ? std::optional<bool>{mFrame->getTrackingParameters()[0].DropTFUponFailure}
               : std::nullopt;
    }
    if (source == ClusterSourceId{1}) {
      return mFrame != nullptr && !mFrame->getTrackingParameters().empty()
               ? std::optional<bool>{mFrame->getTrackingParameters()[0].DropTFUponFailure}
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
  }

  const SurfaceTrackingScratch& getITSScratch() const noexcept { return mFrame->getWorkspace(); }
  const SurfaceTrackingScratch& getMFTScratch() const noexcept { return mFrame->getWorkspace(); }
  gsl::span<const SurfaceId> getITSOrderedSurfaces() const noexcept { return mFrame->getLayout(0).getOrderedSurfaces().first(ITSNLayers); }
  gsl::span<const SurfaceId> getMFTOrderedSurfaces() const noexcept { return mFrame->getLayout(0).getOrderedSurfaces().subspan(ITSNLayers, MFTNLayers); }
  const ITSSharedClusterCompatibility& getITSSharedClusterCompatibility() const noexcept
  {
    return mITSCompatibility;
  }
  const MFTPublicationCompatibility& getMFTPublicationCompatibility() const noexcept
  {
    return mMFTCompatibility;
  }
  TraversalTopologyView getITSLayoutView() const noexcept
  {
    return mFrame != nullptr && mFrame->isConfigured() ? mFrame->getWorkspace().getTraversalWorkspace(0).getTopologyView() : TraversalTopologyView{};
  }
  TraversalTopologyView getMFTLayoutView() const noexcept { return getITSLayoutView(); }

 private:
  TrackerInitialization mConfiguration;
  TimeFrame* mFrame = nullptr;
  std::unique_ptr<Tracker> mTracker;
  std::unique_ptr<TrackerTraits> mTraits;
  std::optional<TrackingResult> mLastResult;
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
  std::shared_ptr<tbb::task_arena> mArena;
};

} // namespace o2::itsmft::tracking::test

#endif // ALICEO2_ITSMFT_TRACKING_TEST_COMBINEDTRACKINGTESTSUPPORT_H_
