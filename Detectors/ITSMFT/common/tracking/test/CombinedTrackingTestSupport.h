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

#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/SurfaceGraphBuilder.h"
#include "ITSMFTTracking/ITSSharedClusterCompatibility.h"
#include "ITSMFTTracking/MFTPublicationCompatibility.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/SurfacePlanTrackingParticipant.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/TrackingParticipant.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"

namespace o2::itsmft::tracking::test
{

// Test-only application-plan fixture. It composes the production concrete
// participants, plans, bindings, and schedule directly. It deliberately owns
// no TimeFrame, publication clock, validity flag, or event-loop decision;
// tests that need those responsibilities keep them in their own workflow
// driver, mirroring CombinedCATrackerDPL.
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
                                                               SurfaceKind::Cylinder, TransitionPolicyTag::CylinderCylinder});
  iteration.bindings.push_back(SurfacePlanBinding::Declaration{ClusterSourceId{1},
                                                               surfaceRangeMask(ITSNLayers, MFTNLayers), mftSurfaces,
                                                               SurfaceKind::Disk, TransitionPolicyTag::DiskDisk});
  configuration.iterations.push_back(std::move(iteration));
  return configuration;
}

class CombinedTrackingParticipantPlan
{
 public:
  CombinedTrackingParticipantPlan(std::vector<TrackingParameters> itsParams, std::vector<TrackingParameters> mftParams)
  {
    if (itsParams.size() != 1 || mftParams.size() != 1) {
      throw std::invalid_argument{"combined test application plan requires one iteration per detector"};
    }

    mConfiguration = makeCombinedConfiguration(itsParams[0], mftParams[0]);
    mITSParticipant = std::make_unique<SurfacePlanTrackingParticipantITS>(ParticipantId{0}, ClusterSourceId{0});
    mMFTParticipant = std::make_unique<SurfacePlanTrackingParticipantMFT>(ParticipantId{1}, ClusterSourceId{1});
    mSchedule = {mITSParticipant.get(), mMFTParticipant.get()};
  }

  CombinedTrackingParticipantPlan(const CombinedTrackingParticipantPlan&) = delete;
  CombinedTrackingParticipantPlan& operator=(const CombinedTrackingParticipantPlan&) = delete;

  void adoptFrame(TimeFrame& frame)
  {
    mFrame = &frame;
    const auto result = mITSParticipant->initialize(frame, mConfiguration);
    if (!result.ok()) {
      throw std::runtime_error{"combined test application plan failed to configure the TimeFrame"};
    }
    mMFTParticipant->adoptConfiguredFrame(frame);
  }
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource>) {}
  void setBz(float bz)
  {
    mITSParticipant->setBz(bz);
    mMFTParticipant->setBz(bz);
  }
  void setNThreads(int n)
  {
    mITSParticipant->setNThreads(n);
    mMFTParticipant->setNThreads(n);
  }

  gsl::span<TrackingParticipant* const> schedule() noexcept { return mSchedule; }
  SurfacePlanTrackingParticipantITS& itsParticipant() noexcept { return *mITSParticipant; }
  SurfacePlanTrackingParticipantMFT& mftParticipant() noexcept { return *mMFTParticipant; }
  const SurfacePlanTrackingParticipantITS& itsParticipant() const noexcept { return *mITSParticipant; }
  const SurfacePlanTrackingParticipantMFT& mftParticipant() const noexcept { return *mMFTParticipant; }

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

  std::array<MultiSourceTimeFrameLoader::AtomicLoadBinding, 2> loadBindings(const ClusterSourceInput& itsSource,
                                                                            const ClusterSourceInput& mftSource) noexcept
  {
    return {MultiSourceTimeFrameLoader::AtomicLoadBinding{itsSource, mITSParticipant->loadTarget()},
            MultiSourceTimeFrameLoader::AtomicLoadBinding{mftSource, mMFTParticipant->loadTarget()}};
  }

  SurfaceCatalogView catalogView() const noexcept { return combinedCatalogView(); }
  std::optional<bool> dropTFUponFailureFor(ClusterSourceId source) const noexcept
  {
    if (source == ClusterSourceId{0}) {
      return mITSParticipant->getDropTFUponFailure();
    }
    if (source == ClusterSourceId{1}) {
      return mMFTParticipant->getDropTFUponFailure();
    }
    return std::nullopt;
  }
  void configureRofTables(const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource)
  {
    mITSParticipant->configureRofTables(itsSource.timing, static_cast<uint32_t>(itsSource.rofs.size()));
    mMFTParticipant->configureRofTables(mftSource.timing, static_cast<uint32_t>(mftSource.rofs.size()));
  }

  const SurfaceTrackingScratch& getITSScratch() const noexcept { return mITSParticipant->getScratch(); }
  const SurfaceTrackingScratch& getMFTScratch() const noexcept { return mMFTParticipant->getScratch(); }
  gsl::span<const SurfaceId> getITSOrderedSurfaces() const noexcept { return mITSParticipant->ownedSurfaces(); }
  gsl::span<const SurfaceId> getMFTOrderedSurfaces() const noexcept { return mMFTParticipant->ownedSurfaces(); }
  const ITSSharedClusterCompatibility& getITSSharedClusterCompatibility() const noexcept
  {
    return *mITSParticipant->getITSSharedClusterCompatibility();
  }
  const MFTPublicationCompatibility& getMFTPublicationCompatibility() const noexcept
  {
    return *mMFTParticipant->getMFTPublicationCompatibility();
  }
  SurfaceGraphView getITSLayoutView() const noexcept { return mFrame != nullptr && mFrame->isConfigured() ? mFrame->getGraph(0).getView() : SurfaceGraphView{}; }
  SurfaceGraphView getMFTLayoutView() const noexcept { return getITSLayoutView(); }

 private:
  TrackerInitialization mConfiguration;
  TimeFrame* mFrame = nullptr;
  std::unique_ptr<SurfacePlanTrackingParticipantITS> mITSParticipant;
  std::unique_ptr<SurfacePlanTrackingParticipantMFT> mMFTParticipant;
  std::array<TrackingParticipant*, 2> mSchedule{};
};

} // namespace o2::itsmft::tracking::test

#endif // ALICEO2_ITSMFT_TRACKING_TEST_COMBINEDTRACKINGTESTSUPPORT_H_
