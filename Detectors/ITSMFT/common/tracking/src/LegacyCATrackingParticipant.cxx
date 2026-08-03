// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/LegacyCATrackingParticipant.h"

#include <utility>

namespace o2::itsmft::tracking
{

template <int NLayers>
LegacyCATrackingParticipant<NLayers>::LegacyCATrackingParticipant(ParticipantId id, std::vector<TrackingParameters> params)
  : mId{id}, mParams(std::move(params)), mTracker(&mTraits), mLoadTarget(mScratch)
{
  mTracker.adoptScratch(mScratch);
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    mTracker.adoptITSSharedClusterCompatibility(static_cast<ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar);
  }
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    mTracker.adoptMFTPublicationCompatibility(static_cast<MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar);
  }
  mTracker.setParameters(mParams);
}

template <int NLayers>
void LegacyCATrackingParticipant<NLayers>::adoptDetectorTraversalBinding(std::unique_ptr<DetectorTraversalBinding> binding)
{
  mBinding = std::move(binding);
  mTraits.adoptDetectorTraversalBinding(mBinding.get());
}

template <int NLayers>
void LegacyCATrackingParticipant<NLayers>::adoptDetectorLayoutSet(const DetectorLayoutSet& plan)
{
  mPlan = &plan;
  mTracker.adoptDetectorLayoutSet(plan);
}

template <int NLayers>
void LegacyCATrackingParticipant<NLayers>::adoptFrame(TimeFrame& frame)
{
  mTracker.adoptFrame(frame);
}

template <int NLayers>
void LegacyCATrackingParticipant<NLayers>::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  mScratch.setMemoryPool(pool);
  mTraits.setMemoryPool(pool);
  mTracker.setMemoryPool(pool);
}

template <int NLayers>
void LegacyCATrackingParticipant<NLayers>::setBz(float bz)
{
  mTracker.setBz(bz);
}

template <int NLayers>
void LegacyCATrackingParticipant<NLayers>::setNThreads(int n)
{
  mTraits.setNThreads(n, mArena);
}

template <int NLayers>
void LegacyCATrackingParticipant<NLayers>::configureRofTables(const ROFTimingConfig& timing, uint32_t nROFsTF)
{
  o2::its::LayerTiming layerTiming{};
  layerTiming.mNROFsTF = nROFsTF;
  layerTiming.mROFLength = timing.rofLength;
  layerTiming.mROFDelay = timing.rofDelay;
  layerTiming.mROFBias = timing.rofBias;
  layerTiming.mROFAddTimeErr = timing.rofAddTimeErr;

  typename ScratchN::ROFOverlapTableN rofTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    rofTable.defineLayer(layer, layerTiming);
  }
  rofTable.init();
  mScratch.setROFOverlapTable(rofTable);

  typename ScratchN::ROFVertexLookupTableN vtxTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    vtxTable.defineLayer(layer, layerTiming);
  }
  vtxTable.init();
  mScratch.setROFVertexLookupTable(vtxTable);

  typename ScratchN::ROFMaskTableN mask{rofTable};
  mask.resetMask();
  for (int layer = 0; layer < NLayers; ++layer) {
    mask.setROFsEnabled(layer, 0, static_cast<int>(nROFsTF), 1);
  }
  mScratch.setMultiplicityCutMask(std::move(mask));

  mScratch.initTrackerTopologies(mParams);
}

template <int NLayers>
void LegacyCATrackingParticipant<NLayers>::clearCompatibility() noexcept
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    static_cast<ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    static_cast<MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }
}

template <int NLayers>
void LegacyCATrackingParticipant<NLayers>::clearPublicationSidecar() noexcept
{
  clearCompatibility();
}

template <int NLayers>
gsl::span<const SurfaceId> LegacyCATrackingParticipant<NLayers>::ownedSurfaces() const noexcept
{
  if (mPlan == nullptr) {
    return {};
  }
  return gsl::span<const SurfaceId>{mPlan->getConfigurationKey().orderedSurfaces};
}

template <int NLayers>
ParticipantTrackingResult LegacyCATrackingParticipant<NLayers>::track(TimeFrame&)
{
  // clustersToTracks() itself only ever *returns* Success or
  // RecoverableDropped (CATracker.h); Structural always escapes as a thrown
  // exception instead, left uncaught here for TrackingEngine::executeEvent
  // ()'s own try/catch to classify -- same division of responsibility the
  // coordinator's process() already relied on before this slice.
  const auto result = mTracker.clustersToTracks();
  if (result.outcome != TrackingOutcome::Success) {
    mTracked = false;
    return {ParticipantOutcome::RecoverableDropped, 0};
  }
  mTracked = true;
  return {ParticipantOutcome::Success, mScratch.getNumberOfTracks()};
}

template <int NLayers>
void LegacyCATrackingParticipant<NLayers>::eventReset(TimeFrame&) noexcept
{
  // Scratch/sidecar only -- never `frame` itself, see TrackingParticipant.h.
  // resetScratch()+the engine's own single TimeFrame::wipe() together
  // reproduce MultiSourceTimeFrameLoader::resetITSAndMFTEvent()'s exact
  // sequencing (LegacyTrackerScratch.h's resetTimeFrameEvent() doc
  // anticipates exactly this: "a future owner with several participating
  // scratches ... must reset every participating scratch first and only
  // then wipe the one shared TimeFrame once").
  mTracked = false;
  mScratch.resetScratch();
  clearCompatibility();
}

template <int NLayers>
std::optional<ParticipantPublicationExport> LegacyCATrackingParticipant<NLayers>::publicationExport() const
{
  if (!mTracked || mPlan == nullptr) {
    return std::nullopt;
  }
  return ParticipantPublicationExport{mId, gsl::span<const SurfaceId>{mPlan->getConfigurationKey().orderedSurfaces}};
}

template <int NLayers>
const ITSSharedClusterCompatibility* LegacyCATrackingParticipant<NLayers>::getITSSharedClusterCompatibility() const noexcept
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    return &static_cast<const ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar;
  }
  return nullptr;
}

template <int NLayers>
const MFTPublicationCompatibility* LegacyCATrackingParticipant<NLayers>::getMFTPublicationCompatibility() const noexcept
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    return &static_cast<const MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar;
  }
  return nullptr;
}

template class LegacyCATrackingParticipant<ITSNLayers>;
template class LegacyCATrackingParticipant<MFTNLayers>;

} // namespace o2::itsmft::tracking
