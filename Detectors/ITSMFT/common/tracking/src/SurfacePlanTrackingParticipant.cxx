// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/SurfacePlanTrackingParticipant.h"

#include <utility>

namespace o2::itsmft::tracking
{

template <int NLayers>
SurfacePlanTrackingParticipant<NLayers>::SurfacePlanTrackingParticipant(ParticipantId id, std::vector<TrackingParameters> params)
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
void SurfacePlanTrackingParticipant<NLayers>::adoptSurfacePlanBinding(std::unique_ptr<SurfacePlanBinding> binding)
{
  mBinding = std::move(binding);
  mTraits.adoptSurfacePlanBinding(mBinding.get());
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::adoptDetectorLayoutSet(const DetectorLayoutSet& plan)
{
  mPlan = &plan;
  // The plan and binding are adopted before this call, so their compact
  // topology counts size the scratch before any loader or tracker access.
  mScratch.adoptPlan(static_cast<std::size_t>(ownedSurfaces().size()),
                     mBinding->getGlobalTransitions().size(),
                     mBinding->getGlobalCells().size());
  mTracker.adoptDetectorLayoutSet(plan);
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::adoptFrame(TimeFrame& frame)
{
  mTracker.adoptFrame(frame);
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  mScratch.setMemoryPool(pool);
  mTraits.setMemoryPool(pool);
  mTracker.setMemoryPool(pool);
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::setBz(float bz)
{
  mTracker.setBz(bz);
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::setNThreads(int n)
{
  mTraits.setNThreads(n, mArena);
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::configureRofTables(const ROFTimingConfig& timing, uint32_t nROFsTF)
{
  o2::its::LayerTiming layerTiming{};
  layerTiming.mNROFsTF = nROFsTF;
  layerTiming.mROFLength = timing.rofLength;
  layerTiming.mROFDelay = timing.rofDelay;
  layerTiming.mROFBias = timing.rofBias;
  layerTiming.mROFAddTimeErr = timing.rofAddTimeErr;

  // M6e2: bound directly on this class's own NLayers, not via a scratch alias:
  // ROFOverlapTableN/etc -- SurfaceTrackingScratch's own storage is now
  // dual-typed (ITS(7)+MFT(10)), so it no longer exposes a single nested
  // alias for these (see SurfaceTrackingScratch.h, TrackingInterface.h's
  // identical fix). setROFOverlapTable()/setROFVertexLookupTable()/
  // setMultiplicityCutMask() are template methods on SurfaceTrackingScratch
  // now, but deduce NLayers from this locally-typed argument with no
  // explicit `<>` at the call site, so no `.template` disambiguator or
  // dispatcher is needed here (unlike the argument-less getters elsewhere).
  o2::its::ROFOverlapTable<NLayers> rofTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    rofTable.defineLayer(layer, layerTiming);
  }
  rofTable.init();
  mScratch.setROFOverlapTable(rofTable);

  o2::its::ROFVertexLookupTable<NLayers> vtxTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    vtxTable.defineLayer(layer, layerTiming);
  }
  vtxTable.init();
  mScratch.setROFVertexLookupTable(vtxTable);

  o2::its::ROFMaskTable<NLayers> mask{rofTable};
  mask.resetMask();
  for (int layer = 0; layer < NLayers; ++layer) {
    mask.setROFsEnabled(layer, 0, static_cast<int>(nROFsTF), 1);
  }
  mScratch.setMultiplicityCutMask(std::move(mask));

  scratchInitTrackerTopologies<NLayers>(mScratch, mParams, static_cast<int>(mScratch.getNOwnedSurfaces()));
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::clearCompatibility() noexcept
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    static_cast<ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    static_cast<MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::clearPublicationSidecar() noexcept
{
  clearCompatibility();
}

template <int NLayers>
gsl::span<const SurfaceId> SurfacePlanTrackingParticipant<NLayers>::ownedSurfaces() const noexcept
{
  if (mBinding == nullptr) {
    return {};
  }
  return mBinding->getOrderedSurfaces();
}

template <int NLayers>
ParticipantTrackingResult SurfacePlanTrackingParticipant<NLayers>::track(TimeFrame& frame)
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
  return {ParticipantOutcome::Success, frame.getCommonTracks().size()};
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::eventReset(TimeFrame&) noexcept
{
  // Scratch/sidecar only -- never `frame` itself, see TrackingParticipant.h.
  // resetScratch()+the engine's own single TimeFrame::wipe() together
  // reproduce the loader's reset sequencing: every participating scratch
  // is cleared before the one shared TimeFrame is wiped.
  mTracked = false;
  mScratch.resetScratch();
  clearCompatibility();
}

template <int NLayers>
std::optional<ParticipantPublicationExport> SurfacePlanTrackingParticipant<NLayers>::publicationExport() const
{
  if (!mTracked || mPlan == nullptr) {
    return std::nullopt;
  }
  return ParticipantPublicationExport{mId, ownedSurfaces()};
}

template <int NLayers>
const ITSSharedClusterCompatibility* SurfacePlanTrackingParticipant<NLayers>::getITSSharedClusterCompatibility() const noexcept
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    return &static_cast<const ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar;
  }
  return nullptr;
}

template <int NLayers>
const MFTPublicationCompatibility* SurfacePlanTrackingParticipant<NLayers>::getMFTPublicationCompatibility() const noexcept
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    return &static_cast<const MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar;
  }
  return nullptr;
}

template class SurfacePlanTrackingParticipant<ITSNLayers>;
template class SurfacePlanTrackingParticipant<MFTNLayers>;

} // namespace o2::itsmft::tracking
