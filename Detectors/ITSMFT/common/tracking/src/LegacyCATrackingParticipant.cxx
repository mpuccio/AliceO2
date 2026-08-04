// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/LegacyCATrackingParticipant.h"

#include <type_traits>
#include <utility>

namespace o2::itsmft::tracking
{

template <int NLayers, typename ScratchT, typename BindingT>
LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::LegacyCATrackingParticipant(ParticipantId id, std::vector<TrackingParameters> params)
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

template <int NLayers, typename ScratchT, typename BindingT>
void LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::adoptDetectorTraversalBinding(std::unique_ptr<BindingT> binding)
{
  mBinding = std::move(binding);
  mTraits.adoptDetectorTraversalBinding(mBinding.get());
}

template <int NLayers, typename ScratchT, typename BindingT>
void LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::adoptDetectorLayoutSet(const DetectorLayoutSet& plan)
{
  mPlan = &plan;
  // M6d: SurfaceTrackingScratch (unlike LegacyTrackerScratch<NLayers>, whose
  // Group A member is a fixed-size std::array<T, NLayers> valid at every
  // index from default construction) needs an explicit plan-adoption step
  // before its own runtime-sized containers are indexable. Both the plan
  // (just adopted above) and the binding (already adopted --
  // ITSMFTLegacyParticipantSet's constructor calls
  // adoptDetectorTraversalBinding() before adoptDetectorLayoutSet(), a
  // precondition for this call, not re-derived here) supply adoptPlan()'s
  // three runtime counts.
  if constexpr (std::is_same_v<ScratchT, SurfaceTrackingScratch>) {
    mScratch.adoptPlan(static_cast<std::size_t>(ownedSurfaces().size()),
                       mBinding->getGlobalTransitions().size(),
                       mBinding->getGlobalCells().size());
  }
  mTracker.adoptDetectorLayoutSet(plan);
}

template <int NLayers, typename ScratchT, typename BindingT>
void LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::adoptFrame(TimeFrame& frame)
{
  mTracker.adoptFrame(frame);
}

template <int NLayers, typename ScratchT, typename BindingT>
void LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  mScratch.setMemoryPool(pool);
  mTraits.setMemoryPool(pool);
  mTracker.setMemoryPool(pool);
}

template <int NLayers, typename ScratchT, typename BindingT>
void LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::setBz(float bz)
{
  mTracker.setBz(bz);
}

template <int NLayers, typename ScratchT, typename BindingT>
void LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::setNThreads(int n)
{
  mTraits.setNThreads(n, mArena);
}

template <int NLayers, typename ScratchT, typename BindingT>
void LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::configureRofTables(const ROFTimingConfig& timing, uint32_t nROFsTF)
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

template <int NLayers, typename ScratchT, typename BindingT>
void LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::clearCompatibility() noexcept
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    static_cast<ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    static_cast<MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }
}

template <int NLayers, typename ScratchT, typename BindingT>
void LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::clearPublicationSidecar() noexcept
{
  clearCompatibility();
}

template <int NLayers, typename ScratchT, typename BindingT>
gsl::span<const SurfaceId> LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::ownedSurfaces() const noexcept
{
  if (mPlan == nullptr) {
    return {};
  }
  return gsl::span<const SurfaceId>{mPlan->getConfigurationKey().orderedSurfaces};
}

template <int NLayers, typename ScratchT, typename BindingT>
ParticipantTrackingResult LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::track(TimeFrame&)
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

template <int NLayers, typename ScratchT, typename BindingT>
void LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::eventReset(TimeFrame&) noexcept
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

template <int NLayers, typename ScratchT, typename BindingT>
std::optional<ParticipantPublicationExport> LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::publicationExport() const
{
  if (!mTracked || mPlan == nullptr) {
    return std::nullopt;
  }
  return ParticipantPublicationExport{mId, gsl::span<const SurfaceId>{mPlan->getConfigurationKey().orderedSurfaces}};
}

template <int NLayers, typename ScratchT, typename BindingT>
const ITSSharedClusterCompatibility* LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::getITSSharedClusterCompatibility() const noexcept
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    return &static_cast<const ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar;
  }
  return nullptr;
}

template <int NLayers, typename ScratchT, typename BindingT>
const MFTPublicationCompatibility* LegacyCATrackingParticipant<NLayers, ScratchT, BindingT>::getMFTPublicationCompatibility() const noexcept
{
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    return &static_cast<const MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar;
  }
  return nullptr;
}

// M6d: ITS unchanged; MFT's sole production instantiation now owns
// SurfaceTrackingScratch/SurfacePlanBinding (LegacyCATrackingParticipant.h's
// own LegacyCATrackingParticipantMFT alias). No default-arg
// LegacyCATrackingParticipant<MFTNLayers> instantiation is needed: nothing
// in production uses it (ITSMFTTrackingInterface<MFTNLayers>, the
// standalone-MFT-workflow's own path, never used this class at all).
template class LegacyCATrackingParticipant<ITSNLayers>;
template class LegacyCATrackingParticipant<MFTNLayers, SurfaceTrackingScratch, SurfacePlanBinding>;

} // namespace o2::itsmft::tracking
