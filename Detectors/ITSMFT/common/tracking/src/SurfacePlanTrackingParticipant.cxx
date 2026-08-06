// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/SurfacePlanTrackingParticipant.h"

#include "ITSMFTTracking/CommonTrackShadow.h"
#include "ITSMFTTracking/DetectorTrackingOperationAdapterSupport.h"
#include "ITStracking/Constants.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace o2::itsmft::tracking
{

template <int NLayers>
SurfacePlanTrackingParticipant<NLayers>::SurfacePlanTrackingParticipant(ParticipantId id, ClusterSourceId source)
  : mId{id}, mSource{source}, mTracker(&mTraits)
{
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    mDetectorPublicationAdapter.adoptITSSharedClusterCompatibility(&static_cast<ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar);
  }
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    mDetectorPublicationAdapter.adoptMFTPublicationCompatibility(&static_cast<MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar);
  }
}

template <int NLayers>
TrackerInitializationResult SurfacePlanTrackingParticipant<NLayers>::initialize(TimeFrame& frame, const TrackerInitialization& configuration)
{
  mTracker.setSource(mSource);
  const auto result = mTracker.initialize(frame, configuration);
  if (!result.ok()) {
    return result;
  }
  adoptConfiguredFrame(frame);
  return result;
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::adoptConfiguredFrame(TimeFrame& frame)
{
  mFrame = &frame;
  mTracker.setSource(mSource);
  mTracker.adoptFrame(frame);
  mTraits.setMemoryPool(frame.getMemoryPool());
  if (!frame.isConfigured()) {
    throw std::runtime_error{"SurfacePlanTrackingParticipant: frame is not configured"};
  }
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

  // Build the frozen fixed-capacity tables at this application-adapter edge.
  // Only their non-owning runtime views cross into the common scratch; no
  // detector-shaped table or layer-count dispatcher is retained there.
  o2::its::ROFOverlapTable<NLayers> rofTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    rofTable.defineLayer(layer, layerTiming);
  }
  rofTable.init();

  o2::its::ROFVertexLookupTable<NLayers> vtxTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    vtxTable.defineLayer(layer, layerTiming);
  }
  vtxTable.init();

  o2::its::ROFMaskTable<NLayers> mask{rofTable};
  mask.resetMask();
  for (int layer = 0; layer < NLayers; ++layer) {
    mask.setROFsEnabled(layer, 0, static_cast<int>(nROFsTF), 1);
  }
  mROFOverlapTable = std::move(rofTable);
  mROFVertexLookupTable = std::move(vtxTable);
  mMultiplicityMask = std::move(mask);
  getScratch().setROFViews(RuntimeROFViews{mROFOverlapTable.getView(), mROFVertexLookupTable.getView(), mMultiplicityMask.getView(), mUPCMask.getView()});
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::clearCompatibility() noexcept
{
  mDetectorPublicationAdapter.reset();
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    static_cast<ITSSharedClusterCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }
  if constexpr (DetId == o2::detectors::DetID::MFT) {
    static_cast<MFTPublicationCompatibilityOwner<NLayers>&>(*this).sidecar.clear();
  }
}

template <int NLayers>
bool SurfacePlanTrackingParticipant<NLayers>::refitSeed(const TrackSeed& seed,
                                                        const TrackingParameters& params,
                                                        float bz,
                                                        SurfaceTrackingScratch& scratch,
                                                        gsl::span<const gsl::span<const SurfaceMeasurement>> layerMeasurements,
                                                        SurfaceCatalogView surfaceCatalog,
                                                        ClusterSourceId expectedSource,
                                                        TrackingCandidate& candidate)
{
  return detail::refitDetectorSeed<DetId>(seed, params, bz, scratch, layerMeasurements, surfaceCatalog, expectedSource, candidate);
}

template <int NLayers>
bool SurfacePlanTrackingParticipant<NLayers>::completeAccepted(gsl::span<const TrackingCandidate> candidates,
                                                               const TrackingParameters& params,
                                                               const SurfaceTrackingScratch& scratch,
                                                               bool final)
{
  return mDetectorPublicationAdapter.completeAccepted(candidates, params, scratch, final);
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::clearPublicationSidecar() noexcept
{
  clearCompatibility();
}

template <int NLayers>
gsl::span<const SurfaceId> SurfacePlanTrackingParticipant<NLayers>::ownedSurfaces() const noexcept
{
  if (mFrame == nullptr) {
    return {};
  }
  const auto* binding = mFrame->getBinding(0, mSource);
  return binding == nullptr ? gsl::span<const SurfaceId>{} : binding->getOrderedSurfaces();
}

template <int NLayers>
ParticipantTrackingResult SurfacePlanTrackingParticipant<NLayers>::track(TimeFrame& frame)
{
  // clustersToTracks() itself only ever *returns* Success or
  // RecoverableDropped (Tracker.h); Structural always escapes as a thrown
  // exception instead, left uncaught here for TrackingEngine::executeEvent
  // ()'s own try/catch to classify -- same division of responsibility the
  // coordinator's process() already relied on before this slice.
  const auto result = mTracker.clustersToTracks(*this);
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
  // Generic event data is reset once by the owning TimeFrame. This hook only
  // clears adapter publication state.
  mTracked = false;
  clearCompatibility();
}

template <int NLayers>
std::optional<ParticipantPublicationExport> SurfacePlanTrackingParticipant<NLayers>::publicationExport() const
{
  if (!mTracked || mFrame == nullptr || !mFrame->isConfigured()) {
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
