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
void SurfacePlanTrackingParticipant<NLayers>::setROFViews(RuntimeROFViews views) noexcept
{
  if (mFrame != nullptr) {
    getScratch().setROFViews(views);
  }
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::clearROFViews() noexcept
{
  if (mFrame != nullptr) {
    getScratch().setROFViews({});
  }
}

template <int NLayers>
void SurfacePlanTrackingParticipant<NLayers>::clearCompatibility() noexcept
{
  if (mDetectorPublicationAdapter != nullptr) {
    mDetectorPublicationAdapter->reset();
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
  if (mDetectorPublicationAdapter == nullptr) {
    return true;
  }
  return mDetectorPublicationAdapter->completeAccepted(candidates, params, scratch, final);
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
  return mDetectorPublicationAdapter == nullptr ? nullptr : mDetectorPublicationAdapter->getITSSharedClusterCompatibility();
}

template <int NLayers>
const MFTPublicationCompatibility* SurfacePlanTrackingParticipant<NLayers>::getMFTPublicationCompatibility() const noexcept
{
  return mDetectorPublicationAdapter == nullptr ? nullptr : mDetectorPublicationAdapter->getMFTPublicationCompatibility();
}

template class SurfacePlanTrackingParticipant<ITSNLayers>;
template class SurfacePlanTrackingParticipant<MFTNLayers>;

} // namespace o2::itsmft::tracking
