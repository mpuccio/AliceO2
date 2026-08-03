// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/ITSMFTLegacyParticipantSet.h"

#include <stdexcept>
#include <utility>

#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"

namespace o2::itsmft::tracking
{

namespace
{

SurfaceCatalogView combinedCatalogView()
{
  return SurfaceCatalogView{kITSMFTCombinedStaticSurfaceCatalog.data(),
                            static_cast<uint32_t>(kITSMFTCombinedStaticSurfaceCatalog.size())};
}

std::vector<SurfaceId> orderedSurfaceRange(uint16_t first, uint16_t count)
{
  std::vector<SurfaceId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

SurfaceMask surfaceRangeMask(uint16_t first, uint16_t count)
{
  SurfaceMask result;
  for (uint16_t i = 0; i < count; ++i) {
    result.set(SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

// Builds the one shared ITS+MFT DetectorLayout -- both detectors' own
// disjoint subgraphs in a single DetectorLayoutBuilder call, exactly as
// testCombinedStaticSurfaceCatalogTopology.cxx /
// testDetectorTraversalBindingOrchestration.cxx already prove -- from each
// detector's own (single-iteration) TrackingParameters. Called exactly once
// by the constructor: this is this set's one authoritative combined-topology
// construction; ownDetectorPlan() below only ever copies its result, never
// rebuilds it.
DetectorLayoutBuildResult buildCombinedLayout(gsl::span<const SurfaceId> itsSurfaces, const TrackingParameters& itsParams,
                                              gsl::span<const SurfaceId> mftSurfaces, const TrackingParameters& mftParams)
{
  DetectorLayoutSubgraph itsSubgraph;
  itsSubgraph.orderedSurfaces.assign(itsSurfaces.begin(), itsSurfaces.end());
  itsSubgraph.maxHoles = itsParams.MaxHoles;
  itsSubgraph.holeSurfaces = positionalSurfaceMask(itsParams.HoleLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));
  itsSubgraph.seedingSurfaces = positionalSurfaceMask(itsParams.StartLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));

  DetectorLayoutSubgraph mftSubgraph;
  mftSubgraph.orderedSurfaces.assign(mftSurfaces.begin(), mftSurfaces.end());
  mftSubgraph.maxHoles = mftParams.MaxHoles;
  mftSubgraph.holeSurfaces = positionalSurfaceMask(mftParams.HoleLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));
  mftSubgraph.seedingSurfaces = positionalSurfaceMask(mftParams.StartLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));

  DetectorLayoutBuilder builder{combinedCatalogView()};
  builder.addSubgraph(std::move(itsSubgraph));
  builder.addSubgraph(std::move(mftSubgraph));
  return builder.build();
}

// Wraps a *copy* of this set's one authoritative combined ITS+MFT
// DetectorLayout (built exactly once, by buildCombinedLayout(), in the
// constructor -- never rebuilt here) into the DetectorLayoutSet `NLayers`
// will actually adopt: its key's orderedSurfaces is that detector's own
// global-id span, so TrackerTraits<NLayers>::initialiseTimeFrame()'s
// legacy-layer/material binding (step 2.5, TrackerTraits.cxx) resolves
// against the right detector, while the layout content itself carries both
// detectors' topology for the adopted DetectorTraversalBinding to scope.
//
// `layouts.push_back(authoritative)` is a plain copy-construction (DetectorLayout/
// SparseTrackingTopology are ordinary copyable value types -- no owning
// resource, no deleted/user-provided special members) of the one already-
// built object, not a second independent construction: DetectorLayoutSet
// has no reference/shared-ownership constructor (it always takes its
// `std::vector<DetectorLayout>` by value), so a passive copy is the only way
// two DetectorLayoutSets can each hold this set's one authoritative
// topology. Both copies are therefore guaranteed byte-identical, including
// every global TransitionId/CellTopologyId's content, by construction --
// there is no code path that could make them diverge.
template <int NLayers>
DetectorLayoutSet ownDetectorPlan(const DetectorLayout& authoritative, gsl::span<const SurfaceId> ownSurfaces,
                                  const TrackingParameters& ownParams)
{
  DetectorLayoutConfigurationKey key;
  key.orderedSurfaces.assign(ownSurfaces.begin(), ownSurfaces.end());
  key.iterations.push_back(DetectorLayoutIterationConfiguration{
    static_cast<uint32_t>(NLayers), ownParams.MaxHoles, ownParams.HoleLayerMask, ownParams.StartLayerMask});
  std::vector<DetectorLayout> layouts;
  layouts.push_back(authoritative);
  return DetectorLayoutSet{std::move(key), combinedCatalogView(), std::move(layouts)};
}

} // namespace

ITSMFTLegacyParticipantSet::ITSMFTLegacyParticipantSet(std::vector<o2::itsmft::TrackingParameters> itsParams,
                                                       std::vector<o2::itsmft::TrackingParameters> mftParams)
  : mITSParticipant(ParticipantId{0}, itsParams), mMFTParticipant(ParticipantId{1}, mftParams)
{
  if (itsParams.size() != 1 || mftParams.size() != 1) {
    throw std::invalid_argument("ITSMFTLegacyParticipantSet requires exactly one TrackingParameters iteration per detector");
  }

  const auto itsSurfaces = orderedSurfaceRange(0, ITSNLayers);
  const auto mftSurfaces = orderedSurfaceRange(ITSNLayers, MFTNLayers);

  // This set's one authoritative combined ITS+MFT topology construction.
  // mITSPlan/mMFTPlan below each get a passive copy of `combinedLayout` (see
  // ownDetectorPlan()'s own doc) -- never a second, independent
  // buildCombinedLayout() call.
  auto combinedBuild = buildCombinedLayout(itsSurfaces, itsParams[0], mftSurfaces, mftParams[0]);
  if (!combinedBuild.ok()) {
    throw std::runtime_error("ITSMFTLegacyParticipantSet: failed to build the shared ITS+MFT DetectorLayout");
  }
  const DetectorLayout& combinedLayout = *combinedBuild.layout;

  mITSPlan.emplace(ownDetectorPlan<ITSNLayers>(combinedLayout, itsSurfaces, itsParams[0]));
  mMFTPlan.emplace(ownDetectorPlan<MFTNLayers>(combinedLayout, mftSurfaces, mftParams[0]));

  auto itsBindingResult = DetectorTraversalBinding::build(mITSPlan->getLayoutView(0), o2::detectors::DetID::ITS, ClusterSourceId{0},
                                                          surfaceRangeMask(0, ITSNLayers), itsSurfaces);
  if (!itsBindingResult.ok()) {
    throw std::runtime_error("ITSMFTLegacyParticipantSet: failed to build the ITS DetectorTraversalBinding");
  }
  auto mftBindingResult = DetectorTraversalBinding::build(mMFTPlan->getLayoutView(0), o2::detectors::DetID::MFT, ClusterSourceId{1},
                                                          surfaceRangeMask(ITSNLayers, MFTNLayers), mftSurfaces);
  if (!mftBindingResult.ok()) {
    throw std::runtime_error("ITSMFTLegacyParticipantSet: failed to build the MFT DetectorTraversalBinding");
  }

  mITSParticipant.adoptDetectorTraversalBinding(std::move(itsBindingResult.binding));
  mMFTParticipant.adoptDetectorTraversalBinding(std::move(mftBindingResult.binding));
  mITSParticipant.adoptDetectorLayoutSet(*mITSPlan);
  mMFTParticipant.adoptDetectorLayoutSet(*mMFTPlan);

  mSchedule = {&mITSParticipant, &mMFTParticipant};
}

void ITSMFTLegacyParticipantSet::adoptFrame(TimeFrame& frame)
{
  mITSParticipant.adoptFrame(frame);
  mMFTParticipant.adoptFrame(frame);
}

void ITSMFTLegacyParticipantSet::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  mITSParticipant.setMemoryPool(pool);
  mMFTParticipant.setMemoryPool(pool);
}

void ITSMFTLegacyParticipantSet::setBz(float bz)
{
  mITSParticipant.setBz(bz);
  mMFTParticipant.setBz(bz);
}

void ITSMFTLegacyParticipantSet::setNThreads(int n)
{
  mITSParticipant.setNThreads(n);
  mMFTParticipant.setNThreads(n);
}

std::optional<LoadSourcesResult> ITSMFTLegacyParticipantSet::validateSources(const ClusterSourceInput& itsSource,
                                                                             const ClusterSourceInput& mftSource) const noexcept
{
  // The fixed ITS=0/MFT=1 source contract lives only here now, not in the
  // generic loadEvent() transaction: a mismatch is synthesized as the exact
  // same LoadSourcesResult loadITSAndMFT()'s equivalent guard used to
  // return, so the coordinator's generic error-handling classifies it
  // identically.
  if (itsSource.id != ClusterSourceId{0} || itsSource.detector != o2::detectors::DetID::ITS) {
    return LoadSourcesResult{MultiSourceLoadError::UnsupportedDetector, itsSource.id};
  }
  if (mftSource.id != ClusterSourceId{1} || mftSource.detector != o2::detectors::DetID::MFT) {
    return LoadSourcesResult{MultiSourceLoadError::UnsupportedDetector, mftSource.id};
  }
  return std::nullopt;
}

std::array<MultiSourceTimeFrameLoader::AtomicLoadBinding, 2> ITSMFTLegacyParticipantSet::loadBindings(
  const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource) noexcept
{
  return {MultiSourceTimeFrameLoader::AtomicLoadBinding{itsSource, mITSParticipant.loadTarget()},
          MultiSourceTimeFrameLoader::AtomicLoadBinding{mftSource, mMFTParticipant.loadTarget()}};
}

SurfaceCatalogView ITSMFTLegacyParticipantSet::catalogView() const noexcept
{
  return combinedCatalogView();
}

std::optional<bool> ITSMFTLegacyParticipantSet::dropTFUponFailureFor(ClusterSourceId source) const noexcept
{
  if (source == ClusterSourceId{0}) {
    return mITSParticipant.getDropTFUponFailure();
  }
  if (source == ClusterSourceId{1}) {
    return mMFTParticipant.getDropTFUponFailure();
  }
  return std::nullopt;
}

void ITSMFTLegacyParticipantSet::configureRofTables(const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource)
{
  mITSParticipant.configureRofTables(itsSource.timing, static_cast<uint32_t>(itsSource.rofs.size()));
  mMFTParticipant.configureRofTables(mftSource.timing, static_cast<uint32_t>(mftSource.rofs.size()));
}

void ITSMFTLegacyParticipantSet::clearPublicationSidecars() noexcept
{
  mITSParticipant.clearPublicationSidecar();
  mMFTParticipant.clearPublicationSidecar();
}

void ITSMFTLegacyParticipantSet::invalidatePublication() noexcept
{
  mITSClock.reset();
  mMFTClock.reset();
  mPublicationValid = false;
}

void ITSMFTLegacyParticipantSet::markPublicationValid() noexcept
{
  mITSClock.emplace(mITSParticipant.getScratch().getROFOverlapTableView().getClockLayer());
  mMFTClock.emplace(mMFTParticipant.getScratch().getROFOverlapTableView().getClockLayer());
  mPublicationValid = true;
}

std::optional<CommonTrackPublicationExport> ITSMFTLegacyParticipantSet::getITSPublicationExport() const
{
  if (!mPublicationValid || !mITSClock) {
    return std::nullopt;
  }
  return CommonTrackPublicationExport{o2::detectors::DetID::ITS, ClusterSourceId{0}, *mITSClock,
                                      gsl::span<const SurfaceId>{mITSPlan->getConfigurationKey().orderedSurfaces}};
}

std::optional<CommonTrackPublicationExport> ITSMFTLegacyParticipantSet::getMFTPublicationExport() const
{
  if (!mPublicationValid || !mMFTClock) {
    return std::nullopt;
  }
  return CommonTrackPublicationExport{o2::detectors::DetID::MFT, ClusterSourceId{1}, *mMFTClock,
                                      gsl::span<const SurfaceId>{mMFTPlan->getConfigurationKey().orderedSurfaces}};
}

} // namespace o2::itsmft::tracking
