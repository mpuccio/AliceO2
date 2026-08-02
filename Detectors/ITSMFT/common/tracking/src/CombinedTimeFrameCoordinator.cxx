// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/CombinedTimeFrameCoordinator.h"

#include <stdexcept>
#include <utility>

#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
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

// Same LayerMask(legacy position)->SurfaceMask(global id) conversion
// DetectorLayoutSet.cxx's buildDetectorLayoutSet() uses internally (private,
// not exported): every set bit is a *position* in `orderedSurfaces`, never a
// numeric comparison against the SurfaceId values it contains.
SurfaceMask positionalSurfaceMask(LayerMask layerMask, gsl::span<const SurfaceId> orderedSurfaces)
{
  SurfaceMask result;
  for (uint32_t position = 0; position < orderedSurfaces.size(); ++position) {
    if (layerMask.has(position)) {
      result.set(orderedSurfaces[position]);
    }
  }
  return result;
}

// Builds the one shared ITS+MFT DetectorLayout -- both detectors' own
// disjoint subgraphs in a single DetectorLayoutBuilder call, exactly as
// testCombinedStaticSurfaceCatalogTopology.cxx /
// testDetectorTraversalBindingOrchestration.cxx already prove -- from each
// detector's own (single-iteration) TrackingParameters. Called twice (once
// per detector's own DetectorLayoutSet, which owns its DetectorLayout by
// value): a pure function of its arguments, so both calls produce
// byte-identical topology content.
DetectorLayoutBuildResult buildCombinedLayout(gsl::span<const SurfaceId> itsSurfaces, const TrackingParameters& itsParams,
                                              gsl::span<const SurfaceId> mftSurfaces, const TrackingParameters& mftParams)
{
  DetectorLayoutSubgraph itsSubgraph;
  itsSubgraph.orderedSurfaces.assign(itsSurfaces.begin(), itsSurfaces.end());
  itsSubgraph.maxHoles = itsParams.MaxHoles;
  itsSubgraph.holeSurfaces = positionalSurfaceMask(itsParams.HoleLayerMask, itsSurfaces);
  itsSubgraph.seedingSurfaces = positionalSurfaceMask(itsParams.StartLayerMask, itsSurfaces);
  itsSubgraph.policyTag = TransitionPolicyTag::CylinderCylinder;

  DetectorLayoutSubgraph mftSubgraph;
  mftSubgraph.orderedSurfaces.assign(mftSurfaces.begin(), mftSurfaces.end());
  mftSubgraph.maxHoles = mftParams.MaxHoles;
  mftSubgraph.holeSurfaces = positionalSurfaceMask(mftParams.HoleLayerMask, mftSurfaces);
  mftSubgraph.seedingSurfaces = positionalSurfaceMask(mftParams.StartLayerMask, mftSurfaces);
  mftSubgraph.policyTag = TransitionPolicyTag::DiskDisk;

  DetectorLayoutBuilder builder{combinedCatalogView()};
  builder.addSubgraph(std::move(itsSubgraph));
  builder.addSubgraph(std::move(mftSubgraph));
  return builder.build();
}

// Wraps one buildCombinedLayout() build (both detectors' subgraphs) into the
// DetectorLayoutSet `NLayers` will actually adopt: its key's orderedSurfaces
// is that detector's own global-id span, so TrackerTraits<NLayers>::
// initialiseTimeFrame()'s legacy-layer/material binding (step 2.5,
// TrackerTraits.cxx) resolves against the right detector, while the layout
// content itself carries both detectors' topology for the adopted
// DetectorTraversalBinding to scope. Each caller supplies its own build
// (buildCombinedLayout() is a pure function of its arguments, so the two
// calls this constructor makes produce byte-identical topology content).
template <int NLayers>
DetectorLayoutSet ownDetectorPlan(DetectorLayoutBuildResult built, gsl::span<const SurfaceId> ownSurfaces,
                                  const TrackingParameters& ownParams, TransitionPolicyTag ownPolicy)
{
  if (!built.ok()) {
    throw std::runtime_error("CombinedTimeFrameCoordinator: failed to build the shared ITS+MFT DetectorLayout");
  }
  DetectorLayoutConfigurationKey key;
  key.orderedSurfaces.assign(ownSurfaces.begin(), ownSurfaces.end());
  key.policyTag = ownPolicy;
  key.iterations.push_back(DetectorLayoutIterationConfiguration{
    static_cast<uint32_t>(NLayers), ownParams.MaxHoles, ownParams.HoleLayerMask, ownParams.StartLayerMask});
  std::vector<DetectorLayout> layouts;
  layouts.push_back(std::move(*built.layout));
  return DetectorLayoutSet{std::move(key), combinedCatalogView(), std::move(layouts)};
}

template <int NLayers>
void configureRofTables(LegacyTrackerScratch<NLayers>& scratch, const ROFTimingConfig& timing, uint32_t nROFsTF,
                        const std::vector<TrackingParameters>& params)
{
  o2::its::LayerTiming layerTiming{};
  layerTiming.mNROFsTF = nROFsTF;
  layerTiming.mROFLength = timing.rofLength;
  layerTiming.mROFDelay = timing.rofDelay;
  layerTiming.mROFBias = timing.rofBias;
  layerTiming.mROFAddTimeErr = timing.rofAddTimeErr;

  typename LegacyTrackerScratch<NLayers>::ROFOverlapTableN rofTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    rofTable.defineLayer(layer, layerTiming);
  }
  rofTable.init();
  scratch.setROFOverlapTable(rofTable);

  typename LegacyTrackerScratch<NLayers>::ROFVertexLookupTableN vtxTable;
  for (int layer = 0; layer < NLayers; ++layer) {
    vtxTable.defineLayer(layer, layerTiming);
  }
  vtxTable.init();
  scratch.setROFVertexLookupTable(vtxTable);

  typename LegacyTrackerScratch<NLayers>::ROFMaskTableN mask{rofTable};
  mask.resetMask();
  for (int layer = 0; layer < NLayers; ++layer) {
    mask.setROFsEnabled(layer, 0, static_cast<int>(nROFsTF), 1);
  }
  scratch.setMultiplicityCutMask(std::move(mask));

  scratch.initTrackerTopologies(params);
}

} // namespace

CombinedTimeFrameCoordinator::CombinedTimeFrameCoordinator(std::vector<o2::itsmft::TrackingParameters> itsParams,
                                                           std::vector<o2::itsmft::TrackingParameters> mftParams)
  : mITSParams(std::move(itsParams)), mMFTParams(std::move(mftParams)), mITSTracker(&mITSTraits), mMFTTracker(&mMFTTraits)
{
  if (mITSParams.size() != 1 || mMFTParams.size() != 1) {
    throw std::invalid_argument("CombinedTimeFrameCoordinator requires exactly one TrackingParameters iteration per detector");
  }

  const auto itsSurfaces = orderedSurfaceRange(0, ITSNLayers);
  const auto mftSurfaces = orderedSurfaceRange(ITSNLayers, MFTNLayers);

  mITSPlan.emplace(ownDetectorPlan<ITSNLayers>(buildCombinedLayout(itsSurfaces, mITSParams[0], mftSurfaces, mMFTParams[0]),
                                               itsSurfaces, mITSParams[0], TransitionPolicyTag::CylinderCylinder));
  mMFTPlan.emplace(ownDetectorPlan<MFTNLayers>(buildCombinedLayout(itsSurfaces, mITSParams[0], mftSurfaces, mMFTParams[0]),
                                               mftSurfaces, mMFTParams[0], TransitionPolicyTag::DiskDisk));

  auto itsBindingResult = DetectorTraversalBinding::build(mITSPlan->getLayoutView(0), o2::detectors::DetID::ITS, ClusterSourceId{0},
                                                          surfaceRangeMask(0, ITSNLayers), itsSurfaces);
  if (!itsBindingResult.ok()) {
    throw std::runtime_error("CombinedTimeFrameCoordinator: failed to build the ITS DetectorTraversalBinding");
  }
  mITSBinding = std::move(itsBindingResult.binding);

  auto mftBindingResult = DetectorTraversalBinding::build(mMFTPlan->getLayoutView(0), o2::detectors::DetID::MFT, ClusterSourceId{1},
                                                          surfaceRangeMask(ITSNLayers, MFTNLayers), mftSurfaces);
  if (!mftBindingResult.ok()) {
    throw std::runtime_error("CombinedTimeFrameCoordinator: failed to build the MFT DetectorTraversalBinding");
  }
  mMFTBinding = std::move(mftBindingResult.binding);

  mITSTraits.adoptDetectorTraversalBinding(mITSBinding.get());
  mMFTTraits.adoptDetectorTraversalBinding(mMFTBinding.get());

  mITSTracker.adoptScratch(mITSScratch);
  mITSTracker.adoptITSSharedClusterCompatibility(mITSCompatibility);
  mITSTracker.adoptDetectorLayoutSet(*mITSPlan);
  mITSTracker.setParameters(mITSParams);

  mMFTTracker.adoptScratch(mMFTScratch);
  mMFTTracker.adoptMFTPublicationCompatibility(mMFTCompatibility);
  mMFTTracker.adoptDetectorLayoutSet(*mMFTPlan);
  mMFTTracker.setParameters(mMFTParams);
}

void CombinedTimeFrameCoordinator::adoptFrame(TimeFrame& frame)
{
  mFrame = &frame;
  mITSTracker.adoptFrame(frame);
  mMFTTracker.adoptFrame(frame);
}

void CombinedTimeFrameCoordinator::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  mMemoryPool = pool;
  if (mFrame != nullptr) {
    mFrame->setMemoryPool(pool);
  }
  mITSScratch.setMemoryPool(pool);
  mMFTScratch.setMemoryPool(pool);
  mITSTraits.setMemoryPool(pool);
  mMFTTraits.setMemoryPool(pool);
  mITSTracker.setMemoryPool(pool);
  mMFTTracker.setMemoryPool(pool);
}

void CombinedTimeFrameCoordinator::setBz(float bz)
{
  if (mFrame != nullptr) {
    mFrame->setBz(bz);
  }
  mITSTracker.setBz(bz);
  mMFTTracker.setBz(bz);
}

void CombinedTimeFrameCoordinator::setNThreads(int n)
{
  mITSTraits.setNThreads(n, mITSArena);
  mMFTTraits.setNThreads(n, mMFTArena);
}

void CombinedTimeFrameCoordinator::resetCombinedEvent() noexcept
{
  MultiSourceTimeFrameLoader::resetITSAndMFTEvent(*mFrame, mITSScratch, mMFTScratch);
  mITSCompatibility.clear();
  mMFTCompatibility.clear();
  mITSClock.reset();
  mMFTClock.reset();
  mPublicationValid = false;
}

CombinedTimeFrameCoordinator::CombinedTrackingResult CombinedTimeFrameCoordinator::process(
  const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource, const o2::InteractionRecord& origin)
{
  if (mFrame == nullptr) {
    throw std::logic_error("CombinedTimeFrameCoordinator::process() called before adoptFrame()");
  }

  mPublicationValid = false;
  mITSClock.reset();
  mMFTClock.reset();

  try {
    const auto loadResult = MultiSourceTimeFrameLoader::loadITSAndMFT(*mFrame, mITSScratch, mMFTScratch, itsSource, mftSource,
                                                                      combinedCatalogView(), origin);
    if (!loadResult.ok()) {
      resetCombinedEvent();
      return {};
    }

    configureRofTables<ITSNLayers>(mITSScratch, itsSource.timing, static_cast<uint32_t>(itsSource.rofs.size()), mITSParams);
    configureRofTables<MFTNLayers>(mMFTScratch, mftSource.timing, static_cast<uint32_t>(mftSource.rofs.size()), mMFTParams);

    // Serial, ITS first: the shared TimeFrame's CommonTrack/TrackClusterIndices
    // storage is append-only (acceptTracks() -> AcceptedTrackShadowPublisher,
    // TrackerTraits.cxx), so this call order alone is what makes accepted
    // CommonTrack publication deterministic ITS-then-MFT.
    const auto itsResult = mITSTracker.clustersToTracks();
    if (itsResult.outcome != TrackingOutcome::Success) {
      resetCombinedEvent();
      return {};
    }

    const auto mftResult = mMFTTracker.clustersToTracks();
    if (mftResult.outcome != TrackingOutcome::Success) {
      resetCombinedEvent();
      return {};
    }
  } catch (...) {
    resetCombinedEvent();
    return {};
  }

  mITSClock.emplace(mITSScratch.getROFOverlapTableView().getClockLayer());
  mMFTClock.emplace(mMFTScratch.getROFOverlapTableView().getClockLayer());
  mPublicationValid = true;

  CombinedTrackingResult result;
  result.outcome = CombinedOutcome::Success;
  result.nITSTracks = mITSScratch.getNumberOfTracks();
  result.nMFTTracks = mMFTScratch.getNumberOfTracks();
  return result;
}

std::optional<CommonTrackPublicationExport> CombinedTimeFrameCoordinator::getITSPublicationExport() const
{
  if (!mPublicationValid || !mITSClock) {
    return std::nullopt;
  }
  return CommonTrackPublicationExport{o2::detectors::DetID::ITS, ClusterSourceId{0}, *mITSClock,
                                      gsl::span<const SurfaceId>{mITSPlan->getConfigurationKey().orderedSurfaces}};
}

std::optional<CommonTrackPublicationExport> CombinedTimeFrameCoordinator::getMFTPublicationExport() const
{
  if (!mPublicationValid || !mMFTClock) {
    return std::nullopt;
  }
  return CommonTrackPublicationExport{o2::detectors::DetID::MFT, ClusterSourceId{1}, *mMFTClock,
                                      gsl::span<const SurfaceId>{mMFTPlan->getConfigurationKey().orderedSurfaces}};
}

} // namespace o2::itsmft::tracking
