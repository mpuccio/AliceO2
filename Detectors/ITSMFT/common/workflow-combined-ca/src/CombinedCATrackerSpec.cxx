// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTCombinedCAWorkflow/CombinedCATrackerSpec.h"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

#include <gsl/span>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/Logger.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTCombinedCAWorkflow/ConfigPreflight.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/CommonTrackOutputAdapter.h"
#include "ITSMFTTracking/DetectorLayoutBuilder.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/ROFTimingUniformity.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/detail/SurfacePlanBinding.h"
#include "ITSMFTTracking/TimeFrameLoadFailure.h"
#include "ITSMFTTracking/TrackingConfigParam.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTTracking/MFTTrackingParam.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

using namespace o2::framework;

namespace o2::itsmft::combined
{

namespace
{

using o2::itsmft::tracking::ClusterSourceId;
using o2::itsmft::tracking::ClusterSourceInput;
using o2::itsmft::tracking::LoadSourcesResult;
using o2::itsmft::tracking::MultiSourceTimeFrameLoader;
using o2::itsmft::tracking::ParticipantOutcome;

// Per-layer DPLAlpideParam<DetId> -> one source-level ROFTimingConfig,
// fatal on non-positive ROF length or non-uniform per-layer staggering --
// same classification ITSMFTTrackingInterface::configureROFLookupTables()
// applies (TrackingInterface.cxx), reused via the same shared
// deriveUniformROFTimingConfig() (ROFTimingUniformity.h). Unlike that
// single-detector path, this workflow does not need to derive nROFsTF here:
// The workflow-owned configureRofTables() (SurfacePlanTrackingParticipant
// .cxx) takes it directly as the workflow's own ClusterSourceInput::rofs
// .size(), the actual ROF count this TF's DPL input already carries.
template <o2::detectors::DetID::ID DetId, int NLayers>
o2::itsmft::tracking::ROFTimingConfig deriveRofTimingConfigOrFatal(const o2::itsmft::TrackingParameters& params)
{
  constexpr const char* detName = DetId == o2::detectors::DetID::ITS ? "ITS" : "MFT";
  const auto& par = o2::itsmft::DPLAlpideParam<DetId>::Instance();
  std::array<o2::its::LayerTiming, NLayers> layerTimings{};
  for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
    const auto rofLengthInBC = par.getROFLengthInBC(iLayer);
    if (rofLengthInBC <= 0) {
      LOGP(fatal, "Combined ITS+MFT tracker: {} per-layer ROF timing configuration has a non-positive ROF length ({}) on layer {}",
           detName, rofLengthInBC, iLayer);
    }
    layerTimings[iLayer].mNROFsTF = 0;
    layerTimings[iLayer].mROFLength = static_cast<uint32_t>(rofLengthInBC);
    layerTimings[iLayer].mROFDelay = static_cast<uint32_t>(par.getROFDelayInBC(iLayer));
    layerTimings[iLayer].mROFBias = static_cast<uint32_t>(par.getROFBiasInBC(iLayer));
    layerTimings[iLayer].mROFAddTimeErr = static_cast<uint32_t>(params.AddTimeError[iLayer]);
  }
  const auto uniform = o2::itsmft::tracking::deriveUniformROFTimingConfig(layerTimings);
  if (!uniform.uniform) {
    LOGP(fatal,
         "Combined ITS+MFT tracker: {} per-layer ROF timing configuration is not uniform; source-level "
         "normalized loading requires equal length/delay/bias/addTimeErr across all {} layers",
         detName, NLayers);
  }
  return uniform.config;
}

o2::InteractionRecord chooseOrigin(gsl::span<const o2::itsmft::ROFRecord> itsRofs, gsl::span<const o2::itsmft::ROFRecord> mftRofs)
{
  if (!itsRofs.empty()) {
    return itsRofs.front().getBCData();
  }
  if (!mftRofs.empty()) {
    return mftRofs.front().getBCData();
  }
  return o2::InteractionRecord{};
}

// The combined ITS/MFT application plan is workflow-owned. It is deliberately
// built here, beside the DPL task that owns its source positions and schedule,
// rather than hidden behind a library-level coordinator.
o2::itsmft::tracking::SurfaceCatalogView combinedCatalogView()
{
  return {o2::itsmft::tracking::kITSMFTCombinedStaticSurfaceCatalog.data(),
          static_cast<uint32_t>(o2::itsmft::tracking::kITSMFTCombinedStaticSurfaceCatalog.size())};
}

std::vector<o2::itsmft::tracking::SurfaceId> orderedSurfaceRange(uint16_t first, uint16_t count)
{
  std::vector<o2::itsmft::tracking::SurfaceId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(o2::itsmft::tracking::SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

o2::itsmft::tracking::SurfaceMask surfaceRangeMask(uint16_t first, uint16_t count)
{
  o2::itsmft::tracking::SurfaceMask result;
  for (uint16_t i = 0; i < count; ++i) {
    result.set(o2::itsmft::tracking::SurfaceId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

o2::itsmft::tracking::DetectorLayoutBuildResult buildCombinedLayout(
  gsl::span<const o2::itsmft::tracking::SurfaceId> itsSurfaces,
  const o2::itsmft::TrackingParameters& itsParams,
  gsl::span<const o2::itsmft::tracking::SurfaceId> mftSurfaces,
  const o2::itsmft::TrackingParameters& mftParams)
{
  o2::itsmft::tracking::DetectorLayoutSubgraph itsSubgraph;
  itsSubgraph.orderedSurfaces.assign(itsSurfaces.begin(), itsSurfaces.end());
  itsSubgraph.maxHoles = itsParams.MaxHoles;
  itsSubgraph.holeSurfaces = o2::itsmft::tracking::positionalSurfaceMask(
    itsParams.HoleLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));
  itsSubgraph.seedingSurfaces = o2::itsmft::tracking::positionalSurfaceMask(
    itsParams.StartLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));

  o2::itsmft::tracking::DetectorLayoutSubgraph mftSubgraph;
  mftSubgraph.orderedSurfaces.assign(mftSurfaces.begin(), mftSurfaces.end());
  mftSubgraph.maxHoles = mftParams.MaxHoles;
  mftSubgraph.holeSurfaces = o2::itsmft::tracking::positionalSurfaceMask(
    mftParams.HoleLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));
  mftSubgraph.seedingSurfaces = o2::itsmft::tracking::positionalSurfaceMask(
    mftParams.StartLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));

  o2::itsmft::tracking::DetectorLayoutBuilder builder{combinedCatalogView()};
  builder.addSubgraph(std::move(itsSubgraph));
  builder.addSubgraph(std::move(mftSubgraph));
  return builder.build();
}

template <int NLayers>
o2::itsmft::tracking::DetectorLayoutSet ownDetectorPlan(
  const o2::itsmft::tracking::DetectorLayout& authoritative,
  gsl::span<const o2::itsmft::tracking::SurfaceId> ownSurfaces,
  const o2::itsmft::TrackingParameters& params)
{
  o2::itsmft::tracking::DetectorLayoutConfigurationKey key;
  key.orderedSurfaces.assign(ownSurfaces.begin(), ownSurfaces.end());
  key.iterations.push_back(o2::itsmft::tracking::DetectorLayoutIterationConfiguration{
    static_cast<uint32_t>(NLayers), params.MaxHoles, params.HoleLayerMask, params.StartLayerMask});
  std::vector<o2::itsmft::tracking::DetectorLayout> layouts;
  layouts.push_back(authoritative);
  return o2::itsmft::tracking::DetectorLayoutSet{std::move(key), combinedCatalogView(), std::move(layouts)};
}

} // namespace

CombinedCATrackerDPL::CombinedCATrackerDPL(std::shared_ptr<o2::base::GRPGeomRequest> gr, bool useMC)
  : mGGCCDBRequest(std::move(gr)), mUseMC(useMC)
{
  mITSDecoder = std::make_unique<o2::itsmft::tracking::GeometryClusterDecoder<o2::detectors::DetID::ITS>>();
  mMFTDecoder = std::make_unique<o2::itsmft::tracking::GeometryClusterDecoder<o2::detectors::DetID::MFT>>();
}

void CombinedCATrackerDPL::init(InitContext&)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
}

void CombinedCATrackerDPL::buildParticipantsOnce()
{
  if (mParticipantsBuilt) {
    return;
  }
  mParticipantsBuilt = true;

  // Sync-only, single iteration per detector -- enforced already by
  // TrackingMode::getTrackingParameters(ITS, Sync) (fatal on any other
  // mode) and by requireSyncTrackingModeOrFatal()/
  // requireDiamondVertexConstraintOrFatal() at the workflow's own preflight
  // (itsmft-combined-ca-tracker-workflow.cxx), before this device was even
  // constructed.
  auto itsParams = o2::itsmft::TrackingMode::getTrackingParameters(o2::detectors::DetID::ITS, o2::itsmft::TrackingMode::Sync);
  auto mftParams = o2::itsmft::TrackingMode::getTrackingParameters(o2::detectors::DetID::MFT, o2::itsmft::TrackingMode::Sync);
  mITSTrackingParams = itsParams[0];
  mMFTTrackingParams = mftParams[0];
  const size_t maxMemory = std::min(itsParams[0].MaxMemory, mftParams[0].MaxMemory);

  const auto itsSurfaces = orderedSurfaceRange(0, o2::itsmft::tracking::ITSNLayers);
  const auto mftSurfaces = orderedSurfaceRange(o2::itsmft::tracking::ITSNLayers, o2::itsmft::tracking::MFTNLayers);
  const auto combinedBuild = buildCombinedLayout(itsSurfaces, itsParams[0], mftSurfaces, mftParams[0]);
  if (!combinedBuild.ok()) {
    throw std::runtime_error{"Combined ITS+MFT tracker: failed to build the application detector layout"};
  }
  const auto& combinedLayout = *combinedBuild.layout;
  mITSPlan.emplace(ownDetectorPlan<o2::itsmft::tracking::ITSNLayers>(combinedLayout, itsSurfaces, itsParams[0]));
  mMFTPlan.emplace(ownDetectorPlan<o2::itsmft::tracking::MFTNLayers>(combinedLayout, mftSurfaces, mftParams[0]));

  mITSParticipant = std::make_unique<o2::itsmft::tracking::SurfacePlanTrackingParticipantITS>(
    o2::itsmft::tracking::ParticipantId{0}, std::move(itsParams));
  mMFTParticipant = std::make_unique<o2::itsmft::tracking::SurfacePlanTrackingParticipantMFT>(
    o2::itsmft::tracking::ParticipantId{1}, std::move(mftParams));

  auto itsBindingResult = o2::itsmft::tracking::SurfacePlanBinding::build(
    mITSPlan->getLayoutView(0), o2::itsmft::tracking::ClusterSourceId{0},
    surfaceRangeMask(0, o2::itsmft::tracking::ITSNLayers), itsSurfaces,
    o2::itsmft::tracking::SurfaceKind::Cylinder, o2::itsmft::tracking::TransitionPolicyTag::CylinderCylinder);
  if (!itsBindingResult.ok()) {
    throw std::runtime_error{"Combined ITS+MFT tracker: failed to build the ITS surface binding"};
  }
  auto mftBindingResult = o2::itsmft::tracking::SurfacePlanBinding::build(
    mMFTPlan->getLayoutView(0), o2::itsmft::tracking::ClusterSourceId{1},
    surfaceRangeMask(o2::itsmft::tracking::ITSNLayers, o2::itsmft::tracking::MFTNLayers), mftSurfaces,
    o2::itsmft::tracking::SurfaceKind::Disk, o2::itsmft::tracking::TransitionPolicyTag::DiskDisk);
  if (!mftBindingResult.ok()) {
    throw std::runtime_error{"Combined ITS+MFT tracker: failed to build the MFT surface binding"};
  }
  mITSParticipant->adoptSurfacePlanBinding(std::move(itsBindingResult.binding));
  mMFTParticipant->adoptSurfacePlanBinding(std::move(mftBindingResult.binding));
  mITSParticipant->adoptDetectorLayoutSet(*mITSPlan);
  mMFTParticipant->adoptDetectorLayoutSet(*mMFTPlan);
  mSchedule = {mITSParticipant.get(), mMFTParticipant.get()};

  mFrame.setMemoryPool(std::make_shared<o2::itsmft::tracking::BoundedMemoryResource>(maxMemory));
  mITSParticipant->adoptFrame(mFrame);
  mMFTParticipant->adoptFrame(mFrame);
  mITSParticipant->setMemoryPool(mFrame.getMemoryPool());
  mMFTParticipant->setMemoryPool(mFrame.getMemoryPool());
  mITSParticipant->setBz(o2::base::Propagator::Instance()->getNominalBz());
  mMFTParticipant->setBz(o2::base::Propagator::Instance()->getNominalBz());

  const int itsNThreads = o2::itsmft::ITSCommonCATrackerParam::Instance().nThreads;
  const int mftNThreads = o2::itsmft::tracking::TrackerParamRef<o2::detectors::DetID::MFT>::get().nThreads;
  const int nThreads = std::max({1, itsNThreads, mftNThreads});
  mITSParticipant->setNThreads(nThreads);
  mMFTParticipant->setNThreads(nThreads);
}

std::optional<LoadSourcesResult> CombinedCATrackerDPL::validateSources(const ClusterSourceInput& itsSource,
                                                                        const ClusterSourceInput& mftSource) const noexcept
{
  if (itsSource.id != ClusterSourceId{0} || itsSource.detector != o2::detectors::DetID::ITS) {
    return LoadSourcesResult{o2::itsmft::tracking::MultiSourceLoadError::UnsupportedDetector, itsSource.id};
  }
  if (mftSource.id != ClusterSourceId{1} || mftSource.detector != o2::detectors::DetID::MFT) {
    return LoadSourcesResult{o2::itsmft::tracking::MultiSourceLoadError::UnsupportedDetector, mftSource.id};
  }
  return std::nullopt;
}

std::array<MultiSourceTimeFrameLoader::AtomicLoadBinding, 2> CombinedCATrackerDPL::loadBindings(
  const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource) noexcept
{
  return {MultiSourceTimeFrameLoader::AtomicLoadBinding{itsSource, mITSParticipant->loadTarget()},
          MultiSourceTimeFrameLoader::AtomicLoadBinding{mftSource, mMFTParticipant->loadTarget()}};
}

o2::itsmft::tracking::SurfaceCatalogView CombinedCATrackerDPL::catalogView() const noexcept
{
  return combinedCatalogView();
}

std::optional<bool> CombinedCATrackerDPL::dropTFUponFailureFor(ClusterSourceId source) const noexcept
{
  if (source == ClusterSourceId{0}) {
    return mITSParticipant->getDropTFUponFailure();
  }
  if (source == ClusterSourceId{1}) {
    return mMFTParticipant->getDropTFUponFailure();
  }
  return std::nullopt;
}

void CombinedCATrackerDPL::configureRofTables(const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource)
{
  mITSParticipant->configureRofTables(itsSource.timing, static_cast<uint32_t>(itsSource.rofs.size()));
  mMFTParticipant->configureRofTables(mftSource.timing, static_cast<uint32_t>(mftSource.rofs.size()));
}

void CombinedCATrackerDPL::clearPublicationSidecars() noexcept
{
  mITSParticipant->clearPublicationSidecar();
  mMFTParticipant->clearPublicationSidecar();
}

void CombinedCATrackerDPL::invalidatePublication() noexcept
{
  mITSClock.reset();
  mMFTClock.reset();
  mPublicationValid = false;
}

void CombinedCATrackerDPL::markPublicationValid() noexcept
{
  mITSClock.emplace(mITSParticipant->getScratch().getROFOverlapTableView<o2::itsmft::tracking::ITSNLayers>().getClockLayer());
  mMFTClock.emplace(mMFTParticipant->getScratch().getROFOverlapTableView<o2::itsmft::tracking::MFTNLayers>().getClockLayer());
  mPublicationValid = true;
}

std::optional<o2::itsmft::tracking::CommonTrackPublicationExport> CombinedCATrackerDPL::getITSPublicationExport() const
{
  if (!mPublicationValid || !mITSClock) {
    return std::nullopt;
  }
  return o2::itsmft::tracking::CommonTrackPublicationExport{o2::detectors::DetID::ITS, ClusterSourceId{0}, *mITSClock,
                                                            getITSOrderedSurfaces()};
}

std::optional<o2::itsmft::tracking::CommonTrackPublicationExport> CombinedCATrackerDPL::getMFTPublicationExport() const
{
  if (!mPublicationValid || !mMFTClock) {
    return std::nullopt;
  }
  return o2::itsmft::tracking::CommonTrackPublicationExport{o2::detectors::DetID::MFT, ClusterSourceId{1}, *mMFTClock,
                                                            getMFTOrderedSurfaces()};
}

ParticipantOutcome CombinedCATrackerDPL::trackFrame(const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource,
                                                    const o2::InteractionRecord& origin)
{
  invalidatePublication();
  // A prior successful trackFrame() call leaves both sidecars sealed
  // (ITSSharedClusterCompatibility::sealFromMarkedTracks(),
  // MFTPublicationCompatibility's per-TF entries): unlike
  // TimeFrame::getCommonTracks()/getTrackClusterIndices(), which
  // commitNormalizedFrame() clears atomically on the *next* load, neither
  // sidecar is cleared by a successful trackFrame() itself. Left sealed,
  // the very next TF's first accepted ITS track would fail
  // AcceptedTrackShadowPublisher<ITSNLayers>::publish()'s already-sealed
  // guard and fatal inside TrackerTraits<NLayers>::acceptTracks() ("CommonTrack
  // shadow construction failed"). Clearing both unconditionally at the top
  // of every trackFrame() call -- success or failure alike -- keeps every
  // TF starting from the same fresh state regardless of how the previous
  // one ended.
  clearPublicationSidecars();

  // The fixed ITS=0/MFT=1 source contract lives in this workflow-owned
  // application composition, never in the generic loader.
  LoadSourcesResult loadResult;
  if (const auto rejected = validateSources(itsSource, mftSource)) {
    loadResult = *rejected;
  } else {
    const auto bindings = loadBindings(itsSource, mftSource);
    loadResult = MultiSourceTimeFrameLoader::loadEvent(mFrame, gsl::span<const MultiSourceTimeFrameLoader::AtomicLoadBinding>{bindings},
                                                       catalogView(), origin);
  }
  if (!loadResult.ok()) {
    // Reuse isRecoverableLoadError() (TimeFrameLoadFailure.h) rather than a
    // parallel taxonomy, then gate it by the *owning* detector's own
    // DropTFUponFailure -- the workflow-owned source mapping carries
    // the fixed ITS/MFT source-position contract. This is a *load*
    // failure: the event was never atomically committed, so
    // TrackingEngine::executeEvent() must never be called on it --
    // resetEvent() alone applies the same all-participant/shared-frame
    // reset contract without ever reaching track().
    const bool errorIsRecoverable = o2::itsmft::tracking::isRecoverableLoadError(loadResult.error, loadResult.timingDetail);
    const auto dropAllowed = dropTFUponFailureFor(loadResult.source);
    const bool sourceRecognized = dropAllowed.has_value();
    if (!sourceRecognized) {
      LOGP(error, "Combined TF load failure reports an unrecognized source id {}; treating as structural", loadResult.source.value());
    }
    const auto outcome = errorIsRecoverable && sourceRecognized && dropAllowed.value_or(false)
                           ? ParticipantOutcome::RecoverableDropped
                           : ParticipantOutcome::Structural;
    LOGP(error, "Combined TF loading failed (source={}, error={}, recoverable={}, dropAllowed={}): outcome={}",
         loadResult.source.value(), static_cast<int>(loadResult.error), errorIsRecoverable, dropAllowed.value_or(false),
         outcome == ParticipantOutcome::RecoverableDropped ? "RecoverableDropped" : "Structural");
    mEngine.resetEvent(mFrame, schedule());
    invalidatePublication();
    return outcome;
  }

  configureRofTables(itsSource, mftSource);

  // The load has committed: executeEvent() may now run. It executes the
  // explicit [ITS, MFT] schedule's track() calls in that exact order into
  // the shared TimeFrame (accepted CommonTracks therefore append
  // ITS-then-MFT), and on any non-Success outcome or exception already
  // applies the same whole-event reset the load-failure branch above
  // applies -- see TrackingEngine::executeEvent()'s own doc.
  const auto eventResult = mEngine.executeEvent(mFrame, schedule());
  if (eventResult.outcome != ParticipantOutcome::Success) {
    LOGP(error, "Combined TF tracking failed via the delegated engine (outcome={})",
         eventResult.outcome == ParticipantOutcome::RecoverableDropped ? "RecoverableDropped" : "Structural");
    // executeEvent() already reset every participant and wiped the shared
    // TimeFrame; only the publication/timing bridge state remains to
    // invalidate.
    invalidatePublication();
    return eventResult.outcome;
  }

  markPublicationValid();
  return ParticipantOutcome::Success;
}

void CombinedCATrackerDPL::run(ProcessingContext& pc)
{
  updateTimeDependentParams(pc);
  buildParticipantsOnce();

  auto itsRofsInput = pc.inputs().get<const std::vector<o2::itsmft::ROFRecord>>("ROframesITS");
  auto mftRofsInput = pc.inputs().get<const std::vector<o2::itsmft::ROFRecord>>("ROframesMFT");
  auto itsCompClusters = pc.inputs().get<const std::vector<o2::itsmft::CompClusterExt>>("compClustersITS");
  auto mftCompClusters = pc.inputs().get<const std::vector<o2::itsmft::CompClusterExt>>("compClustersMFT");
  gsl::span<const unsigned char> itsPatterns = pc.inputs().get<gsl::span<unsigned char>>("patternsITS");
  gsl::span<const unsigned char> mftPatterns = pc.inputs().get<gsl::span<unsigned char>>("patternsMFT");

  const dataformats::MCTruthContainer<MCCompLabel>* itsLabels = nullptr;
  const dataformats::MCTruthContainer<MCCompLabel>* mftLabels = nullptr;
  if (mUseMC) {
    if (pc.inputs().getPos("labelsITS") >= 0) {
      itsLabels = pc.inputs().get<const dataformats::MCTruthContainer<MCCompLabel>*>("labelsITS").release();
    }
    if (pc.inputs().getPos("labelsMFT") >= 0) {
      mftLabels = pc.inputs().get<const dataformats::MCTruthContainer<MCCompLabel>*>("labelsMFT").release();
    }
  }

  LOGP(info, "Combined ITS+MFT CA input pulled {} ITS clusters ({} ROFs), {} MFT clusters ({} ROFs)",
       itsCompClusters.size(), itsRofsInput.size(), mftCompClusters.size(), mftRofsInput.size());

  const gsl::span<const o2::itsmft::ROFRecord> itsRofs(itsRofsInput.data(), itsRofsInput.size());
  const gsl::span<const o2::itsmft::ROFRecord> mftRofs(mftRofsInput.data(), mftRofsInput.size());

  // Direct field-by-field construction: each source input is workflow-owned
  // and built fresh for this TF.
  ClusterSourceInput itsSource{};
  itsSource.id = ClusterSourceId{0};
  itsSource.detector = o2::detectors::DetID::ITS;
  itsSource.clusters = gsl::span<const o2::itsmft::CompClusterExt>(itsCompClusters.data(), itsCompClusters.size());
  itsSource.patterns = itsPatterns;
  itsSource.rofs = itsRofs;
  itsSource.dictionary = mITSDict;
  itsSource.labels = itsLabels;
  // layerToSurface built directly from the workflow-owned plan's always-valid
  // ordered-surface getter -- never re-derived by hand as a literal
  // ITS=0..6/MFT=7..16 offset.
  itsSource.layerToSurface = getITSOrderedSurfaces();
  itsSource.timing = deriveRofTimingConfigOrFatal<o2::detectors::DetID::ITS, o2::itsmft::tracking::ITSNLayers>(mITSTrackingParams);
  itsSource.decoder = mITSDecoder.get();

  ClusterSourceInput mftSource{};
  mftSource.id = ClusterSourceId{1};
  mftSource.detector = o2::detectors::DetID::MFT;
  mftSource.clusters = gsl::span<const o2::itsmft::CompClusterExt>(mftCompClusters.data(), mftCompClusters.size());
  mftSource.patterns = mftPatterns;
  mftSource.rofs = mftRofs;
  mftSource.dictionary = mMFTDict;
  mftSource.labels = mftLabels;
  mftSource.layerToSurface = getMFTOrderedSurfaces();
  mftSource.timing = deriveRofTimingConfigOrFatal<o2::detectors::DetID::MFT, o2::itsmft::tracking::MFTNLayers>(mMFTTrackingParams);
  mftSource.decoder = mMFTDecoder.get();

  const auto origin = chooseOrigin(itsRofs, mftRofs);
  const auto outcome = trackFrame(itsSource, mftSource, origin);

  if (outcome == ParticipantOutcome::RecoverableDropped) {
    LOGP(error,
         "Combined ITS+MFT CA tracking recoverably dropped this TF ({} ITS ROFs/{} ITS clusters, {} MFT "
         "ROFs/{} MFT clusters); publishing nothing and continuing with the next TimeFrame",
         itsRofs.size(), itsCompClusters.size(), mftRofs.size(), mftCompClusters.size());
    return;
  }
  if (outcome == ParticipantOutcome::Structural) {
    throw std::runtime_error{"Combined ITS+MFT CA tracking hit a structural failure"};
  }

  const auto itsExport = getITSPublicationExport();
  const auto mftExport = getMFTPublicationExport();
  if (!itsExport || !mftExport) {
    throw std::runtime_error{"Combined ITS+MFT CommonTrack output publication context is unavailable after a successful trackFrame()"};
  }
  const o2::itsmft::tracking::CommonTrackPublicationContext itsContext{
    itsExport->detector, itsExport->source, itsRofs, itsExport->clock, itsExport->orderedSurfaces};
  const o2::itsmft::tracking::CommonTrackPublicationContext mftContext{
    mftExport->detector, mftExport->source, mftRofs, mftExport->clock, mftExport->orderedSurfaces};

  // Stage both detectors' outputs completely before any pc.outputs() call:
  // if either staging fails, neither stream is ever requested.
  o2::itsmft::tracking::CommonTrackOutputAdapterError itsError = o2::itsmft::tracking::CommonTrackOutputAdapterError::None;
  const auto stagedITS = o2::itsmft::tracking::stageITSCommonTrackOutput(
    mFrame, itsContext, getITSSharedClusterCompatibility(), mUseMC, itsError);
  o2::itsmft::tracking::CommonTrackOutputAdapterError mftError = o2::itsmft::tracking::CommonTrackOutputAdapterError::None;
  const auto stagedMFT = o2::itsmft::tracking::stageMFTCommonTrackOutput(
    mFrame, mftContext, getMFTPublicationCompatibility(), mUseMC, mftError);
  if (!stagedITS || !stagedMFT) {
    throw std::runtime_error{"Combined ITS+MFT CommonTrack output staging failed"};
  }

  auto& itsTrackROFs = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{"ITS", "ITSTrackROF", 0},
                                                                             stagedITS->trackROFs.begin(), stagedITS->trackROFs.end());
  auto& itsTracks = pc.outputs().make<std::vector<o2::its::TrackITS>>(Output{"ITS", "TRACKS", 0});
  itsTracks.assign(stagedITS->tracks.begin(), stagedITS->tracks.end());
  auto& itsClusIdx = pc.outputs().make<std::vector<int>>(Output{"ITS", "TRACKCLSID", 0});
  itsClusIdx.assign(stagedITS->clusterIndices.begin(), stagedITS->clusterIndices.end());
  LOGP(info, "Combined ITS+MFT CA pushed {} ITS tracks in {} ROFs", itsTracks.size(), itsTrackROFs.size());
  if (mUseMC) {
    pc.outputs().snapshot(Output{"ITS", "TRACKSMCTR", 0}, stagedITS->labels);
  }

  auto& mftTrackROFs = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{"MFT", "MFTTrackROF", 0},
                                                                             stagedMFT->trackROFs.begin(), stagedMFT->trackROFs.end());
  auto& mftTracks = pc.outputs().make<std::vector<o2::mft::TrackMFT>>(Output{"MFT", "TRACKS", 0});
  mftTracks.assign(stagedMFT->tracks.begin(), stagedMFT->tracks.end());
  auto& mftClusIdx = pc.outputs().make<std::vector<int>>(Output{"MFT", "TRACKCLSID", 0});
  mftClusIdx.assign(stagedMFT->clusterIndices.begin(), stagedMFT->clusterIndices.end());
  auto& mftSeedPatterns = pc.outputs().make<std::vector<uint16_t>>(Output{"MFT", "TRACKSEEDPAT", 0});
  mftSeedPatterns.assign(stagedMFT->seedPatterns.begin(), stagedMFT->seedPatterns.end());
  LOGP(info, "Combined ITS+MFT CA pushed {} MFT tracks in {} ROFs", mftTracks.size(), mftTrackROFs.size());
  if (mUseMC) {
    pc.outputs().snapshot(Output{"MFT", "TRACKSMCTR", 0}, stagedMFT->labels);
  }
}

void CombinedCATrackerDPL::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool geometryFilled = false;
  if (!geometryFilled) {
    geometryFilled = true;
    pc.inputs().get<o2::itsmft::TopologyDictionary*>("itscldict");
    pc.inputs().get<o2::itsmft::TopologyDictionary*>("mftcldict");
    pc.inputs().get<o2::its::GeometryTGeo*>("itsTGeo");
    pc.inputs().get<o2::mft::GeometryTGeo*>("mftTGeo");
    o2::its::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                                                o2::math_utils::TransformType::T2GRot,
                                                                                o2::math_utils::TransformType::T2G));
    o2::mft::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                                                o2::math_utils::TransformType::T2GRot,
                                                                                o2::math_utils::TransformType::T2G,
                                                                                o2::math_utils::TransformType::L2G));
  }
  if (!mContinuousReadoutChecked && o2::base::GRPGeomHelper::instance().getGRPECS() != nullptr) {
    mContinuousReadoutChecked = true;
    const bool continuous = o2::base::GRPGeomHelper::instance().getGRPECS()->isDetContinuousReadOut(o2::detectors::DetID::MFT);
    requireContinuousMFTReadoutOrFatal(continuous);
  }
}

void CombinedCATrackerDPL::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    LOG(info) << "Combined ITS+MFT CA input ITS cluster dictionary updated";
    mITSDict = static_cast<const o2::itsmft::TopologyDictionary*>(obj);
    return;
  }
  if (matcher == ConcreteDataMatcher("MFT", "CLUSDICT", 0)) {
    LOG(info) << "Combined ITS+MFT CA input MFT cluster dictionary updated";
    mMFTDict = static_cast<const o2::itsmft::TopologyDictionary*>(obj);
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "GEOMTGEO", 0)) {
    LOG(info) << "Combined ITS+MFT CA input ITS GeometryTGeo loaded from CCDB";
    o2::its::GeometryTGeo::adopt(static_cast<o2::its::GeometryTGeo*>(obj));
    o2::its::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                                                o2::math_utils::TransformType::T2GRot,
                                                                                o2::math_utils::TransformType::T2G));
    return;
  }
  if (matcher == ConcreteDataMatcher("MFT", "GEOMTGEO", 0)) {
    LOG(info) << "Combined ITS+MFT CA input MFT GeometryTGeo loaded from CCDB";
    o2::mft::GeometryTGeo::adopt(static_cast<o2::mft::GeometryTGeo*>(obj));
    o2::mft::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                                                o2::math_utils::TransformType::T2GRot,
                                                                                o2::math_utils::TransformType::T2G,
                                                                                o2::math_utils::TransformType::L2G));
    return;
  }
}

DataProcessorSpec getCombinedCATrackerSpec(bool useMC)
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("compClustersITS", "ITS", "COMPCLUSTERS", 0, Lifetime::Timeframe);
  inputs.emplace_back("patternsITS", "ITS", "PATTERNS", 0, Lifetime::Timeframe);
  inputs.emplace_back("ROframesITS", "ITS", "CLUSTERSROF", 0, Lifetime::Timeframe);
  inputs.emplace_back("itscldict", "ITS", "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec("ITS/Calib/ClusterDictionary"));
  inputs.emplace_back("compClustersMFT", "MFT", "COMPCLUSTERS", 0, Lifetime::Timeframe);
  inputs.emplace_back("patternsMFT", "MFT", "PATTERNS", 0, Lifetime::Timeframe);
  inputs.emplace_back("ROframesMFT", "MFT", "CLUSTERSROF", 0, Lifetime::Timeframe);
  inputs.emplace_back("mftcldict", "MFT", "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec("MFT/Calib/ClusterDictionary"));
  if (useMC) {
    inputs.emplace_back("labelsITS", "ITS", "CLUSTERSMCTR", 0, Lifetime::Timeframe);
    inputs.emplace_back("labelsMFT", "MFT", "CLUSTERSMCTR", 0, Lifetime::Timeframe);
  }
  // Deliberately no IRFramesITS input -- see ConfigPreflight.h's
  // requireNoMFTIRFrameConfigOrFatal()/requireContinuousMFTReadoutOrFatal().

  // GeomRequest::None: clusters entering this tracker are already aligned,
  // so this workflow never requests the full aligned global geometry
  // (GeomRequest::Aligned -> CCDB "GLO/Config/GeometryAligned"). It always
  // requests each detector's own nominal/ideal GeometryTGeo snapshot
  // directly instead -- the same object both single-detector opt-in
  // workflows already default to (their own useGeom=false path) -- since
  // GeometryClusterDecoder<DetId>::decode() only ever consults
  // getMatrixT2L()/getMatrixL2G() from it (IOUtils.cxx), never an
  // alignment-adjusted matrix.
  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false, true, false, true, true,
                                                              o2::base::GRPGeomRequest::None, inputs, true);
  ggRequest->addInput({"itsTGeo", "ITS", "GEOMTGEO", 0, Lifetime::Condition, ccdbParamSpec("ITS/Config/Geometry")}, inputs);
  ggRequest->addInput({"mftTGeo", "MFT", "GEOMTGEO", 0, Lifetime::Condition, ccdbParamSpec("MFT/Config/Geometry")}, inputs);

  std::vector<OutputSpec> outputs;
  outputs.emplace_back("ITS", "TRACKS", 0, Lifetime::Timeframe);
  outputs.emplace_back("ITS", "TRACKCLSID", 0, Lifetime::Timeframe);
  outputs.emplace_back("ITS", "ITSTrackROF", 0, Lifetime::Timeframe);
  outputs.emplace_back("MFT", "TRACKS", 0, Lifetime::Timeframe);
  outputs.emplace_back("MFT", "MFTTrackROF", 0, Lifetime::Timeframe);
  outputs.emplace_back("MFT", "TRACKCLSID", 0, Lifetime::Timeframe);
  outputs.emplace_back("MFT", "TRACKSEEDPAT", 0, Lifetime::Timeframe);
  if (useMC) {
    outputs.emplace_back("ITS", "TRACKSMCTR", 0, Lifetime::Timeframe);
    outputs.emplace_back("MFT", "TRACKSMCTR", 0, Lifetime::Timeframe);
  }

  return DataProcessorSpec{
    "itsmft-combined-ca-tracker",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<CombinedCATrackerDPL>(ggRequest, useMC)},
    Options{}};
}

} // namespace o2::itsmft::combined
