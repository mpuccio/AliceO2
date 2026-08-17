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
#include "ITSMFTTracking/GenericTrackOutputAdapter.h"
#include "ITSMFTTracking/detail/DetectorRefitSupport.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/SurfaceTiming.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
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
using o2::itsmft::tracking::loadTimeFrameSources;
using o2::itsmft::tracking::TrackingOutcome;

template <int NLayers, o2::itsmft::tracking::SurfaceKind Kind>
bool completePublication(o2::itsmft::tracking::DetectorPublicationAdapter<NLayers>& publication,
                         const o2::itsmft::tracking::TimeFrame& frame,
                         const o2::itsmft::tracking::TrackingResult& result)
{
  const auto& parameters = frame.getTrackingParameters();
  const auto& scratch = frame.getWorkspace();
  for (std::size_t iteration = 0; iteration < parameters.size(); ++iteration) {
    const auto& candidates = scratch.getTraversalWorkspace(iteration).acceptedTracks;
    if (iteration >= result.acceptedTrackCounts.size() || result.acceptedTrackCounts[iteration] != candidates.size()) {
      return false;
    }
    std::vector<o2::itsmft::tracking::TrackingCandidate> selected;
    for (std::size_t index = 0; index < result.acceptedTrackCounts[iteration]; ++index) {
      if (candidates[index].track.innerState.kind == Kind) {
        selected.push_back(candidates[index]);
      }
    }
    if (!publication.completeAccepted(selected, parameters[iteration], scratch, iteration + 1 == parameters.size())) {
      return false;
    }
  }
  return true;
}

// Derive one source-level timing configuration from per-layer parameters;
// reject non-positive lengths and non-uniform staggering.
template <o2::detectors::DetID::ID DetId, int NLayers>
o2::itsmft::tracking::ROFTimingConfig deriveRofTimingConfigOrFatal(const o2::itsmft::TrackingParameters& params,
                                                                   std::size_t surfaceOffset)
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
    layerTimings[iLayer].mROFAddTimeErr = static_cast<uint32_t>(params.AddTimeError[surfaceOffset + iLayer]);
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

// The workflow owns the combined application plan and source schedule.
o2::itsmft::tracking::SurfaceCatalogView combinedCatalogView()
{
  return {o2::itsmft::tracking::kITSMFTCombinedStaticSurfaceCatalog.data(),
          static_cast<uint32_t>(o2::itsmft::tracking::kITSMFTCombinedStaticSurfaceCatalog.size())};
}

std::vector<o2::itsmft::tracking::LayerId> orderedSurfaceRange(uint16_t first, uint16_t count)
{
  std::vector<o2::itsmft::tracking::LayerId> result;
  result.reserve(count);
  for (uint16_t i = 0; i < count; ++i) {
    result.push_back(o2::itsmft::tracking::LayerId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

o2::itsmft::tracking::SurfaceMask surfaceRangeMask(uint16_t first, uint16_t count)
{
  o2::itsmft::tracking::SurfaceMask result;
  for (uint16_t i = 0; i < count; ++i) {
    result.set(o2::itsmft::tracking::LayerId{static_cast<uint16_t>(first + i)});
  }
  return result;
}

o2::itsmft::tracking::SurfaceLayoutDefinition combinedLayoutDefinition(
  gsl::span<const o2::itsmft::tracking::LayerId> itsSurfaces,
  const o2::itsmft::TrackingParameters& itsParams,
  gsl::span<const o2::itsmft::tracking::LayerId> mftSurfaces,
  const o2::itsmft::TrackingParameters& mftParams)
{
  o2::itsmft::tracking::SurfaceLayoutDefinition definition;
  definition.orderedSurfaces.reserve(itsSurfaces.size() + mftSurfaces.size());
  definition.orderedSurfaces.insert(definition.orderedSurfaces.end(), itsSurfaces.begin(), itsSurfaces.end());
  definition.orderedSurfaces.insert(definition.orderedSurfaces.end(), mftSurfaces.begin(), mftSurfaces.end());

  definition.componentOffsets = {0, static_cast<uint16_t>(itsSurfaces.size())};

  // MaxHoles is a graph-wide scalar. The combined workflow's scalar baseline
  // is the ITS Sync configuration; detector-local masks are projected into
  // the same global surface rank space below.
  definition.maxHoles = itsParams.MaxHoles;
  definition.holeSurfaces = o2::itsmft::tracking::positionalSurfaceMask(
    itsParams.HoleLayerMask, itsSurfaces, static_cast<uint32_t>(itsSurfaces.size()));
  definition.holeSurfaces |= o2::itsmft::tracking::positionalSurfaceMask(
    mftParams.HoleLayerMask, mftSurfaces, static_cast<uint32_t>(mftSurfaces.size()));
  definition.seedingSurfaces = surfaceRangeMask(0, static_cast<uint16_t>(definition.orderedSurfaces.size()));
  return definition;
}

template <typename T>
void appendSurfaceValues(std::vector<T>& destination, const std::vector<T>& first, const std::vector<T>& second)
{
  destination.clear();
  destination.reserve(first.size() + second.size());
  destination.insert(destination.end(), first.begin(), first.end());
  destination.insert(destination.end(), second.begin(), second.end());
}

o2::itsmft::TrackingParameters combinedTrackingParameters(const o2::itsmft::TrackingParameters& itsParams,
                                                          const o2::itsmft::TrackingParameters& mftParams)
{
  auto result = itsParams;
  result.NLayers = o2::itsmft::tracking::ITSNLayers + o2::itsmft::tracking::MFTNLayers;

  // These are the surface-indexed fields consumed by the common tracker.
  // Preserve the actual detector values and concatenate them in global rank
  // order (ITS [0, 7), MFT [7, 17)).
  appendSurfaceValues(result.AddTimeError, itsParams.AddTimeError, mftParams.AddTimeError);
  appendSurfaceValues(result.LayerZ, itsParams.LayerZ, mftParams.LayerZ);
  appendSurfaceValues(result.LayerRadii, itsParams.LayerRadii, mftParams.LayerRadii);
  appendSurfaceValues(result.LayerxX0, itsParams.LayerxX0, mftParams.LayerxX0);
  appendSurfaceValues(result.LayerResolution, itsParams.LayerResolution, mftParams.LayerResolution);
  appendSurfaceValues(result.SystError2Row, itsParams.SystError2Row, mftParams.SystError2Row);
  appendSurfaceValues(result.SystError2Col, itsParams.SystError2Col, mftParams.SystError2Col);

  // ITS uses LayerZ as its effective column extent; MFT supplies an explicit
  // global-X extent. Materialize that fallback before concatenating so the
  // resulting 17-layer configuration is complete for the index-table setup.
  const auto& itsColExtent = itsParams.LayerColHalfExtent.empty() ? itsParams.LayerZ : itsParams.LayerColHalfExtent;
  const auto& mftColExtent = mftParams.LayerColHalfExtent.empty() ? mftParams.LayerZ : mftParams.LayerColHalfExtent;
  appendSurfaceValues(result.LayerColHalfExtent, itsColExtent, mftColExtent);

  // StartLayerMask is intentionally independent of either detector's local
  // defaults: all 17 combined surfaces are valid road starts.
  result.StartLayerMask = (1u << result.NLayers) - 1u;
  result.HoleLayerMask = o2::itsmft::tracking::LayerMask{
    itsParams.HoleLayerMask.value() |
    (mftParams.HoleLayerMask.value() << o2::itsmft::tracking::ITSNLayers)};
  return result;
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

  // Preflight guarantees Sync mode and one iteration per detector.
  auto itsParams = o2::itsmft::TrackingMode::getTrackingParameters(o2::detectors::DetID::ITS, o2::itsmft::TrackingMode::Sync);
  auto mftParams = o2::itsmft::TrackingMode::getTrackingParameters(o2::detectors::DetID::MFT, o2::itsmft::TrackingMode::Sync);
  const auto itsSurfaces = orderedSurfaceRange(0, o2::itsmft::tracking::ITSNLayers);
  const auto mftSurfaces = orderedSurfaceRange(o2::itsmft::tracking::ITSNLayers, o2::itsmft::tracking::MFTNLayers);
  const auto combinedParams = combinedTrackingParameters(itsParams[0], mftParams[0]);

  mTraits = std::make_unique<o2::itsmft::tracking::TrackerTraits>();
  mTracker = std::make_unique<o2::itsmft::tracking::Tracker>(
    &o2::itsmft::tracking::detail::refitSurfaceSeed);

  o2::itsmft::tracking::TrackerInitialization configuration;
  configuration.catalog = combinedCatalogView();
  configuration.memoryPool = std::make_shared<o2::itsmft::tracking::BoundedMemoryResource>(
    std::min(itsParams[0].MaxMemory, mftParams[0].MaxMemory));
  o2::itsmft::tracking::TrackerIterationConfiguration iteration;
  iteration.layout = combinedLayoutDefinition(itsSurfaces, itsParams[0], mftSurfaces, mftParams[0]);
  iteration.parameters = combinedParams;
  configuration.iterations.push_back(std::move(iteration));

  const auto initialization = mTracker->initialize(mFrame, configuration);
  if (!initialization.ok()) {
    throw std::runtime_error{"Combined ITS+MFT tracker: failed to commit static TimeFrame configuration"};
  }
  mITSPublicationAdapter.adoptITSSharedClusterCompatibility(&mITSCompatibility);
  mMFTPublicationAdapter.adoptMFTPublicationCompatibility(&mMFTCompatibility);
  mFrame.setBz(o2::base::Propagator::Instance()->getNominalBz());

  const int itsNThreads = o2::itsmft::ITSCommonCATrackerParam::Instance().nThreads;
  const int mftNThreads = o2::itsmft::tracking::TrackerParamRef<o2::detectors::DetID::MFT>::get().nThreads;
  const int nThreads = std::max({1, itsNThreads, mftNThreads});
  std::shared_ptr<tbb::task_arena> arena;
  mTraits->setNThreads(nThreads, arena);
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

o2::itsmft::tracking::SurfaceCatalogView CombinedCATrackerDPL::catalogView() const noexcept
{
  return combinedCatalogView();
}

std::optional<bool> CombinedCATrackerDPL::dropTFUponFailureFor(ClusterSourceId source) const noexcept
{
  if (source == ClusterSourceId{0} || source == ClusterSourceId{1}) {
    const auto& parameters = mFrame.getTrackingParameters();
    return parameters.empty() ? std::optional<bool>{} : std::optional<bool>{parameters[0].DropTFUponFailure};
  }
  return std::nullopt;
}

void CombinedCATrackerDPL::configureRofTables(const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource)
{
  auto configure = [](auto& overlap, auto& vertex, auto& mask, const auto& timing, uint32_t nROFs, int layers) {
    using LayerTiming = o2::its::LayerTiming;
    LayerTiming layerTiming{};
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
  configure(mITSROFOverlapTable, mITSROFVertexLookupTable, mITSMultiplicityMask, itsSource.timing, static_cast<uint32_t>(itsSource.rofs.size()), o2::itsmft::tracking::ITSNLayers);
  configure(mMFTROFOverlapTable, mMFTROFVertexLookupTable, mMFTMultiplicityMask, mftSource.timing, static_cast<uint32_t>(mftSource.rofs.size()), o2::itsmft::tracking::MFTNLayers);
}

void CombinedCATrackerDPL::clearPublicationSidecars() noexcept
{
  mITSPublicationAdapter.reset();
  mMFTPublicationAdapter.reset();
}

void CombinedCATrackerDPL::clearRofViews() noexcept
{
  mFrame.getWorkspace().setROFViews({});
}

void CombinedCATrackerDPL::invalidatePublication() noexcept
{
  mITSClock.reset();
  mMFTClock.reset();
  mPublicationValid = false;
  clearRofViews();
}

void CombinedCATrackerDPL::markPublicationValid() noexcept
{
  mITSClock.emplace(mITSROFOverlapTable.getView().getClockLayer());
  mMFTClock.emplace(mMFTROFOverlapTable.getView().getClockLayer());
  mPublicationValid = true;
}

std::optional<o2::itsmft::tracking::GenericTrackPublicationExport> CombinedCATrackerDPL::getITSPublicationExport() const
{
  if (!mPublicationValid || !mITSClock) {
    return std::nullopt;
  }
  return o2::itsmft::tracking::GenericTrackPublicationExport{o2::detectors::DetID::ITS, ClusterSourceId{0}, *mITSClock,
                                                             getITSOrderedSurfaces()};
}

std::optional<o2::itsmft::tracking::GenericTrackPublicationExport> CombinedCATrackerDPL::getMFTPublicationExport() const
{
  if (!mPublicationValid || !mMFTClock) {
    return std::nullopt;
  }
  return o2::itsmft::tracking::GenericTrackPublicationExport{o2::detectors::DetID::MFT, ClusterSourceId{1}, *mMFTClock,
                                                             getMFTOrderedSurfaces()};
}

TrackingOutcome CombinedCATrackerDPL::trackFrame(const ClusterSourceInput& itsSource, const ClusterSourceInput& mftSource,
                                                 const o2::InteractionRecord& origin)
{
  invalidatePublication();
  // Sidecars are sealed on success and are not cleared by the next frame's
  // normalized-frame commit; clear them before each new transaction.
  clearPublicationSidecars();

  configureRofTables(itsSource, mftSource);
  auto itsInput = itsSource;
  auto mftInput = mftSource;
  itsInput.rofViews = {mITSROFOverlapTable.getView(), mITSROFVertexLookupTable.getView(),
                       mITSMultiplicityMask.getView(), mITSUPCMask.getView()};
  mftInput.rofViews = {mMFTROFOverlapTable.getView(), mMFTROFVertexLookupTable.getView(),
                       mMFTMultiplicityMask.getView(), mMFTUPCMask.getView()};

  // The fixed ITS=0/MFT=1 source contract lives in this workflow-owned
  // application composition, never in the generic loader.
  LoadSourcesResult loadResult;
  if (const auto rejected = validateSources(itsSource, mftSource)) {
    loadResult = *rejected;
  } else {
    const std::array<ClusterSourceInput, 2> sources{itsInput, mftInput};
    loadResult = loadTimeFrameSources(mFrame, gsl::span<const ClusterSourceInput>{sources}, catalogView(), origin);
  }
  if (!loadResult.ok()) {
    // Apply the shared recoverable-load taxonomy and the owning detector's
    // DropTFUponFailure policy. No tracker runs after an uncommitted load.
    const bool errorIsRecoverable = o2::itsmft::tracking::isRecoverableLoadError(loadResult.error, loadResult.timingDetail);
    const auto dropAllowed = dropTFUponFailureFor(loadResult.source);
    const bool sourceRecognized = dropAllowed.has_value();
    if (!sourceRecognized) {
      LOGP(error, "Combined TF load failure reports an unrecognized source id {}; treating as structural", loadResult.source.value());
    }
    const auto outcome = errorIsRecoverable && sourceRecognized && dropAllowed.value_or(false)
                           ? TrackingOutcome::RecoverableDropped
                           : TrackingOutcome::Structural;
    LOGP(error, "Combined TF loading failed (source={}, error={}, recoverable={}, dropAllowed={}): outcome={}",
         loadResult.source.value(), static_cast<int>(loadResult.error), errorIsRecoverable, dropAllowed.value_or(false),
         outcome == TrackingOutcome::RecoverableDropped ? "RecoverableDropped" : "Structural");
    mITSPublicationAdapter.reset();
    mMFTPublicationAdapter.reset();
    mFrame.resetTimeFrame();
    invalidatePublication();
    return outcome;
  }

  // Loading commits the normalized multi-source event and the one global
  // workspace before the tracker executes the disconnected components.
  try {
    const auto result = mTracker->run(mFrame, *mTraits);
    if (result.outcome != TrackingOutcome::Success) {
      mITSPublicationAdapter.reset();
      mMFTPublicationAdapter.reset();
      invalidatePublication();
      return result.outcome;
    }
    if (!completePublication<o2::itsmft::tracking::ITSNLayers, o2::itsmft::tracking::SurfaceKind::Cylinder>(
          mITSPublicationAdapter, mFrame, result)) {
      mITSPublicationAdapter.reset();
      throw std::runtime_error{"failed to seal ITS tracking compatibility"};
    }
    if (!completePublication<o2::itsmft::tracking::MFTNLayers, o2::itsmft::tracking::SurfaceKind::Disk>(
          mMFTPublicationAdapter, mFrame, result)) {
      mMFTPublicationAdapter.reset();
      throw std::runtime_error{"failed to seal MFT tracking compatibility"};
    }
  } catch (...) {
    mITSPublicationAdapter.reset();
    mMFTPublicationAdapter.reset();
    invalidatePublication();
    throw;
  }

  if (!mFrame.isConfigured()) {
    invalidatePublication();
    return TrackingOutcome::Structural;
  }

  markPublicationValid();
  return TrackingOutcome::Success;
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

  // Build fresh workflow-owned source inputs for this TF.
  ClusterSourceInput itsSource{};
  itsSource.id = ClusterSourceId{0};
  itsSource.detector = o2::detectors::DetID::ITS;
  itsSource.clusters = gsl::span<const o2::itsmft::CompClusterExt>(itsCompClusters.data(), itsCompClusters.size());
  itsSource.patterns = itsPatterns;
  itsSource.rofs = itsRofs;
  itsSource.dictionary = mITSDict;
  itsSource.labels = itsLabels;
  // Use the plan's ordered surfaces rather than re-deriving detector offsets.
  itsSource.layerToSurface = getITSOrderedSurfaces();
  itsSource.timing = deriveRofTimingConfigOrFatal<o2::detectors::DetID::ITS, o2::itsmft::tracking::ITSNLayers>(
    *mFrame.getTrackingParameters(0), 0);
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
  mftSource.timing = deriveRofTimingConfigOrFatal<o2::detectors::DetID::MFT, o2::itsmft::tracking::MFTNLayers>(
    *mFrame.getTrackingParameters(0), o2::itsmft::tracking::ITSNLayers);
  mftSource.decoder = mMFTDecoder.get();

  const auto origin = chooseOrigin(itsRofs, mftRofs);
  const auto resetCount = mFrame.getEventResetCount();
  try {
    const auto outcome = trackFrame(itsSource, mftSource, origin);

    if (outcome == TrackingOutcome::RecoverableDropped) {
      LOGP(error,
           "Combined ITS+MFT CA tracking recoverably dropped this TF ({} ITS ROFs/{} ITS clusters, {} MFT "
           "ROFs/{} MFT clusters); publishing nothing and continuing with the next TimeFrame",
           itsRofs.size(), itsCompClusters.size(), mftRofs.size(), mftCompClusters.size());
      return;
    }
    if (outcome == TrackingOutcome::Structural) {
      throw std::runtime_error{"Combined ITS+MFT CA tracking hit a structural failure"};
    }

    const auto itsExport = getITSPublicationExport();
    const auto mftExport = getMFTPublicationExport();
    if (!itsExport || !mftExport) {
      throw std::runtime_error{"Combined ITS+MFT GenericTrack output publication context is unavailable after a successful trackFrame()"};
    }
    const o2::itsmft::tracking::GenericTrackPublicationContext itsContext{
      itsExport->detector, itsExport->source, itsRofs, itsExport->clock, itsExport->orderedSurfaces};
    const o2::itsmft::tracking::GenericTrackPublicationContext mftContext{
      mftExport->detector, mftExport->source, mftRofs, mftExport->clock, mftExport->orderedSurfaces};

    // Stage both outputs before requesting either output stream.
    o2::itsmft::tracking::GenericTrackOutputAdapterError itsError = o2::itsmft::tracking::GenericTrackOutputAdapterError::None;
    const auto stagedITS = o2::itsmft::tracking::stageITSGenericTrackOutput(
      mFrame, itsContext, getITSSharedClusterCompatibility(), mUseMC, itsError);
    o2::itsmft::tracking::GenericTrackOutputAdapterError mftError = o2::itsmft::tracking::GenericTrackOutputAdapterError::None;
    const auto stagedMFT = o2::itsmft::tracking::stageMFTGenericTrackOutput(
      mFrame, mftContext, getMFTPublicationCompatibility(), mUseMC, mftError);
    if (!stagedITS || !stagedMFT) {
      throw std::runtime_error{"Combined ITS+MFT GenericTrack output staging failed"};
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
    mITSPublicationAdapter.reset();
    mMFTPublicationAdapter.reset();
    mFrame.resetTimeFrame();
    invalidatePublication();
  } catch (...) {
    if (mFrame.getEventResetCount() == resetCount) {
      mITSPublicationAdapter.reset();
      mMFTPublicationAdapter.reset();
      mFrame.resetTimeFrame();
    }
    invalidatePublication();
    throw;
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

  // Clusters are already aligned; request each detector's nominal GeometryTGeo
  // snapshot used by GeometryClusterDecoder.
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
