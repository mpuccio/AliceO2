// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file CATrackerSpec.h
/// \brief ITS common-CA tracker DPL device with tracker-only outputs.

#ifndef O2_ITS_CA_WORKFLOW_CATRACKERSPEC_H_
#define O2_ITS_CA_WORKFLOW_CATRACKERSPEC_H_

#include <cassert>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

#include <gsl/span>

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITSMFTTracking/GenericTrackOutputAdapter.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ClusterDecoding.h"
#include "ITSMFTTracking/detail/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/IOUtils.h"
#include "ITSMFTTracking/ITSMFTDetectorDefinitions.h"
#include "ITSMFTTracking/ROFViews.h"
#include "ITSMFTTracking/ROFLookupTables.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::its::ca
{

/// Failure-contract publication decision, factored out as a pure function so
/// the publish-vs-skip contract is testable without a DPL ProcessingContext.
enum class CATrackerPublicationAction {
  PublishInactiveEmpty, ///< tracker not configured/active: publish empty outputs
  PublishActiveResult,  ///< active tracking returned a non-dropped result (including valid-empty input)
  SkipDroppedTimeFrame, ///< active tracking recoverably dropped this TF: publish nothing
};

CATrackerPublicationAction decideCATrackerPublicationAction(bool trackerActive, o2::itsmft::tracking::TrackingOutcome outcome) noexcept;

/// Converts one o2::its::TrackITSExt (the common tracker's internal barrel
/// track representation, see ITSMFTTracking/Cell.h's CATrackTypeHelper<7>)
/// into the persisted o2::its::TrackITS, appending this track's flattened
/// external cluster indices onto `clusterIndices` and setting the
/// RangeRefComp cluster range (o2::its::TrackITS::setFirstClusterEntry(),
/// via the CA algorithm's own pre-set cluster count) accordingly.
///
/// Copies cluster indices in decreasing layer order into `clusterIndices`,
/// preserving the persisted outer-to-inner order.
///
/// Cluster sizes have already been captured on the track before local indices
/// are rewritten to external indices; do not query the TimeFrame here.
///
/// Takes `track` by value so RangeRefComp mutation cannot alter TimeFrame-owned
/// state before the output copy is appended.
///
/// The vectors are templated to accept DPL's polymorphic-allocator outputs.
template <typename ClusterIdxVec, typename TracksVec>
void convertTrackITSExtToTrackITS(o2::its::TrackITSExt track,
                                  ClusterIdxVec& clusterIndices,
                                  TracksVec& tracks)
{
  track.setFirstClusterEntry(static_cast<int>(clusterIndices.size()));
  const int nclExpected = track.getNumberOfClusters();
  int nclFound = 0;
  for (int layer = o2::its::TrackITSExt::MaxClusters; layer--;) {
    const int clid = track.getClusterIndex(layer);
    if (clid < 0) {
      continue;
    }
    clusterIndices.push_back(clid);
    ++nclFound;
  }
  assert(nclFound == nclExpected);
  (void)nclExpected;
  (void)nclFound;
  tracks.push_back(track);
}

/// ITS common-CA tracker DPL task. Owns the TimeFrame and composes the
/// workflow input/timing/publication edge with Tracker.
class CATrackerDPL : public o2::framework::Task
{
 public:
  CATrackerDPL(std::shared_ptr<o2::base::GRPGeomRequest> gr, bool useMC, o2::itsmft::TrackingMode::Type trMode);
  ~CATrackerDPL() override = default;

  void init(framework::InitContext& ic) final;
  void run(framework::ProcessingContext& pc) final;
  void finaliseCCDB(framework::ConcreteDataMatcher& matcher, void* obj) final;

 private:
  void updateTimeDependentParams(framework::ProcessingContext& pc);
  void configureROFViews(gsl::span<const o2::itsmft::ROFRecord> rofs);
  void invalidatePublication() noexcept;
  void initialiseTracking();
  o2::itsmft::tracking::TrackingOutcome processTimeFrame(
    gsl::span<const o2::itsmft::ROFRecord> rofs,
    gsl::span<const o2::itsmft::CompClusterExt> clusters,
    gsl::span<const unsigned char> patterns,
    const o2::dataformats::MCTruthContainer<MCCompLabel>* labels);
  void resetTimeFrame() noexcept;
  bool isActive() const noexcept { return mTracker != nullptr && mTracker->isConfiguredFor(mFrame); }
  const o2::itsmft::tracking::TimeFrameScratch& getScratch() const noexcept
  {
    return mFrame.getScratch();
  }
  o2::itsmft::tracking::TimeFrameScratch& getScratch() noexcept
  {
    return mFrame.getScratch();
  }

  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  bool mUseMC = false;
  bool mTrackingInitialised = false;
  o2::itsmft::TrackingMode::Type mTrackingMode = o2::itsmft::TrackingMode::Unset;
  o2::itsmft::tracking::TimeFrame mFrame;
  std::vector<std::vector<uint32_t>> mExternalIndicesBySurface;
  std::vector<std::vector<uint32_t>> mClusterSizesBySurface;
  std::unique_ptr<o2::itsmft::tracking::TrackerTraits> mTrackerTraits;
  std::unique_ptr<o2::itsmft::tracking::Tracker> mTracker;
  std::unique_ptr<o2::itsmft::tracking::ClusterDecoder> mClusterDecoder;
  const o2::itsmft::TopologyDictionary* mDictionary = nullptr;
  o2::itsmft::tracking::DetectorPublicationAdapter<o2::itsmft::tracking::ITSNLayers> mPublicationAdapter;
  o2::itsmft::tracking::ITSSharedClusterCompatibility mCompatibility;
  o2::its::ROFOverlapTable<o2::itsmft::tracking::ITSNLayers> mROFOverlapTable;
  o2::its::ROFVertexLookupTable<o2::itsmft::tracking::ITSNLayers> mROFVertexLookupTable;
  o2::its::ROFMaskTable<o2::itsmft::tracking::ITSNLayers> mMultiplicityMask;
  o2::its::ROFMaskTable<o2::itsmft::tracking::ITSNLayers> mUPCMask;
  std::optional<o2::itsmft::tracking::ClockTimingPublicationView> mPublicationClock;
};

o2::framework::DataProcessorSpec getCATrackerSpec(bool useMC, bool useGeom, o2::itsmft::TrackingMode::Type trMode);

} // namespace o2::its::ca

#endif // O2_ITS_CA_WORKFLOW_CATRACKERSPEC_H_
