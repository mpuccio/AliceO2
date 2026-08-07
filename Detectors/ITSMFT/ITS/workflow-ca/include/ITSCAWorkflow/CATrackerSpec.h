// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
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
/// \brief Opt-in ITS common-CA tracker DPL device (Gate 3 workflow-onboarding
///        Slice 2). Modeled on MFTWorkflow/CATrackerSpec.h where contracts
///        match; publishes tracker-only outputs (no VERTICES/VERTICESROF).

#ifndef O2_ITS_CA_WORKFLOW_CATRACKERSPEC_H_
#define O2_ITS_CA_WORKFLOW_CATRACKERSPEC_H_

#include <cassert>
#include <memory>
#include <optional>
#include <vector>

#include <gsl/span>

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ClockTimingPublicationView.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/DetectorPublicationAdapter.h"
#include "ITSMFTTracking/TrackingOperationAdapter.h"
#include "ITSMFTTracking/Tracker.h"
#include "ITSMFTTracking/TrackerTraits.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/TimeFrameLoadFailure.h"
#include "ITSMFTTracking/StaticDetectorCatalogs.h"
#include "ITSMFTTracking/ROFViews.h"
#include "ITStracking/ROFLookupTables.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::its::ca
{

/// Gate 3 failure-contract publication decision, mirroring
/// o2::mft::decideCATrackerPublicationAction() exactly (same three-way
/// outcome, same meaning) -- factored out as a pure function so the
/// publish-vs-skip contract can be exercised without a DPL ProcessingContext.
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
/// Byte-for-byte the same per-track logic as the legacy
/// ITSWorkflow/TrackerSpec.cxx::run() loop: cluster indices are read
/// per-layer (`track.getClusterIndex(layer)`, -1 if that layer has no hit)
/// in decreasing layer order and each valid one is pushed onto
/// `clusterIndices`; the resulting flattened order is
/// outer-to-inner-most-populated-first, not physical layer order -- this
/// matches every existing ITS track consumer's expectation (they all read
/// through RangeRefComp + TRACKCLSID, never by assuming physical layer
/// order).
///
/// Unlike the legacy loop, this does NOT re-derive the packed cluster size
/// here: by the time a track reaches `tf.getTracks()`,
/// `track.getClusterIndex(layer)` has already been rewritten by
/// Tracker<NLayers>::rectifyClusterIndices() (CATracker.cxx) from this
/// layer's own local cluster identity -- the domain mClusterSize is keyed
/// by -- to the external/global one, so the local identity needed to
/// address mClusterSize is gone by the time this function runs.
/// rectifyClusterIndices() captures `tf.getClusterSize(layer, localIndex)`
/// onto the track itself (via TrackITS::setClusterSize()) while the local
/// index is still available, so this function only has to carry that
/// already-correct per-layer size through the by-value copy into the output
/// TrackITS -- it must not query the TimeFrame with the (by-then external)
/// `clid` here, which is exactly the local-vs-external index confusion this
/// function used to have.
///
/// Takes `track` by value: the RangeRefComp mutation happens on the
/// caller's own copy of the source TrackITSExt before it is sliced into the
/// output TrackITS, so the caller's TimeFrame-owned track is never mutated
/// by this call.
///
/// `clusterIndices`/`tracks` are templated (not plain std::vector<>&): DPL's
/// pc.outputs().make<std::vector<T>>() actually returns a pmr-allocator
/// vector (std::vector<T, polymorphic_allocator<T>>), not std::vector<T>, so
/// a fixed std::vector<T>& parameter would reject the real call site in
/// CATrackerSpec.cxx while only ever being exercised with plain
/// std::vector<T> in tests -- genericity here is required for the real
/// caller, not speculative.
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

/// ITS common-CA tracker DPL task. Owns the standalone TimeFrame and
/// composes the workflow-owned input/timing/publication edge with Tracker.
/// Frozen legacy o2::its::TrackerDPL
/// (ITSWorkflow/TrackerSpec.h) and o2-its-reco-workflow are untouched by
/// this class; it lives in an isolated library/executable
/// (o2-its-ca-tracker-workflow) with no link-graph overlap with ITSWorkflow.
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
  void resetEvent() noexcept;
  bool isActive() const noexcept { return mFrame.isConfigured() && !mFrame.getTrackingParameters().empty(); }
  const o2::itsmft::tracking::SurfaceTrackingScratch& getScratch() const noexcept
  {
    return mFrame.getWorkspace(o2::itsmft::tracking::ClusterSourceId{0});
  }
  o2::itsmft::tracking::SurfaceTrackingScratch& getScratch() noexcept
  {
    return mFrame.getWorkspace(o2::itsmft::tracking::ClusterSourceId{0});
  }

  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  bool mUseMC = false;
  bool mTrackingInitialised = false;
  o2::itsmft::TrackingMode::Type mTrackingMode = o2::itsmft::TrackingMode::Unset;
  o2::itsmft::tracking::TimeFrame mFrame;
  std::unique_ptr<o2::itsmft::tracking::TrackerTraits> mTrackerTraits;
  std::unique_ptr<o2::itsmft::tracking::Tracker> mTracker;
  std::unique_ptr<o2::itsmft::tracking::TrackingOperationAdapter> mOperationAdapter;
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
