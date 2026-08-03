// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file testITSMFTLegacyParticipantSet.cxx
/// \brief M2c (GenericTrackingEngineMigration.md; ADR 0007) focused tests for
/// the one explicit ITS/MFT-named application-layer participant set:
///  - global-plan authority: both detector legs' DetectorLayoutViews derive
///    from the one authoritative combined ITS+MFT topology this set builds
///    exactly once in its constructor;
///  - source-qualified binding order: validateSources()/loadBindings()
///    enforce and preserve the fixed ITS=0/MFT=1 contract;
///  - explicit schedule order: schedule() is [ITS, MFT], never re-derived
///    from participant construction order;
///  - publication-export lifetimes: getITSPublicationExport()/
///    getMFTPublicationExport() are engaged only between
///    markPublicationValid() and the next invalidatePublication().
///
/// These exercise ITSMFTLegacyParticipantSet directly, not through a
/// composed load/track/publish pipeline -- testCombinedTrackingComposition
/// .cxx keeps the full load/track/publish composition-level coverage
/// (ported from the pre-M3 combined-coordinator test suite this milestone
/// deleted).

#define BOOST_TEST_MODULE ITSMFT ITSMFTLegacyParticipantSet
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/ITSMFTLegacyParticipantSet.h"
#include "ITSMFTTracking/MultiSourceTimeFrameLoader.h"
#include "ITSMFTTracking/TimeFrame.h"

using namespace o2::itsmft;
using namespace o2::itsmft::tracking;

namespace
{

TrackingParameters makeItsParams()
{
  TrackingParameters p;
  resetDetectorDefaults(p, o2::detectors::DetID::ITS);
  return p;
}

TrackingParameters makeMftParams()
{
  TrackingParameters p;
  resetDetectorDefaults(p, o2::detectors::DetID::MFT);
  return p;
}

ITSMFTLegacyParticipantSet makeSet()
{
  return ITSMFTLegacyParticipantSet{std::vector<TrackingParameters>{makeItsParams()},
                                    std::vector<TrackingParameters>{makeMftParams()}};
}

// A minimal, structurally valid ClusterSourceInput: only the fields
// validateSources()/loadBindings()/configureRofTables() actually read
// (id, detector, timing, rofs) -- no decoder/dictionary/clusters, since
// none of the tests below drive an actual atomic load or track() call.
ClusterSourceInput minimalSource(ClusterSourceId id, o2::detectors::DetID::ID det)
{
  ClusterSourceInput source{};
  source.id = id;
  source.detector = det;
  source.timing = ROFTimingConfig{40, 0, 0, 0};
  source.rofs = std::vector<ROFRecord>{ROFRecord{{100, 5}, 0, 0, 0}};
  return source;
}

} // namespace

BOOST_AUTO_TEST_CASE(GlobalPlanAuthorityBothDetectorLayoutsAreByteIdentical)
{
  // Single authoritative combined topology (M2c carries forward the same
  // proof the pre-M3 combined-coordinator test suite already established at
  // the composition level): the set builds its one shared ITS+MFT
  // DetectorLayout exactly once in its constructor, and both
  // DetectorLayoutSets/DetectorTraversalBindings only ever receive a
  // passive copy of that one built object.
  auto participants = makeSet();
  const auto itsView = participants.getITSLayoutView();
  const auto mftView = participants.getMFTLayoutView();

  BOOST_REQUIRE_EQUAL(itsView.nSurfaces, mftView.nSurfaces);
  BOOST_REQUIRE_EQUAL(itsView.nSurfaces, static_cast<uint32_t>(ITSNLayers + MFTNLayers));
  BOOST_CHECK(itsView.surfaces == mftView.surfaces); // pointer identity: same static catalog storage
  BOOST_CHECK(itsView.cylinderSurfaces == mftView.cylinderSurfaces);
  BOOST_CHECK(itsView.diskSurfaces == mftView.diskSurfaces);

  const auto& itsTopo = itsView.topology;
  const auto& mftTopo = mftView.topology;
  BOOST_REQUIRE_EQUAL(itsTopo.nTransitions, mftTopo.nTransitions);
  BOOST_REQUIRE_EQUAL(itsTopo.nCells, mftTopo.nCells);
  for (uint32_t t = 0; t < itsTopo.nTransitions; ++t) {
    const auto& a = itsTopo.getTransition(TransitionId{static_cast<uint16_t>(t)});
    const auto& b = mftTopo.getTransition(TransitionId{static_cast<uint16_t>(t)});
    BOOST_CHECK(a.from == b.from);
    BOOST_CHECK(a.to == b.to);
    BOOST_CHECK(a.policyTag == b.policyTag);
    BOOST_CHECK(a.skippedSurfaces == b.skippedSurfaces);
  }

  // Ordered-surface spans agree with the fixed ITS 0..6 / MFT 7..16 global
  // offsets, valid immediately after construction.
  BOOST_REQUIRE_EQUAL(participants.getITSOrderedSurfaces().size(), static_cast<size_t>(ITSNLayers));
  BOOST_REQUIRE_EQUAL(participants.getMFTOrderedSurfaces().size(), static_cast<size_t>(MFTNLayers));
  BOOST_CHECK(participants.getITSOrderedSurfaces()[0] == SurfaceId{0});
  BOOST_CHECK(participants.getMFTOrderedSurfaces()[0] == SurfaceId{ITSNLayers});
}

BOOST_AUTO_TEST_CASE(ExplicitScheduleOrderIsITSThenMFTNeverReDerived)
{
  // ADR 0007 decision 6: the schedule is explicit plan data, not an
  // inherent consequence of participant declaration/construction order.
  auto participants = makeSet();
  const auto schedule = participants.schedule();
  BOOST_REQUIRE_EQUAL(schedule.size(), 2u);
  BOOST_REQUIRE(schedule[0] != nullptr);
  BOOST_REQUIRE(schedule[1] != nullptr);
  BOOST_CHECK(schedule[0]->id() == ParticipantId{0});
  BOOST_CHECK(schedule[1]->id() == ParticipantId{1});
  BOOST_CHECK_EQUAL(schedule[0]->ownedSurfaces().size(), static_cast<size_t>(ITSNLayers));
  BOOST_CHECK_EQUAL(schedule[1]->ownedSurfaces().size(), static_cast<size_t>(MFTNLayers));
  BOOST_CHECK(schedule[0]->ownedSurfaces()[0] == SurfaceId{0});
  BOOST_CHECK(schedule[1]->ownedSurfaces()[0] == SurfaceId{ITSNLayers});
}

BOOST_AUTO_TEST_CASE(SourceQualifiedBindingOrderAcceptsFixedContractAndOrdersITSThenMFT)
{
  auto participants = makeSet();
  const auto itsSource = minimalSource(ClusterSourceId{0}, o2::detectors::DetID::ITS);
  const auto mftSource = minimalSource(ClusterSourceId{1}, o2::detectors::DetID::MFT);

  BOOST_CHECK(!participants.validateSources(itsSource, mftSource).has_value());

  const auto bindings = participants.loadBindings(itsSource, mftSource);
  BOOST_REQUIRE_EQUAL(bindings.size(), 2u);
  BOOST_CHECK(bindings[0].source.id == ClusterSourceId{0});
  BOOST_CHECK(bindings[0].source.detector == o2::detectors::DetID::ITS);
  BOOST_CHECK(bindings[1].source.id == ClusterSourceId{1});
  BOOST_CHECK(bindings[1].source.detector == o2::detectors::DetID::MFT);

  // This set's combined catalog view is the fixed 17-surface ITS+MFT static
  // catalog every binding above is scoped against.
  BOOST_CHECK_EQUAL(participants.catalogView().nSurfaces, static_cast<uint32_t>(ITSNLayers + MFTNLayers));
}

BOOST_AUTO_TEST_CASE(SourceQualifiedBindingOrderRejectsSwappedOrMismatchedIds)
{
  auto participants = makeSet();

  // Swapped ids/detectors: neither slot matches its own fixed position.
  const auto swappedIts = minimalSource(ClusterSourceId{1}, o2::detectors::DetID::MFT);
  const auto swappedMft = minimalSource(ClusterSourceId{0}, o2::detectors::DetID::ITS);
  const auto swappedRejection = participants.validateSources(swappedIts, swappedMft);
  BOOST_REQUIRE(swappedRejection.has_value());
  BOOST_CHECK(swappedRejection->error == MultiSourceLoadError::UnsupportedDetector);
  BOOST_CHECK(swappedRejection->source == swappedIts.id);

  // A correct ITS slot but an unrecognized MFT id/detector.
  const auto itsSource = minimalSource(ClusterSourceId{0}, o2::detectors::DetID::ITS);
  const auto badMftSource = minimalSource(ClusterSourceId{5}, o2::detectors::DetID::MFT);
  const auto mftRejection = participants.validateSources(itsSource, badMftSource);
  BOOST_REQUIRE(mftRejection.has_value());
  BOOST_CHECK(mftRejection->error == MultiSourceLoadError::UnsupportedDetector);
  BOOST_CHECK(mftRejection->source == badMftSource.id);

  // An unrecognized load-failure source is reported as such by
  // dropTFUponFailureFor() too -- never silently attributed to either
  // detector's own DropTFUponFailure.
  BOOST_CHECK(!participants.dropTFUponFailureFor(ClusterSourceId{5}).has_value());
  BOOST_CHECK(participants.dropTFUponFailureFor(ClusterSourceId{0}).has_value());
  BOOST_CHECK(participants.dropTFUponFailureFor(ClusterSourceId{1}).has_value());
}

BOOST_AUTO_TEST_CASE(PublicationExportLifetimesEngagedOnlyBetweenMarkAndInvalidate)
{
  auto participants = makeSet();
  TimeFrame frame;
  participants.adoptFrame(frame);

  // Never valid before any markPublicationValid() call.
  BOOST_CHECK(!participants.getITSPublicationExport().has_value());
  BOOST_CHECK(!participants.getMFTPublicationExport().has_value());

  // configureRofTables() is the only precondition markPublicationValid()
  // needs (it reads each participant's own ROF overlap table clock layer) --
  // no atomic load or track() call is required for this lifetime proof.
  const auto itsSource = minimalSource(ClusterSourceId{0}, o2::detectors::DetID::ITS);
  const auto mftSource = minimalSource(ClusterSourceId{1}, o2::detectors::DetID::MFT);
  participants.configureRofTables(itsSource, mftSource);

  participants.markPublicationValid();
  const auto itsExport = participants.getITSPublicationExport();
  const auto mftExport = participants.getMFTPublicationExport();
  BOOST_REQUIRE(itsExport.has_value());
  BOOST_REQUIRE(mftExport.has_value());
  BOOST_CHECK(itsExport->detector == o2::detectors::DetID::ITS);
  BOOST_CHECK(itsExport->source == ClusterSourceId{0});
  BOOST_CHECK_EQUAL(itsExport->orderedSurfaces.size(), static_cast<size_t>(ITSNLayers));
  BOOST_CHECK(mftExport->detector == o2::detectors::DetID::MFT);
  BOOST_CHECK(mftExport->source == ClusterSourceId{1});
  BOOST_CHECK_EQUAL(mftExport->orderedSurfaces.size(), static_cast<size_t>(MFTNLayers));

  // invalidatePublication() disengages both, without needing another
  // eventReset()/wipe() -- the same narrow scope the combined DPL task's
  // own trackFrame() composition relies on (CombinedCATrackerSpec.cxx).
  participants.invalidatePublication();
  BOOST_CHECK(!participants.getITSPublicationExport().has_value());
  BOOST_CHECK(!participants.getMFTPublicationExport().has_value());

  // clearPublicationSidecars() is a distinct, narrower operation (detector
  // compatibility sidecars only) -- exercised here just to prove it never
  // engages or disengages the publication/timing bridge itself.
  participants.markPublicationValid();
  BOOST_REQUIRE(participants.getITSPublicationExport().has_value());
  participants.clearPublicationSidecars();
  BOOST_CHECK(participants.getITSPublicationExport().has_value());
}

BOOST_AUTO_TEST_CASE(ConstructorRejectsMultiIterationParameters)
{
  // The only shape this set's single shared combined layout can represent
  // is exactly one TrackingParameters iteration per detector -- ported from
  // the pre-M3 combined-coordinator test suite this milestone deleted,
  // since this constructor validation now lives only here.
  std::vector<TrackingParameters> twoIterations(2);
  BOOST_CHECK_THROW((ITSMFTLegacyParticipantSet{twoIterations, {makeMftParams()}}), std::invalid_argument);
  BOOST_CHECK_THROW((ITSMFTLegacyParticipantSet{{makeItsParams()}, twoIterations}), std::invalid_argument);
}
