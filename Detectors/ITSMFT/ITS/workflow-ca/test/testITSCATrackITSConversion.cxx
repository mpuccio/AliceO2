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

// Gate 3 workflow-onboarding Slice 2: focused, independent tests for
// convertTrackITSExtToTrackITS() (CATrackerSpec.h) -- the TrackITSExt ->
// TrackITS conversion used by the opt-in ITS common-CA tracker's output
// stage. Deliberately does not construct a TimeFrame: exercises the per-track
// flattening/RangeRefComp logic directly against a hand-built TrackITSExt.

#define BOOST_TEST_MODULE ITSMFT ITSCATrackITSConversion
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include "ITSCAWorkflow/CATrackerSpec.h"

using namespace o2::its::ca;

namespace
{
o2::its::TrackITSExt makeTrack(int nClusters)
{
  o2::its::TrackITSExt track;
  track.setNumberOfClusters(nClusters);
  return track;
}
} // namespace

BOOST_AUTO_TEST_CASE(NoHitsProducesEmptyClusterRange)
{
  auto track = makeTrack(0);
  std::vector<int> clusterIndices;
  std::vector<o2::its::TrackITS> tracks;

  convertTrackITSExtToTrackITS(track, clusterIndices, tracks);

  BOOST_REQUIRE_EQUAL(tracks.size(), 1u);
  BOOST_CHECK(clusterIndices.empty());
  BOOST_CHECK_EQUAL(tracks[0].getFirstClusterEntry(), 0);
  BOOST_CHECK_EQUAL(tracks[0].getNumberOfClusters(), 0);
}

BOOST_AUTO_TEST_CASE(HitsAreFlattenedInDecreasingLayerOrder)
{
  // setExternalClusterIndex(..., true) bumps the cluster count by one per
  // call (see TrackITS.h), so start from an empty track rather than
  // pre-setting the count.
  auto track = makeTrack(0);
  // hits on layers 0, 3, 6 with distinct external cluster indices
  track.setExternalClusterIndex(0, 100, true);
  track.setExternalClusterIndex(3, 103, true);
  track.setExternalClusterIndex(6, 106, true);
  // Per-layer cluster sizes are set upstream, by
  // Tracker<NLayers>::rectifyClusterIndices(), before the cluster index is
  // rewritten to the external identity above -- convertTrackITSExtToTrackITS
  // must carry these through unchanged, not re-derive them from the
  // (already external) cluster index.
  track.setClusterSize(0, 2);
  track.setClusterSize(3, 5);
  track.setClusterSize(6, 9);

  std::vector<int> clusterIndices;
  std::vector<o2::its::TrackITS> tracks;

  convertTrackITSExtToTrackITS(track, clusterIndices, tracks);

  BOOST_REQUIRE_EQUAL(tracks.size(), 1u);
  // Decreasing layer order: layer 6 first, then 3, then 0.
  BOOST_REQUIRE_EQUAL(clusterIndices.size(), 3u);
  BOOST_CHECK_EQUAL(clusterIndices[0], 106);
  BOOST_CHECK_EQUAL(clusterIndices[1], 103);
  BOOST_CHECK_EQUAL(clusterIndices[2], 100);

  BOOST_CHECK_EQUAL(tracks[0].getFirstClusterEntry(), 0);
  BOOST_CHECK_EQUAL(tracks[0].getNumberOfClusters(), 3);
  BOOST_CHECK_EQUAL(tracks[0].getClusterSize(6), 9);
  BOOST_CHECK_EQUAL(tracks[0].getClusterSize(3), 5);
  BOOST_CHECK_EQUAL(tracks[0].getClusterSize(0), 2);
}

BOOST_AUTO_TEST_CASE(FirstClusterEntryOffsetsIntoSharedFlattenedArray)
{
  // Two tracks converted in sequence must append into the same
  // clusterIndices array, with the second track's firstClusterEntry offset
  // by the first track's cluster count -- this is the multi-track
  // RangeRefComp contract the writer/downstream readers depend on.
  auto trackA = makeTrack(0);
  trackA.setExternalClusterIndex(0, 10, true);
  trackA.setExternalClusterIndex(1, 11, true);

  auto trackB = makeTrack(0);
  trackB.setExternalClusterIndex(2, 22, true);

  std::vector<int> clusterIndices;
  std::vector<o2::its::TrackITS> tracks;

  convertTrackITSExtToTrackITS(trackA, clusterIndices, tracks);
  convertTrackITSExtToTrackITS(trackB, clusterIndices, tracks);

  BOOST_REQUIRE_EQUAL(tracks.size(), 2u);
  BOOST_REQUIRE_EQUAL(clusterIndices.size(), 3u);

  BOOST_CHECK_EQUAL(tracks[0].getFirstClusterEntry(), 0);
  BOOST_CHECK_EQUAL(tracks[0].getNumberOfClusters(), 2);
  BOOST_CHECK_EQUAL(tracks[1].getFirstClusterEntry(), 2);
  BOOST_CHECK_EQUAL(tracks[1].getNumberOfClusters(), 1);
  BOOST_CHECK_EQUAL(clusterIndices[tracks[1].getFirstClusterEntry()], 22);
}

BOOST_AUTO_TEST_CASE(SourceTrackIsNotMutatedThroughCallerReference)
{
  // convertTrackITSExtToTrackITS takes `track` by value, so a caller passing
  // a TimeFrame-owned track by (const) reference must observe it unchanged
  // afterwards.
  auto original = makeTrack(0);
  original.setExternalClusterIndex(4, 40, true);
  const auto originalFirstEntry = original.getFirstClusterEntry();
  const auto originalNCl = original.getNumberOfClusters();

  std::vector<int> clusterIndices;
  std::vector<o2::its::TrackITS> tracks;
  // Bias clusterIndices with unrelated content so setFirstClusterEntry
  // inside the conversion would visibly diverge from originalFirstEntry if
  // it mutated the caller's object.
  clusterIndices.push_back(-999);

  convertTrackITSExtToTrackITS(original, clusterIndices, tracks);

  BOOST_CHECK_EQUAL(original.getFirstClusterEntry(), originalFirstEntry);
  BOOST_CHECK_EQUAL(original.getNumberOfClusters(), originalNCl);
  BOOST_CHECK_EQUAL(tracks[0].getFirstClusterEntry(), 1); // offset by the pre-existing sentinel entry
}

BOOST_AUTO_TEST_CASE(ClusterSizeSurvivesLargeNonMonotonicExternalIndicesUnchanged)
{
  // Regression test for the local-vs-external cluster-index bug this
  // function used to have: it used to re-derive each cluster's size via a
  // caller-supplied `clusterSizeAt(layer, clid)` callback, with `clid`
  // already being the external/global cluster identity by the time a track
  // reaches this function (Tracker<NLayers>::rectifyClusterIndices() in
  // CATracker.cxx has already overwritten the layer-local identity in
  // place) -- exactly the domain mismatch that produced wrong or
  // out-of-bounds reads whenever a TimeFrame's per-layer mClusterSize
  // vector was indexed with that external id instead of a layer-local one.
  //
  // convertTrackITSExtToTrackITS no longer takes a TimeFrame/callback at
  // all: it only carries through whatever size rectifyClusterIndices()
  // already set on the track (via TrackITS::setClusterSize()) while the
  // layer-local identity was still available. This test proves that by
  // construction: the external cluster indices below are deliberately
  // large and non-monotonic in layer order (the exact "future multi-source
  // identity" scenario the fix must not assume away), and the emitted
  // sizes must match the sizes set on the input track exactly, regardless
  // of those external index values -- there is no TimeFrame here for a
  // wrong index to alias into.
  auto track = makeTrack(0);
  track.setExternalClusterIndex(0, 500000, true); // large, out of any plausible per-layer range
  track.setExternalClusterIndex(2, 7, true);      // small, but out of layer order (non-monotonic)
  track.setExternalClusterIndex(5, 123456, true);
  track.setClusterSize(0, 3);
  track.setClusterSize(2, 11);
  track.setClusterSize(5, 15); // clamped to TrackITS::setClusterSize's 15-max, still distinctive

  std::vector<int> clusterIndices;
  std::vector<o2::its::TrackITS> tracks;

  convertTrackITSExtToTrackITS(track, clusterIndices, tracks);

  BOOST_REQUIRE_EQUAL(tracks.size(), 1u);
  BOOST_REQUIRE_EQUAL(clusterIndices.size(), 3u);
  // Decreasing layer order: layer 5, then 2, then 0.
  BOOST_CHECK_EQUAL(clusterIndices[0], 123456);
  BOOST_CHECK_EQUAL(clusterIndices[1], 7);
  BOOST_CHECK_EQUAL(clusterIndices[2], 500000);

  BOOST_CHECK_EQUAL(tracks[0].getClusterSize(5), 15);
  BOOST_CHECK_EQUAL(tracks[0].getClusterSize(2), 11);
  BOOST_CHECK_EQUAL(tracks[0].getClusterSize(0), 3);
}
