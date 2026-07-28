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

// Gate 3 workflow-onboarding Slice 2: focused, independent tests for
// convertTrackITSExtToTrackITS() (CATrackerSpec.h) -- the TrackITSExt ->
// TrackITS conversion used by the opt-in ITS common-CA tracker's output
// stage. Deliberately does not construct a TimeFrame: exercises the per-track
// flattening/RangeRefComp logic directly against a hand-built TrackITSExt.

#define BOOST_TEST_MODULE ITSMFT ITSCATrackITSConversion
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <map>
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

  convertTrackITSExtToTrackITS(track, clusterIndices, tracks, [](int, int) { return 0; });

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

  std::vector<int> clusterIndices;
  std::vector<o2::its::TrackITS> tracks;
  const std::map<int, int> sizeByExtIdx{{100, 2}, {103, 5}, {106, 9}};

  convertTrackITSExtToTrackITS(track, clusterIndices, tracks,
                               [&sizeByExtIdx](int /*layer*/, int extIdx) { return sizeByExtIdx.at(extIdx); });

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
  auto clusterSizeAt = [](int, int) { return 1; };

  convertTrackITSExtToTrackITS(trackA, clusterIndices, tracks, clusterSizeAt);
  convertTrackITSExtToTrackITS(trackB, clusterIndices, tracks, clusterSizeAt);

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

  convertTrackITSExtToTrackITS(original, clusterIndices, tracks, [](int, int) { return 3; });

  BOOST_CHECK_EQUAL(original.getFirstClusterEntry(), originalFirstEntry);
  BOOST_CHECK_EQUAL(original.getNumberOfClusters(), originalNCl);
  BOOST_CHECK_EQUAL(tracks[0].getFirstClusterEntry(), 1); // offset by the pre-existing sentinel entry
}
