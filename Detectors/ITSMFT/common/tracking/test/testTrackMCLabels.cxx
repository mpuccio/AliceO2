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

#define BOOST_TEST_MODULE ITSMFT track MC labels
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <array>
#include <initializer_list>
#include <memory>
#include <vector>

#include <gsl/span>

#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/GenericTrack.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "TraversalTestSupport.h"

namespace
{
using o2::MCCompLabel;
using namespace o2::itsmft::tracking;

MCCompLabel label(int track, bool fake = false, int event = 0, int source = 0)
{
  return {track, event, source, fake};
}

void checkIdentityAndFake(const MCCompLabel& actual, const MCCompLabel& expected, bool fake)
{
  BOOST_CHECK(actual == expected);
  BOOST_CHECK_EQUAL(actual.getTrackEventSourceID(), expected.getTrackEventSourceID());
  BOOST_CHECK_EQUAL(actual.isFake(), fake);
}

class LabelFixture
{
 public:
  LabelFixture()
  {
    for (uint16_t layer = 0; layer < mSurfaces.size(); ++layer) {
      mSurfaces[layer] = {layer, static_cast<uint8_t>(o2::detectors::DetID::ITS), SurfaceKind::Cylinder};
    }
    BOOST_REQUIRE(mFrame.configure(DetectorLayout{mSurfaces, makeDetectorLayout()}, 0, 0,
                                   std::make_shared<BoundedMemoryResource>()));
    mFrame.setHasMCInformation(true);
  }

  TrackClusterReference addCluster(std::initializer_list<MCCompLabel> labels)
  {
    BOOST_REQUIRE_LT(mNextLayer, mSurfaces.size());
    const LayerId layer{static_cast<uint16_t>(mNextLayer++)};
    mFrame.addMeasurement(layer, GlobalMeasurement{}, SurfaceMeasurement{},
                          gsl::span<const MCCompLabel>{labels.begin(), labels.size()});
    return {layer, 0, 0};
  }

  void addTrack(std::initializer_list<TrackClusterReference> references)
  {
    GenericTrack track;
    track.firstClusterRef = static_cast<uint32_t>(mFrame.getTrackClusterIndices().size());
    mFrame.getTrackClusterIndices().insert(mFrame.getTrackClusterIndices().end(), references.begin(), references.end());
    track.clusterRefEnd = static_cast<uint32_t>(mFrame.getTrackClusterIndices().size());
    mFrame.getGenericTracks().push_back(track);
  }

  void compute() { TrackerTestAccess::computeTracksMClabels(mTracker, mFrame); }

  TimeFrame& frame() { return mFrame; }

 private:
  std::array<SurfaceDescriptor, 8> mSurfaces{};
  std::size_t mNextLayer{0};
  TimeFrame mFrame;
  Tracker mTracker;
};
} // namespace

BOOST_AUTO_TEST_CASE(labels_are_stored_in_generic_track_order)
{
  LabelFixture fixture;
  const auto a = label(1);
  const auto b = label(2);
  const auto a0 = fixture.addCluster({a});
  const auto a1 = fixture.addCluster({a});
  const auto b0 = fixture.addCluster({b});
  fixture.addTrack({a0, a1});
  fixture.addTrack({b0});

  fixture.compute();

  const auto& labels = fixture.frame().getTrackLabels();
  BOOST_REQUIRE_EQUAL(labels.size(), 2u);
  checkIdentityAndFake(labels[0], a, false);
  checkIdentityAndFake(labels[1], b, false);
}

BOOST_AUTO_TEST_CASE(one_stray_or_empty_cluster_makes_track_fake)
{
  LabelFixture fixture;
  const auto a = label(1);
  const auto b = label(2);
  const auto a0 = fixture.addCluster({a});
  const auto a1 = fixture.addCluster({a});
  const auto stray = fixture.addCluster({b});
  const auto empty = fixture.addCluster({});
  fixture.addTrack({a0, a1, stray});
  fixture.addTrack({a0, a1, empty});

  fixture.compute();

  const auto& labels = fixture.frame().getTrackLabels();
  BOOST_REQUIRE_EQUAL(labels.size(), 2u);
  checkIdentityAndFake(labels[0], a, true);
  checkIdentityAndFake(labels[1], a, true);
}

BOOST_AUTO_TEST_CASE(all_candidates_on_each_cluster_are_counted_once)
{
  LabelFixture fixture;
  const auto a = label(1);
  const auto b = label(2);
  const auto aOnly = fixture.addCluster({a});
  const auto both = fixture.addCluster({a, b});
  const auto duplicateB = fixture.addCluster({b, b, b});
  const auto bOnly = fixture.addCluster({b});
  fixture.addTrack({aOnly, both, duplicateB, bOnly});

  fixture.compute();

  const auto& labels = fixture.frame().getTrackLabels();
  BOOST_REQUIRE_EQUAL(labels.size(), 1u);
  // B occurs on three clusters and A on two. The old ITS loop failed to add B
  // as a candidate on `both`, while duplicate labels could overcount a cluster.
  checkIdentityAndFake(labels[0], b, true);
}

BOOST_AUTO_TEST_CASE(ties_use_the_first_encountered_identity)
{
  LabelFixture fixture;
  const auto a = label(1);
  const auto b = label(2);
  const auto bCluster = fixture.addCluster({b});
  const auto aCluster = fixture.addCluster({a});
  fixture.addTrack({bCluster, aCluster});

  fixture.compute();

  const auto& labels = fixture.frame().getTrackLabels();
  BOOST_REQUIRE_EQUAL(labels.size(), 1u);
  checkIdentityAndFake(labels[0], b, true);
}

BOOST_AUTO_TEST_CASE(no_labels_returns_fake_unset_and_no_mc_clears_the_sidecar)
{
  LabelFixture fixture;
  const auto empty = fixture.addCluster({});
  fixture.addTrack({empty});
  fixture.addTrack({});

  fixture.compute();

  const auto& labels = fixture.frame().getTrackLabels();
  BOOST_REQUIRE_EQUAL(labels.size(), 2u);
  for (const auto& result : labels) {
    BOOST_CHECK(result.isFake());
    BOOST_CHECK(!result.isSet());
  }

  fixture.frame().setHasMCInformation(false);
  fixture.compute();
  BOOST_CHECK(fixture.frame().getTrackLabels().empty());
}

BOOST_AUTO_TEST_CASE(timeframe_reset_clears_tracks_labels_and_references_together)
{
  LabelFixture fixture;
  const auto a = label(1);
  const auto cluster = fixture.addCluster({a});
  fixture.addTrack({cluster});
  fixture.compute();
  BOOST_REQUIRE_EQUAL(fixture.frame().getGenericTracks().size(), 1u);
  BOOST_REQUIRE_EQUAL(fixture.frame().getTrackLabels().size(), 1u);
  BOOST_REQUIRE_EQUAL(fixture.frame().getTrackClusterIndices().size(), 1u);

  fixture.frame().resetTimeFrame();

  BOOST_CHECK(fixture.frame().getGenericTracks().empty());
  BOOST_CHECK(fixture.frame().getTrackLabels().empty());
  BOOST_CHECK(fixture.frame().getTrackClusterIndices().empty());
}
