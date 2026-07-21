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

#define BOOST_TEST_MODULE ITSMFT MC label accumulator
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <gsl/gsl>

#include "ITSMFTTracking/MCLabelAccumulator.h"

namespace
{
using o2::MCCompLabel;
using o2::itsmft::tracking::MCLabelAccumulator;

MCCompLabel label(int track, bool fake = false, int event = 0, int source = 0)
{
  return {track, event, source, fake};
}

MCCompLabel reduce(std::initializer_list<MCCompLabel> labels)
{
  MCLabelAccumulator accumulator;
  accumulator.addCluster({labels.begin(), labels.size()});
  return accumulator.finalize();
}

void checkIdentityAndFake(const MCCompLabel& actual, const MCCompLabel& expected, bool fake)
{
  BOOST_CHECK(actual == expected);
  BOOST_CHECK_EQUAL(actual.getTrackEventSourceID(), expected.getTrackEventSourceID());
  BOOST_CHECK_EQUAL(actual.isFake(), fake);
}
} // namespace

BOOST_AUTO_TEST_CASE(one_identity_on_every_cluster_is_correct)
{
  MCLabelAccumulator accumulator;
  const auto a = label(1);
  accumulator.addCluster({&a, 1});
  accumulator.addCluster({&a, 1});
  accumulator.addCluster({&a, 1});
  checkIdentityAndFake(accumulator.finalize(), a, false);
}

BOOST_AUTO_TEST_CASE(disagreeing_or_empty_cluster_makes_winner_fake)
{
  MCLabelAccumulator disagreement;
  const auto a = label(1);
  const auto b = label(2);
  disagreement.addCluster({&a, 1});
  disagreement.addCluster({&a, 1});
  disagreement.addCluster({&b, 1});
  checkIdentityAndFake(disagreement.finalize(), a, true);

  MCLabelAccumulator empty;
  empty.addCluster({&a, 1});
  empty.addCluster({});
  empty.addCluster({&a, 1});
  checkIdentityAndFake(empty.finalize(), a, true);
}

BOOST_AUTO_TEST_CASE(no_attached_clusters_or_no_labels_returns_default_fake_unset)
{
  MCLabelAccumulator noClusters;
  const auto noClustersResult = noClusters.finalize();
  BOOST_CHECK(noClustersResult.isFake());
  BOOST_CHECK(!noClustersResult.isSet());

  MCLabelAccumulator noLabels;
  noLabels.addCluster({});
  noLabels.addCluster({});
  const auto noLabelsResult = noLabels.finalize();
  BOOST_CHECK(noLabelsResult.isFake());
  BOOST_CHECK(!noLabelsResult.isSet());
}

BOOST_AUTO_TEST_CASE(existing_and_new_candidate_on_same_cluster_are_both_processed)
{
  MCLabelAccumulator accumulator;
  const auto a = label(1);
  const auto b = label(2);
  accumulator.addCluster({&a, 1});
  const MCCompLabel both[] = {a, b};
  accumulator.addCluster(both);
  accumulator.addCluster({&b, 1});
  accumulator.addCluster({&b, 1});

  // The old cluster-global-found implementation dropped B from `both`,
  // leaving A and B tied at two clusters and incorrectly selecting A.
  // B must instead win with three of four attached clusters and be fake.
  checkIdentityAndFake(accumulator.finalize(), b, true);
}

BOOST_AUTO_TEST_CASE(duplicates_within_one_cluster_count_once)
{
  MCLabelAccumulator accumulator;
  const auto a = label(1);
  const MCCompLabel duplicate[] = {a, a, a};
  accumulator.addCluster(duplicate);
  accumulator.addCluster({});
  checkIdentityAndFake(accumulator.finalize(), a, true);
}

BOOST_AUTO_TEST_CASE(fake_variants_share_identity_and_first_representative)
{
  const auto correct = label(1);
  const auto fake = label(1, true);

  MCLabelAccumulator firstCorrect;
  const MCCompLabel correctThenFake[] = {correct, fake};
  firstCorrect.addCluster(correctThenFake);
  firstCorrect.addCluster({&fake, 1});
  checkIdentityAndFake(firstCorrect.finalize(), correct, false);

  MCLabelAccumulator firstFake;
  const MCCompLabel fakeThenCorrect[] = {fake, correct};
  firstFake.addCluster(fakeThenCorrect);
  firstFake.addCluster({&correct, 1});
  checkIdentityAndFake(firstFake.finalize(), fake, true);
}

BOOST_AUTO_TEST_CASE(tie_uses_first_occurrence_and_later_higher_count_wins)
{
  const auto a = label(1);
  const auto b = label(2);

  MCLabelAccumulator tie;
  tie.addCluster({&a, 1});
  tie.addCluster({&b, 1});
  checkIdentityAndFake(tie.finalize(), a, true);

  MCLabelAccumulator higher;
  higher.addCluster({&a, 1});
  higher.addCluster({&b, 1});
  higher.addCluster({&b, 1});
  checkIdentityAndFake(higher.finalize(), b, true);
}

BOOST_AUTO_TEST_CASE(noise_qed_unset_and_input_fake_are_not_filtered)
{
  const MCCompLabel noise{true};
  const auto qed = label(3, false, 0, 99);
  const MCCompLabel unset{};
  const auto inputFake = label(4, true);

  checkIdentityAndFake(reduce({noise}), noise, noise.isFake());
  BOOST_CHECK(reduce({noise}).isNoise());
  checkIdentityAndFake(reduce({qed}), qed, false);
  BOOST_CHECK(reduce({qed}).isNoise());
  const auto unsetResult = reduce({unset});
  BOOST_CHECK(unsetResult == unset);
  BOOST_CHECK(unsetResult.isFake());
  BOOST_CHECK(!unsetResult.isSet());
  checkIdentityAndFake(reduce({inputFake}), inputFake, true);
}

BOOST_AUTO_TEST_CASE(reduction_is_deterministic_and_order_controls_only_documented_contract)
{
  const auto a = label(1);
  const auto b = label(2);
  const auto fakeA = label(1, true);

  MCLabelAccumulator first;
  MCLabelAccumulator second;
  const MCCompLabel sequence[][2] = {{a, b}, {fakeA, b}, {a, b}};
  for (const auto& cluster : sequence) {
    first.addCluster(cluster);
    second.addCluster(cluster);
  }
  BOOST_CHECK_EQUAL(first.finalize().getRawValue(), second.finalize().getRawValue());

  MCLabelAccumulator representativeOrder;
  const MCCompLabel fakeFirst[] = {fakeA, a};
  representativeOrder.addCluster(fakeFirst);
  checkIdentityAndFake(representativeOrder.finalize(), fakeA, true);

  MCLabelAccumulator tieOrder;
  tieOrder.addCluster({&b, 1});
  tieOrder.addCluster({&a, 1});
  checkIdentityAndFake(tieOrder.finalize(), b, true);
}
