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

/// \file testStageComparison.cxx
/// \brief compareStageOutputs reports each kind of difference as what it is.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK StageComparison
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Study/StageComparison.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "SerialReferences.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstring>
#include <cstdint>
#include <memory>
#include <vector>

namespace o2::gtrk::study
{
namespace
{

using Output = StageOutput<TrackletT>;
using Comparison = StageComparison<TrackletT>;
using Elements = std::vector<std::vector<Tracklet>>;

constexpr std::size_t kSourceColumn = 0U;
constexpr std::size_t kTargetColumn = 1U;
constexpr std::size_t kLambdaColumn = 3U;
constexpr std::size_t kPhiColumn = 4U;

/// Four elements, the second empty, the others with one to three records.
Elements sample()
{
  const auto tracklet = [](std::uint32_t source, std::uint32_t target, Scalar lambda) {
    return Tracklet{.source = source, .target = target, .link = 1U, .lambda = lambda, .phi = 0.5F};
  };
  return {{tracklet(0U, 10U, 0.1F), tracklet(0U, 11U, 0.2F)},
          {},
          {tracklet(2U, 12U, -0.3F)},
          {tracklet(3U, 13U, 0.4F), tracklet(3U, 14U, 0.5F), tracklet(3U, 15U, 0.6F)}};
}

std::unique_ptr<Output> output(const Elements& elements)
{
  auto result = std::make_unique<Output>();
  result->setMemoryPool(std::make_shared<BoundedMemoryResource>());
  result->prepare(static_cast<Output::Index>(elements.size()));
  const auto generate = [&](Output::Index element, auto& emitter) {
    for (const auto& record : elements[element]) {
      emitter.emit(record);
    }
  };
  BOOST_REQUIRE(runSerial(*result, generate) == StageStatus::kComplete);
  return result;
}

} // namespace

/// A cut decided differently: one element has a record more.
BOOST_AUTO_TEST_CASE(AChangedCountIsReportedPerElement)
{
  auto changed = sample();
  changed[1].push_back(Tracklet{.source = 1U, .target = 9U, .link = 0U, .lambda = 0.F, .phi = 0.F});
  const auto comparison = compareStageOutputs(*output(sample()), *output(changed));
  BOOST_TEST(!comparison.identical());
  BOOST_TEST(comparison.elementsWithOtherCount == 1U);
  BOOST_TEST(comparison.recordsDiffering == 0U);
  BOOST_REQUIRE(comparison.first.size() == 1U);
  BOOST_TEST(comparison.first[0].element == 1U);
  BOOST_TEST(comparison.first[0].column == Comparison::kCount);
  BOOST_TEST(comparison.first[0].reference == 0.);
  BOOST_TEST(comparison.first[0].other == 1.);
  // The later elements, shifted by the extra record, still compare as equal.
  BOOST_TEST(comparison.elementsReordered == 0U);
}

/// A value in its last bit: counted in its column, with the difference.
BOOST_AUTO_TEST_CASE(AChangedFieldIsReportedInItsColumn)
{
  auto changed = sample();
  changed[3][1].lambda = std::nextafter(changed[3][1].lambda, 1.F);
  changed[0][0].target = 99U;
  const auto comparison = compareStageOutputs(*output(sample()), *output(changed));
  BOOST_TEST(!comparison.identical());
  BOOST_TEST(comparison.elementsWithOtherCount == 0U);
  BOOST_TEST(comparison.elementsReordered == 0U);
  BOOST_TEST(comparison.recordsDiffering == 2U);
  BOOST_TEST(comparison.differingByColumn[kLambdaColumn] == 1U);
  BOOST_TEST(comparison.differingByColumn[kTargetColumn] == 1U);
  BOOST_TEST(comparison.differingByColumn[kSourceColumn] == 0U);
  BOOST_TEST(comparison.differingByColumn[kPhiColumn] == 0U);
  BOOST_TEST(comparison.maxAbsDifference[kLambdaColumn] > 0.);
  BOOST_TEST(comparison.maxAbsDifference[kLambdaColumn] < 1.e-6);
  BOOST_TEST(comparison.maxAbsDifference[kTargetColumn] == 89.);
  BOOST_REQUIRE(comparison.first.size() == 2U);
  BOOST_TEST(comparison.first[0].element == 0U);
  BOOST_TEST(comparison.first[0].ordinal == 0U);
  BOOST_TEST(comparison.first[0].column == kTargetColumn);
  BOOST_TEST(comparison.first[1].element == 3U);
  BOOST_TEST(comparison.first[1].ordinal == 1U);
  BOOST_TEST(comparison.first[1].column == kLambdaColumn);
}

/// The same records of an element in another order: its order differs, not
/// its content.
BOOST_AUTO_TEST_CASE(AMovedRecordIsAReordering)
{
  auto moved = sample();
  std::swap(moved[3][0], moved[3][2]);
  const auto comparison = compareStageOutputs(*output(sample()), *output(moved));
  BOOST_TEST(!comparison.identical());
  BOOST_TEST(comparison.elementsReordered == 1U);
  BOOST_TEST(comparison.recordsDiffering == 0U);
  BOOST_TEST(comparison.elementsWithOtherCount == 0U);
}

BOOST_AUTO_TEST_CASE(OtherElementCountsAreNotCompared)
{
  auto shorter = sample();
  shorter.pop_back();
  const auto comparison = compareStageOutputs(*output(sample()), *output(shorter));
  BOOST_TEST(!comparison.sameElements);
  BOOST_TEST(comparison.referenceElements == 4U);
  BOOST_TEST(comparison.otherElements == 3U);
  BOOST_TEST(!comparison.identical());
}

namespace
{

constexpr std::size_t kRoadRows = 17U;
using Records = StageOutput<SelectionRecordOf<kRoadRows>::RecordT>;
using SelectionRecord = Records::Record;

SelectionRecord selectionRecord(std::uint32_t first, std::uint8_t count, float chi2)
{
  SelectionRecord record;
  for (std::uint32_t i = 0; i < kRoadRows; ++i) {
    record.clusters.index[i] = i < count ? first + i : ClusterList<kRoadRows>::kNoCluster;
  }
  record.count = count;
  record.chi2 = chi2;
  record.kinks = 0U;
  return record;
}

std::unique_ptr<Records> recordsOf(const std::vector<SelectionRecord>& records)
{
  auto output = std::make_unique<Records>();
  output->setMemoryPool(std::make_shared<BoundedMemoryResource>());
  output->prepare(static_cast<Records::Index>(records.size()));
  BOOST_REQUIRE(runSerial(*output, [&](Records::Index element, auto& emitter) { emitter.emit(records[element]); }) == StageStatus::kComplete);
  return output;
}

} // namespace

/// A column wider than 64 bits -- a road's 17 rows -- is compared by its
/// bytes; it has no numeric difference.
BOOST_AUTO_TEST_CASE(AWideColumnIsComparedByItsBytes)
{
  const std::vector<SelectionRecord> records{selectionRecord(10U, 8U, 1.5F), selectionRecord(40U, 11U, 2.5F), selectionRecord(70U, 17U, 0.5F)};
  auto changed = records;
  changed[1].clusters.index[9] += 1U;
  const auto comparison = compareStageOutputs(*recordsOf(records), *recordsOf(changed));
  BOOST_TEST(!comparison.identical());
  BOOST_TEST(comparison.recordsDiffering == 1U);
  BOOST_TEST(comparison.differingByColumn[0] == 1U);
  BOOST_TEST(comparison.differingByColumn[1] == 0U);
  BOOST_TEST(comparison.differingByColumn[2] == 0U);
  BOOST_TEST(comparison.maxAbsDifference[0] == 0.);
  BOOST_REQUIRE(comparison.first.size() == 1U);
  BOOST_TEST(comparison.first[0].element == 1U);
  BOOST_TEST(comparison.first[0].column == 0U);
  BOOST_TEST(compareStageOutputs(*recordsOf(records), *recordsOf(records)).identical());
}

/// Column sets without elements -- an iteration's gathered records -- row by
/// row.
BOOST_AUTO_TEST_CASE(ColumnsCompareRowByRow)
{
  const std::vector<SelectionRecord> records{selectionRecord(1U, 4U, 1.F), selectionRecord(9U, 5U, 3.F)};
  auto changed = records;
  changed[1].chi2 = 3.25F;
  const auto reference = recordsOf(records);
  const auto other = recordsOf(changed);
  const auto same = compareColumns<SelectionRecordOf<kRoadRows>::RecordT>(reference->view(), reference->size(), reference->view(), reference->size());
  BOOST_TEST(same.identical());
  const auto comparison = compareColumns<SelectionRecordOf<kRoadRows>::RecordT>(reference->view(), reference->size(), other->view(), other->size());
  BOOST_TEST(comparison.recordsDiffering == 1U);
  BOOST_TEST(comparison.differingByColumn[2] == 1U);
  BOOST_TEST(comparison.maxAbsDifference[2] == 0.25);
  BOOST_REQUIRE(comparison.first.size() == 1U);
  BOOST_TEST(comparison.first[0].element == 1U);
  const auto shorter = compareColumns<SelectionRecordOf<kRoadRows>::RecordT>(reference->view(), reference->size(), other->view(), 1U);
  BOOST_TEST(!shorter.sameElements);
}

/// Tracks row by row: a covariance entry of a state makes its record differ,
/// in the state's column alone.
BOOST_AUTO_TEST_CASE(TracksCompareTheirStates)
{
  using Tracks = BoundedColumns<TrackOf<kRoadRows>::RecordT>;
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto make = [&]() {
    auto tracks = std::make_unique<Tracks>();
    tracks->setMemoryResource(pool.get());
    tracks->resize(2U);
    const auto view = tracks->view();
    for (std::size_t row = 0; row < 2U; ++row) {
      for (auto* state : {&view.state[row], &view.outwardState[row]}) {
        *state = BoundTrackState{};
        state->par[0] = 1.F + static_cast<Scalar>(row);
        state->cov[3] = 0.5F;
        state->surface = 7U;
      }
      view.timePS[row] = 0U;
      view.widthPS[row] = 0U;
      view.chi2[row] = 0.F;
      view.outwardChi2[row] = 0.F;
      view.ndf[row] = 0;
      view.status[row] = PropagationStatus::kOk;
      view.count[row] = 0U;
      view.iteration[row] = 0U;
      view.addedOutward[row] = 0U;
      view.addedInward[row] = 0U;
      view.kinks[row] = 0U;
      view.kinkCluster[row] = 0U;
      view.clusters[row] = {};
    }
    return tracks;
  };
  const auto reference = make();
  const auto other = make();
  BOOST_TEST(compareColumns<TrackOf<kRoadRows>::RecordT>(std::as_const(*reference).view(), 2U, std::as_const(*other).view(), 2U).identical());

  other->view().state[1].cov[3] = 0.75F;
  const auto differing = compareColumns<TrackOf<kRoadRows>::RecordT>(std::as_const(*reference).view(), 2U, std::as_const(*other).view(), 2U);
  constexpr auto kState = static_cast<std::size_t>(TrackOf<kRoadRows>::RecordT<memlayout::value>::Column::state);
  BOOST_TEST(differing.recordsDiffering == 1U);
  BOOST_TEST(differing.differingByColumn[kState] == 1U);
  BOOST_TEST(std::ranges::count_if(differing.differingByColumn, [](std::size_t count) { return count != 0U; }) == 1);
  BOOST_REQUIRE(differing.first.size() == 1U);
  BOOST_TEST(differing.first[0].element == 1U);
}

} // namespace o2::gtrk::study
