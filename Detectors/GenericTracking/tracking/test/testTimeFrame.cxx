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

/// \file testTimeFrame.cxx
/// \brief Block-aligned slices over borrowed, grid-ordered measurements.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TimeFrame
#include <boost/test/unit_test.hpp>

#include "ToyBarrel3.h"
#include "GenericTracking/Tracking/TimeFrame.h"

#include "GridTestSupport.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>

namespace o2::gtrk
{

namespace
{

using Frame = TimeFrame<detectors::ToyBarrel3>;
using ToyGridGeometry = GridGeometry<detectors::ToyBarrel3>;

ToyGridGeometry gridGeometry()
{
  GridParams params;
  params.blockWidth = 10'000;
  return test::cappedGrid<detectors::ToyBarrel3>(30., test::withToyFrames(params), SearchWindowParams{}, VertexConstraint{}, 1, 1);
}

/// Rows with only their box starts [ns, in ps] and identities set, as a frame
/// would read them.
template <std::size_t N>
MeasurementStore timeColumns(const std::array<double, N>& times)
{
  MeasurementStore columns;
  resizeMeasurementColumns(columns, N);
  for (std::size_t row = 0; row < N; ++row) {
    columns.timePS[row] = static_cast<std::uint32_t>(std::llround(times[row] * 1000.));
    columns.clusterIndex[row] = static_cast<std::uint32_t>(row);
  }
  return columns;
}

} // namespace

BOOST_AUTO_TEST_CASE(HalfOpenCoresAndWideHalos)
{
  const std::array<double, 6> times{0., 10., 5., 10., 10., 20.};
  const auto columns = timeColumns(times);
  const auto geometry = gridGeometry();
  const std::array<Frame::Index, 10> offsets{0, 1, 2, 2, 3, 4, 4, 4, 5, 6};
  const Frame frame{measurementView(columns), geometry, offsets};

  BOOST_TEST(frame.size() == 6U);
  BOOST_TEST(frame.grid().layerRows(1).size() == 2U);
  BOOST_TEST(frame.grid().layerRows(1).begin == 2U);
  BOOST_TEST(&frame.measurement(3).timePS == &columns.timePS[3]);

  const Frame::TimeSlice slice{frame, {.begin = 1, .end = 2}, 5'000};
  BOOST_TEST(slice.coreBlocks.begin == 1);
  BOOST_TEST(slice.coreBlocks.end == 2);
  BOOST_TEST(slice.haloBlocks.begin == 0);
  BOOST_TEST(slice.haloBlocks.end == 3);
  BOOST_TEST(slice.core[0].begin == 1U);
  BOOST_TEST(slice.core[0].end == 2U);
  BOOST_TEST(slice.core[1].begin == 3U);
  BOOST_TEST(slice.core[1].end == 4U);
  BOOST_TEST(slice.core[2].begin == 4U);
  BOOST_TEST(slice.core[2].end == 5U);
  BOOST_TEST(slice.halo[1].begin == 2U);
  BOOST_TEST(slice.halo[1].end == 4U);
  BOOST_TEST(slice.halo[2].end == 6U);

  std::array<unsigned, 6> coreVisits{};
  for (std::int32_t block = 0; block < geometry.blockCount(); ++block) {
    const Frame::TimeSlice part{frame, {.begin = block, .end = block + 1}, 0};
    for (const auto range : part.core) {
      for (auto index = range.begin; index < range.end; ++index) {
        ++coreVisits[index];
      }
    }
  }
  for (const auto visits : coreVisits) {
    BOOST_TEST(visits == 1U);
  }

  unsigned haloVisits = 0;
  for (std::int32_t block = 0; block < geometry.blockCount(); ++block) {
    const Frame::TimeSlice part{frame, {.begin = block, .end = block + 1}, 12'000};
    if (part.halo[0].begin <= 1U && 1U < part.halo[0].end) {
      ++haloVisits;
    }
  }
  BOOST_TEST(haloVisits > 2U);
}

BOOST_AUTO_TEST_CASE(EmptyLayersAndEmptyTimeframe)
{
  const std::array<double, 3> times{0., 10., 15.};
  const auto columns = timeColumns(times);
  const auto geometry = gridGeometry();
  const std::array<Frame::Index, 10> offsets{0, 1, 2, 2, 2, 2, 2, 2, 3, 3};
  const Frame frame{measurementView(columns), geometry, offsets};
  BOOST_TEST(frame.grid().layerRows(1).empty());

  const Frame::TimeSlice slice{frame, {.begin = 1, .end = 2}, 10'000};
  BOOST_TEST(slice.core[1].begin == 2U);
  BOOST_TEST(slice.core[1].end == 2U);
  BOOST_TEST(slice.halo[1].begin == 2U);
  BOOST_TEST(slice.halo[1].end == 2U);

  const MeasurementView emptyView{MeasurementT<memlayout::const_pointer>{}};
  const std::array<Frame::Index, 10> emptyOffsets{};
  const Frame empty{emptyView, geometry, emptyOffsets};
  BOOST_TEST(empty.size() == 0U);
  const Frame::TimeSlice emptySlice{empty, {.begin = 0, .end = 1}, 20'000};
  for (const auto range : emptySlice.core) {
    BOOST_TEST(range.empty());
  }
  for (const auto range : emptySlice.halo) {
    BOOST_TEST(range.empty());
  }
}

/// A frame given its row count is the frame that reads it from the last
/// offset: what a frame over device columns relies on.
BOOST_AUTO_TEST_CASE(TheRowCountCanBeGivenInsteadOfRead)
{
  const std::array<double, 6> times{0., 10., 5., 10., 10., 20.};
  const auto columns = timeColumns(times);
  const auto geometry = gridGeometry();
  const std::array<Frame::Index, 10> offsets{0, 1, 2, 2, 3, 4, 4, 4, 5, 6};
  const Frame read{measurementView(columns), geometry, offsets};
  const Frame given{measurementView(columns), geometry, offsets, 6U};
  BOOST_TEST(given.size() == read.size());
  BOOST_TEST(given.grid().offsets().data() == offsets.data());
  BOOST_TEST(given.grid().offsets().size() == offsets.size());
  for (std::size_t layer = 0; layer < Frame::kLayerCount; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    BOOST_TEST(given.grid().layerRows(id).begin == read.grid().layerRows(id).begin);
    BOOST_TEST(given.grid().layerRows(id).end == read.grid().layerRows(id).end);
  }
}

} // namespace o2::gtrk
