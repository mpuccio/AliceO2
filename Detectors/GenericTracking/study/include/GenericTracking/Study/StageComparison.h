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

/// \file StageComparison.h
/// \brief Two outputs of one stage compared so that a difference is located,
///        not only detected.
///
/// The outputs are canonical: element e's records
/// at offset[e] + k. Two runs of the same stage over the same inputs -- the
/// host and a device, two thread counts -- must agree byte for byte. When
/// they do not, the comparison says which kind of difference it is, element
/// by element:
///
///   another count   the element emitted a different number of records: a cut
///                   decided differently, or an input that differs
///   reordered       the same records in another order: the output's order,
///                   not its content, differs
///   records differ  equal counts, other values: counted per column, with the
///                   largest difference of each arithmetic column's values
///
/// Records are compared as their exact bytes, whatever the column types (a
/// track state has no padding, BoundTrackState.h). compareColumns does the
/// same row by row for column sets without elements. Each keeps the first few
/// differences as (element, ordinal, column, both values). Host code.

#ifndef GTRK_STUDY_STAGECOMPARISON_H_
#define GTRK_STUDY_STAGECOMPARISON_H_

#include "GenericTracking/Tracking/StageOutput.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <vector>

namespace o2::gtrk::study
{

template <template <template <class> class> class S>
struct StageComparison {
  static constexpr std::size_t kColumns = static_cast<std::size_t>(S<memlayout::value>::Column::ColumnCount);
  /// The column of a difference in an element's record count.
  static constexpr std::size_t kCount = kColumns;

  struct Difference {
    std::uint32_t element = 0U;
    std::uint32_t ordinal = 0U; ///< within the element; for kCount, 0
    std::size_t column = 0U;    ///< memlayout's column order, or kCount
    double reference = 0.;
    double other = 0.;
  };

  bool sameElements = true; ///< the same number of elements; nothing else is compared otherwise
  std::size_t referenceElements = 0U;
  std::size_t otherElements = 0U;
  std::size_t referenceRecords = 0U;
  std::size_t otherRecords = 0U;
  std::size_t elementsWithOtherCount = 0U;
  std::size_t elementsReordered = 0U;
  std::size_t recordsDiffering = 0U; ///< in elements of equal count, not merely reordered
  std::array<std::size_t, kColumns> differingByColumn{};
  std::array<double, kColumns> maxAbsDifference{};
  std::vector<Difference> first;

  [[nodiscard]] bool identical() const noexcept { return sameElements && elementsWithOtherCount == 0U && elementsReordered == 0U && recordsDiffering == 0U; }
};

namespace stage_comparison_detail
{

/// A record's columns laid end to end, as their exact bytes: any column type,
/// arithmetic or not (a ClusterList of 17 rows is 68 bytes).
template <template <template <class> class> class S>
struct RecordLayout {
  static constexpr std::size_t kColumns = StageComparison<S>::kColumns;
  static constexpr std::size_t kBytes = S<memlayout::pointer>{}.apply([](auto*... column) { return (sizeof(*column) + ... + std::size_t{0}); });
  static constexpr std::array<std::size_t, kColumns + 1U> kOffsets = S<memlayout::pointer>{}.apply([](auto*... column) {
    std::array<std::size_t, kColumns + 1U> offsets{};
    std::size_t index = 0U;
    ((offsets[index + 1U] = offsets[index] + sizeof(*column), ++index), ...);
    return offsets;
  });
};

template <template <template <class> class> class S>
struct Record {
  using Layout = RecordLayout<S>;
  std::array<std::byte, Layout::kBytes> bytes{};
  std::array<double, Layout::kColumns> values{}; ///< arithmetic columns only; 0 for the others

  [[nodiscard]] bool operator==(const Record& other) const { return bytes == other.bytes; }
  [[nodiscard]] bool sameColumn(const Record& other, std::size_t column) const
  {
    const auto begin = Layout::kOffsets[column];
    return std::memcmp(bytes.data() + begin, other.bytes.data() + begin, Layout::kOffsets[column + 1U] - begin) == 0;
  }
};

template <template <template <class> class> class S, typename View>
[[nodiscard]] Record<S> recordAt(const View& view, std::size_t at)
{
  Record<S> record;
  std::size_t column = 0U;
  view.apply([&](const auto&... columns) {
    const auto take = [&](const auto& value) {
      using Value = std::remove_cvref_t<decltype(value)>;
      std::byte* const to = record.bytes.data() + RecordLayout<S>::kOffsets[column];
      std::memcpy(to, &value, sizeof(Value));
      if constexpr (std::is_arithmetic_v<Value>) {
        record.values[column] = static_cast<double>(value);
      }
      ++column;
    };
    (take(columns[at]), ...);
  });
  return record;
}

/// Count how \a right differs from \a left, column by column, into \a
/// result.
template <template <template <class> class> class S>
void noteRecord(const Record<S>& left, const Record<S>& right, std::uint32_t element, std::uint32_t ordinal, std::size_t keep, StageComparison<S>& result)
{
  ++result.recordsDiffering;
  for (std::size_t column = 0; column < StageComparison<S>::kColumns; ++column) {
    if (left.sameColumn(right, column)) {
      continue;
    }
    ++result.differingByColumn[column];
    const double difference = std::abs(left.values[column] - right.values[column]);
    result.maxAbsDifference[column] = std::max(result.maxAbsDifference[column], difference);
    if (result.first.size() < keep) {
      result.first.push_back({.element = element, .ordinal = ordinal, .column = column, .reference = left.values[column], .other = right.values[column]});
    }
  }
}

} // namespace stage_comparison_detail

/// Compare two column sets row by row: \a other's \a otherRows rows with
/// \a reference's \a referenceRows. Each row is an element of one record, so
/// a difference is (row, 0, column, both values); rows of another number are
/// not compared. The views are any memlayout column views of `S`.
template <template <template <class> class> class S, typename View>
[[nodiscard]] StageComparison<S> compareColumns(const View& reference, std::size_t referenceRows, const View& other, std::size_t otherRows, std::size_t keep = 8U)
{
  StageComparison<S> result;
  result.referenceElements = referenceRows;
  result.otherElements = otherRows;
  result.referenceRecords = referenceRows;
  result.otherRecords = otherRows;
  if (referenceRows != otherRows) {
    result.sameElements = false;
    return result;
  }
  for (std::size_t row = 0; row < referenceRows; ++row) {
    const auto left = stage_comparison_detail::recordAt<S>(reference, row);
    const auto right = stage_comparison_detail::recordAt<S>(other, row);
    if (!(left == right)) {
      stage_comparison_detail::noteRecord(left, right, static_cast<std::uint32_t>(row), 0U, keep, result);
    }
  }
  return result;
}

/// Compare \a other with \a reference, element by element, keeping the first
/// \a keep differences.
template <template <template <class> class> class S>
[[nodiscard]] StageComparison<S> compareStageOutputs(const StageOutput<S>& reference, const StageOutput<S>& other, std::size_t keep = 8U)
{
  using Comparison = StageComparison<S>;
  using Record = stage_comparison_detail::Record<S>;
  Comparison result;
  result.referenceElements = reference.elementCount();
  result.otherElements = other.elementCount();
  result.referenceRecords = reference.size();
  result.otherRecords = other.size();
  if (reference.elementCount() != other.elementCount()) {
    result.sameElements = false;
    return result;
  }
  const auto referenceOffsets = reference.offsets();
  const auto otherOffsets = other.offsets();
  const auto referenceView = reference.view();
  const auto otherView = other.view();
  std::vector<Record> left;
  std::vector<Record> right;
  for (std::uint32_t element = 0; element < reference.elementCount(); ++element) {
    const auto referenceBegin = referenceOffsets[element];
    const auto otherBegin = otherOffsets[element];
    const auto count = referenceOffsets[element + 1U] - referenceBegin;
    const auto otherCount = otherOffsets[element + 1U] - otherBegin;
    if (count != otherCount) {
      ++result.elementsWithOtherCount;
      if (result.first.size() < keep) {
        result.first.push_back({.element = element, .ordinal = 0U, .column = Comparison::kCount, .reference = static_cast<double>(count), .other = static_cast<double>(otherCount)});
      }
      continue;
    }
    left.clear();
    right.clear();
    bool equal = true;
    for (std::uint32_t ordinal = 0; ordinal < count; ++ordinal) {
      left.push_back(stage_comparison_detail::recordAt<S>(referenceView, referenceBegin + ordinal));
      right.push_back(stage_comparison_detail::recordAt<S>(otherView, otherBegin + ordinal));
      equal = equal && left.back() == right.back();
    }
    if (equal) {
      continue;
    }
    auto sortedLeft = left;
    auto sortedRight = right;
    std::ranges::sort(sortedLeft, {}, &Record::bytes);
    std::ranges::sort(sortedRight, {}, &Record::bytes);
    if (sortedLeft == sortedRight) {
      ++result.elementsReordered;
      continue;
    }
    for (std::uint32_t ordinal = 0; ordinal < count; ++ordinal) {
      if (!(left[ordinal] == right[ordinal])) {
        stage_comparison_detail::noteRecord(left[ordinal], right[ordinal], element, ordinal, keep, result);
      }
    }
  }
  return result;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_STAGECOMPARISON_H_
