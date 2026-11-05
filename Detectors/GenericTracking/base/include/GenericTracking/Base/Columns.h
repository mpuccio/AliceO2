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

/// \file Columns.h
/// \brief Owned SoA columns of any memlayout struct, from a bounded pool.
///
/// Every artefact the tracker produces -- tracklets, triplets, and whatever
/// follows -- is declared once as `S<C>` with MEMLAYOUT_APPLY. This owns its
/// columns as bounded vectors and hands out the pointer views kernels write
/// and read. Each operation is written over `S::apply`, so a new artefact
/// needs no code here. Host-side: a device backend owns device allocations
/// and builds the same views over them.

#ifndef GTRK_BASE_COLUMNS_H_
#define GTRK_BASE_COLUMNS_H_

#include "GenericTracking/Base/BoundedAllocator.h"

#include <memlayout/wrapper.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory_resource>

namespace o2::gtrk
{

template <template <template <class> class> class S>
class BoundedColumns
{
 public:
  using Record = S<memlayout::value>;
  using View = memlayout::Wrapper<S, memlayout::pointer, memlayout::Layout::soa>;
  using ConstView = memlayout::Wrapper<S, memlayout::const_pointer, memlayout::Layout::soa>;

  /// Reseat every column on \a resource, releasing what it held.
  void setMemoryResource(std::pmr::memory_resource* resource)
  {
    mColumns.apply([resource](auto&... column) { (deepVectorClear(column, resource), ...); });
  }

  /// Every column to \a size live elements; capacity is kept when shrinking.
  /// All or nothing: if a column's allocation throws, every column goes back
  /// to its size before -- which allocates nothing -- and the exception goes
  /// on. Columns of different sizes would let `size()`, which reads the first,
  /// promise rows a later column does not have.
  void resize(std::size_t size)
  {
    const std::size_t before = this->size();
    try {
      mColumns.apply([size](auto&... column) { (column.resize(size), ...); });
    } catch (...) {
      mColumns.apply([before](auto&... column) { (column.resize(before), ...); });
      throw;
    }
  }

  /// Every column at \a size rows for a pass that writes them all before it
  /// reads one: `resizeBoundedVector` per column -- in its storage where that
  /// holds them, otherwise rebuilt at exactly \a size, the old storage released
  /// first, so that a growing slice neither copies rows it overwrites nor
  /// doubles the capacity. If an allocation
  /// throws, every column is left empty, and the exception goes on.
  void resizeRewritten(std::size_t size)
  {
    try {
      mColumns.apply([size](auto&... column) { (resizeBoundedVector(column, size), ...); });
    } catch (...) {
      mColumns.apply([](auto&... column) { (column.clear(), ...); });
      throw;
    }
  }

  [[nodiscard]] std::size_t size() const
  {
    return mColumns.apply([](const auto& first, const auto&... /*rest*/) { return first.size(); });
  }

  [[nodiscard]] View view()
  {
    return mColumns.apply([](auto&... column) { return View{S<memlayout::pointer>{column.data()...}}; });
  }

  [[nodiscard]] ConstView view() const
  {
    return mColumns.apply([](const auto&... column) { return ConstView{S<memlayout::const_pointer>{column.data()...}}; });
  }

  void write(std::size_t index, const Record& record)
  {
    assert(index < size());
    mColumns.apply([&](auto&... column) { record.apply([&](const auto&... value) { ((column[index] = value), ...); }); });
  }

  [[nodiscard]] Record read(std::size_t index) const
  {
    assert(index < size());
    return mColumns.apply([index](const auto&... column) { return Record{column[index]...}; });
  }

  /// Copy \a count rows starting at \a from into \a destination at \a to.
  void copyTo(std::size_t from, BoundedColumns& destination, std::size_t to, std::size_t count) const
  {
    assert(from + count <= size() && to + count <= destination.size());
    mColumns.apply([&](const auto&... source) {
      destination.mColumns.apply([&](auto&... target) {
        (std::copy_n(source.begin() + static_cast<std::ptrdiff_t>(from), count, target.begin() + static_cast<std::ptrdiff_t>(to)), ...);
      });
    });
  }

 private:
  S<BoundedVector> mColumns;
};

} // namespace o2::gtrk

#endif // GTRK_BASE_COLUMNS_H_
