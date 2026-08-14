// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3, version 3, copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEGRAPH_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEGRAPH_H_

#include <cstddef>
#include <cstdint>
#include <type_traits>

#ifndef GPUCA_GPUCODE
#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include <gsl/gsl>
#endif

#include "GPUCommonDef.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/IdTypes.h"
#include "ITSMFTTracking/SurfaceMask.h"

namespace o2::itsmft::tracking
{

struct SurfaceLink {
  SurfaceId from{};
  SurfaceId to{};
  SurfaceMask skippedSurfaces{};
  uint16_t flags{0};
};

struct SurfaceCellTopology {
  LinkId firstLink{};
  LinkId secondLink{};
  SurfaceMask hitSurfaces{};
};

struct TopologyRange {
  uint32_t firstEntry{0};
  uint32_t entries{0};

  GPUhdi() constexpr uint32_t getFirstEntry() const noexcept { return firstEntry; }
  GPUhdi() constexpr uint32_t getEntries() const noexcept { return entries; }
  GPUhdi() constexpr uint32_t getEntriesBound() const noexcept { return firstEntry + entries; }
};

// Device-facing graph representation containing the surface catalog,
// traversal order, sparse adjacency, and seed mask in one POD view.
// All pointers are borrowed from one immutable SurfaceGraph owner.
struct SurfaceGraphView {
  const SurfaceDescriptor* surfaces{nullptr};
  uint32_t nSurfaces{0};
  const SurfaceId* orderedSurfaces{nullptr};
  uint32_t nOrderedSurfaces{0};
  SurfaceMask cylinderSurfaces{};
  SurfaceMask diskSurfaces{};
  const SurfaceLink* links{nullptr};
  const SurfaceCellTopology* cells{nullptr};
  const uint32_t* cellsByFirstLinkOffsets{nullptr};
  const CellTopologyId* cellsByFirstLink{nullptr};
  SurfaceMask seedingSurfaces{};
  uint32_t nLinks{0};
  uint32_t nCells{0};
  const uint8_t* surfaceIndicesById{nullptr};

  GPUhdi() uint32_t getSurfaceIndex(SurfaceId id) const
  {
    if (!id.isValid() || id.value() >= MaxLayoutSurfaces) {
      return nSurfaces;
    }
    if (surfaceIndicesById != nullptr) {
      return surfaceIndicesById[id.value()];
    }
    return id.value();
  }
  GPUhdi() const SurfaceDescriptor& getSurface(SurfaceId id) const { return surfaces[getSurfaceIndex(id)]; }
  GPUhdi() SurfaceCatalogView getSurfaceCatalogView() const noexcept { return SurfaceCatalogView{surfaces, nSurfaces, surfaceIndicesById}; }
  GPUhdi() const SurfaceLink& getLink(LinkId id) const { return links[id.value()]; }
  GPUhdi() const SurfaceCellTopology& getCell(CellTopologyId id) const { return cells[id.value()]; }
  GPUhdi() TopologyRange getCellsStartingWithLink(LinkId link) const
  {
    const uint32_t id = link.value();
    return TopologyRange{cellsByFirstLinkOffsets[id], cellsByFirstLinkOffsets[id + 1] - cellsByFirstLinkOffsets[id]};
  }
};

static_assert(std::is_standard_layout_v<SurfaceLink> && std::is_trivially_copyable_v<SurfaceLink>);
static_assert(std::is_standard_layout_v<SurfaceCellTopology> && std::is_trivially_copyable_v<SurfaceCellTopology>);
static_assert(std::is_standard_layout_v<TopologyRange> && std::is_trivially_copyable_v<TopologyRange>);
static_assert(std::is_standard_layout_v<SurfaceGraphView> && std::is_trivially_copyable_v<SurfaceGraphView>);
static_assert(sizeof(SurfaceLink) == 12);
static_assert(sizeof(SurfaceCellTopology) == 8);
static_assert(offsetof(SurfaceLink, from) == 0);
static_assert(offsetof(SurfaceLink, to) == 2);
static_assert(offsetof(SurfaceLink, skippedSurfaces) == 4);
static_assert(offsetof(SurfaceLink, flags) == 8);

#ifndef GPUCA_GPUCODE

inline std::pair<SurfaceMask, SurfaceMask> computeSurfaceKindMasks(gsl::span<const SurfaceDescriptor> surfaces) noexcept
{
  SurfaceMask cylinderSurfaces;
  SurfaceMask diskSurfaces;
  for (const auto& surface : surfaces) {
    (surface.kind == SurfaceKind::Cylinder ? cylinderSurfaces : diskSurfaces).set(surface.id);
  }
  return {cylinderSurfaces, diskSurfaces};
}

enum class SurfaceGraphError : uint8_t {
  None,
  TooManySurfaces,
  NonDenseSurfaceIds,
  SurfaceCountMismatch,
  NotFinalized,
  MixedSurfaceLink,
  SurfaceKindMismatch
};

enum class SurfaceGraphTopologyError : uint8_t {
  None,
  InvalidSurfaceCount,
  InvalidSurface,
  SelfLink,
  DuplicateLink,
  TooManyLinks,
  InvalidLink,
  DisconnectedLinks,
  RepeatedSurface,
  DuplicateCell,
  TooManyCells,
  AlreadyFinalized,
  NotFinalized
};

// Owns one complete surface graph that is immutable after finalization.
class SurfaceGraph
{
 public:
  explicit SurfaceGraph(uint32_t surfaceCount, SurfaceMask seedingSurfaces = {})
    : mSurfaces(surfaceCount), mSurfaceCount{surfaceCount}, mSeedingSurfaces{seedingSurfaces}
  {
    for (uint32_t i = 0; i < surfaceCount; ++i) {
      mSurfaces[i].id = SurfaceId{static_cast<uint16_t>(i)};
    }
    initializeSurfaceIndices();
    if (surfaceCount > MaxLayoutSurfaces) {
      mTopologyError = SurfaceGraphTopologyError::InvalidSurfaceCount;
    } else if (!seedingSurfaces.isSubsetOf(activeSurfaceMask())) {
      mTopologyError = SurfaceGraphTopologyError::InvalidSurface;
    }
  }

  explicit SurfaceGraph(gsl::span<const SurfaceDescriptor> surfaces, SurfaceMask seedingSurfaces = {})
    : mSurfaces{surfaces.begin(), surfaces.end()}, mSurfaceCount{static_cast<uint32_t>(surfaces.size())}, mSeedingSurfaces{seedingSurfaces}
  {
    initializeSurfaceIndices();
    if (mSurfaceCount > MaxLayoutSurfaces || !seedingSurfaces.isSubsetOf(activeSurfaceMask())) {
      mTopologyError = SurfaceGraphTopologyError::InvalidSurface;
    }
  }

  SurfaceGraph(gsl::span<const SurfaceDescriptor> surfaces, SurfaceGraph topology)
    : mSurfaces{surfaces.begin(), surfaces.end()},
      mSurfaceCount{static_cast<uint32_t>(surfaces.size())},
      mOrderedSurfaces{std::move(topology.mOrderedSurfaces)},
      mSeedingSurfaces{topology.mSeedingSurfaces},
      mLinks{std::move(topology.mLinks)},
      mCells{std::move(topology.mCells)},
      mCellsByFirstLinkOffsets{std::move(topology.mCellsByFirstLinkOffsets)},
      mCellsByFirstLink{std::move(topology.mCellsByFirstLink)},
      mTopologyError{topology.mTopologyError},
      mFinalized{topology.mFinalized}
  {
    initializeSurfaceIndices();
    validate();
  }

  SurfaceGraph(SurfaceGraph&&) noexcept = default;
  SurfaceGraph& operator=(SurfaceGraph&&) noexcept = default;
  SurfaceGraph(const SurfaceGraph&) = default;
  SurfaceGraph& operator=(const SurfaceGraph&) = default;

  void setOrderedSurfaces(std::vector<SurfaceId> ordered) { mOrderedSurfaces = std::move(ordered); }

  LinkId addLink(SurfaceLink link)
  {
    if (!canModify()) {
      return LinkId::invalid();
    }
    if (!isSurfaceInGraph(link.from) || !isSurfaceInGraph(link.to) ||
        !link.skippedSurfaces.isSubsetOf(activeSurfaceMask()) ||
        link.skippedSurfaces.has(link.from) || link.skippedSurfaces.has(link.to)) {
      mTopologyError = SurfaceGraphTopologyError::InvalidSurface;
      return LinkId::invalid();
    }
    if (link.from == link.to) {
      mTopologyError = SurfaceGraphTopologyError::SelfLink;
      return LinkId::invalid();
    }
    const auto duplicate = std::find_if(mLinks.begin(), mLinks.end(), [&](const auto& existing) {
      return existing.from == link.from && existing.to == link.to;
    });
    if (duplicate != mLinks.end()) {
      mTopologyError = SurfaceGraphTopologyError::DuplicateLink;
      return LinkId::invalid();
    }
    if (mLinks.size() >= MaxLayoutLinks) {
      mTopologyError = SurfaceGraphTopologyError::TooManyLinks;
      return LinkId::invalid();
    }
    const auto id = LinkId{static_cast<uint16_t>(mLinks.size())};
    mLinks.push_back(link);
    return id;
  }

  CellTopologyId addCell(LinkId first, LinkId second)
  {
    if (!canModify()) {
      return CellTopologyId::invalid();
    }
    if (!isLinkInGraph(first) || !isLinkInGraph(second)) {
      mTopologyError = SurfaceGraphTopologyError::InvalidLink;
      return CellTopologyId::invalid();
    }
    const auto& firstLink = mLinks[first.value()];
    const auto& secondLink = mLinks[second.value()];
    if (firstLink.to != secondLink.from) {
      mTopologyError = SurfaceGraphTopologyError::DisconnectedLinks;
      return CellTopologyId::invalid();
    }
    if (firstLink.from == secondLink.to) {
      mTopologyError = SurfaceGraphTopologyError::RepeatedSurface;
      return CellTopologyId::invalid();
    }
    const auto duplicate = std::find_if(mCells.begin(), mCells.end(), [&](const auto& existing) {
      return existing.firstLink == first && existing.secondLink == second;
    });
    if (duplicate != mCells.end()) {
      mTopologyError = SurfaceGraphTopologyError::DuplicateCell;
      return CellTopologyId::invalid();
    }
    if (mCells.size() >= MaxLayoutCellTopologies) {
      mTopologyError = SurfaceGraphTopologyError::TooManyCells;
      return CellTopologyId::invalid();
    }
    SurfaceMask hits;
    hits.set(firstLink.from);
    hits.set(firstLink.to);
    hits.set(secondLink.to);
    const auto id = CellTopologyId{static_cast<uint16_t>(mCells.size())};
    mCells.push_back(SurfaceCellTopology{first, second, hits});
    return id;
  }

  bool finalize()
  {
    if (mFinalized) {
      return true;
    }
    if (mTopologyError != SurfaceGraphTopologyError::None) {
      return false;
    }
    if (mOrderedSurfaces.empty()) {
      mOrderedSurfaces.resize(mSurfaceCount);
      for (uint32_t i = 0; i < mSurfaceCount; ++i) {
        mOrderedSurfaces[i] = SurfaceId{static_cast<uint16_t>(i)};
      }
    }
    mCellsByFirstLinkOffsets.assign(mLinks.size() + 1, 0);
    for (const auto& cell : mCells) {
      ++mCellsByFirstLinkOffsets[cell.firstLink.value() + 1];
    }
    for (size_t i = 1; i < mCellsByFirstLinkOffsets.size(); ++i) {
      mCellsByFirstLinkOffsets[i] += mCellsByFirstLinkOffsets[i - 1];
    }
    mCellsByFirstLink.resize(mCells.size());
    auto cursor = mCellsByFirstLinkOffsets;
    for (uint32_t cell = 0; cell < mCells.size(); ++cell) {
      const auto link = mCells[cell].firstLink.value();
      mCellsByFirstLink[cursor[link]++] = CellTopologyId{static_cast<uint16_t>(cell)};
    }
    mFinalized = true;
    validate();
    return valid();
  }

  bool valid() const noexcept { return mError == SurfaceGraphError::None && mTopologyError == SurfaceGraphTopologyError::None && mFinalized; }
  SurfaceGraphError getError() const noexcept { return mError; }
  SurfaceGraphTopologyError getTopologyError() const noexcept { return mTopologyError; }
  bool isFinalized() const noexcept { return mFinalized; }
  uint32_t getSurfaceCount() const noexcept { return mSurfaceCount; }
  const std::vector<SurfaceDescriptor>& getSurfaces() const noexcept { return mSurfaces; }
  const std::vector<SurfaceId>& getOrderedSurfaces() const noexcept { return mOrderedSurfaces; }
  const auto& getLinks() const noexcept { return mLinks; }
  const auto& getCells() const noexcept { return mCells; }
  const auto& getCellsByFirstLinkOffsets() const noexcept { return mCellsByFirstLinkOffsets; }
  const auto& getCellsByFirstLink() const noexcept { return mCellsByFirstLink; }
  SurfaceCatalogView getSurfaceCatalog() const noexcept { return SurfaceCatalogView{mSurfaces.data(), mSurfaceCount, mSurfaceIndicesById.data()}; }
  SurfaceGraphView getView() const noexcept
  {
    if (!valid()) {
      return {};
    }
    return SurfaceGraphView{mSurfaces.data(), mSurfaceCount, mOrderedSurfaces.data(), static_cast<uint32_t>(mOrderedSurfaces.size()),
                            mCylinderSurfaces, mDiskSurfaces, mLinks.data(), mCells.data(),
                            mCellsByFirstLinkOffsets.data(), mCellsByFirstLink.data(), mSeedingSurfaces,
                            static_cast<uint32_t>(mLinks.size()), static_cast<uint32_t>(mCells.size()), mSurfaceIndicesById.data()};
  }

 private:
  bool canModify()
  {
    if (mFinalized) {
      mTopologyError = SurfaceGraphTopologyError::AlreadyFinalized;
      return false;
    }
    return mTopologyError == SurfaceGraphTopologyError::None;
  }
  bool isSurfaceInGraph(SurfaceId id) const noexcept
  {
    return id.isValid() && id.value() < MaxLayoutSurfaces && mSurfaceIndicesById[id.value()] != InvalidSurfaceIndex;
  }
  bool isLinkInGraph(LinkId id) const noexcept { return id.isValid() && id.value() < mLinks.size(); }
  SurfaceMask activeSurfaceMask() const noexcept
  {
    SurfaceMask result;
    for (uint32_t i = 0; i < mSurfaceCount; ++i) {
      result.set(mSurfaces[i].id);
    }
    return result;
  }
  static constexpr uint8_t InvalidSurfaceIndex = 0xff;
  void initializeSurfaceIndices() noexcept
  {
    mSurfaceIndicesById.fill(InvalidSurfaceIndex);
    for (uint32_t i = 0; i < mSurfaces.size(); ++i) {
      const auto id = mSurfaces[i].id.value();
      if (id < MaxLayoutSurfaces && mSurfaceIndicesById[id] == InvalidSurfaceIndex) {
        mSurfaceIndicesById[id] = static_cast<uint8_t>(i);
      }
    }
  }
  void validate()
  {
    if (mSurfaceCount > MaxLayoutSurfaces) {
      mError = SurfaceGraphError::TooManySurfaces;
      return;
    }
    if (mSurfaces.size() != mSurfaceCount) {
      mError = SurfaceGraphError::SurfaceCountMismatch;
      return;
    }
    initializeSurfaceIndices();
    for (uint32_t i = 0; i < mSurfaces.size(); ++i) {
      if (!mSurfaces[i].id.isValid() || mSurfaces[i].id.value() >= MaxLayoutSurfaces ||
          mSurfaceIndicesById[mSurfaces[i].id.value()] != i) {
        mError = SurfaceGraphError::NonDenseSurfaceIds;
        return;
      }
    }
    if (!mOrderedSurfaces.empty()) {
      SurfaceMask orderedMask;
      for (const auto id : mOrderedSurfaces) {
        if (!isSurfaceInGraph(id) || orderedMask.has(id)) {
          mError = SurfaceGraphError::NonDenseSurfaceIds;
          return;
        }
        orderedMask.set(id);
      }
    }
    mCylinderSurfaces = {};
    mDiskSurfaces = {};
    for (const auto& surface : mSurfaces) {
      (surface.kind == SurfaceKind::Cylinder ? mCylinderSurfaces : mDiskSurfaces).set(surface.id);
    }
  }

  std::vector<SurfaceDescriptor> mSurfaces;
  uint32_t mSurfaceCount{0};
  std::vector<SurfaceId> mOrderedSurfaces;
  SurfaceMask mSeedingSurfaces{};
  SurfaceMask mCylinderSurfaces{};
  SurfaceMask mDiskSurfaces{};
  std::vector<SurfaceLink> mLinks;
  std::vector<SurfaceCellTopology> mCells;
  std::vector<uint32_t> mCellsByFirstLinkOffsets;
  std::vector<CellTopologyId> mCellsByFirstLink;
  std::array<uint8_t, MaxLayoutSurfaces> mSurfaceIndicesById{};
  SurfaceGraphError mError{SurfaceGraphError::None};
  SurfaceGraphTopologyError mTopologyError{SurfaceGraphTopologyError::None};
  bool mFinalized{false};
};

#endif // GPUCA_GPUCODE

} // namespace o2::itsmft::tracking

#endif
