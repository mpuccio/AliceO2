// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACELAYOUT_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACELAYOUT_H_

#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

#include <gsl/span>

#include "ITSMFTTracking/SurfaceDescriptor.h"
#include "ITSMFTTracking/SurfaceMask.h"

namespace o2::itsmft::tracking
{

enum class SurfaceLayoutError : uint8_t {
  None,
  EmptyCatalog,
  EmptyOrder,
  TooManySurfaces,
  InvalidSurface,
  DuplicateSurface,
  InvalidComponentBoundary,
  NegativeMaxHoles,
  HoleSurfacesOutsideLayout,
  SeedingSurfacesOutsideLayout
};

struct SurfaceLayoutDefinition {
  std::vector<LayerId> orderedSurfaces;
  // First position of each component. Position zero is always required.
  std::vector<uint16_t> componentOffsets{0};
  int maxHoles{0};
  SurfaceMask holeSurfaces{};
  SurfaceMask seedingSurfaces{};
};

inline SurfaceLayoutDefinition makeSurfaceLayoutChain(gsl::span<const LayerId> orderedSurfaces,
                                                      int maxHoles = 0,
                                                      SurfaceMask holeSurfaces = {},
                                                      SurfaceMask seedingSurfaces = {})
{
  SurfaceLayoutDefinition definition;
  definition.orderedSurfaces.assign(orderedSurfaces.begin(), orderedSurfaces.end());
  definition.maxHoles = maxHoles;
  definition.holeSurfaces = holeSurfaces;
  definition.seedingSurfaces = seedingSurfaces;
  return definition;
}

// Immutable layout configuration. Runtime-expanded topology belongs to a
// future TraversalWorkspace; this type intentionally owns no edges, paths,
// adjacency, schedules, or mutable pass state.
class SurfaceLayout
{
 public:
  SurfaceLayout() = default;
  SurfaceLayout(gsl::span<const SurfaceDescriptor> catalog, SurfaceLayoutDefinition definition)
    : mCatalog{catalog.begin(), catalog.end()}, mOrderedSurfaces{std::move(definition.orderedSurfaces)}, mComponentOffsets{std::move(definition.componentOffsets)}, mMaxHoles{definition.maxHoles}, mHoleSurfaces{definition.holeSurfaces}, mSeedingSurfaces{definition.seedingSurfaces}
  {
    validate();
  }

  bool valid() const noexcept { return mError == SurfaceLayoutError::None; }
  SurfaceLayoutError getError() const noexcept { return mError; }
  gsl::span<const LayerId> getOrderedSurfaces() const noexcept { return mOrderedSurfaces; }
  gsl::span<const uint16_t> getComponentOffsets() const noexcept { return mComponentOffsets; }
  int getMaxHoles() const noexcept { return mMaxHoles; }
  SurfaceMask getHoleSurfaces() const noexcept { return mHoleSurfaces; }
  SurfaceMask getSeedingSurfaces() const noexcept { return mSeedingSurfaces; }
  SurfaceCatalogView getSurfaceCatalog() const noexcept { return {mCatalog.data(), static_cast<uint32_t>(mCatalog.size()), mSurfaceIndexById.data()}; }

  bool sameComponent(uint16_t first, uint16_t second) const noexcept
  {
    if (first >= mOrderedSurfaces.size() || second >= mOrderedSurfaces.size()) {
      return false;
    }
    const auto component = [this](uint16_t position) {
      return std::upper_bound(mComponentOffsets.begin(), mComponentOffsets.end(), position) - mComponentOffsets.begin();
    };
    return component(first) == component(second);
  }

 private:
  void validate() noexcept
  {
    mSurfaceIndexById.fill(0xff);
    if (mCatalog.empty()) {
      mError = SurfaceLayoutError::EmptyCatalog;
      return;
    }
    if (mCatalog.size() > MaxLayoutSurfaces) {
      mError = SurfaceLayoutError::TooManySurfaces;
      return;
    }
    for (uint32_t index = 0; index < mCatalog.size(); ++index) {
      const auto id = mCatalog[index].id;
      if (!id.isValid() || id.value() >= MaxLayoutSurfaces) {
        mError = SurfaceLayoutError::InvalidSurface;
        return;
      }
      if (mSurfaceIndexById[id.value()] != 0xff) {
        mError = SurfaceLayoutError::DuplicateSurface;
        return;
      }
      mSurfaceIndexById[id.value()] = static_cast<uint8_t>(index);
    }
    if (mOrderedSurfaces.empty()) {
      mError = SurfaceLayoutError::EmptyOrder;
      return;
    }
    SurfaceMask layoutSurfaces;
    for (const auto id : mOrderedSurfaces) {
      if (!id.isValid() || id.value() >= MaxLayoutSurfaces || mSurfaceIndexById[id.value()] == 0xff) {
        mError = SurfaceLayoutError::InvalidSurface;
        return;
      }
      if (layoutSurfaces.has(id)) {
        mError = SurfaceLayoutError::DuplicateSurface;
        return;
      }
      layoutSurfaces.set(id);
    }
    if (mComponentOffsets.empty() || mComponentOffsets.front() != 0 || mComponentOffsets.back() >= mOrderedSurfaces.size() ||
        !std::is_sorted(mComponentOffsets.begin(), mComponentOffsets.end()) ||
        std::adjacent_find(mComponentOffsets.begin(), mComponentOffsets.end()) != mComponentOffsets.end()) {
      mError = SurfaceLayoutError::InvalidComponentBoundary;
      return;
    }
    if (mMaxHoles < 0) {
      mError = SurfaceLayoutError::NegativeMaxHoles;
      return;
    }
    if (!mHoleSurfaces.isSubsetOf(layoutSurfaces)) {
      mError = SurfaceLayoutError::HoleSurfacesOutsideLayout;
      return;
    }
    if (!mSeedingSurfaces.isSubsetOf(layoutSurfaces)) {
      mError = SurfaceLayoutError::SeedingSurfacesOutsideLayout;
      return;
    }
    mError = SurfaceLayoutError::None;
  }

  std::vector<SurfaceDescriptor> mCatalog;
  std::vector<LayerId> mOrderedSurfaces;
  std::vector<uint16_t> mComponentOffsets;
  int mMaxHoles{0};
  SurfaceMask mHoleSurfaces{};
  SurfaceMask mSeedingSurfaces{};
  std::array<uint8_t, MaxLayoutSurfaces> mSurfaceIndexById{};
  SurfaceLayoutError mError{SurfaceLayoutError::EmptyCatalog};
};

} // namespace o2::itsmft::tracking

#endif
