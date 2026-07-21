// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACESPEC_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACESPEC_H_

#include <array>
#include <cstddef>
#include <limits>
#include <tuple>
#include <type_traits>

#include "ITSMFTTracking/StaticSurfaceDescriptor.h"

namespace o2::itsmft::tracking
{

namespace detail
{
template <typename T>
struct IsStaticSurfaceArray : std::false_type {
};

template <std::size_t N>
struct IsStaticSurfaceArray<std::array<StaticSurfaceDescriptor, N>> : std::true_type {
};

constexpr bool isFinite(float value) noexcept
{
  constexpr auto maximum = std::numeric_limits<float>::max();
  return value == value && value >= -maximum && value <= maximum;
}

constexpr bool isEnabled(SurfaceKind kind) noexcept
{
  return kind == SurfaceKind::Cylinder || kind == SurfaceKind::Disk;
}

constexpr bool isEnabled(SurfaceIndexingFamily family) noexcept
{
  return family == SurfaceIndexingFamily::CylindricalPhiZ || family == SurfaceIndexingFamily::CartesianXY;
}

template <typename Spec>
consteval bool hasStaticSurfaceArray()
{
  if constexpr (!requires { Spec::surfaces; }) {
    return false;
  } else {
    using Array = std::remove_cv_t<decltype(Spec::surfaces)>;
    if constexpr (!IsStaticSurfaceArray<Array>::value) {
      return false;
    } else {
      // Being usable as a non-type template argument proves static lifetime
      // and constant initialization of the canonical array.
      return requires { typename std::integral_constant<const StaticSurfaceDescriptor*, Spec::surfaces.data()>; };
    }
  }
}
} // namespace detail

template <typename Spec>
concept SurfaceSpec = detail::hasStaticSurfaceArray<Spec>();

template <SurfaceSpec Spec>
inline constexpr std::size_t SurfaceCount = std::tuple_size_v<std::remove_cv_t<decltype(Spec::surfaces)>>;

template <typename Spec>
consteval bool validateSurfaceSpec()
{
  if constexpr (!SurfaceSpec<Spec>) {
    return false;
  } else {
    constexpr auto count = SurfaceCount<Spec>;
    if constexpr (count > MaxLayoutSurfaces) {
      return false;
    }

    for (std::size_t i = 0; i < count; ++i) {
      const auto& surface = Spec::surfaces[i];
      if (!surface.id.isValid() || surface.id.value() != i || !detail::isEnabled(surface.kind) ||
          !detail::isFinite(surface.nominalReferenceCoordinate) || !detail::isEnabled(surface.indexingFamily)) {
        return false;
      }

      const auto& acceptance = surface.nominalTrackingAcceptance;
      const bool acceptanceMatches =
        (surface.kind == SurfaceKind::Cylinder && acceptance.kind == SurfaceAcceptanceKind::CylinderZ) ||
        (surface.kind == SurfaceKind::Disk && acceptance.kind == SurfaceAcceptanceKind::DiskRadius);
      if (!acceptanceMatches || !detail::isFinite(acceptance.min) || !detail::isFinite(acceptance.max) ||
          acceptance.min > acceptance.max || (surface.kind == SurfaceKind::Disk && acceptance.min < 0.f)) {
        return false;
      }

      if (!detail::isFinite(surface.material.xOverX0) || surface.material.xOverX0 <= 0.f) {
        return false;
      }

      for (std::size_t other = i + 1; other < count; ++other) {
        if (surface.identity == Spec::surfaces[other].identity) {
          return false;
        }
      }
    }

    // Detector-local indices form an independent dense [0, N) range for
    // every arbitrary detector ID represented in the catalogue.
    for (std::size_t i = 0; i < count; ++i) {
      const auto detectorId = Spec::surfaces[i].identity.detectorId;
      std::size_t detectorCount = 0;
      for (const auto& surface : Spec::surfaces) {
        detectorCount += surface.identity.detectorId == detectorId;
      }
      for (std::size_t expected = 0; expected < detectorCount; ++expected) {
        bool found = false;
        for (const auto& surface : Spec::surfaces) {
          found = found || (surface.identity.detectorId == detectorId &&
                            surface.identity.detectorSurfaceIndex == expected);
        }
        if (!found) {
          return false;
        }
      }
    }
    return true;
  }
}

template <SurfaceSpec A, SurfaceSpec B>
inline constexpr bool SurfaceSpecsCanBeConcatenated = SurfaceCount<A> + SurfaceCount<B> <= MaxLayoutSurfaces;

namespace detail
{
template <SurfaceSpec A, SurfaceSpec B>
consteval auto concatenateAndRebase()
{
  std::array<StaticSurfaceDescriptor, SurfaceCount<A> + SurfaceCount<B>> result{};
  std::size_t output = 0;
  for (const auto& surface : A::surfaces) {
    result[output] = surface;
    result[output].id = SurfaceId{static_cast<uint16_t>(output)};
    ++output;
  }
  for (const auto& surface : B::surfaces) {
    result[output] = surface;
    result[output].id = SurfaceId{static_cast<uint16_t>(output)};
    ++output;
  }
  return result;
}
} // namespace detail

template <SurfaceSpec A, SurfaceSpec B>
struct ConcatenatedSurfaceSpec {
  static_assert(SurfaceSpecsCanBeConcatenated<A, B>, "a SurfaceSpec catalogue cannot contain more than 32 surfaces");
  inline static constexpr auto surfaces = detail::concatenateAndRebase<A, B>();
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACESPEC_H_ */
