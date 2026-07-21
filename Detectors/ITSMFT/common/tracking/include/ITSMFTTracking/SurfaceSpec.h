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

template <std::size_t N>
consteval bool validateSurfaceArray(const std::array<StaticSurfaceDescriptor, N>& surfaces)
{
  if constexpr (N > MaxLayoutSurfaces) {
    return false;
  }

  for (std::size_t i = 0; i < N; ++i) {
    const auto& surface = surfaces[i];
    if (!surface.id.isValid() || surface.id.value() != i || !isEnabled(surface.kind) ||
        !isFinite(surface.nominalReferenceCoordinate) ||
        (surface.kind == SurfaceKind::Cylinder && surface.nominalReferenceCoordinate <= 0.f) ||
        !isEnabled(surface.indexingFamily)) {
      return false;
    }

    const auto& acceptance = surface.nominalTrackingAcceptance;
    const bool acceptanceMatches =
      (surface.kind == SurfaceKind::Cylinder && acceptance.kind == SurfaceAcceptanceKind::CylinderZ) ||
      (surface.kind == SurfaceKind::Disk && acceptance.kind == SurfaceAcceptanceKind::DiskRadius);
    if (!acceptanceMatches || !isFinite(acceptance.min) || !isFinite(acceptance.max) || acceptance.min > acceptance.max ||
        (surface.kind == SurfaceKind::Disk && acceptance.min < 0.f)) {
      return false;
    }

    if (!isFinite(surface.material.xOverX0) || surface.material.xOverX0 <= 0.f) {
      return false;
    }

    for (std::size_t other = i + 1; other < N; ++other) {
      if (surface.identity == surfaces[other].identity) {
        return false;
      }
    }
  }

  // Detector-local indices form an independent dense [0, N) range for every
  // arbitrary detector ID represented in the catalogue.
  for (std::size_t i = 0; i < N; ++i) {
    const auto detectorId = surfaces[i].identity.detectorId;
    std::size_t detectorCount = 0;
    for (const auto& surface : surfaces) {
      detectorCount += surface.identity.detectorId == detectorId;
    }
    for (std::size_t expected = 0; expected < detectorCount; ++expected) {
      bool found = false;
      for (const auto& surface : surfaces) {
        found = found ||
                (surface.identity.detectorId == detectorId && surface.identity.detectorSurfaceIndex == expected);
      }
      if (!found) {
        return false;
      }
    }
  }
  return true;
}

template <typename Spec>
consteval bool validateSurfaceSpecDefinition()
{
  return validateSurfaceArray(Spec::surfaces);
}
} // namespace detail

template <typename Spec>
// Structural requirement only: canonical inline-static-array shape and
// lifetime. Catalogue contents are deliberately not accepted by this concept.
concept SurfaceSpecDefinition = detail::hasStaticSurfaceArray<Spec>();

template <typename Spec>
// A consumer-facing SurfaceSpec has both the required definition shape and a
// fully validated catalogue.
concept SurfaceSpec = SurfaceSpecDefinition<Spec> && detail::validateSurfaceSpecDefinition<Spec>();

template <SurfaceSpec Spec>
inline constexpr std::size_t SurfaceCount = std::tuple_size_v<std::remove_cv_t<decltype(Spec::surfaces)>>;

template <typename Spec>
consteval bool validateSurfaceSpec()
{
  if constexpr (!SurfaceSpecDefinition<Spec>) {
    return false;
  } else {
    return detail::validateSurfaceSpecDefinition<Spec>();
  }
}

namespace detail
{
template <SurfaceSpecDefinition A, SurfaceSpecDefinition B>
consteval auto concatenateAndRebase()
{
  constexpr auto countA = std::tuple_size_v<std::remove_cv_t<decltype(A::surfaces)>>;
  constexpr auto countB = std::tuple_size_v<std::remove_cv_t<decltype(B::surfaces)>>;
  std::array<StaticSurfaceDescriptor, countA + countB> result{};
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

template <SurfaceSpecDefinition A, SurfaceSpecDefinition B>
consteval bool surfaceSpecsCanBeConcatenated()
{
  if constexpr (!SurfaceSpec<A> || !SurfaceSpec<B>) {
    return false;
  } else if constexpr (SurfaceCount<A> + SurfaceCount<B> > MaxLayoutSurfaces) {
    return false;
  } else {
    return validateSurfaceArray(concatenateAndRebase<A, B>());
  }
}
} // namespace detail

template <SurfaceSpecDefinition A, SurfaceSpecDefinition B>
inline constexpr bool SurfaceSpecsCanBeConcatenated = detail::surfaceSpecsCanBeConcatenated<A, B>();

template <SurfaceSpec A, SurfaceSpec B>
struct ConcatenatedSurfaceSpec {
  static_assert(SurfaceSpecsCanBeConcatenated<A, B>, "SurfaceSpecs cannot be concatenated into a valid catalogue");
  inline static constexpr auto surfaces = detail::concatenateAndRebase<A, B>();
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_SURFACESPEC_H_ */
