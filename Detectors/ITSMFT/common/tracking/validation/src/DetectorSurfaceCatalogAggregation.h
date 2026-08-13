// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORSURFACECATALOGAGGREGATION_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORSURFACECATALOGAGGREGATION_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace o2::itsmft::tracking::detail
{

struct GeometryPoint {
  double x{0.};
  double y{0.};
  double z{0.};
};

struct LocalActiveArea {
  double xMin{0.};
  double xMax{0.};
  double zMin{0.};
  double zMax{0.};
};

enum class SurfaceReferenceCoordinate : uint8_t {
  MeanRadius,
  MeanZ
};

enum class SurfaceGeometryAggregationError : uint8_t {
  None,
  SurfaceOutOfRange,
  EmptySurface,
  NonFiniteCoordinate,
  InvalidDerivedGeometry
};

struct AggregatedSurfaceGeometry {
  float referenceCoordinate{0.f};
  float radialMin{0.f};
  float radialMax{0.f};
};

struct SurfaceGeometryAggregationResult {
  std::vector<AggregatedSurfaceGeometry> surfaces{};
  SurfaceGeometryAggregationError error{SurfaceGeometryAggregationError::None};

  bool ok() const noexcept { return error == SurfaceGeometryAggregationError::None; }
};

using SurfaceForChip = std::function<int(size_t)>;
using LocalToGlobal = std::function<GeometryPoint(size_t, const GeometryPoint&)>;

/// Aggregate chip active areas into beam-axis bounds for detector surfaces.
/// ITS and MFT providers supply the layer lookup and local-to-global transform.
SurfaceGeometryAggregationResult aggregateSurfaceGeometry(size_t chipCount,
                                                          size_t surfaceCount,
                                                          const LocalActiveArea& activeArea,
                                                          SurfaceReferenceCoordinate referenceCoordinate,
                                                          const SurfaceForChip& surfaceForChip,
                                                          const LocalToGlobal& localToGlobal);

} // namespace o2::itsmft::tracking::detail

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTORSURFACECATALOGAGGREGATION_H_ */
