// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_DETECTORSURFACECATALOGPROVIDER_H_
#define ALICEO2_ITSMFT_TRACKING_DETECTORSURFACECATALOGPROVIDER_H_

#include <cstdint>

#ifndef GPUCA_GPUCODE
#include <vector>

#include "ITSMFTTracking/SurfaceDescriptor.h"

namespace o2::itsmft::tracking
{

enum class DetectorSurfaceCatalogError : uint8_t {
  None,
  Unavailable,
  GeometryFailure
};

struct DetectorSurfaceCatalogResult {
  std::vector<SurfaceDescriptor> catalog{};
  DetectorSurfaceCatalogError error{DetectorSurfaceCatalogError::None};

  bool ok() const noexcept { return error == DetectorSurfaceCatalogError::None; }
};

/// Host geometry boundary used to construct the canonical, complete surface
/// catalog. Detector geometry types deliberately do not cross this interface.
class DetectorSurfaceCatalogProvider
{
 public:
  virtual ~DetectorSurfaceCatalogProvider() = default;
  virtual DetectorSurfaceCatalogResult buildCatalog() const = 0;
};

} // namespace o2::itsmft::tracking
#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_DETECTORSURFACECATALOGPROVIDER_H_ */
