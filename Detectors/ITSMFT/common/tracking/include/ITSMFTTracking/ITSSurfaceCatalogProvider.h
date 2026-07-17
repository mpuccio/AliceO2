// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_ITSSURFACECATALOGPROVIDER_H_
#define ALICEO2_ITSMFT_TRACKING_ITSSURFACECATALOGPROVIDER_H_

#ifndef GPUCA_GPUCODE

#include "ITSMFTTracking/DetectorSurfaceCatalogProvider.h"

namespace o2::itsmft::tracking
{

class ITSSurfaceCatalogProvider final : public DetectorSurfaceCatalogProvider
{
 public:
  DetectorSurfaceCatalogResult buildCatalog(const DetectorSurfaceCatalogRequest& request) const final;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE

#endif /* ALICEO2_ITSMFT_TRACKING_ITSSURFACECATALOGPROVIDER_H_ */
