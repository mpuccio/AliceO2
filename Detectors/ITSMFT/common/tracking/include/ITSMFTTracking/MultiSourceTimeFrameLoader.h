// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

// L6: the loader is a non-owning direct component. It stages every source and
// installs one complete event into a configured TimeFrame; it does not track,
// publish, own raw ROFs, or make workflow decisions.

#ifndef ALICEO2_ITSMFT_TRACKING_MULTISOURCETIMEFRAMELOADER_H_
#define ALICEO2_ITSMFT_TRACKING_MULTISOURCETIMEFRAMELOADER_H_

#include <gsl/gsl>

#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TrackingConfigParam.h"

namespace o2::itsmft::tracking
{

class MultiSourceTimeFrameLoader
{
 public:
  // All normalized data and source workspaces are staged locally. A failure
  // leaves the frame's previous event and configuration untouched.
  static LoadSourcesResult load(TimeFrame& frame, gsl::span<const ClusterSourceInput> sources,
                                SurfaceCatalogView catalog, const o2::InteractionRecord& origin);
};

} // namespace o2::itsmft::tracking

#endif
