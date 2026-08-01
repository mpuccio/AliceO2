// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#ifndef ALICEO2_ITSMFT_TRACKING_MULTISOURCETIMEFRAMELOADER_H_
#define ALICEO2_ITSMFT_TRACKING_MULTISOURCETIMEFRAMELOADER_H_

#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/LegacyTrackerScratch.h"

namespace o2::itsmft::tracking
{

// Dormant Gate 4 B3.2 owner-level bridge. This is loading only: it creates
// neither a mixed-detector topology nor any tracker/cell/road/refit/track
// state. The two source IDs must be dense (ITS=0, MFT=1); their mappings are
// global SurfaceIds in the one supplied catalog.
class MultiSourceTimeFrameLoader
{
 public:
  static LoadSourcesResult loadITSAndMFT(TimeFrame& frame,
                                         LegacyTrackerScratchITS& itsScratch,
                                         LegacyTrackerScratchMFT& mftScratch,
                                         const ClusterSourceInput& itsSource,
                                         const ClusterSourceInput& mftSource,
                                         SurfaceCatalogView catalog,
                                         const o2::InteractionRecord& origin);

  // Full shared-event reset: scratch state first (so no legacy cache can
  // outlive normalized measurements), then the common owner exactly once.
  static void resetITSAndMFTEvent(TimeFrame& frame,
                                  LegacyTrackerScratchITS& itsScratch,
                                  LegacyTrackerScratchMFT& mftScratch) noexcept;

 private:
  template <int NLayers>
  static void prepareStagingScratch(LegacyTrackerScratch<NLayers>& staged, const LegacyTrackerScratch<NLayers>& live);
  template <int NLayers>
  static void commitNormalizedBackfill(LegacyTrackerScratch<NLayers>& live, LegacyTrackerScratch<NLayers>& staged) noexcept;
};

} // namespace o2::itsmft::tracking

#endif
