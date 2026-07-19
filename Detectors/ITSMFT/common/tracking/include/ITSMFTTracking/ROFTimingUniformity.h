// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
///
/// \file ROFTimingUniformity.h
/// \brief Host-only detail: collapse per-layer ROF timing into one source-level ROFTimingConfig
///

#ifndef ALICEO2_ITSMFT_TRACKING_ROFTIMINGUNIFORMITY_H_
#define ALICEO2_ITSMFT_TRACKING_ROFTIMINGUNIFORMITY_H_

// Host loading-boundary detail only (ITSMFTTrackingInterface::loadTimeFrame()
// internals). Not a general public API: DPLAlpideParam-derived LayerTiming
// construction, and every caller of this function, are host-only already.
#ifndef GPUCA_GPUCODE

#include <gsl/gsl>

#include "ITSMFTTracking/SurfaceTiming.h"
#include "ITStracking/ROFLookupTables.h"

namespace o2::itsmft::tracking
{

struct UniformROFTimingResult {
  ROFTimingConfig config{};
  bool uniform{false};
};

// TimeFrame::loadNormalizedSource() takes one ROFTimingConfig per source,
// but DPLAlpideParam supports genuine per-layer staggering (independent
// roFrameLayerLengthInBC/roFrameLayerBiasInBC/roFrameLayerDelayInBC
// overrides per layer; both ITS and MFT default every override to zero, so
// all layers resolve to the shared global value out of the box, but this is
// a real, supported production knob, not a theoretical one). Collapsing to
// one source-level config without checking would silently discard a
// legitimately divergent configuration. `uniform` is false whenever any
// layer's length/delay/bias/addTimeErr disagrees with layer 0's; `config` is
// only meaningful when `uniform` is true.
inline UniformROFTimingResult deriveUniformROFTimingConfig(gsl::span<const o2::its::LayerTiming> perLayer) noexcept
{
  if (perLayer.empty()) {
    return {};
  }
  const auto& ref = perLayer[0];
  for (const auto& lt : perLayer) {
    if (lt.mROFLength != ref.mROFLength || lt.mROFDelay != ref.mROFDelay ||
        lt.mROFBias != ref.mROFBias || lt.mROFAddTimeErr != ref.mROFAddTimeErr) {
      return {};
    }
  }
  return {ROFTimingConfig{static_cast<TFBC>(ref.mROFLength), static_cast<TFBC>(ref.mROFDelay),
                          static_cast<TFBC>(ref.mROFBias), static_cast<TFBC>(ref.mROFAddTimeErr)},
          true};
}

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE
#endif /* ALICEO2_ITSMFT_TRACKING_ROFTIMINGUNIFORMITY_H_ */
