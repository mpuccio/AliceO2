// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_MULTISOURCELOADING_H_
#define ALICEO2_ITSMFT_TRACKING_MULTISOURCELOADING_H_

#include <cstdint>
#include <limits>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/MultiSourceFrame.h"

namespace o2::itsmft::tracking
{

enum class MultiSourceLoadError : uint8_t {
  None,
  NonDenseSourceIds,
  DuplicateSourceId,
  UnsupportedDetector,
  MissingDecoder,
  InvalidROFRange,
  InvalidLayerMapping,
  DetectorSurfaceMismatch,
  TimingError
};

struct LoadSourcesResult {
  MultiSourceLoadError error{MultiSourceLoadError::None};
  ClusterSourceId source{};
  uint32_t rof{std::numeric_limits<uint32_t>::max()};
  uint32_t clusterIndex{std::numeric_limits<uint32_t>::max()};

  bool ok() const noexcept { return error == MultiSourceLoadError::None; }
};

// Loads one or more detector-qualified sources into `frame`. Geometry
// lookup, pattern consumption, covariance calculation and systematic-error
// application occur exactly once per cluster. Transactional: on any
// validation failure `frame` is left completely unmodified (whatever it
// held before the call) and the returned result identifies the first
// offending source/ROF/cluster. `origin` is the single explicit
// InteractionRecord chosen for the loaded frame; every source ROF is
// converted to a signed TF-relative BC interval against it.
LoadSourcesResult loadSources(MultiSourceFrame& frame,
                              const DetectorLayoutView& layout,
                              gsl::span<const ClusterSourceInput> sources,
                              const o2::InteractionRecord& origin);

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_MULTISOURCELOADING_H_ */
