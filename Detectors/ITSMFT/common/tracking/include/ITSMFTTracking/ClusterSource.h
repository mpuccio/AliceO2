// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_CLUSTERSOURCE_H_
#define ALICEO2_ITSMFT_TRACKING_CLUSTERSOURCE_H_

#include <gsl/gsl>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/SurfaceId.h"
#include "ITSMFTTracking/SurfaceTiming.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

namespace o2::itsmft::tracking
{

// Non-owning, host-side description of one detector-qualified input stream
// (Architecture.md section 5, "cluster source"). Source IDs must be unique,
// dense, TimeFrame-local and stable for the loaded frame; the all-ones value
// is reserved as invalid (see ClusterSourceId). Multiple sources from the
// same detector are legal, as are multiple sources targeting the same
// surface. The source retains ownership of all detector-specific inputs;
// loadSources() copies out only normalized measurements and timing.
struct ClusterSourceInput {
  ClusterSourceId id{};
  o2::detectors::DetID::ID detector{o2::detectors::DetID::ITS};
  gsl::span<const o2::itsmft::CompClusterExt> clusters{};
  gsl::span<const unsigned char> patterns{};
  gsl::span<const o2::itsmft::ROFRecord> rofs{};
  const o2::itsmft::TopologyDictionary* dictionary{nullptr};
  const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels{nullptr};
  // Detector-local layer (as discovered by the decoder) -> global SurfaceId.
  gsl::span<const SurfaceId> layerToSurface{};
  ROFTimingConfig timing{};
  const ClusterDecoder* decoder{nullptr};
  bool applySysErrors{true};
};

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_CLUSTERSOURCE_H_ */
