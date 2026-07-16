// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_LEGACYTIMEFRAMEBRIDGE_H_
#define ALICEO2_ITSMFT_TRACKING_LEGACYTIMEFRAMEBRIDGE_H_

#include <cstdint>
#include <vector>

#include <gsl/gsl>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSMFTTracking/ClusterDecoder.h"
#include "ITSMFTTracking/DetectorLayout.h"
#include "ITSMFTTracking/MultiSourceFrame.h"
#include "ITSMFTTracking/MultiSourceLoading.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

// Gate 1 compatibility bridge (Architecture.md / AgentCoordination.md):
// exposes the same raw cluster streams that the existing single-detector
// loaders -- o2::its::TimeFrame<NLayers>::loadROFrameData
// (Detectors/ITSMFT/ITS/tracking) and o2::mft::ioutils::loadROFrameData /
// o2::mft::ROframe (Detectors/ITSMFT/MFT/tracking) -- already consume,
// through the normalized MultiSourceFrame owner via loadSources(). This is
// purely additive: it does not touch either production loader, any
// workflow, TrackerTraits/CA orchestration, GPU kernels, output conversion,
// tracking parameters or vertexing ownership. It performs no decoding of its
// own -- the caller-supplied ClusterDecoder (e.g. the existing
// GeometryClusterDecoder<DetId> production adapter) remains the single place
// each compact cluster is decoded, exactly as for any other loadSources()
// caller.
namespace o2::itsmft::tracking::bridge
{

// A trivial disconnected single-detector layout: `nLayers` surfaces, one per
// detector-local layer, in layer order, with no transitions. CA topology is
// out of scope for a loading-only bridge -- only DetectorLayout's role in
// loadSources() (surface count/kind/detector validation) is needed here.
DetectorLayout makeSingleDetectorLayout(o2::detectors::DetID::ID detector, uint16_t nLayers);

// Identity detector-local-layer -> global SurfaceId mapping matching the
// surface numbering produced by makeSingleDetectorLayout() for the same
// `nLayers`.
std::vector<SurfaceId> identityLayerToSurface(uint16_t nLayers);

// Loads one legacy single-detector cluster stream -- the same raw
// ROFRecord/CompClusterExt/pattern/dictionary/label inputs the existing
// single-detector loaders already consume -- into `frame` through
// loadSources(), via the caller-supplied `decoder` (production callers pass
// the existing ITSGeometryClusterDecoder/MFTGeometryClusterDecoder; tests may
// pass a host-only decoder, exactly like any other loadSources() caller).
// `applySysErrors` defaults to false because neither
// o2::its::TimeFrame<NLayers>::loadROFrameData nor
// o2::mft::ioutils::loadROFrameData applies systematic errors while loading
// clusters (see AgentCoordination.md decision log / Gate 1 report); passing
// `true` is only meaningful for callers that intentionally want the
// GeometryClusterDecoder sys-error convention instead of legacy parity.
// Transactional and single-decode, like loadSources() itself.
LoadSourcesResult loadLegacySource(MultiSourceFrame& frame,
                                   const DetectorLayoutView& layout,
                                   o2::detectors::DetID::ID detector,
                                   gsl::span<const SurfaceId> layerToSurface,
                                   ClusterSourceId sourceId,
                                   gsl::span<const o2::itsmft::ROFRecord> rofs,
                                   gsl::span<const o2::itsmft::CompClusterExt> clusters,
                                   gsl::span<const unsigned char> patterns,
                                   const o2::itsmft::TopologyDictionary* dictionary,
                                   const o2::dataformats::MCTruthContainer<o2::MCCompLabel>* labels,
                                   const ClusterDecoder& decoder,
                                   const o2::InteractionRecord& origin,
                                   ROFTimingConfig timing = {},
                                   bool applySysErrors = false);

} // namespace o2::itsmft::tracking::bridge

#endif /* ALICEO2_ITSMFT_TRACKING_LEGACYTIMEFRAMEBRIDGE_H_ */
