// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "ITSMFTTracking/LegacyTimeFrameBridge.h"

#include "ITSMFTTracking/ClusterSource.h"

namespace o2::itsmft::tracking::bridge
{

DetectorLayout makeSingleDetectorLayout(o2::detectors::DetID::ID detector, uint16_t nLayers)
{
  SparseTrackingTopology topology{nLayers};
  topology.finalize();
  std::vector<SurfaceDescriptor> surfaces;
  surfaces.reserve(nLayers);
  const auto kind = (detector == o2::detectors::DetID::MFT) ? SurfaceKind::Disk : SurfaceKind::Cylinder;
  for (uint16_t i = 0; i < nLayers; ++i) {
    surfaces.push_back(SurfaceDescriptor{SurfaceId{i}, i, static_cast<uint8_t>(detector), kind});
  }
  return DetectorLayout{std::move(surfaces), std::move(topology)};
}

std::vector<SurfaceId> identityLayerToSurface(uint16_t nLayers)
{
  std::vector<SurfaceId> mapping;
  mapping.reserve(nLayers);
  for (uint16_t i = 0; i < nLayers; ++i) {
    mapping.push_back(SurfaceId{i});
  }
  return mapping;
}

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
                                   ROFTimingConfig timing,
                                   bool applySysErrors)
{
  ClusterSourceInput src;
  src.id = sourceId;
  src.detector = detector;
  src.clusters = clusters;
  src.patterns = patterns;
  src.rofs = rofs;
  src.dictionary = dictionary;
  src.labels = labels;
  src.layerToSurface = layerToSurface;
  src.timing = timing;
  src.decoder = &decoder;
  src.applySysErrors = applySysErrors;
  return loadSources(frame, layout, gsl::span<const ClusterSourceInput>(&src, 1), origin);
}

} // namespace o2::itsmft::tracking::bridge
