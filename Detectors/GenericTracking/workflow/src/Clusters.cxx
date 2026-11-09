// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Alice/Placement.h"

#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITSMFTTracking/IOUtils.h"

#include <gsl/span>

#include <cstdint>

namespace o2::gtrk::alice
{

std::vector<LocalCluster> localClusters(std::span<const o2::itsmft::CompClusterExt> clusters, std::span<const unsigned char> patterns, const o2::itsmft::TopologyDictionary& dictionary)
{
  std::vector<LocalCluster> result;
  result.reserve(clusters.size());
  const gsl::span<const unsigned char> stream{patterns.data(), patterns.size()};
  auto pattern = stream.begin();
  for (const auto& cluster : clusters) {
    const auto read = pattern;
    const auto decoded = o2::itsmft::ioutils::extractClusterData<float>(cluster, pattern, &dictionary);
    LocalCluster local{
      .chip = cluster.getSensorID(),
      .x = decoded.coordinates.X(),
      .z = decoded.coordinates.Z(),
      .sigma2X = decoded.sig2Row,
      .sigma2Z = decoded.sig2Col,
      .nPixels = decoded.nPixels,
    };
    const auto spans = [&local](const o2::itsmft::ClusterPattern& patt) {
      local.rowSpan = static_cast<std::uint8_t>(patt.getRowSpan());
      local.colSpan = static_cast<std::uint8_t>(patt.getColumnSpan());
    };
    // The pattern's spans: a group's and an invalid ID's from the stream,
    // read again where the decode read it; a common topology's from the
    // dictionary.
    if (pattern != read) {
      auto again = read;
      spans(o2::itsmft::ClusterPattern{again});
    } else {
      spans(dictionary.getPattern(cluster.getPatternID()));
    }
    result.push_back(local);
  }
  return result;
}

std::vector<Cluster> toClusters(std::span<const LocalCluster> local, std::span<const o2::itsmft::ROFRecord> rofs, ReadoutFrames frames, std::uint32_t tfFirstOrbit) { return detail::toClusters(local, rofs, frames, tfFirstOrbit, kActiveCentreX); }

ReadoutFrames readoutFrames(std::span<const o2::itsmft::ROFRecord> rofs) { return detail::readoutFrames(rofs); }

} // namespace o2::gtrk::alice
