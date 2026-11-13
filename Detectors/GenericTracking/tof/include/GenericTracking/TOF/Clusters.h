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

#ifndef GTRK_TOF_CLUSTERS_H_
#define GTRK_TOF_CLUSTERS_H_
#include "GenericTracking/TOF/Matcher.h"
#include <span>
#include <vector>
namespace o2::iotof
{
class Cluster;
class GeometryTGeo;
} // namespace o2::iotof
namespace o2::gtrk::tof
{
/// The only bridge to the temporary cluster API. Drops clusters whose pattern
/// is unavailable; preserves indices in the original array for every retained hit.
std::vector<Hit> convertClusters(std::span<const o2::iotof::Cluster> clusters, const o2::iotof::GeometryTGeo& geometry, float sigmaTimeNS);
} // namespace o2::gtrk::tof
#endif
