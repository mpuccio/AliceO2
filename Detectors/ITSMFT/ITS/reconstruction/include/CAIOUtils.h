// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAIOUtils.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CAEVENTLOADER_H_
#define TRACKINGITSU_INCLUDE_CAEVENTLOADER_H_

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "CAEvent.h"
#include "CALabel.h"
#include "CARoad.h"

namespace CAIOUtils {
std::vector<CAEvent> loadEventData(const std::string&);
std::vector<std::unordered_map<int, CALabel>> loadLabels(const int, const std::string&);
void writeRoadsReport(std::ofstream&, std::ofstream&, std::ofstream&, const std::vector<std::vector<CARoad>>&,
    const std::unordered_map<int, CALabel>&);
}

#endif /* TRACKINGITSU_INCLUDE_CAEVENTLOADER_H_ */
