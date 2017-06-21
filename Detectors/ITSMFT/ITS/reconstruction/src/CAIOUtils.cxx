// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAIOUtils.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "CAIOUtils.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <utility>

#include "CAConstants.h"

namespace {

constexpr int PrimaryVertexLayerId = -1;
constexpr int EventLabelsSeparator = -1;
}

std::vector<CAEvent> CAIOUtils::loadEventData(const std::string& fileName)
{
  std::vector<CAEvent> events;
  std::ifstream inputStream;
  std::string line, unusedVariable;
  int layerId, monteCarlo;
  int clusterId = -1;
  float xCoordinate, yCoordinate, zCoordinate, alphaAngle;

  inputStream.open(fileName);

  while (std::getline(inputStream, line)) {

    std::istringstream inputStringStream(line);

    if (inputStringStream >> layerId >> xCoordinate >> yCoordinate >> zCoordinate) {

      if (layerId == PrimaryVertexLayerId) {

        if (clusterId != 0) {

          events.emplace_back(events.size());
        }

        events.back().addPrimaryVertex(xCoordinate, yCoordinate, zCoordinate);
        clusterId = 0;

      } else {

        if (inputStringStream >> unusedVariable >> unusedVariable >> unusedVariable >> alphaAngle >> monteCarlo) {

          events.back().pushClusterToLayer(layerId, clusterId, xCoordinate, yCoordinate, zCoordinate, alphaAngle,
              monteCarlo);
          ++clusterId;
        }
      }
    }
  }

  return events;
}

std::vector<std::unordered_map<int, CALabel>> CAIOUtils::loadLabels(const int eventsNum, const std::string& fileName)
{
  std::vector<std::unordered_map<int, CALabel>> labelsMap;
  std::unordered_map<int, CALabel> currentEventLabelsMap;
  std::ifstream inputStream;
  std::string line;

  int monteCarloId, pdgCode, numberOfClusters;
  float transverseMomentum, phiCoordinate, pseudorapidity;

  labelsMap.reserve(eventsNum);

  inputStream.open(fileName);
  std::getline(inputStream, line);

  while (std::getline(inputStream, line)) {

    std::istringstream inputStringStream(line);

    if (inputStringStream >> monteCarloId) {

      if (monteCarloId == EventLabelsSeparator) {

        labelsMap.emplace_back(currentEventLabelsMap);
        currentEventLabelsMap.clear();

      } else {

        if (inputStringStream >> transverseMomentum >> phiCoordinate >> pseudorapidity >> pdgCode >> numberOfClusters) {

          if (std::abs(pdgCode) == CAConstants::PDGCodes::PionCode && numberOfClusters == 7) {

            currentEventLabelsMap.emplace(std::piecewise_construct, std::forward_as_tuple(monteCarloId),
                std::forward_as_tuple(monteCarloId, transverseMomentum, phiCoordinate, pseudorapidity, pdgCode,
                    numberOfClusters));
          }
        }
      }
    }
  }

  labelsMap.emplace_back(currentEventLabelsMap);

  return labelsMap;
}

void CAIOUtils::writeRoadsReport(std::ofstream& correctRoadsOutputStream, std::ofstream& duplicateRoadsOutputStream,
    std::ofstream& fakeRoadsOutputStream, const std::vector<std::vector<CARoad>>& roads,
    const std::unordered_map<int, CALabel>& labelsMap)
{
  const int numVertices = roads.size();
  std::unordered_set<int> foundMonteCarloIds;

  correctRoadsOutputStream << EventLabelsSeparator << std::endl;
  fakeRoadsOutputStream << EventLabelsSeparator << std::endl;

  for (int iVertex = 0; iVertex < numVertices; ++iVertex) {

    const std::vector<CARoad>& currentVertexRoads = roads[iVertex];
    const int numRoads = currentVertexRoads.size();

    for (int iRoad = 0; iRoad < numRoads; ++iRoad) {

      const CARoad& currentRoad = currentVertexRoads[iRoad];
      const int currentRoadLabel = currentRoad.getLabel();

      if (!labelsMap.count(currentRoadLabel)) {

        continue;
      }

      const CALabel& currentLabel = labelsMap.at(currentRoadLabel);

      if (currentRoad.isFakeRoad()) {

        fakeRoadsOutputStream << currentLabel << std::endl;

      } else {

        if (foundMonteCarloIds.count(currentLabel.monteCarloId)) {

          duplicateRoadsOutputStream << currentLabel << std::endl;

        } else {

          correctRoadsOutputStream << currentLabel << std::endl;
          foundMonteCarloIds.emplace(currentLabel.monteCarloId);
        }
      }
    }
  }
}
