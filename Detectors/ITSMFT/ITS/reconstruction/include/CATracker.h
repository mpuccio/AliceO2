// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CATracker.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CATRACKER_H_
#define TRACKINGITSU_INCLUDE_CATRACKER_H_

#include <iostream>
#include <vector>

#include "CARoad.h"

class CAEvent;
struct CAPrimaryVertexContext;

class CATracker final
{
  public:
    explicit CATracker(const CAEvent&);

    CATracker(const CATracker&) = delete;
    CATracker &operator=(const CATracker&) = delete;

    std::vector<std::vector<CARoad>> clustersToTracks();
    std::vector<std::vector<CARoad>> clustersToTracksVerbose();
    std::vector<std::vector<CARoad>> clustersToTracksMemoryBenchmark(std::ofstream&);

  protected:
    void computeTracklets(CAPrimaryVertexContext&);
    void computeCells(CAPrimaryVertexContext&);
    void findCellsNeighbours(CAPrimaryVertexContext&);
    void findTracks(CAPrimaryVertexContext&);
    void traverseCellsTree(CAPrimaryVertexContext&, const int, const int);
    void computeMontecarloLabels(CAPrimaryVertexContext&);

  private:
    const CAEvent& mEvent;
    std::vector<int> mUsedClustersTable;
};

#endif /* TRACKINGITSU_INCLUDE_CATRACKER_H_ */
