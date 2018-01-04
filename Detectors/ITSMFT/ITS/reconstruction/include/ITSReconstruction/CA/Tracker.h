// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file Tracker.h
/// \brief 
///

#ifndef TRACKINGITSU_INCLUDE_TRACKER_H_
#define TRACKINGITSU_INCLUDE_TRACKER_H_

#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include "ITSReconstruction/CA/Definitions.h"
#include "ITSReconstruction/CA/Event.h"
#include "ITSReconstruction/CA/MathUtils.h"
#include "ITSReconstruction/CA/PrimaryVertexContext.h"
#include "ITSReconstruction/CA/Road.h"

namespace o2
{
namespace ITS
{
namespace CA
{

template<bool IsGPU>
class TrackerTraits
{
  public:
    void computeLayerTracklets(PrimaryVertexContext&);
    void computeLayerCells(PrimaryVertexContext&);

  protected:
    ~TrackerTraits() = default;
};

template<bool IsGPU>
class Tracker: private TrackerTraits<IsGPU>
{
  private:
    typedef TrackerTraits<IsGPU> Trait;
  public:
    Tracker();

    Tracker(const Tracker&) = delete;
    Tracker &operator=(const Tracker&) = delete;

    void setBz(float bz);
    float getBz() const;

    std::vector<std::vector<Road>> clustersToTracks(const Event&);
    std::vector<std::vector<Road>> clustersToTracksVerbose(const Event&);
    std::vector<std::vector<Road>> clustersToTracksMemoryBenchmark(const Event&, std::ofstream&);
    std::vector<std::vector<Road>> clustersToTracksTimeBenchmark(const Event&, std::ofstream&);

  private:
    Base::Track::TrackParCov buildTrackSeed(const Cluster& cluster1, const Cluster& cluster2,
        const Cluster& cluster3, const TrackingFrameInfo& tf3);
    void computeTracklets();
    void computeCells();
    void findCellsNeighbours();
    void findRoads();
    void findTracks(const Event& ev);
    void traverseCellsTree(const int, const int);
    void computeMontecarloLabels();

    float evaluateTask(void (Tracker<IsGPU>::*)(void), const char*);
    float evaluateTask(void (Tracker<IsGPU>::*)(void), const char*, std::ostream&);

    PrimaryVertexContext mPrimaryVertexContext;
    float                mBz = 0.5f;
};

template<bool IsGPU>
float Tracker<IsGPU>::getBz() const
{
  return mBz;
}

template<bool IsGPU>
void Tracker<IsGPU>::setBz(float bz)
{
  mBz = bz;
}

template<> void TrackerTraits<TRACKINGITSU_GPU_MODE>::computeLayerTracklets(PrimaryVertexContext&);
template<> void TrackerTraits<TRACKINGITSU_GPU_MODE>::computeLayerCells(PrimaryVertexContext&);

}
}
}

#endif /* TRACKINGITSU_INCLUDE_TRACKER_H_ */
