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
/// \file TimeFrame.h
/// \brief
///

#ifndef TRACKINGITSU_INCLUDE_TIMEFRAME_H_
#define TRACKINGITSU_INCLUDE_TIMEFRAME_H_

#include <array>
#include <vector>
#include <utility>
#include <cassert>
#include <gsl/gsl>

#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITStracking/Cell.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Configuration.h"
#include "ITStracking/Constants.h"
#include "ITStracking/Definitions.h"
#include "ITStracking/Road.h"
#include "ITStracking/Tracklet.h"

#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

namespace o2
{

namespace itsmft
{
class Cluster;
class CompClusterExt;
class TopologyDictionary;
}

namespace its
{

class TimeFrame final
{
 public:
  const float3& getPrimaryVertex(const int) const;
  gsl::span<const float3> getPrimaryVertices(int tf) const;
  gsl::span<const float3> getPrimaryVertices(int romin, int romax) const;
  int getPrimaryVerticesNum(int rofID = -1) const;
  void addPrimaryVertices(const std::vector<std::pair<float3, int>>& vertices);
  int loadROFrameData(const o2::itsmft::ROFRecord& rof, gsl::span<const itsmft::Cluster> clusters,
                      const dataformats::MCTruthContainer<MCCompLabel>* mcLabels = nullptr);
  int loadROFrameData(const o2::itsmft::ROFRecord& rof, gsl::span<const itsmft::CompClusterExt> clusters, gsl::span<const unsigned char>::iterator& pattIt,
                      const itsmft::TopologyDictionary& dict, const dataformats::MCTruthContainer<MCCompLabel>* mcLabels = nullptr);
  int getTotalClusters() const;
  bool empty() const;

  int   getSortedIndex(int rof, int layer, int i) const;
  int   getNrof() const;
  float getBeamX() const;
  float getBeamY() const;
  gsl::span<Cluster> getClustersOnLayer(int rofId, int layerId);
  index_table_t& getIndexTables(int tf);
  const std::vector<TrackingFrameInfo>& getTrackingFrameInfoOnLayer(int layerId) const;

  const TrackingFrameInfo& getClusterTrackingFrameInfo(int layerId, const Cluster& cl) const;
  const MCCompLabel& getClusterLabels(int layerId, const Cluster& cl) const;
  const MCCompLabel& getClusterLabels(int layerId, const int clId) const;
  int getClusterExternalIndex(int layerId, const int clId) const;

  bool hasMCinformation() const;
  void initialise(const int iteration, const MemoryParameters& memParam);

  bool isClusterUsed(int layer, int clusterId) const;
  void markUsedCluster(int layer, int clusterId);

  std::array<std::vector<Tracklet>, constants::its::TrackletsPerRoad>& getTracklets();
  std::array<std::vector<int>, constants::its::CellsPerRoad>& getTrackletsLookupTable();

  std::array<std::vector<Cluster>, constants::its::LayersNumber>& getClusters();
  std::array<std::vector<Cluster>, constants::its::LayersNumber>& getUnsortedClusters();
  std::array<std::vector<Cell>, constants::its::CellsPerRoad>& getCells();
  std::array<std::vector<int>, constants::its::CellsPerRoad - 1>& getCellsLookupTable();
  std::array<std::vector<std::vector<int>>, constants::its::CellsPerRoad - 1>& getCellsNeighbours();
  std::vector<Road>& getRoads();

  void initialiseRoadLabels();
  void setRoadLabel(int i, const unsigned long long& lab, bool fake);
  const unsigned long long& getRoadLabel(int i) const;
  bool isRoadFake(int i) const;

  void clear();

 private:
  template <typename... T>
  void addClusterToLayer(int layer, T&&... args);
  template <typename... T>
  void addTrackingFrameInfoToLayer(int layer, T&&... args);
  void addClusterLabelToLayer(int layer, const MCCompLabel label);
  void addClusterExternalIndexToLayer(int layer, const int idx);
  gsl::span<const Cluster> getUnsortedClustersOnLayer(int rofId, int layerId);

  int                                                                      mNrof = 0;
  int                                                                      mBeamPosWeight = 0;
  float                                                                    mBeamPos[2] = {0.f,0.f};
  std::vector<int>                                                         mROframesPV = {0};
  std::array<std::vector<int>, constants::its::LayersNumber>               mROframesClusters = {std::vector<int>(1,0), std::vector<int>(1,0), std::vector<int>(1,0), std::vector<int>(1,0), std::vector<int>(1,0), std::vector<int>(1,0), std::vector<int>(1,0)};
  std::vector<float3>                                                      mPrimaryVertices;
  std::array<std::vector<Cluster>, constants::its::LayersNumber>           mClusters;
  std::array<std::vector<Cluster>, constants::its::LayersNumber>           mUnsortedClusters;
  std::array<std::vector<bool>, constants::its::LayersNumber>              mUsedClusters;
  std::array<std::vector<TrackingFrameInfo>, constants::its::LayersNumber> mTrackingFrameInfo;
  std::array<std::vector<MCCompLabel>, constants::its::LayersNumber>       mClusterLabels;
  std::array<std::vector<int>, constants::its::LayersNumber>               mClusterExternalIndices;
  std::array<std::vector<Cell>, constants::its::CellsPerRoad>              mCells;
  std::array<std::vector<int>, constants::its::CellsPerRoad - 1>           mCellsLookupTable;
  std::array<std::vector<std::vector<int>>, constants::its::CellsPerRoad - 1> mCellsNeighbours;
  std::vector<Road>                                                        mRoads;

  std::vector<index_table_t>                                               mIndexTables;
  std::array<std::vector<Tracklet>, constants::its::TrackletsPerRoad>      mTracklets;
  std::array<std::vector<int>, constants::its::CellsPerRoad>               mTrackletsLookupTable;

  std::vector<std::pair<unsigned long long, bool>>                         mRoadLabels;
};


inline const float3& TimeFrame::getPrimaryVertex(const int vertexIndex) const { return mPrimaryVertices[vertexIndex]; }

inline gsl::span<const float3> TimeFrame::getPrimaryVertices(int tf) const
{ 
  const int start = tf > 0 ? tf - 1 : 0;
  const int stop = tf >= mNrof - 1 ? mNrof : tf + 2;
  return {&mPrimaryVertices[start], mROframesPV[stop] - mROframesPV[start]};
}

inline gsl::span<const float3> TimeFrame::getPrimaryVertices(int romin, int romax) const 
{
  return {&mPrimaryVertices[romin], mROframesPV[romax + 1] - mROframesPV[romin]};
}

inline int TimeFrame::getPrimaryVerticesNum(int rofID) const
{
  return rofID < 0 ? mPrimaryVertices.size() : mROframesPV[rofID + 1] - mROframesPV[rofID];
}

inline bool TimeFrame::empty() const { return getTotalClusters() == 0; }

inline int TimeFrame::getSortedIndex(int rof, int layer, int index) const { return mROframesClusters[layer][rof] + index; }

inline int TimeFrame::getNrof() const { return mNrof; };

inline float TimeFrame::getBeamX() const { return mBeamPos[0]; }

inline float TimeFrame::getBeamY() const { return mBeamPos[1]; }

inline gsl::span<Cluster> TimeFrame::getClustersOnLayer(int rofId, int layerId)
{
  gsl::span<Cluster> retSpan = {&mClusters[layerId][mROframesClusters[layerId][rofId]], mROframesClusters[layerId][rofId + 1] - mROframesClusters[layerId][rofId]};
  if (mROframesClusters[layerId][rofId]+retSpan.size() > mClusters[layerId].size())
    std::cout << "Catastrofic problem in getClustersOnLayer" << std::endl;
  return retSpan;
}

inline gsl::span<const Cluster> TimeFrame::getUnsortedClustersOnLayer(int rofId, int layerId)
{
  return {&mUnsortedClusters[layerId][mROframesClusters[layerId][rofId]], mROframesClusters[layerId][rofId + 1] - mROframesClusters[layerId][rofId]};
}

inline const std::vector<TrackingFrameInfo>& TimeFrame::getTrackingFrameInfoOnLayer(int layerId) const
{
  return mTrackingFrameInfo[layerId];
}

inline const TrackingFrameInfo& TimeFrame::getClusterTrackingFrameInfo(int layerId, const Cluster& cl) const
{
  return mTrackingFrameInfo[layerId][cl.clusterId];
}

inline const MCCompLabel& TimeFrame::getClusterLabels(int layerId, const Cluster& cl) const
{
  return mClusterLabels[layerId][cl.clusterId];
}

inline const MCCompLabel& TimeFrame::getClusterLabels(int layerId, const int clId) const
{
  return mClusterLabels[layerId][clId];
}

inline int TimeFrame::getClusterExternalIndex(int layerId, const int clId) const
{
  return mClusterExternalIndices[layerId][clId];
}

inline index_table_t& TimeFrame::getIndexTables(int tf)
{
  return mIndexTables[tf];
}

template <typename... T>
void TimeFrame::addClusterToLayer(int layer, T&&... values)
{
  mUnsortedClusters[layer].emplace_back(std::forward<T>(values)...);
}

template <typename... T>
void TimeFrame::addTrackingFrameInfoToLayer(int layer, T&&... values)
{
  mTrackingFrameInfo[layer].emplace_back(std::forward<T>(values)...);
}

inline void TimeFrame::addClusterLabelToLayer(int layer, const MCCompLabel label) { mClusterLabels[layer].emplace_back(label); }

inline void TimeFrame::addClusterExternalIndexToLayer(int layer, const int idx)
{
  mClusterExternalIndices[layer].push_back(idx);
}

inline void TimeFrame::clear()
{
  for (int iL = 0; iL < constants::its::LayersNumber; ++iL) {
    mClusters[iL].clear();
    mTrackingFrameInfo[iL].clear();
    mClusterLabels[iL].clear();
    mClusterExternalIndices[iL].clear();
  }
  mPrimaryVertices.clear();
}

inline bool TimeFrame::hasMCinformation() const
{
  for (const auto& vect : mClusterLabels) {
    if (!vect.empty()) {
      return true;
    }
  }
  return false;
}

inline bool TimeFrame::isClusterUsed(int layer, int clusterId) const
{
  return mUsedClusters[layer][clusterId];
}

inline void TimeFrame::markUsedCluster(int layer, int clusterId) { mUsedClusters[layer][clusterId] = true; }

inline std::array<std::vector<Tracklet>, constants::its::TrackletsPerRoad>& TimeFrame::getTracklets()
{
  return mTracklets;
}

inline std::array<std::vector<int>, constants::its::CellsPerRoad>& TimeFrame::getTrackletsLookupTable()
{
  return mTrackletsLookupTable;
}

inline void TimeFrame::initialiseRoadLabels()
{
  mRoadLabels.clear();
  mRoadLabels.resize(mRoads.size());
}

inline void TimeFrame::setRoadLabel(int i, const unsigned long long& lab, bool fake)
{
  mRoadLabels[i].first = lab;
  mRoadLabels[i].second = fake;
}

inline const unsigned long long& TimeFrame::getRoadLabel(int i) const
{
  return mRoadLabels[i].first;
}

inline bool TimeFrame::isRoadFake(int i) const
{
  return mRoadLabels[i].second;
}

inline std::array<std::vector<Cluster>, constants::its::LayersNumber>& TimeFrame::getClusters()
{
  return mClusters;
}

inline std::array<std::vector<Cluster>, constants::its::LayersNumber>& TimeFrame::getUnsortedClusters()
{
  return mUnsortedClusters;
}


inline std::array<std::vector<Cell>, constants::its::CellsPerRoad>& TimeFrame::getCells() { return mCells; }

inline std::array<std::vector<int>, constants::its::CellsPerRoad - 1>& TimeFrame::getCellsLookupTable()
{
  return mCellsLookupTable;
}

inline std::array<std::vector<std::vector<int>>, constants::its::CellsPerRoad - 1>&
  TimeFrame::getCellsNeighbours()
{
  return mCellsNeighbours;
}

inline std::vector<Road>& TimeFrame::getRoads() { return mRoads; }

} // namespace its
} // namespace o2

#endif /* TRACKINGITSU_INCLUDE_TimeFrame_H_ */
