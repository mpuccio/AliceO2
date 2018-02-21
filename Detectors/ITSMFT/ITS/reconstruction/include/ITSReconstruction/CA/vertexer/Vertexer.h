// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Vertexer.h
/// \brief
/// \author matteo.concas@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_VERTEXER_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_VERTEXER_H_

#include <vector>
#include <array>

#include "ITSReconstruction/CA/Constants.h"
#include "ITSReconstruction/CA/Definitions.h"
#include "ITSReconstruction/CA/Event.h"

namespace o2
{
namespace ITS
{
namespace CA
{
class Cluster;
class Event;

class Vertexer
{
public:
  explicit Vertexer(const Event& event);
  virtual ~Vertexer();
  Vertexer(const Vertexer&) = delete;
  Vertexer& operator=(const Vertexer&) = delete;

  void initialize(const float zCut, const float phiCut);
  void computeTriplets();
  void checkTriplets();
  void FindVertices();
  const std::vector<std::pair<int, int>> selectClusters(
    const std::array<int, Constants::IndexTable::ZBins * Constants::IndexTable::PhiBins + 1> &indexTable,
    const std::array<int, 4> &selectedBinsRect);

protected:
  bool mVertexerInitialized;
  float mDeltaRadii10, mDeltaRadii21;
  float mZCut, mPhiCut;
  int mPhiSpan, mZSpan;
  std::array<float, 3> mAverageClustersRadii; 
  std::array<float, Constants::ITS::LayersNumber> mITSRadii;
  float mZBinSize;
  std::vector<std::pair<int, int>> mClustersToProcessInner;
  std::vector<std::pair<int, int>> mClustersToProcessOuter;
  Event mEvent;
  std::vector<std::array<float, 3>> mVertices;
  std::array<std::vector<Cluster>, Constants::ITS::LayersNumberVertexer> mClusters;
  std::array<std::array<int, Constants::IndexTable::ZBins * Constants::IndexTable::PhiBins + 1>,
            Constants::ITS::LayersNumberVertexer> mIndexTables;
  std::vector<std::array<int, 3>> mTriplets;
};

}
}
}

#endif /* O2_ITSMFT_RECONSTRUCTION_CA_VERTEXER_H_ */ 