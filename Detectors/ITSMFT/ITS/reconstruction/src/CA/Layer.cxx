// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Layer.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#include "ITSReconstruction/CA/Layer.h"

#include "ITSReconstruction/CA/Constants.h"

namespace o2
{
namespace ITS
{
namespace CA
{

Layer::Layer() : mLayerIndex{ Constants::ITS::UnusedIndex }
{
  // Nothing to do
}

Layer::Layer(const int layerIndex) : mLayerIndex{ layerIndex }
{
  // Nothing to do
}

void Layer::addCluster(const int clusterId, const float cluster_xyz[3])
{
  mClusters.emplace_back(clusterId, cluster_xyz[0], cluster_xyz[1], cluster_xyz[2], mLayerIndex);
}

}
}
}
