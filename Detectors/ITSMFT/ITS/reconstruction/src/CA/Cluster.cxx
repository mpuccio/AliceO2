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
/// \file Cluster.cxx
/// \brief
///

#include "ITSReconstruction/CA/Cluster.h"

#include "ITSReconstruction/CA/IndexTableUtils.h"
#include "ITSReconstruction/CA/MathUtils.h"

namespace o2
{
namespace ITS
{
namespace CA
{

Cluster::Cluster(const int clusterId, const int layerIndex, const float xCoordinate, const float yCoordinate,
    const float zCoordinate, const float alphaAngle, const int monteCarloId)
    : xCoordinate { xCoordinate }, yCoordinate { yCoordinate }, zCoordinate { zCoordinate }, phiCoordinate { 0 }, rCoordinate {
        0 }, clusterId { clusterId }, alphaAngle { alphaAngle }, monteCarloId { monteCarloId }, indexTableBinIndex { 0 }
{
  // Nothing to do
}

Cluster::Cluster(const int layerIndex, const float3 &primaryVertex, const Cluster& other)
    : xCoordinate { other.xCoordinate }, yCoordinate { other.yCoordinate }, zCoordinate { other.zCoordinate }, phiCoordinate {
        MathUtils::getNormalizedPhiCoordinate(
            MathUtils::calculatePhiCoordinate(xCoordinate - primaryVertex.x, yCoordinate - primaryVertex.y)) }, rCoordinate {
        MathUtils::calculateRCoordinate(xCoordinate - primaryVertex.x, yCoordinate - primaryVertex.y) }, clusterId {
        other.clusterId }, alphaAngle { other.alphaAngle }, monteCarloId { other.monteCarloId }, indexTableBinIndex {
        IndexTableUtils::getBinIndex(IndexTableUtils::getZBinIndex(layerIndex, zCoordinate),
            IndexTableUtils::getPhiBinIndex(phiCoordinate)) }
{
  // Nothing to do
}

}
}
}
