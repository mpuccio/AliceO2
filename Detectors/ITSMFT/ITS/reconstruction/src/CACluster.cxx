// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CACluster.cxx
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch
#include <CACluster.h>

#include <CAIndexTableUtils.h>
#include <CAMathUtils.h>

CACluster::CACluster(const int clusterId, const float xCoordinate, const float yCoordinate, const float zCoordinate,
    const float alphaAngle, const int monteCarloId)
    : clusterId { clusterId }, xCoordinate { xCoordinate }, yCoordinate { yCoordinate }, zCoordinate { zCoordinate }, alphaAngle {
        alphaAngle }, monteCarloId { monteCarloId }, phiCoordinate { 0 }, rCoordinate { 0 }, indexTableBinIndex { 0 }
{
}

CACluster::CACluster(const int layerIndex, const std::array<float, 3> &primaryVertex, const CACluster& other)
    : clusterId { other.clusterId }, xCoordinate { other.xCoordinate }, yCoordinate { other.yCoordinate }, zCoordinate {
        other.zCoordinate }, alphaAngle { other.alphaAngle }, monteCarloId { other.monteCarloId }, phiCoordinate {
        CAMathUtils::getNormalizedPhiCoordinate(
            CAMathUtils::calculatePhiCoordinate(xCoordinate - primaryVertex[0], yCoordinate - primaryVertex[1])) }, rCoordinate {
        CAMathUtils::calculateRCoordinate(xCoordinate - primaryVertex[0], yCoordinate - primaryVertex[1]) }, indexTableBinIndex {
        CAIndexTableUtils::getBinIndex(CAIndexTableUtils::getZBinIndex(layerIndex, zCoordinate),
            CAIndexTableUtils::getPhiBinIndex(phiCoordinate)) }
{
}
