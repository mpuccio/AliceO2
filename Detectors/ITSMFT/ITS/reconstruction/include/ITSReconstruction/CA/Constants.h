// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Constants.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_CONSTANTS_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_CONSTANTS_H_

#include <array>

#include "DetectorsBase/Constants.h"
#include "ITSReconstruction/CA/ArrayUtils.h"

namespace o2
{
namespace ITS
{
namespace CA
{
namespace Constants
{

namespace ITS
{
constexpr int LayersNumber{ 7 };
constexpr int TrackletsPerRoad{ LayersNumber - 1 };
constexpr int CellsPerRoad{ LayersNumber - 2 };
constexpr int UnusedIndex{ -1 };
constexpr float Resolution{ 0.0005f };

constexpr std::array<float, LayersNumber> LayersZCoordinate{ { 16.333f, 16.333f, 16.333f, 42.140f, 42.140f, 73.745f,
                                                               73.745f } };
constexpr std::array<float, LayersNumber> LayersRCoordinate{ { 2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f,
                                                               34.388f, 39.3329f } };
}

namespace Thresholds
{
constexpr std::array<float, ITS::TrackletsPerRoad> TrackletMaxDeltaZ{ { 0.1f, 0.1f, 0.3f, 0.3f, 0.3f, 0.3f } };
constexpr float CellMaxDeltaTanLambda{ 0.025f };
constexpr std::array<float, ITS::CellsPerRoad> CellMaxDeltaZ{ { 0.2f, 0.4f, 0.5f, 0.6f, 3.0f } };
constexpr std::array<float, ITS::CellsPerRoad> CellMaxDCA{ { 0.05f, 0.04f, 0.05f, 0.2f, 0.4f } };
constexpr float CellMaxDeltaPhiThreshold{ 0.14f };
constexpr float ZCoordinateCut{ 0.5f };
constexpr float PhiCoordinateCut{ 0.3f };
constexpr std::array<float, ITS::CellsPerRoad - 1> NeighbourCellMaxNormalVectorsDelta{ { 0.002f, 0.009f, 0.002f,
                                                                                         0.005f } };
constexpr std::array<float, ITS::CellsPerRoad - 1> NeighbourCellMaxCurvaturesDelta{ { 0.008f, 0.0025f, 0.003f,
                                                                                      0.0035f } };
constexpr int CellsMinLevel{ 5 };
}

namespace IndexTable
{
constexpr int ZBins{ 20 };
constexpr int PhiBins{ 20 };
constexpr float InversePhiBinSize{ Constants::IndexTable::PhiBins / Base::Constants::k2PI };
constexpr float getInverseBinSize(const int layerIndex) { return 0.5f * ZBins / ITS::LayersZCoordinate[layerIndex]; }

constexpr std::array<float, ITS::LayersNumber> InverseZBinSize{ ArrayUtils::fillArray<float, ITS::LayersNumber>(
  getInverseBinSize) };
}

namespace Memory
{
constexpr std::array<float, ITS::TrackletsPerRoad> TrackletsMemoryCoefficients{
  { 0.0016353f, 0.0013627f, 0.000984f, 0.00078135f, 0.00057934f, 0.00052217f }
};
constexpr std::array<float, ITS::CellsPerRoad> CellsMemoryCoefficients{ { 2.3208e-08f, 2.104e-08f, 1.6432e-08f,
                                                                          1.2412e-08f, 1.3543e-08f } };
}

namespace PDGCodes
{
constexpr int PionCode{ 211 };
}
}

}
}
}

#endif /* O2_ITSMFT_RECONSTRUCTION_CA_CONSTANTS_H_ */
