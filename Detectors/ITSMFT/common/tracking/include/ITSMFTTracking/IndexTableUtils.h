// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file IndexTableUtils.h
/// \brief Shared index-table utilities for ITS (phi-z) and MFT (x-y)
///

#ifndef ALICEO2_ITSMFT_TRACKING_INDEXTABLEUTILS_H_
#define ALICEO2_ITSMFT_TRACKING_INDEXTABLEUTILS_H_

#include <array>

#include "CommonConstants/MathConstants.h"
#include "GPUCommonMath.h"
#include "GPUCommonDef.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITStracking/Cluster.h"
#include "MFTTracking/Constants.h"

namespace o2::itsmft
{

enum class IndexTableCoordType : uint8_t { PhiZ, XY };

namespace index_table_utils
{
GPUhdi() float getNormalizedPhi(float phi)
{
  phi -= o2::constants::math::TwoPI * o2::gpu::GPUCommonMath::Floor(phi * (1.f / o2::constants::math::TwoPI));
  return phi;
}
} // namespace index_table_utils

/// Row/column LUT helper (ITS: row=phi, col=z; MFT: row=y, col=x).
template <int nLayers>
class IndexTableUtils
{
 public:
  /// Configure LUT geometry. ITS (PhiZ): row = phi [0, TwoPI), col = z; MFT (XY): row = y, col = x.
  void setIndexTableParams(IndexTableCoordType coordType, int nRowBins, int nColBins,
                           float rowMin, float rowMax,
                           const std::array<float, nLayers>& layerColHalfExtent)
  {
    mCoordType = coordType;
    mRowOrigin = (coordType == IndexTableCoordType::PhiZ) ? 0.f : rowMin;
    mRowCoordinateSpan = rowMax - rowMin;
    mInverseRowBinSize = (mRowCoordinateSpan > 0.f) ? static_cast<float>(nRowBins) / mRowCoordinateSpan : 0.f;
    mNcolBins = nColBins;
    mNrowBins = nRowBins;
    for (int iLayer{0}; iLayer < nLayers; ++iLayer) {
      mLayerColHalfExtent[iLayer] = layerColHalfExtent[iLayer];
      mInverseColBinSize[iLayer] = 0.5f * nColBins / layerColHalfExtent[iLayer];
    }
  }

  /// Fill LUT geometry from any struct exposing RowBins, ColBins and LayerZ (ITS phi-z).
  template <class T>
  void setTrackingParameters(const T& params)
  {
    setIndexTableParams(IndexTableCoordType::PhiZ, params.RowBins, params.ColBins,
                        0.f, o2::constants::math::TwoPI, layerColHalfExtentFrom(params));
  }

  /// Fill LUT geometry for MFT (row = global y, col = global x).
  template <class T>
  void setTrackingParametersXY(const T& params, float rowMin, float rowMax)
  {
    setIndexTableParams(IndexTableCoordType::XY, params.RowBins, params.ColBins,
                        rowMin, rowMax, layerColHalfExtentFrom(params));
  }

  GPUhdi() float getInverseColCoordinate(const int layerIndex) const
  {
    return 0.5f * mNcolBins / mLayerColHalfExtent[layerIndex];
  }

  GPUhdi() int getColBinIndex(const int layerIndex, const float colCoordinate) const
  {
    return (colCoordinate + mLayerColHalfExtent[layerIndex]) * mInverseColBinSize[layerIndex];
  }

  GPUhdi() int getRowBinIndex(const float rowCoordinate) const
  {
    if (mCoordType == IndexTableCoordType::PhiZ) {
      return rowCoordinate * mInverseRowBinSize;
    }
    return (rowCoordinate - mRowOrigin) * mInverseRowBinSize;
  }

  GPUhdi() int getBinIndex(const int colIndex, const int rowIndex) const
  {
    return o2::gpu::GPUCommonMath::Min(rowIndex * mNcolBins + colIndex, (mNcolBins * mNrowBins) - 1);
  }

  GPUhdi() int countRowSelectedBins(const int* indexTable, const int rowBinIndex,
                                    const int minColBinIndex, const int maxColBinIndex) const
  {
    const int firstBinIndex{getBinIndex(minColBinIndex, rowBinIndex)};
    const int maxBinIndex{firstBinIndex + maxColBinIndex - minColBinIndex + 1};

    return indexTable[maxBinIndex] - indexTable[firstBinIndex];
  }

  void print() const;

  GPUhdi() int getNcolBins() const { return mNcolBins; }
  GPUhdi() int getNrowBins() const { return mNrowBins; }
  GPUhdi() float getLayerColHalfExtent(int i) const { return mLayerColHalfExtent[i]; }
  GPUhdi() void setNcolBins(const int colBins) { mNcolBins = colBins; }
  GPUhdi() void setNrowBins(const int rowBins) { mNrowBins = rowBins; }
  GPUhdi() IndexTableCoordType getCoordType() const { return mCoordType; }

 private:
  template <class T>
  static std::array<float, nLayers> layerColHalfExtentFrom(const T& params)
  {
    std::array<float, nLayers> extents{};
    if constexpr (requires { params.LayerColHalfExtent; }) {
      const auto& colExtents = params.LayerColHalfExtent.empty() ? params.LayerZ : params.LayerColHalfExtent;
      for (int iLayer{0}; iLayer < nLayers && iLayer < static_cast<int>(colExtents.size()); ++iLayer) {
        extents[iLayer] = colExtents[iLayer];
      }
    } else {
      for (int iLayer{0}; iLayer < nLayers && iLayer < static_cast<int>(params.LayerZ.size()); ++iLayer) {
        extents[iLayer] = params.LayerZ[iLayer];
      }
    }
    return extents;
  }

  int mNcolBins = 0;
  int mNrowBins = 0;
  float mInverseRowBinSize = 0.f;
  float mRowOrigin = 0.f;
  float mRowCoordinateSpan = o2::constants::math::TwoPI;
  IndexTableCoordType mCoordType{IndexTableCoordType::PhiZ};
  std::array<float, nLayers> mLayerColHalfExtent{};
  std::array<float, nLayers> mInverseColBinSize{};
};

template <int nLayers>
void IndexTableUtils<nLayers>::print() const
{
  printf("NcolBins: %d, NrowBins: %d, InverseRowBinSize: %f\n", mNcolBins, mNrowBins, mInverseRowBinSize);
  for (int iLayer{0}; iLayer < nLayers; ++iLayer) {
    printf("Layer %d: ColHalfExtent: %f, InverseColBinSize: %f\n", iLayer, mLayerColHalfExtent[iLayer], mInverseColBinSize[iLayer]);
  }
}

/// ITS: row = phi, col = z.
template <int nLayers>
GPUhdi() int4 getBinsPhiZ(float phi, const int layerIndex,
                          float z1, float z2, float maxDeltaCol, float maxDeltaRow,
                          const IndexTableUtils<nLayers>& utils)
{
  const float colRangeMin = o2::gpu::GPUCommonMath::Min(z1, z2) - maxDeltaCol;
  const float rowRangeMin = (maxDeltaRow > o2::constants::math::PI) ? 0.f : phi - maxDeltaRow;
  const float colRangeMax = o2::gpu::GPUCommonMath::Max(z1, z2) + maxDeltaCol;
  const float rowRangeMax = (maxDeltaRow > o2::constants::math::PI) ? o2::constants::math::TwoPI : phi + maxDeltaRow;

  if (colRangeMax < -utils.getLayerColHalfExtent(layerIndex) ||
      colRangeMin > utils.getLayerColHalfExtent(layerIndex) || colRangeMin > colRangeMax) {
    return int4{-1, -1, -1, -1};
  }

  return int4{o2::gpu::GPUCommonMath::Max(0, utils.getColBinIndex(layerIndex, colRangeMin)),
              utils.getRowBinIndex(index_table_utils::getNormalizedPhi(rowRangeMin)),
              o2::gpu::GPUCommonMath::Min(utils.getNcolBins() - 1, utils.getColBinIndex(layerIndex, colRangeMax)),
              utils.getRowBinIndex(index_table_utils::getNormalizedPhi(rowRangeMax))};
}

/// MFT: row = y, col = x.
template <int nLayers>
GPUhdi() int4 getBinsXY(float x, float y, const int layerIndex,
                        float x1, float x2, float maxDeltaCol, float maxDeltaRow,
                        const IndexTableUtils<nLayers>& utils)
{
  const float colRangeMin = o2::gpu::GPUCommonMath::Min(x1, x2) - maxDeltaCol;
  const float rowRangeMin = y - maxDeltaRow;
  const float colRangeMax = o2::gpu::GPUCommonMath::Max(x1, x2) + maxDeltaCol;
  const float rowRangeMax = y + maxDeltaRow;

  const float colHalf = utils.getLayerColHalfExtent(layerIndex);
  if (colRangeMax < -colHalf || colRangeMin > colHalf || colRangeMin > colRangeMax) {
    return int4{-1, -1, -1, -1};
  }

  return int4{o2::gpu::GPUCommonMath::Max(0, utils.getColBinIndex(layerIndex, colRangeMin)),
              o2::gpu::GPUCommonMath::Max(0, utils.getRowBinIndex(rowRangeMin)),
              o2::gpu::GPUCommonMath::Min(utils.getNcolBins() - 1, utils.getColBinIndex(layerIndex, colRangeMax)),
              utils.getRowBinIndex(rowRangeMax)};
}

/// MFT cone projection from one half-layer to another (used for x-y index-table lookup only).
template <int nLayers>
GPUhdi() void mftConeProject(const o2::its::Cluster& cluster, int fromLayer, int toLayer, float& xProj, float& yProj)
{
  const auto& layerZ = o2::mft::constants::mft::LayerZCoordinate();
  const auto& invLayerZ = o2::mft::constants::mft::InverseLayerZCoordinate();
  const float scale = (layerZ[toLayer] - layerZ[fromLayer]) * invLayerZ[fromLayer];
  xProj = cluster.xCoordinate * (1.f + scale);
  yProj = cluster.yCoordinate * (1.f + scale);
}

/// Cluster-driven LUT window: phi-z for ITS, cone-projected x-y for MFT.
template <int nLayers>
GPUhdi() int4 getBinsRectCluster(const o2::its::Cluster& cluster, int fromLayer, int toLayer,
                                 float colRangeMin, float colRangeMax, float maxDeltaCol, float maxDeltaRow,
                                 const IndexTableUtils<nLayers>& utils)
{
  if (utils.getCoordType() == IndexTableCoordType::XY) {
    float xProj = 0.f;
    float yProj = 0.f;
    mftConeProject<nLayers>(cluster, fromLayer, toLayer, xProj, yProj);
    return getBinsXY(xProj, yProj, toLayer, xProj, xProj, maxDeltaCol, maxDeltaRow, utils);
  }
  return getBinsPhiZ(cluster.phi, toLayer, colRangeMin, colRangeMax, maxDeltaCol, maxDeltaRow, utils);
}

} // namespace o2::itsmft

#endif /* ALICEO2_ITSMFT_TRACKING_INDEXTABLEUTILS_H_ */
