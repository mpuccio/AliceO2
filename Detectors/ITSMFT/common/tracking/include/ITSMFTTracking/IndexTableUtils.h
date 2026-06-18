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
                        float x1, float x2, float y1, float y2,
                        float maxDeltaCol, float maxDeltaRow,
                        const IndexTableUtils<nLayers>& utils)
{
  const float colRangeMin = o2::gpu::GPUCommonMath::Min(x1, x2) - maxDeltaCol;
  const float rowRangeMin = o2::gpu::GPUCommonMath::Min(y1, y2) - maxDeltaRow;
  const float colRangeMax = o2::gpu::GPUCommonMath::Max(x1, x2) + maxDeltaCol;
  const float rowRangeMax = o2::gpu::GPUCommonMath::Max(y1, y2) + maxDeltaRow;

  const float colHalf = utils.getLayerColHalfExtent(layerIndex);
  if (colRangeMax < -colHalf || colRangeMin > colHalf || colRangeMin > colRangeMax) {
    return int4{-1, -1, -1, -1};
  }

  return int4{o2::gpu::GPUCommonMath::Max(0, utils.getColBinIndex(layerIndex, colRangeMin)),
              o2::gpu::GPUCommonMath::Max(0, utils.getRowBinIndex(rowRangeMin)),
              o2::gpu::GPUCommonMath::Min(utils.getNcolBins() - 1, utils.getColBinIndex(layerIndex, colRangeMax)),
              utils.getRowBinIndex(rowRangeMax)};
}

/// MFT: extrapolate cluster to toLayer on the line through the primary vertex.
template <int nLayers>
GPUhdi() void mftConeProject(const o2::its::Cluster& cluster, int fromLayer, int toLayer,
                             float pvX, float pvY, float pvZ, float& xProj, float& yProj)
{
  const auto& layerZ = o2::mft::constants::mft::LayerZCoordinate();
  const float zFrom = layerZ[fromLayer];
  const float zTo = layerZ[toLayer];
  const float dz0 = zFrom - pvZ;
  if (o2::gpu::GPUCommonMath::Abs(dz0) < 1e-6f) {
    xProj = cluster.xCoordinate;
    yProj = cluster.yCoordinate;
    return;
  }
  const float w = (zTo - pvZ) / dz0;
  xProj = pvX + w * (cluster.xCoordinate - pvX);
  yProj = pvY + w * (cluster.yCoordinate - pvY);
}

/// MFT LUT window around a precomputed (x, y) projection on toLayer.
template <int nLayers>
GPUhdi() int4 getBinsRectClusterAtProj(float xProj, float yProj, int toLayer,
                                       float colRangeMin, float colRangeMax, float maxDeltaCol, float maxDeltaRow,
                                       const IndexTableUtils<nLayers>& utils)
{
  const float rProj = o2::gpu::GPUCommonMath::Hypot(xProj, yProj);
  float x1 = xProj;
  float x2 = xProj;
  float y1 = yProj;
  float y2 = yProj;
  if (rProj > 0.f) {
    const float invRProj = 1.f / rProj;
    x1 = colRangeMin * xProj * invRProj;
    x2 = colRangeMax * xProj * invRProj;
    y1 = colRangeMin * yProj * invRProj;
    y2 = colRangeMax * yProj * invRProj;
  }
  return getBinsXY(xProj, yProj, toLayer, x1, x2, y1, y2, maxDeltaCol, maxDeltaRow, utils);
}

/// Cluster-driven LUT window: phi-z for ITS, projected x-y for MFT.
/// ITS: colRangeMin/Max = z window; MFT: colRangeMin/Max = rMin/rMax at toLayer from diamond z spread.
template <int nLayers>
GPUhdi() int4 getBinsRectCluster(const o2::its::Cluster& cluster, int fromLayer, int toLayer,
                                 float colRangeMin, float colRangeMax, float maxDeltaCol, float maxDeltaRow,
                                 const IndexTableUtils<nLayers>& utils,
                                 float pvX = 0.f, float pvY = 0.f, float pvZ = 0.f)
{
  if (utils.getCoordType() == IndexTableCoordType::XY) {
    float xProj = 0.f;
    float yProj = 0.f;
    mftConeProject<nLayers>(cluster, fromLayer, toLayer, pvX, pvY, pvZ, xProj, yProj);
    return getBinsRectClusterAtProj<nLayers>(xProj, yProj, toLayer, colRangeMin, colRangeMax, maxDeltaCol, maxDeltaRow, utils);
  }
  return getBinsPhiZ(cluster.phi, toLayer, colRangeMin, colRangeMax, maxDeltaCol, maxDeltaRow, utils);
}

} // namespace o2::itsmft

#endif /* ALICEO2_ITSMFT_TRACKING_INDEXTABLEUTILS_H_ */
