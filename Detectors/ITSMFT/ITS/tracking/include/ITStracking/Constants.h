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
/// \file Constants.h
/// \brief
///

#ifndef TRACKINGITSU_INCLUDE_CONSTANTS_H_
#define TRACKINGITSU_INCLUDE_CONSTANTS_H_

#ifndef GPUCA_GPUCODE_DEVICE
#include <climits>
#include <vector>
#endif

#include "ITStracking/Definitions.h"

namespace o2
{
namespace its
{

namespace constants
{

constexpr bool DoTimeBenchmarks = true;

namespace math
{
constexpr float Pi{3.14159265359f};
constexpr float TwoPi{2.0f * Pi};
constexpr float FloatMinThreshold{1e-20f};
} // namespace math

namespace its
{
constexpr int LayersNumberVertexer{3};
constexpr int ClustersPerCell{3};
constexpr int UnusedIndex{-1};
constexpr float Resolution{0.0005f};

GPU_HOST_DEVICE constexpr GPUArray<float, 3> VertexerHistogramVolume()
{
  return GPUArray<float, 3>{{1.98, 1.98, 40.f}};
}
} // namespace its

namespace pdgcodes
{
constexpr int PionCode{211};
}
} // namespace constants
} // namespace its
} // namespace o2

#endif /* TRACKINGITSU_INCLUDE_CONSTANTS_H_ */
