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

/// \file KinkTrack.h
/// \brief
///

#ifndef _ALICEO2_KINKTRACK_
#define _ALICEO2_KINKTRACK_

#include <array>
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h" //(?)

using GIndex = o2::dataformats::VtxTrackIndex;

namespace o2
{
namespace dataformats
{

struct KinkTrack {
  unsigned int mITSRef = -1;
  int mNLayers = -1;
  GIndex mTrackIdx;
  float mITSClusSize;

  float mChi2Vertex = -1;
  float mChi2Match = -1;

  std::array<float, 3> mMotherP = {0, 0, 0};
  std::array<float, 3> mDaughterP = {0, 0, 0};
  std::array<float, 3> mDecayVtx;
  std::array<float, 2> mMasses; // Hypertriton and Sigma
};

} // namespace dataformats
} // namespace o2

#endif // _ALICEO2_KINKTRACK_