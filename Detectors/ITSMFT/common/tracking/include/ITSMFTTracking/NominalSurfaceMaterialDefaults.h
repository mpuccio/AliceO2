// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_NOMINALSURFACEMATERIALDEFAULTS_H_
#define ALICEO2_ITSMFT_TRACKING_NOMINALSURFACEMATERIALDEFAULTS_H_

#include <array>

#include "ITSMFTTracking/TrackingConfigParam.h"

// Common source for detector-default nominal xOverX0 per tracking surface,
// used by TrackingParameters and the ITS/MFT geometry catalogs. The MFT value
// is an intentional common-CA nominal model.
namespace o2::itsmft::tracking
{

static_assert(MFTNLayers % 2 == 0, "MFTNLayers must be even to derive MFTDisks");
/// Common-CA half-disk count, derived from MFTNLayers rather than
/// o2::mft::constants::mft::DisksNumber (numerically identical: 10 / 2 == 5).
inline constexpr int MFTDisks = MFTNLayers / 2;

/// Existing ITS detector-default per-layer xOverX0.
inline constexpr std::array<float, ITSNLayers> kNominalITSLayerX0{
  5.e-3f, 5.e-3f, 5.e-3f, 1.e-2f, 1.e-2f, 1.e-2f, 1.e-2f};

/// Existing MFT average material budget within acceptance, split evenly across
/// every tracking surface.
inline constexpr float kMFTNominalRadLength = 0.042f;

constexpr std::array<float, MFTNLayers> makeNominalMFTLayerX0()
{
  std::array<float, MFTNLayers> values{};
  for (auto& value : values) {
    value = kMFTNominalRadLength / static_cast<float>(MFTDisks);
  }
  return values;
}

inline constexpr std::array<float, MFTNLayers> kNominalMFTLayerX0 = makeNominalMFTLayerX0();

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_NOMINALSURFACEMATERIALDEFAULTS_H_ */
