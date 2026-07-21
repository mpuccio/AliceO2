// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#ifndef ALICEO2_ITSMFT_TRACKING_SURFACEKINEMATICSTATEADAPTERS_H_
#define ALICEO2_ITSMFT_TRACKING_SURFACEKINEMATICSTATEADAPTERS_H_

#include <cmath>

#include "ITSMFTTracking/SurfaceKinematicState.h"
#include "ReconstructionDataFormats/TrackFwd.h"
#include "ReconstructionDataFormats/Track.h"

namespace o2::itsmft::tracking
{

// Host-only compatibility adapters. They only copy same-family representations;
// they never propagate or convert between barrel and forward parameterizations.
inline bool fromBarrelTrackParCov(const o2::track::TrackParCovF& source, SurfaceKinematicState& destination) noexcept
{
  SurfaceKinematicState scratch{};
  scratch.referenceCoordinate = source.getX();
  scratch.alpha = source.getAlpha();
  for (uint8_t i = 0; i < 5; ++i) {
    scratch.parameters[i] = source.getParam(i);
  }
  for (uint8_t i = 0; i < 15; ++i) {
    scratch.covariance[i] = source.getCov()[i];
  }
  scratch.family = StateFamily::Barrel;
  scratch.absCharge = static_cast<uint8_t>(source.getAbsCharge());
  scratch.pid = source.getPID();
  destination = scratch;
  return true;
}

inline bool toBarrelTrackParCov(const SurfaceKinematicState& source, o2::track::TrackParCovF& destination) noexcept
{
  if (source.family != StateFamily::Barrel) {
    return false;
  }
  o2::track::TrackParCovF::params_t parameters{};
  o2::track::TrackParCovF::covMat_t covariance{};
  for (uint8_t i = 0; i < 5; ++i) {
    parameters[i] = source.parameters[i];
  }
  for (uint8_t i = 0; i < 15; ++i) {
    covariance[i] = source.covariance[i];
  }
  const o2::track::TrackParCovF scratch{source.referenceCoordinate, source.alpha, parameters, covariance, source.absCharge, source.pid};
  destination = scratch;
  return true;
}

inline bool fromForwardTrackParCov(const o2::track::TrackParCovFwd& source, SurfaceKinematicState& destination) noexcept
{
  SurfaceKinematicState scratch{};
  const auto& covariance = source.getCovariances();
  const double parameters[] = {source.getX(), source.getY(), source.getPhi(), source.getTanl(), source.getInvQPt()};
  if (!std::isfinite(source.getZ())) {
    return false;
  }
  for (uint8_t i = 0; i < 5; ++i) {
    if (!std::isfinite(parameters[i])) {
      return false;
    }
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      if (!std::isfinite(covariance(row, column))) {
        return false;
      }
    }
  }
  scratch.referenceCoordinate = static_cast<float>(source.getZ());
  scratch.alpha = 0.f;
  for (uint8_t i = 0; i < 5; ++i) {
    scratch.parameters[i] = static_cast<float>(parameters[i]);
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      scratch.covariance[packedCovarianceIndex(row, column)] = static_cast<float>(covariance(row, column));
    }
  }
  scratch.family = StateFamily::Forward;
  scratch.absCharge = 1;
  scratch.pid = o2::track::PID::Pion;
  destination = scratch;
  return true;
}

inline bool toForwardTrackParCov(const SurfaceKinematicState& source, o2::track::TrackParCovFwd& destination) noexcept
{
  if (source.family != StateFamily::Forward) {
    return false;
  }
  o2::track::SMatrix5 parameters{};
  o2::track::SMatrix55Sym covariance{};
  for (uint8_t i = 0; i < 5; ++i) {
    parameters(i) = static_cast<double>(source.parameters[i]);
  }
  for (uint8_t row = 0; row < 5; ++row) {
    for (uint8_t column = 0; column <= row; ++column) {
      covariance(row, column) = static_cast<double>(source.covariance[packedCovarianceIndex(row, column)]);
    }
  }
  const o2::track::TrackParCovFwd scratch{static_cast<double>(source.referenceCoordinate), parameters, covariance, 0.};
  destination = scratch;
  return true;
}

} // namespace o2::itsmft::tracking

#endif // ALICEO2_ITSMFT_TRACKING_SURFACEKINEMATICSTATEADAPTERS_H_
