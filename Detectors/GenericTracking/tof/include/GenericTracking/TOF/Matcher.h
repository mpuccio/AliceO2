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

#ifndef GTRK_TOF_MATCHER_H_
#define GTRK_TOF_MATCHER_H_

#include "GenericTracking/TOF/Match.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include <functional>
#include <span>
#include <vector>

namespace o2::gtrk::tof
{
enum class Layer : std::uint8_t { Inner,
                                  Outer };

struct Hit {
  TrackingSurface surface;
  Timestamp time;
  float loc0 = 0.F;
  float loc1 = 0.F;
  float varLoc0 = 0.F;
  float varLoc1 = 0.F;
  float sigmaTimeNS = 0.020F;
  std::int32_t index = -1; ///< original cluster index, not the converted-hit row
  Layer layer = Layer::Inner;
};

struct Track {
  BoundTrackState state;
  BoundTrackState outwardState;
  TrackingSurface surface;
  TrackingSurface outwardSurface;
  Timestamp boxStart;
  double boxWidthNS = 0.;
  bool valid = false;
};

struct Prediction {
  float spatialChi2 = kUnset;
  std::array<double, kHypotheses> flightNS{};
  // Derivatives with respect to the five parameters of Track::state, with
  // the source frame held fixed. Shared errors cancel in outer-inner timing.
  std::array<std::array<double, kNBoundParams>, kHypotheses> derivative{};
  std::array<bool, kHypotheses> valid{};
};

struct Parameters {
  double maxSpatialChi2 = 25.;
  double maxTimeChi2 = 25.;
  double timeMarginNS = 0.1; ///< enlargement of the coarse tracking interval
  double maxStepCM = 2.;
  double maxPathCM = 1500.;
};

using Predictor = std::function<Prediction(const Track&, const Hit&)>;

/// Outer first; time-guided inner for outer-matched tracks second; remaining
/// tracks against unused inner clusters last. Deterministic greedy arbitration
/// within each pass, ordered by score, input track row, original cluster index.
/// Exhaustive reference implementation: no approximate spatial preselection.
std::vector<Match> match(std::span<const Track> tracks, std::span<const Hit> hits, const Predictor& predict, const Parameters& parameters = {});

} // namespace o2::gtrk::tof
#endif
