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

#include "GenericTracking/TOF/Matcher.h"
#include "CommonConstants/LHCConstants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <unordered_set>

namespace o2::gtrk::tof
{
double differenceNS(Timestamp a, Timestamp b)
{
  const double bc = a.bc >= b.bc ? static_cast<double>(a.bc - b.bc) : -static_cast<double>(b.bc - a.bc);
  return bc * o2::constants::lhc::LHCBunchSpacingNS + static_cast<double>(a.fineNS) - b.fineNS;
}

Timestamp shifted(Timestamp time, double deltaNS)
{
  constexpr double spacing = o2::constants::lhc::LHCBunchSpacingNS;
  const double value = static_cast<double>(time.fineNS) + deltaNS;
  const double whole = std::floor(value / spacing);
  if (!std::isfinite(whole) || whole <= static_cast<double>(std::numeric_limits<std::int64_t>::min()) || whole >= static_cast<double>(std::numeric_limits<std::int64_t>::max())) {
    throw std::out_of_range("TOF time shift is not representable");
  }
  const auto shift = static_cast<std::int64_t>(whole);
  if ((shift < 0 && static_cast<std::uint64_t>(-shift) > time.bc) ||
      (shift > 0 && static_cast<std::uint64_t>(shift) > std::numeric_limits<std::uint64_t>::max() - time.bc)) {
    throw std::out_of_range("TOF time shift overflows the BC counter");
  }
  time.bc = shift < 0 ? time.bc - static_cast<std::uint64_t>(-shift) : time.bc + static_cast<std::uint64_t>(shift);
  time.fineNS = static_cast<float>(value - whole * spacing);
  if (time.fineNS >= spacing) {
    if (time.bc == std::numeric_limits<std::uint64_t>::max()) {
      throw std::out_of_range("TOF time rounding overflows the BC counter");
    }
    ++time.bc;
    time.fineNS = 0.F;
  }
  return time;
}

namespace
{
using Derivative = std::array<double, kNBoundParams>;
double covariance(const Track& track, const Derivative& a, const Derivative& b)
{
  double value = 0.;
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j < kNBoundParams; ++j) {
      value += a[i] * track.state.covAt(i, j) * b[j];
    }
  }
  return value;
}

struct Candidate {
  std::size_t track = 0;
  std::size_t hit = 0;
  double score = 0.;
  Prediction prediction;
  std::array<TimeHypothesis, kHypotheses> hypotheses{};
};

bool inBox(const Track& track, const Hit& hit, double flight, double margin)
{
  // Work relative to the measured timestamp; no absolute time is a double.
  const double time = differenceNS(hit.time, track.boxStart) - flight;
  return time >= -margin && time < track.boxWidthNS + margin;
}
} // namespace

std::vector<Match> match(std::span<const Track> tracks, std::span<const Hit> hits, const Predictor& predict, const Parameters& parameters)
{
  std::vector<std::size_t> rows(hits.size());
  std::iota(rows.begin(), rows.end(), 0);
  return matchSelected(tracks, hits, predict, [&](std::size_t) { return std::span<const std::size_t>(rows); }, parameters);
}

std::vector<Match> matchSelected(std::span<const Track> tracks, std::span<const Hit> hits, const Predictor& predict, const CandidateRows& rows, const Parameters& parameters)
{
  if (!(parameters.maxSpatialChi2 > 0.) || !(parameters.maxTimeChi2 > 0.) || !(parameters.timeMarginNS >= 0.)) {
    throw std::invalid_argument("invalid TOF matching cuts");
  }
  std::unordered_set<std::int32_t> identities;
  for (const auto& hit : hits) {
    if (hit.index < 0 || !identities.insert(hit.index).second || !(hit.sigmaTimeNS > 0.F) ||
        !std::isfinite(hit.sigmaTimeNS) || !(hit.time.fineNS >= 0.F) || !(hit.time.fineNS < o2::constants::lhc::LHCBunchSpacingNS)) {
      throw std::invalid_argument("invalid or duplicate TOF cluster identity/timing");
    }
  }
  std::vector<Match> result(tracks.size());
  std::vector<bool> used(hits.size(), false);
  std::vector<std::size_t> outerRow(tracks.size(), hits.size());
  std::vector<Prediction> outerPrediction(tracks.size());

  for (int pass = 0; pass < 3; ++pass) {
    std::vector<Candidate> candidates;
    for (std::size_t t = 0; t < tracks.size(); ++t) {
      const auto& track = tracks[t];
      if (!track.valid || !(track.boxWidthNS > 0.) || (pass == 1 && result[t].outer < 0) || (pass == 2 && result[t].outer >= 0)) {
        continue;
      }
      std::size_t previous = 0;
      bool first = true;
      for (const auto h : rows(t)) {
        if (!first && h <= previous) {
          throw std::invalid_argument("TOF candidate rows must be strictly increasing");
        }
        previous = h;
        first = false;
        if (h >= hits.size()) {
          throw std::out_of_range("TOF candidate hit row outside the input");
        }
        const auto& hit = hits[h];
        if (used[h] || hit.layer != (pass == 0 ? Layer::Outer : Layer::Inner)) {
          continue;
        }
        Candidate candidate{.track = t, .hit = h, .prediction = predict(track, hit)};
        const auto& prediction = candidate.prediction;
        if (!(prediction.spatialChi2 >= 0.F && prediction.spatialChi2 <= parameters.maxSpatialChi2)) {
          continue;
        }
        double bestTimeChi2 = std::numeric_limits<double>::infinity();
        bool accepted = false;
        for (std::size_t mass = 0; mass < kHypotheses; ++mass) {
          const double flight = prediction.flightNS[mass];
          if (!prediction.valid[mass] || !std::isfinite(flight) || !inBox(track, hit, flight, parameters.timeMarginNS)) {
            continue;
          }
          auto& hypothesis = candidate.hypotheses[mass];
          const auto& derivative = prediction.derivative[mass];
          const double variance = hit.sigmaTimeNS * hit.sigmaTimeNS + covariance(track, derivative, derivative);
          if (!(variance > 0.) || !std::isfinite(variance)) {
            continue;
          }
          try {
            hypothesis.time = shifted(hit.time, -flight);
          } catch (const std::out_of_range&) {
            continue; // Reference time before the start of the BC counter.
          }
          hypothesis.sigmaNS = static_cast<float>(std::sqrt(variance));
          hypothesis.valid = true;
          if (pass == 0) {
            hypothesis.outerFlightNS = static_cast<float>(flight);
          } else {
            hypothesis.innerFlightNS = static_cast<float>(flight);
          }
          if (pass == 1) {
            const auto& outer = result[t].hypotheses[mass];
            if (!outer.valid) {
              hypothesis.valid = false;
              continue;
            }
            const auto& outerHit = hits[outerRow[t]];
            const auto& other = outerPrediction[t];
            Derivative delta;
            for (int i = 0; i < kNBoundParams; ++i) {
              delta[i] = derivative[i] - other.derivative[mass][i];
            }
            // Evaluate the interval directly at the hits, avoiding rounding
            // two separately corrected reference timestamps before gating.
            const double residual = differenceNS(outerHit.time, hit.time) - (other.flightNS[mass] - flight);
            const double pairVariance = outerHit.sigmaTimeNS * outerHit.sigmaTimeNS + hit.sigmaTimeNS * hit.sigmaTimeNS + covariance(track, delta, delta);
            const double chi2 = residual * residual / pairVariance;
            if (!(pairVariance > 0.) || !(chi2 <= parameters.maxTimeChi2)) {
              hypothesis.valid = false;
              continue;
            }
            // Correlated BLUE: the same track's flight-time errors enter
            // both measurements. Their difference sets the inner weight.
            const double weight = (outerHit.sigmaTimeNS * outerHit.sigmaTimeNS - covariance(track, other.derivative[mass], delta)) / pairVariance;
            Derivative combined;
            for (int i = 0; i < kNBoundParams; ++i) {
              combined[i] = weight * derivative[i] + (1. - weight) * other.derivative[mass][i];
            }
            const double combinedVariance = weight * weight * hit.sigmaTimeNS * hit.sigmaTimeNS +
                                            (1. - weight) * (1. - weight) * outerHit.sigmaTimeNS * outerHit.sigmaTimeNS + covariance(track, combined, combined);
            if (!(combinedVariance > 0.) || !std::isfinite(combinedVariance)) {
              hypothesis.valid = false;
              continue;
            }
            // Anchor the combined estimate to the raw outer timestamp.
            try {
              hypothesis.time = shifted(outerHit.time, -other.flightNS[mass] - weight * residual);
            } catch (const std::out_of_range&) {
              hypothesis.valid = false;
              continue;
            }
            hypothesis.sigmaNS = static_cast<float>(std::sqrt(combinedVariance));
            hypothesis.outerFlightNS = static_cast<float>(other.flightNS[mass]);
            hypothesis.pairChi2 = static_cast<float>(chi2);
            bestTimeChi2 = std::min(bestTimeChi2, chi2);
          }
          accepted = accepted || hypothesis.valid;
        }
        if (!accepted) {
          continue;
        }
        candidate.score = prediction.spatialChi2 + (pass == 1 ? bestTimeChi2 : 0.);
        if (pass == 0) {
          ++result[t].outerCandidates;
        } else {
          ++result[t].innerCandidates;
        }
        candidates.push_back(candidate);
      }
    }
    std::sort(candidates.begin(), candidates.end(), [&](const Candidate& a, const Candidate& b) {
      return std::tuple{a.score, a.track, hits[a.hit].index} < std::tuple{b.score, b.track, hits[b.hit].index};
    });
    for (const auto& candidate : candidates) {
      auto& output = result[candidate.track];
      if (used[candidate.hit] || (pass == 0 ? output.outer >= 0 : output.inner >= 0)) {
        continue;
      }
      used[candidate.hit] = true;
      output.hypotheses = candidate.hypotheses;
      if (pass == 0) {
        output.outer = hits[candidate.hit].index;
        output.outerSpatialChi2 = candidate.prediction.spatialChi2;
        outerRow[candidate.track] = candidate.hit;
        outerPrediction[candidate.track] = candidate.prediction;
      } else {
        output.inner = hits[candidate.hit].index;
        output.innerSpatialChi2 = candidate.prediction.spatialChi2;
      }
    }
  }
  return result;
}
} // namespace o2::gtrk::tof
