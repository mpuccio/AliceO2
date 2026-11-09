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

/// \file Validation.h
/// \brief The point-only overloads of the triplet and scattering functions,
///        and the stage rules, that the studies call.
///
/// What a stage has at hand -- the kink terms, the scattering -- computed from
/// the points alone, the link and chi2 rules the stages apply, and the road
/// fit stage's tests on one road. Host code; the tests' own overloads are
/// test/Oracles.h's.

#ifndef GTRK_STUDY_VALIDATION_H_
#define GTRK_STUDY_VALIDATION_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/NeighbourCut.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackFitter.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/TripletCut.h"
#include "GenericTracking/Tracking/TripletFit.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>

namespace o2::gtrk::study
{

// ---- the triplet and neighbour fits from the points alone ----

/// The triplet's kink variance with theta from the first to the third
/// measurement.
template <DetectorLayout Det>
[[nodiscard]] inline Scalar tripletScattering(const ScatteringTable<Det>& table, const PairCut& cut, LayerID middleLayer, const PairPoint& first, const PairPoint& middle, const PairPoint& third)
{
  return tripletScattering(table, cut, middleLayer, first, middle, third, detail::referenceSinTheta(first, third));
}

/// The fit factor from the points and a given angular variance, the kink
/// terms computed here.
[[nodiscard]] inline bool makeTripletFitFactor(const PairPoint& first, const PairPoint& middle, const PairPoint& third, Scalar angularVariance, TripletFitFactor& result)
{
  KinkTerms terms;
  return makeKinkTerms(first, middle, third, terms) && makeTripletFitFactor(terms, first, middle, third, angularVariance, result);
}

/// Whether the iteration takes \a link: it is among its source's out-links,
/// which hold the links `takenLinks` gives: the skip rule, the layer mask and
/// the path length.
template <DetectorLayout Det>
[[nodiscard]] inline bool linkTaken(const IterationParameters<Det>& parameters, std::size_t link)
{
  const auto source = static_cast<std::size_t>(parameters.linkSources[link]);
  for (auto at = parameters.outLinkOffsets[source]; at < parameters.outLinkOffsets[source + 1U]; ++at) {
    if (parameters.outLinkIndices[at] == link) {
      return true;
    }
  }
  return false;
}

/// The road fit's chi2 bound for \a clusters clusters and \a kinks kinks.
template <DetectorLayout Det>
[[nodiscard]] inline Scalar roadChi2Cut(const IterationParameters<Det>& parameters, std::uint32_t clusters, std::uint32_t kinks = 0U)
{
  return parameters.roadChi2MaxByNdf[(2U * clusters) - kNBoundParams - (2U * kinks)];
}

/// Why the road fit stage refuses a road, in the order of its tests, or that
/// it takes it.
enum class FitRefusal : std::uint8_t { kTime,
                                       kChi2,
                                       kPropagation,
                                       kPtFloor,
                                       kBeamLine,
                                       kTaken };
inline constexpr std::size_t kFitRefusals = 6U;
inline constexpr std::array<std::string_view, kFitRefusals> kFitRefusalNames{"time", "chi2", "propagation", "pT floor", "beam line", "taken"};

/// The road fit (`o2::gtrk::fitRoad`) with the host's choice at run time,
/// where the iteration has the allowance (`takesLastStepAllowance`), for the
/// studies that repeat it; the stages choose at compile time.
template <typename Parameters, typename Prop, DetectorLayout Det, typename Material, typename Chain>
[[nodiscard]] RoadFitResult fitRoad(const Parameters& parameters, const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n)
{
  return takesLastStepAllowance(parameters) ? o2::gtrk::fitRoad<true>(parameters, prop, geometry, material, measurements, n) : o2::gtrk::fitRoad<false>(parameters, prop, geometry, material, measurements, n);
}

/// The road fit stage's tests on the road of \a rows of \a frame, innermost
/// first, in its order (RoadFit.h, `fitRows`): its boxes share no time, the
/// outward fit's chi2 or propagation (the run-time `fitRoad`, the last-step
/// allowance as the stage's), the pT floor for its length, and with the beam
/// constraint the beam line. The road as given: `iterationInnerTrim`'s
/// shorter fits of a refused road are not replayed.
template <DetectorLayout Det, typename Material>
[[nodiscard]] FitRefusal roadFitRefusal(const IterationParameters<Det>& parameters, const Propagator& propagator, const TrackingGeometry<Det>& tracking, const Material& material, const TimeFrame<Det>& frame, std::span<const std::uint32_t> rows)
{
  const auto columns = frame.measurements();
  auto common = timeBoxAt(columns, rows[0]);
  for (std::size_t i = 1; i < rows.size(); ++i) {
    const auto box = timeBoxAt(columns, rows[i]);
    if (!common.isCompatible(box)) {
      return FitRefusal::kTime;
    }
    common += box;
  }
  const auto count = static_cast<std::uint32_t>(rows.size());
  const auto fit = fitRoad(parameters, propagator, tracking, material, RowChain<Det>{.frame = &frame, .rows = rows.data()}, static_cast<int>(count));
  if (fit.status == RoadFitStatus::kRejected) {
    return FitRefusal::kChi2;
  }
  if (fit.status != RoadFitStatus::kAccepted) {
    return FitRefusal::kPropagation;
  }
  if (!clearsPtFloor(parameters, count, fit)) {
    return FitRefusal::kPtFloor;
  }
  if (Scalar beamChi2 = 0.F; parameters.beamConstraint && !clearsBeamLine(propagator, parameters, frame, tracking, material, rows.data(), count, fit.state, beamChi2)) {
    return FitRefusal::kBeamLine;
  }
  return FitRefusal::kTaken;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_VALIDATION_H_
