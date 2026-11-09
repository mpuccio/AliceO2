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

/// \file MCTruth.h
/// \brief The study's MC label rules: a primary's first outward pass, a
///        particle's origin, its collision, a triplet's common particle.

#ifndef GTRK_STUDY_MCTRUTH_H_
#define GTRK_STUDY_MCTRUTH_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Study/TrackComparison.h"

#include "SimulationDataFormat/DigitizationContext.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCUtils.h"

#include <TMCProcess.h>
#include <numbers>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace o2::gtrk::study
{

using workflow::carriesParticle;

/// Whether a primary's measurement lies on its first outward helix pass.
template <DetectorLayout Det>
[[nodiscard]] bool onFirstOutwardPass(const o2::MCTrack& mc, LayerID layer, Scalar phi, Scalar second, double bzKiloGauss)
{
  constexpr double kPhiTolerance = 0.3;
  constexpr double kZTolerance = 2.;
  if (!Detector<Det>::isBarrelLayer(layer)) {
    return true;
  }
  const double pt = mc.GetPt();
  const double radius = curvatureRadius(pt, bzKiloGauss);
  const double layerRadius = layerDistance<Det>(layer);
  if (!(2. * radius > layerRadius)) {
    return false;
  }
  const double deflection = std::asin(layerRadius / (2. * radius));
  const double cotTheta = mc.GetStartVertexMomentumZ() / pt;
  const double predictedZ = mc.Vz() + (cotTheta * 2. * radius * deflection);
  const double zTolerance = std::max(kZTolerance, std::numbers::pi_v<double> * radius * std::abs(cotTheta));
  if (std::abs(static_cast<double>(second) - predictedZ) > zTolerance) {
    return false;
  }
  const double phi0 = std::atan2(mc.GetStartVertexMomentumY(), mc.GetStartVertexMomentumX());
  const auto offset = [phi](double predicted) { return std::abs(std::remainder(static_cast<double>(phi) - predicted, 2 * std::numbers::pi_v<double>)); };
  return std::min(offset(phi0 + deflection), offset(phi0 - deflection)) < kPhiTolerance;
}

/// Where \a particle comes from, given whether it is a physical primary:
/// otherwise a decay's daughter if the generator made it -- feed-down it
/// decayed -- or transport made it in a decay, and the material's if neither.
[[nodiscard]] inline study::ParticleOrigin particleOrigin(const o2::MCTrack& particle, bool physicalPrimary)
{
  if (physicalPrimary) {
    return study::ParticleOrigin::kPhysicalPrimary;
  }
  if (particle.isPrimary() || particle.getProcess() == kPDecay) {
    return study::ParticleOrigin::kDecay;
  }
  return study::ParticleOrigin::kMaterial;
}

/// The same for \a particle, one of \a event's, by O2's physical-primary
/// test, which takes a transported particle's mother for granted: without
/// one, it is not a physical primary.
[[nodiscard]] inline study::ParticleOrigin particleOrigin(const o2::MCTrack& particle, const std::vector<o2::MCTrack>& event)
{
  using o2::mcutils::MCTrackNavigator;
  const bool navigable = particle.isPrimary() || MCTrackNavigator::getMother(particle, event) != nullptr;
  return particleOrigin(particle, navigable && MCTrackNavigator::isPhysicalPrimary(particle, event));
}

/// The collision of a packed particle identity (`study::TruthRow`'s
/// particles): its source and event, the bits above the track's.
[[nodiscard]] inline std::uint64_t collisionOf(std::uint64_t particle) { return particle >> o2::MCCompLabel::nbitsTrackID; }

/// Calls \a body for each particle of the first of \a held's rows that all
/// the other rows carry, in that first row's label order; whether there was one.
template <typename Body>
bool forEachCommonParticle(std::span<const TruthRow> rows, std::span<const std::uint32_t> held, const Body& body)
{
  if (held.empty()) {
    return false;
  }
  bool any = false;
  const auto& first = rows[held.front()];
  for (std::size_t at = 0; at < first.count; ++at) {
    const auto id = first.particles[at];
    if (std::all_of(held.begin() + 1, held.end(), [&](std::uint32_t row) { return rows[row].carries(id); })) {
      body(id);
      any = true;
    }
  }
  return any;
}

/// The particle a triplet's three rows \a triple of \a rows share: the first
/// of the middle row's labels the other two also carry; none if no label is
/// on all three.
[[nodiscard]] inline std::optional<std::uint64_t> commonParticle(std::span<const study::TruthRow> rows, const std::array<std::uint32_t, 3>& triple)
{
  std::optional<std::uint64_t> result;
  const std::array ordered{triple[1], triple[0], triple[2]};
  forEachCommonParticle(rows, ordered, [&](std::uint64_t particle) {
    if (!result) {
      result = particle;
    }
  });
  return result;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_MCTRUTH_H_
