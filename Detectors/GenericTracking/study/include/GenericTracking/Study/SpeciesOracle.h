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

/// \file SpeciesOracle.h
/// \brief The MC oracle refits: a good track's rows
///        fitted again with its particle's true mass, the scattering's and
///        the energy loss's apart, for the `oracle` QA tree.
///
/// The fit takes one mass hypothesis for both (`Propagator::massGeV`). An
/// `OraclePropagator` carries the energy loss's own (`lossMassGeV`);
/// `stepPropagator` hands it on whole, so the fit's steps see it. The
/// tracking never instantiates it.

#ifndef GTRK_STUDY_SPECIESORACLE_H_
#define GTRK_STUDY_SPECIESORACLE_H_

#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/TrackFitter.h"

#include "CommonConstants/PhysicsConstants.h"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <optional>

namespace o2::gtrk::study
{

/// A propagator whose energy loss takes \a lossMass, its scattering
/// `massGeV`.
struct OraclePropagator : Propagator {
  Scalar lossMass = static_cast<Scalar>(o2::constants::physics::MassPionCharged);
};

[[nodiscard]] inline Scalar lossMassGeV(const OraclePropagator& prop) { return prop.lossMass; }

/// The oracle's propagator for every step: itself, as `stepPropagator`
/// gives a `Propagator`.
template <typename Geometry, typename Chain>
[[nodiscard]] const OraclePropagator& stepPropagator(const OraclePropagator& prop, const Geometry& /*geometry*/, const Chain& /*measurements*/, int /*from*/, int /*to*/)
{
  return prop; // NOLINT(bugprone-return-const-ref-from-parameter)
}

/// The refits, each without kinks: the pion's two passes (the control); the
/// true mass's, for both, for the scattering alone, for the energy loss
/// alone; the true mass's inward pass alone from the track's own outward
/// state (the pion's), as the tracking's refit. Then the electron
/// hypothesis for every species, the inward pass from the same outward state.
enum class OracleVariant : std::uint8_t { kPion,
                                          kTrue,
                                          kTrueScattering,
                                          kTrueLoss,
                                          kTrueInward,
                                          kElectronInward,
                                          kCount };
inline constexpr std::size_t kOracleVariants = static_cast<std::size_t>(OracleVariant::kCount);

/// The mass of \a pdg's species [GeV]: e, mu, pi, K, p; none for another.
[[nodiscard]] inline std::optional<Scalar> oracleMass(std::int32_t pdg)
{
  namespace physics = o2::constants::physics;
  switch (std::abs(pdg)) {
    case 11:
      return static_cast<Scalar>(physics::MassElectron);
    case 13:
      return static_cast<Scalar>(physics::MassMuon);
    case 211:
      return static_cast<Scalar>(physics::MassPionCharged);
    case 321:
      return static_cast<Scalar>(physics::MassKaonCharged);
    case 2212:
      return static_cast<Scalar>(physics::MassProton);
    default:
      return std::nullopt;
  }
}

/// \a variant's refit of the \a n rows of \a chain, with \a base's field and
/// tolerances and \a mass, the particle's; \a outward is the track's outward
/// state for the inward-only variants.
template <DetectorLayout Det, typename Material, typename Chain>
[[nodiscard]] TrackFitResult oracleFit(OracleVariant variant, Scalar mass, const Propagator& base, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& chain, int n, const BoundTrackState& outward)
{
  const auto pion = static_cast<Scalar>(o2::constants::physics::MassPionCharged);
  const auto electron = static_cast<Scalar>(o2::constants::physics::MassElectron);
  const auto with = [&](Scalar scattering, Scalar loss) {
    OraclePropagator prop;
    static_cast<Propagator&>(prop) = base;
    prop.massGeV = scattering;
    prop.lossMass = loss;
    return prop;
  };
  switch (variant) {
    case OracleVariant::kPion:
      return fitMeasurements(with(pion, pion), geometry, material, chain, n);
    case OracleVariant::kTrue:
      return fitMeasurements(with(mass, mass), geometry, material, chain, n);
    case OracleVariant::kTrueScattering:
      return fitMeasurements(with(mass, pion), geometry, material, chain, n);
    case OracleVariant::kTrueLoss:
      return fitMeasurements(with(pion, mass), geometry, material, chain, n);
    case OracleVariant::kTrueInward:
      return fitInward(with(mass, mass), geometry, material, chain, n, outward);
    case OracleVariant::kElectronInward:
    case OracleVariant::kCount:
      break;
  }
  return fitInward(with(electron, electron), geometry, material, chain, n, outward);
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_SPECIESORACLE_H_
