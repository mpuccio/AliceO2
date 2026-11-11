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

/// \file TrackingQA.h
/// \brief One tracker run's MC truth QA as trees, binned only by the plotting
///        macro.
///
/// Per timeframe, one record per findable particle -- by TrackComparison.h's
/// `isFindable`, its `findable`, or at the same pT floor by its MC hits'
/// outward leg on `minOutwardLayers` layers or more (`outwardLayers`) -- and
/// one per track:
///
///   particles: the truth's pT, eta and phi, its origin, `foundBy`, a bit
///              per iteration with a good track of it (`goodParticle`, the
///              text report's label rule), and `seededBy`, a bit per
///              iteration with a road of its own before the road fit
///              (StageFunnel.h); its PDG codes, decay geometry, layers,
///              the funnel's verdict and the other tracker's
///              (`QAParticle`)
///   tracks:    the fit's pT, eta and phi, its iteration, fake and clone;
///              fit status, generated pT/eta, origin and relative inverse-pT
///              residual for matched tracks, without requiring hit files
///
/// Every efficiency follows from the two; for iteration i, `bit` = 1 << i:
///
///   exclusive:   found by i                       (foundBy & bit) of all
///   cumulative:  found by some iteration <= i     (foundBy & (2 bit - 1)) of all
///   conditional: found by i, not found before i   (foundBy & bit) of !(foundBy & (bit - 1))
///
/// and the same for seeding with `seededBy`, the conditional one's condition
/// still "not found before i".
///
/// A fake is a track no particle is on every measurement of; a clone is a good
/// track of a particle with a good track before it, in iteration and then
/// track order.

#ifndef GTRK_STUDY_TRACKINGQA_H_
#define GTRK_STUDY_TRACKINGQA_H_

#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Study/TrackComparison.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <span>
#include <unordered_map>
#include <vector>

namespace o2::utils
{
class TreeStreamRedirector;
}

namespace o2::gtrk::study
{

/// The eta range the QA's axes span [low, high).
struct EtaRange {
  double low = 0.;
  double high = 0.;
};

/// The eta range a layout's silicon reaches from the origin: the extreme
/// eta of the envelopes' corners in r and z, rounded outwards to 0.5. The
/// comparison's eta edges are report bins, not an acceptance.
[[nodiscard]] EtaRange etaReach(std::span<const LayerEnvelope> envelopes);

/// \a phi [rad] in [0, 2 pi), the one wrapping of every phi the QA keeps:
/// in double for any number of turns, then float (O2's `bringTo02Pi` works in
/// float).
[[nodiscard]] inline float wrapPhi(double phi)
{
  constexpr double kTwoPi = 2. * std::numbers::pi;
  const auto wrapped = static_cast<float>(phi - (kTwoPi * std::floor(phi / kTwoPi)));
  return wrapped < static_cast<float>(kTwoPi) ? wrapped : 0.F;
}

/// A track as the QA reads it besides its rows: its fitted pT, eta and phi
/// at its innermost measurement, the iteration that found it, and its fit's
/// chi2 and degrees of freedom.
struct TrackKinematics {
  float pt = 0.F; ///< [GeV]
  float eta = 0.F;
  float phi = 0.F; ///< [rad]
  std::uint8_t iteration = 0U;
  float chi2 = 0.F;
  std::int32_t ndf = 0;
  bool fitValid = true; ///< the refit's propagation status is kOk
};

/// Iterations a particle's `foundBy` holds.
inline constexpr int kQAMaxIterations = std::numeric_limits<std::uint64_t>::digits;

/// The value of `QAParticle::lostAt` and `pairLoss` where nothing was
/// observed: no stage funnel, or a particle it does not hold.
inline constexpr std::uint8_t kQANotObserved = std::numeric_limits<std::uint8_t>::max();

/// What the sink knows of a findable particle beyond the truth table: its
/// decay geometry -- where it was
/// produced and how it passes the beam line, and its mother's direction
/// there -- and its stage funnel's verdict.
struct ParticleExtras {
  float vr = 0.F; ///< the production point's transverse distance from the beam line [cm]
  /// The truth helix's transverse distance of closest approach to the beam
  /// line, signed: the beam line's distance from the circle's centre less the
  /// radius; a straight line's, signed by the side, for a neutral or no field
  /// [cm]. NaN without the particle's kinematics.
  float d0 = std::numeric_limits<float>::quiet_NaN();
  float z0 = std::numeric_limits<float>::quiet_NaN(); ///< its z there less the beam spot's [cm]
  /// The angle between its momentum at production and its mother's at the
  /// same point [rad]; NaN without a mother.
  float openingAngle = std::numeric_limits<float>::quiet_NaN();
  /// How far the production point lies from the mother's circle,
  /// transversely [cm]: the check on the mother's helix; NaN without a
  /// mother, 0 for a neutral one's straight line through it.
  float motherMiss = std::numeric_limits<float>::quiet_NaN();
  /// The stage funnel's first stage without a chain of its own as long as the
  /// minimum, a FunnelStage (StageFunnel.h), kFound if found.
  std::uint8_t lostAt = kQANotObserved;
  /// Of its adjacent-layer pairs no iteration kept, the test most fail first,
  /// a PairLoss (StageFunnel.h), over the iterations whose pT_min it clears;
  /// kQANotObserved if none was lost.
  std::uint8_t pairLoss = kQANotObserved;
  /// Lost at the roads: why its chain of neighbours gave no road of its own,
  /// a RoadLoss (StageFunnel.h), in the last iteration that saw one;
  /// kQANotObserved otherwise.
  std::uint8_t roadLoss = kQANotObserved;
  /// Lost at the road fit: why its longest road of its own was no record, a
  /// FitLoss (StageFunnel.h), in the last iteration that fitted one, and that
  /// road's layers, a bit each; kQANotObserved and 0 otherwise.
  std::uint8_t fitLoss = kQANotObserved;
  std::uint64_t fitRoadLayers = 0U;
  /// The same particles: every reason any iteration gave, a bit per RoadLoss
  /// and per FitLoss; 0 otherwise.
  std::uint8_t roadLosses = 0U;
  std::uint8_t fitLosses = 0U;
  /// Lost at the roads to a chain raised from below it: the road fit of its
  /// own part, a FitLoss, kOther if the fit takes it; kQANotObserved
  /// otherwise.
  std::uint8_t raisedFit = kQANotObserved;
  std::uint8_t outwardLayers = 0U;      ///< the layers of its MC hits' outward leg (`outwardLayers`)
  std::uint64_t expectedLayerMask = 0U; ///< full-length reference: detector layers, or ALICE3 outward MC traversal
};

/// A findable particle of the `particles` tree.
struct QAParticle {
  float pt = 0.F; ///< [GeV]
  float eta = 0.F;
  float phi = 0.F;             ///< [rad], [0, 2 pi)
  std::uint8_t origin = 0U;    ///< ParticleOrigin
  std::uint64_t foundBy = 0U;  ///< a bit per iteration with a good track of it
  std::uint64_t seededBy = 0U; ///< a bit per iteration with a road of its own
  std::int32_t pdg = 0;
  std::int32_t motherPdg = 0; ///< 0 without a mother
  float vz = 0.F;             ///< the production point's z [cm]
  std::uint64_t layers = 0U;  ///< the layers it has measurements on, a bit each
  bool foundByO2 = false;     ///< a good track of it among the other tracker's
  bool findable = true;       ///< by `isFindable`; otherwise by its outward leg alone
  ParticleExtras extras;      ///< its decay geometry and funnel verdict, the tree's branches in between
};

/// One crossing of a sensor by a particle, from its MC hit [cm].
struct HitCrossing {
  double time = 0.; ///< any clock along the particle's path
  std::array<double, 3> position{};
  LayerID layer = 0;
  std::array<double, 3> momentum{}; ///< [GeV]
};

/// A particle's outward leg: its \a crossings in time order up to the first
/// one nearer its \a production point than the one before -- a curler's turn
/// back.
[[nodiscard]] std::vector<HitCrossing> outwardLeg(std::vector<HitCrossing> crossings, const std::array<double, 3>& production);

/// The layers of an outward \a leg (`outwardLeg`), a layer counted once
/// however often it is crossed.
[[nodiscard]] std::uint8_t outwardLayers(std::span<const HitCrossing> leg);

/// A particle's helix as the MC kinematics give it at its production point.
struct TruthHelix {
  std::array<double, 3> position{}; ///< [cm]
  std::array<double, 3> momentum{}; ///< [GeV]
  double charge = 0.;               ///< [e]
};

/// Where a helix passes an axis parallel to z: its transverse distance of
/// closest approach, signed as `QAParticle::d0`, and its z there [cm]; NaN
/// both for kinematics or a field not finite, or no transverse momentum.
struct HelixImpact {
  double d0 = 0.;
  double z = 0.;
};

/// \a helix in the uniform field \a bzKiloGauss against the axis through
/// (\a x, \a y): the circle's centre lies on the side the Lorentz force
/// points to, and of the closest points the nearer pass to the production
/// point is taken, the turn from it in (-pi, pi]. A neutral, or no field: the
/// straight line.
[[nodiscard]] HelixImpact helixImpact(const TruthHelix& helix, double bzKiloGauss, double x, double y);

/// The unit direction of \a helix at the transverse point (\a x, \a y): the
/// tangent to its circle there in the sense of motion, the point taken by
/// its direction from the centre; the same direction as at production for a
/// neutral, no field or no transverse momentum. NaN for kinematics, a field or
/// a point not finite, no momentum, or the point at the centre.
[[nodiscard]] std::array<double, 3> helixDirectionAt(const TruthHelix& helix, double bzKiloGauss, double x, double y);

/// The transverse distance of (\a x, \a y) from \a helix's circle, or its
/// straight line [cm], a line along z without transverse momentum; NaN for
/// kinematics, a field or a point not finite, or no momentum.
[[nodiscard]] double helixMiss(const TruthHelix& helix, double bzKiloGauss, double x, double y);

/// Each particle's extras, by particle id.
using ParticleExtrasTable = std::unordered_map<std::uint64_t, ParticleExtras>;

/// A track of the `tracks` tree.
struct QATrack {
  float pt = 0.F; ///< [GeV]
  float eta = 0.F;
  float phi = 0.F; ///< [rad], [0, 2 pi)
  std::uint8_t iteration = 0U;
  bool fake = false;
  bool clone = false;
  float chi2 = 0.F; ///< the fit's; TMath::Prob(chi2, ndf) is its probability
  std::int32_t ndf = 0;
  bool fitValid = true;
  std::uint64_t expectedLayerMask = 0U;                  ///< zero when a full-length reference is unavailable
  std::uint64_t layerMask = 0U;                          ///< distinct reconstructed detector layers
  std::uint64_t mcParticle = kNoTrackParticle;           ///< source/event/track identity, local to this timeframe
  float ptGen = std::numeric_limits<float>::quiet_NaN(); ///< generated momentum at production [GeV]
  float etaGen = std::numeric_limits<float>::quiet_NaN();
  std::uint8_t origin = 0xFFU;                                   ///< ParticleOrigin; unknown for unlabelled tracks
  float invPtResidual = std::numeric_limits<float>::quiet_NaN(); ///< (1/pT - 1/pTgen) * pTgen
};

/// A good track against its particle's MC crossing of its innermost
/// measurement's sensor, on the plane its state is bound to
/// (TrackingGeometry::planeOf): an entry of the `resolution` tree.
/// Residuals are fitted less true, in
/// the state's units (cm, rad, 1/GeV); sigmas the fit's. pT directly: its
/// relative residual, and its sigma from the theta-q/p block of the
/// covariance, J C J^T.
struct QAResolution {
  float pt = 0.F; ///< the truth's [GeV]
  float eta = 0.F;
  float phi = 0.F;        ///< [rad], [0, 2 pi)
  std::int8_t charge = 0; ///< [e]
  std::int32_t pdg = 0;
  std::uint8_t origin = 0U; ///< ParticleOrigin
  std::uint8_t iteration = 0U;
  std::uint8_t kinks = 0U;
  std::uint8_t clusters = 0U;
  std::uint8_t layer = 0U;         ///< the innermost measurement's, where the state is
  std::array<float, 5> residual{}; ///< loc0, loc1, phi, theta, q/p
  std::array<float, 5> sigma{};
  float fitPt = 0.F;      ///< [GeV]
  float ptResidual = 0.F; ///< (fitPt - pt) / pt
  float ptSigma = 0.F;    ///< of fitPt [GeV]
  /// The closest approach to the MC event's vertex, d0 and z0 [cm]:
  /// on the plane through the vertex
  /// along the beam, normal to the track's transverse direction there; the
  /// fit's sigma of each (J C J^T), and the scattering's between the vertex
  /// and the innermost measurement; NaN where the propagation failed.
  std::array<float, 2> dca{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  std::array<float, 2> dcaSigma{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  std::array<float, 2> dcaScattering{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  /// The particle's production point's distance from that vertex [cm]: 0
  /// for one made at the collision, a heavy-flavour decay's daughter's
  /// flight otherwise, whose d0 is its own and not a residual.
  float productionOffset = std::numeric_limits<float>::quiet_NaN();
};

/// A good track's oracle refit (SpeciesOracle.h) against
/// its particle's MC crossing of its innermost measurement's sensor: an
/// entry of the `oracle` tree, one per track and variant, as `QAResolution`'s
/// residuals and sigmas. \a momentum is the crossing's.
struct QAOracle {
  std::int32_t pdg = 0;
  std::int8_t charge = 0; ///< [e]
  std::uint8_t origin = 0U;
  float eta = 0.F; ///< the truth's at production
  std::uint8_t layer = 0U;
  std::uint8_t kinks = 0U; ///< the tracking's own refit's
  std::uint8_t clusters = 0U;
  float productionOffset = std::numeric_limits<float>::quiet_NaN();
  float momentum = 0.F; ///< at the crossing [GeV]
  /// At the crossing of the outermost measurement's sensor [GeV]: the
  /// momentum lost inside the fitted chain, 1 - pOuter / momentum.
  float pOuter = std::numeric_limits<float>::quiet_NaN();
  std::uint8_t variant = 0U; ///< OracleVariant
  std::uint8_t status = 0U;  ///< the refit's PropagationStatus; the rest only for kOk
  float chi2 = 0.F;
  std::int16_t ndf = 0;
  std::array<float, 5> residual{};
  std::array<float, 5> sigma{};
  float ptResidual = 0.F; ///< fitted pT over the crossing's, less 1
};

/// A findable particle's momentum along its outward leg: an entry of the
/// `momenta` tree, with --species-oracle, for every electron and every 16th
/// pion. \a pAtLayer: |p| entering the
/// layer's first crossing on the outward leg (`outwardLeg`), NaN
/// without one [GeV].
struct QAMomenta {
  std::int32_t pdg = 0;
  std::uint8_t origin = 0U;
  float pt = 0.F; ///< at production [GeV]
  float eta = 0.F;
  float p0 = 0.F; ///< |p| at production [GeV]
  bool found = false;
  std::uint8_t lostAt = 0U;   ///< QAParticle's
  std::uint8_t pairLoss = 0U; ///< QAParticle's
  std::uint64_t layers = 0U;  ///< the layers with a measurement of it
  /// Its fakes -- tracks most of whose measurements carry it, not all: the
  /// fewest wrong measurements among them, 0 without one; with one wrong, its
  /// layer (0xFF otherwise).
  std::uint8_t fakeWrong = 0U;
  std::uint8_t fakeWrongLayer = 0xFFU;
  std::vector<float> pAtLayer;
};

/// \a entries into \a output's `momenta` tree.
void streamMomenta(o2::utils::TreeStreamRedirector& output, std::span<const QAMomenta> entries);

/// \a entries into \a output's `oracle` tree.
void streamOracle(o2::utils::TreeStreamRedirector& output, std::span<const QAOracle> entries);

/// One timeframe's records.
struct TrackingQA {
  std::vector<QAParticle> particles;
  std::vector<QATrack> tracks; ///< in the order of the tracks given
};

/// Each particle's `seededBy`, by particle id (StageFunnel::seededBy).
using SeededBy = std::unordered_map<std::uint64_t, std::uint64_t>;

/// One timeframe: its rows' truth, its particles and its tracks, as rows
/// and as kinematics, track by track, the seeded, each particle's extras,
/// and the other tracker's tracks, if any. Throws std::invalid_argument for a
/// track's iteration outside \a iterations, or more than kQAMaxIterations.
[[nodiscard]] TrackingQA trackingQA(std::span<const TruthRow> rows, const ParticleTable& truth, const TrackRows& tracks, std::span<const TrackKinematics> kinematics, const SeededBy& seeded, const ParticleExtrasTable& extras, const TrackRows* theirs, int iterations, const ComparisonSettings& settings);

/// \a qa's records into \a output's `particles` and `tracks` trees.
void stream(o2::utils::TreeStreamRedirector& output, const TrackingQA& qa);

/// Write the same track schema for a comparison tracker, without duplicating
/// the generated-particle tree. Truth and clone rules are those of trackingQA().
void streamTracks(o2::utils::TreeStreamRedirector& output, std::span<const QATrack> tracks, const char* tree = "tracks");

/// \a entries into \a output's `resolution` tree.
void streamResolution(o2::utils::TreeStreamRedirector& output, std::span<const QAResolution> entries);

/// The run's one entry of the `run` tree: what the macro's axes need,
/// \a settings' findable rules (`minLayers`, `minPtGeV`, `minOutwardLayers`)
/// for its definitions, and the run's MC \a events. A nonzero fullLengthMask
/// requires all detector layers; zero denotes per-particle MC traversal masks.
void streamRun(o2::utils::TreeStreamRedirector& output, int iterations, EtaRange eta, const ComparisonSettings& settings, std::size_t events, std::uint64_t fullLengthMask = 0U);

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_TRACKINGQA_H_
