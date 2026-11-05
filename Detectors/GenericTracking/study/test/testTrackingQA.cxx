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

/// \file testTrackingQA.cxx
/// \brief The QA's records per timeframe (TrackingQA.h), against the text
///        score (TrackComparison.h).

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK tracking QA
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Study/TrackingQA.h"
#include "GenericTracking/Study/TrackComparison.h"

#include "CommonConstants/MathConstants.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <numbers>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

namespace
{

TruthRow row(LayerID layer, std::initializer_list<std::uint64_t> particles)
{
  TruthRow result;
  result.layer = layer;
  for (const auto particle : particles) {
    result.particles[result.count++] = particle;
  }
  return result;
}

/// Three layers, two adjacent ones findable, two iterations.
///   particle 1: a primary on every layer, found by iteration 0 and again by 1
///   particle 2: a decay daughter on every layer, found by iteration 1 only
///   particle 3: a primary on layers 0 and 2 only, not findable
///   particle 4: from the material on layers 1 and 2, findable, not found
///   particle 5: a primary on every layer below the pT floor, not findable
struct Toy {
  std::vector<TruthRow> rows{row(0, {1}), row(1, {1}), row(2, {1}), row(0, {2}), row(1, {2}), row(2, {2}), row(0, {3}),
                             row(2, {3}), row(1, {4}), row(2, {4}), row(0, {5}), row(1, {5}), row(2, {5})};
  ParticleTable truth{{1, {.pt = 1.F, .primary = true, .eta = 0.2F, .phi = 1.F, .origin = ParticleOrigin::kPhysicalPrimary}},
                      {2, {.pt = 0.5F, .eta = -0.3F, .phi = -1.F, .origin = ParticleOrigin::kDecay}},
                      {3, {.pt = 2.F, .primary = true, .eta = 0.4F, .origin = ParticleOrigin::kPhysicalPrimary}},
                      {4, {.pt = 0.3F, .eta = 0.9F, .phi = 3.F, .origin = ParticleOrigin::kMaterial}},
                      {5, {.pt = 0.05F, .primary = true, .eta = 0.F, .origin = ParticleOrigin::kPhysicalPrimary}}};
  ComparisonSettings settings{.nLayers = 3, .minLayers = 2, .minPtGeV = 0.1, .chainLinks = {{.from = 0, .to = 1}, {.from = 1, .to = 2}}};
  // Particle 1 in iteration 1, particle 2, a fake of particles 4 and 2, and
  // particle 1 in iteration 0.
  TrackRows tracks{{0, 1, 2}, {3, 4, 5}, {8, 5}, {0, 1, 2}};
  std::vector<TrackKinematics> kinematics{{.pt = 1.F, .eta = 0.2F, .phi = 1.F, .iteration = 1},
                                          {.pt = 0.5F, .eta = -0.3F, .phi = -1.F, .iteration = 1},
                                          {.pt = 0.4F, .eta = 0.5F, .phi = 2.F, .iteration = 0},
                                          {.pt = 1.F, .eta = 0.2F, .phi = 1.F, .iteration = 0, .chi2 = 3.5F, .ndf = 4}};
  // Particle 1 seeded by both, 2 by both, 4 by iteration 0; 3 is not findable.
  SeededBy seeded{{1, 0b11U}, {2, 0b11U}, {3, 0b01U}, {4, 0b01U}};
};

/// \a helix carried along its path by the Lorentz force alone, in \a steps
/// steps of \a step [cm] of transverse path, backward for a negative one:
/// the field's equations of motion, midpoint rule, with no circle in them.
TruthHelix integrate(TruthHelix helix, double bzKiloGauss, double step, int steps)
{
  const double pt = std::hypot(helix.momentum[0], helix.momentum[1]);
  // du/ds = (q Bz c / pT) (uy, -ux): F = q v x B with B along z.
  const double k = helix.charge * bzKiloGauss * 0.299792458e-3 / pt;
  double ux = helix.momentum[0] / pt;
  double uy = helix.momentum[1] / pt;
  for (int i = 0; i < steps; ++i) {
    const double mx = ux + (0.5 * step * k * uy);
    const double my = uy - (0.5 * step * k * ux);
    helix.position[0] += step * mx;
    helix.position[1] += step * my;
    helix.position[2] += step * helix.momentum[2] / pt;
    const double nx = ux + (step * k * my);
    const double ny = uy - (step * k * mx);
    const double norm = std::hypot(nx, ny);
    ux = nx / norm;
    uy = ny / norm;
  }
  helix.momentum[0] = ux * pt;
  helix.momentum[1] = uy * pt;
  return helix;
}

const QAParticle* particleOf(const TrackingQA& qa, float pt)
{
  const auto found = std::ranges::find(qa.particles, pt, &QAParticle::pt);
  return found == qa.particles.end() ? nullptr : &*found;
}

} // namespace

BOOST_AUTO_TEST_CASE(TheFindableAndTheIterationsThatFoundAndSeededThem)
{
  const Toy toy;
  const auto qa = trackingQA(toy.rows, toy.truth, toy.tracks, toy.kinematics, toy.seeded, {}, nullptr, 2, toy.settings);
  BOOST_TEST(qa.particles.size() == 3U); // 1, 2 and 4
  const auto* first = particleOf(qa, 1.F);
  const auto* decay = particleOf(qa, 0.5F);
  const auto* material = particleOf(qa, 0.3F);
  BOOST_REQUIRE(first != nullptr);
  BOOST_REQUIRE(decay != nullptr);
  BOOST_REQUIRE(material != nullptr);
  BOOST_TEST(first->foundBy == 0b11U);
  BOOST_TEST(first->origin == static_cast<std::uint8_t>(ParticleOrigin::kPhysicalPrimary));
  BOOST_TEST(decay->foundBy == 0b10U);
  BOOST_TEST(decay->origin == static_cast<std::uint8_t>(ParticleOrigin::kDecay));
  BOOST_TEST(material->foundBy == 0U);
  BOOST_TEST(first->seededBy == 0b11U);
  BOOST_TEST(decay->seededBy == 0b11U);
  BOOST_TEST(material->seededBy == 0b01U);
  BOOST_TEST(material->origin == static_cast<std::uint8_t>(ParticleOrigin::kMaterial));
}

BOOST_AUTO_TEST_CASE(FakesAndClonesInIterationOrder)
{
  const Toy toy;
  const auto qa = trackingQA(toy.rows, toy.truth, toy.tracks, toy.kinematics, toy.seeded, {}, nullptr, 2, toy.settings);
  BOOST_REQUIRE(qa.tracks.size() == 4U);
  // In the order given; particle 1's track of iteration 1 is the clone,
  // although it comes first.
  BOOST_TEST((!qa.tracks[0].fake && qa.tracks[0].clone));
  BOOST_TEST((!qa.tracks[1].fake && !qa.tracks[1].clone));
  BOOST_TEST((qa.tracks[2].fake && !qa.tracks[2].clone));
  BOOST_TEST((!qa.tracks[3].fake && !qa.tracks[3].clone));
  BOOST_TEST(qa.tracks[2].iteration == 0U);
  BOOST_TEST(qa.tracks[3].chi2 == 3.5F);
  BOOST_TEST(qa.tracks[3].ndf == 4);
}

BOOST_AUTO_TEST_CASE(EveryPhiInZeroToTwoPi)
{
  const Toy toy;
  const auto qa = trackingQA(toy.rows, toy.truth, toy.tracks, toy.kinematics, toy.seeded, {}, nullptr, 2, toy.settings);
  BOOST_TEST(particleOf(qa, 0.5F)->phi == wrapPhi(-1.));
  BOOST_TEST(qa.tracks[1].phi == wrapPhi(-1.));
  constexpr double kTwoPi = 2. * std::numbers::pi;
  BOOST_TEST(wrapPhi(-1.) == static_cast<float>(kTwoPi - 1.));
  BOOST_TEST(wrapPhi(kTwoPi) == 0.F);
  BOOST_TEST(wrapPhi(std::nextafter(kTwoPi, 0.)) == 0.F); // 2 pi as a float
  BOOST_TEST(wrapPhi(-std::numbers::pi) == static_cast<float>(std::numbers::pi));
  BOOST_TEST(wrapPhi(7.) == static_cast<float>(7. - kTwoPi));
}

BOOST_AUTO_TEST_CASE(TheTotalsAreTheTextScore)
{
  const Toy toy;
  const auto qa = trackingQA(toy.rows, toy.truth, toy.tracks, toy.kinematics, toy.seeded, {}, nullptr, 2, toy.settings);
  const auto score = scoreTracks(toy.rows, toy.truth, toy.tracks, toy.settings);
  for (std::size_t origin = 0; origin < score.originFindable.size(); ++origin) {
    const auto of = [&](const QAParticle& p) { return p.origin == origin; };
    const auto findable = std::ranges::count_if(qa.particles, of);
    const auto found = std::ranges::count_if(qa.particles, [&](const QAParticle& p) { return of(p) && p.foundBy != 0U; });
    BOOST_TEST(static_cast<std::uint64_t>(findable) == score.originFindable[origin].findable);
    BOOST_TEST(static_cast<std::uint64_t>(found) == score.originFindable[origin].found);
  }
  // The primaries' two classes split the findable primaries.
  BOOST_TEST(score.primaryRegion[0].findable + score.primaryRegion[1].findable == score.originFindable[0].findable);
  BOOST_TEST(score.primaryRegion[0].found + score.primaryRegion[1].found == score.originFindable[0].found);
  BOOST_TEST(qa.tracks.size() == score.tracks);
  BOOST_TEST(static_cast<std::uint64_t>(std::ranges::count_if(qa.tracks, &QATrack::fake)) == score.fakes);
}

BOOST_AUTO_TEST_CASE(TheEtaReachOfEnvelopes)
{
  // A barrel layer at r 2.3 cm, |z| < 13.5 cm, and a disc at z 77 cm over r 3.4-15 cm.
  const std::array<LayerEnvelope, 2> envelopes{LayerEnvelope{.position = 2.3, .rMin = 2.3, .rMax = 2.4, .zMin = -13.5, .zMax = 13.5},
                                               LayerEnvelope{.position = 77., .rMin = 3.4, .rMax = 15., .zMin = 76.9, .zMax = 77.1}};
  const auto reach = etaReach(envelopes);
  BOOST_TEST(reach.low == -2.5); // asinh(-13.5 / 2.3) = -2.47
  BOOST_TEST(reach.high == 4.0); // asinh(77.1 / 3.4) = 3.81
}

BOOST_AUTO_TEST_CASE(TheExtrasLayersAndTheOtherTrackersParticles)
{
  const Toy toy;
  const ParticleExtrasTable extras{{2, {.vr = 1.5F, .d0 = 0.25F, .z0 = -0.5F, .openingAngle = 0.1F, .motherMiss = 0.01F, .lostAt = 6U, .raisedFit = 5U}},
                                   {4, {.vr = 3.F, .lostAt = 4U, .pairLoss = 6U, .fitLoss = 1U, .fitRoadLayers = 0b110U, .fitLosses = 0b1010U}}};
  const TrackRows theirs{{3, 4, 5}, {8, 9, 0}}; // particle 2, and a fake
  const auto qa = trackingQA(toy.rows, toy.truth, toy.tracks, toy.kinematics, toy.seeded, extras, &theirs, 2, toy.settings);
  const auto* first = particleOf(qa, 1.F);
  const auto* decay = particleOf(qa, 0.5F);
  const auto* material = particleOf(qa, 0.3F);
  BOOST_REQUIRE((first != nullptr && decay != nullptr && material != nullptr));
  BOOST_TEST(first->layers == 0b111U);
  BOOST_TEST(material->layers == 0b110U);
  BOOST_TEST(decay->extras.d0 == 0.25F);

  BOOST_TEST(decay->foundByO2);
  BOOST_TEST(!first->foundByO2);
  BOOST_TEST(!material->foundByO2);
  // Without extras: nothing observed, the geometry NaN.
  BOOST_TEST(first->extras.lostAt == kQANotObserved);
  BOOST_TEST(first->extras.pairLoss == kQANotObserved);
  BOOST_TEST(std::isnan(first->extras.d0));
  BOOST_TEST(std::isnan(first->extras.openingAngle));
}

BOOST_AUTO_TEST_CASE(TheOutwardLegEndsWhereTheCurlerTurnsBack)
{
  // Out through layers 0, 1 (twice: an overlap) and 2, back through 1, then on
  // layer 3 inside layer 2's radius: the return's layers do not count. Given
  // out of time order.
  const std::vector<HitCrossing> crossings{{.time = 5., .position = {1.5, 0., 0.}, .layer = 3}, {.time = 1., .position = {1., 0., 0.}, .layer = 0}, {.time = 4., .position = {0., 2., 0.}, .layer = 1}, {.time = 2., .position = {2., 0., 0.}, .layer = 1}, {.time = 2.5, .position = {2.01, 0.1, 0.}, .layer = 1}, {.time = 3., .position = {3., 0., 0.}, .layer = 2}};
  BOOST_TEST(outwardLayers(outwardLeg(crossings, {0., 0., 0.})) == 3U);
  // Seen from a production point beyond them, the leg ends at its first crossing.
  BOOST_TEST(outwardLayers(outwardLeg(crossings, {10., 0., 0.})) == 1U);
}

BOOST_AUTO_TEST_CASE(TheOutwardLegMakesAParticleFindableAboveTheFloor)
{
  const Toy toy;
  // Particle 3 (no chain) and particle 5 (below the floor) on the outward rule's layers.
  const ParticleExtrasTable extras{{3, {.outwardLayers = 7U}}, {5, {.outwardLayers = 7U}}, {4, {.outwardLayers = 2U}}};
  const auto qa = trackingQA(toy.rows, toy.truth, toy.tracks, toy.kinematics, toy.seeded, extras, nullptr, 2, toy.settings);
  BOOST_TEST(qa.particles.size() == 4U); // 1, 2, 4 and 3
  const auto* outward = particleOf(qa, 2.F);
  const auto* material = particleOf(qa, 0.3F);
  BOOST_REQUIRE((outward != nullptr && material != nullptr));
  BOOST_TEST(!outward->findable);
  BOOST_TEST(outward->extras.outwardLayers == 7U);
  BOOST_TEST(material->findable);
  BOOST_TEST(material->extras.outwardLayers == 2U);
  BOOST_TEST(particleOf(qa, 0.05F) == nullptr);
}

BOOST_AUTO_TEST_CASE(AParticleFromTheBeamLineMissesItByNothing)
{
  for (const double charge : {-1., 1.}) {
    for (const double bz : {-5., 5.}) {
      const TruthHelix helix{.position = {0.03, -0.02, 2.5}, .momentum = {-0.2, 0.15, 0.4}, .charge = charge};
      const auto impact = helixImpact(helix, bz, 0.03, -0.02);
      BOOST_TEST(std::abs(impact.d0) < 1e-12);
      BOOST_TEST(std::abs(impact.z - 2.5) < 1e-12);
      const auto direction = helixDirectionAt(helix, bz, 0.03, -0.02);
      BOOST_TEST(std::abs(direction[0] - (-0.2 / std::hypot(0.2, 0.15, 0.4))) < 1e-12);
      BOOST_TEST(std::abs(direction[2] - (0.4 / std::hypot(0.2, 0.15, 0.4))) < 1e-12);
    }
  }
}

BOOST_AUTO_TEST_CASE(ADaughterOnItsMothersHelixPointsBackAlongIt)
{
  // From the beam line, 10 cm of transverse path by the equations of motion:
  // the helix there passes the beam line where it started, and the helix at
  // the start has there the direction the motion gives.
  for (const double charge : {-1., 1.}) {
    for (const double bz : {-5., 5.}) {
      const TruthHelix start{.position = {0., 0., 0.5}, .momentum = {0.3, 0.1, 0.2}, .charge = charge};
      const auto later = integrate(start, bz, 1e-3, 10000);
      const auto impact = helixImpact(later, bz, 0., 0.);
      BOOST_TEST(std::abs(impact.d0) < 1e-5);
      BOOST_TEST(std::abs(impact.z - 0.5) < 1e-5);
      BOOST_TEST(helixMiss(start, bz, later.position[0], later.position[1]) < 1e-6);
      const auto direction = helixDirectionAt(start, bz, later.position[0], later.position[1]);
      const double p = std::hypot(later.momentum[0], later.momentum[1], later.momentum[2]);
      for (std::size_t i = 0; i < 3U; ++i) {
        BOOST_TEST(std::abs(direction[i] - (later.momentum[i] / p)) < 1e-6);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ADisplacedHelixMissesTheBeamLineByItsCircle)
{
  // A circle about (3, 4), 5 cm from the beam line: radius 2 cm misses it by
  // +3, 7 cm by -2. The closest point's z from the motion run backward.
  for (const double charge : {-1., 1.}) {
    for (const double bz : {-5., 5.}) {
      for (const double radius : {2., 7.}) {
        const double pt = radius * std::abs(charge * bz) * 0.299792458e-3;
        const double angle = 0.3;
        // Clockwise for q Bz > 0: the direction is the radius turned by -90 degrees.
        const double sense = charge * bz > 0. ? -1. : 1.;
        const TruthHelix helix{.position = {3. + (radius * std::cos(angle)), 4. + (radius * std::sin(angle)), 1.},
                               .momentum = {-sense * pt * std::sin(angle), sense * pt * std::cos(angle), 0.5 * pt},
                               .charge = charge};
        // The motion keeps to that circle: the setup's sense is the field's.
        const auto ahead = integrate(helix, bz, 1e-3, 3000);
        BOOST_TEST(std::abs(std::hypot(ahead.position[0] - 3., ahead.position[1] - 4.) - radius) < 1e-6);
        const auto impact = helixImpact(helix, bz, 0., 0.);
        // The radius from O2's B2C, a float: 1e-7 of it.
        BOOST_TEST(std::abs(impact.d0 - (5. - radius)) < 1e-5);
        // The nearer pass: half a turn back and half a turn ahead, step by step.
        double best = std::hypot(helix.position[0], helix.position[1]);
        double bestZ = helix.position[2];
        const int halfTurn = static_cast<int>(std::numbers::pi * radius / 1e-3);
        for (const double step : {-1e-3, 1e-3}) {
          auto walker = helix;
          for (int i = 0; i < halfTurn; ++i) {
            walker = integrate(walker, bz, step, 1);
            if (const double distance = std::hypot(walker.position[0], walker.position[1]); distance < best) {
              best = distance;
              bestZ = walker.position[2];
            }
          }
        }
        BOOST_TEST(std::abs(impact.z - bestZ) < 1e-3);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(HalfATurnAwayIsTheTurnBehind)
{
  // The unit circle about the origin, clockwise (q Bz > 0) from (0, 1, 0) with
  // pz = pT; the axis at (0, -1) is half a turn away either way: the
  // convention's turn is +pi, behind, so z = -pi.
  const double bz = 5.;
  const double pt = std::abs(static_cast<double>(o2::constants::math::B2C) * bz); // radius exactly 1
  const TruthHelix helix{.position = {0., 1., 0.}, .momentum = {pt, 0., pt}, .charge = 1.};
  const auto impact = helixImpact(helix, bz, 0., -1.);
  BOOST_TEST(std::abs(impact.d0) < 1e-12);
  BOOST_TEST(std::abs(impact.z + std::numbers::pi) < 1e-12);
}

BOOST_AUTO_TEST_CASE(KinematicsNoHelixFollowsAreNaN)
{
  // A NaN momentum, a NaN field, no momentum at all.
  const TruthHelix invalid{.position = {1., 2., 3.}, .momentum = {std::nan(""), 0., 1.}, .charge = 1.};
  const TruthHelix fine{.position = {1., 2., 3.}, .momentum = {0.3, 0., 1.}, .charge = 1.};
  const TruthHelix atRest{.position = {1., 2., 3.}, .momentum = {0., 0., 0.}, .charge = 1.};
  for (const auto& [helix, bz] : {std::pair{invalid, -5.}, std::pair{fine, std::nan("")}, std::pair{atRest, -5.}}) {
    const auto impact = helixImpact(helix, bz, 0., 0.);
    BOOST_TEST(std::isnan(impact.d0));
    BOOST_TEST(std::isnan(impact.z));
    BOOST_TEST(std::isnan(helixMiss(helix, bz, 0., 0.)));
    BOOST_TEST(std::isnan(helixDirectionAt(helix, bz, 0., 0.)[0]));
  }
  // No transverse momentum: a line along z, with no closest point to an axis.
  const TruthHelix alongZ{.position = {1., 2., 3.}, .momentum = {0., 0., 1.}, .charge = 1.};
  BOOST_TEST(std::isnan(helixImpact(alongZ, -5., 0., 0.).d0));
  BOOST_TEST(std::abs(helixMiss(alongZ, -5., 0., 0.) - std::sqrt(5.)) < 1e-12);
  BOOST_TEST(std::abs(helixDirectionAt(alongZ, -5., 0., 0.)[2] - 1.) < 1e-12);
  // The point at the circle's centre has no tangent.
  const auto r = curvatureRadius(0.3, -5.);
  const auto atCentre = helixDirectionAt(TruthHelix{.position = {0., 0., 0.}, .momentum = {0.3, 0., 0.}, .charge = 1.}, -5., 0., r);
  BOOST_TEST(std::isnan(atCentre[0]));
}

BOOST_AUTO_TEST_CASE(ANeutralsLineMissesByItsDistance)
{
  const TruthHelix neutral{.position = {1., 0., 2.}, .momentum = {0., 0.5, 0.5}, .charge = 0.};
  const auto impact = helixImpact(neutral, -5., 0., 0.);
  BOOST_TEST(std::abs(std::abs(impact.d0) - 1.) < 1e-12);
  BOOST_TEST(std::abs(impact.z - 2.) < 1e-12);
  BOOST_TEST(std::abs(helixMiss(neutral, -5., 1., 3.)) < 1e-12);
  const auto direction = helixDirectionAt(neutral, -5., 1., 3.);
  BOOST_TEST(std::abs(direction[1] - std::numbers::sqrt2 / 2.) < 1e-12);
}

BOOST_AUTO_TEST_CASE(AnIterationOutsideTheConfiguredIsRefused)
{
  Toy toy;
  toy.kinematics[0].iteration = 2;
  BOOST_CHECK_THROW((void)trackingQA(toy.rows, toy.truth, toy.tracks, toy.kinematics, toy.seeded, {}, nullptr, 2, toy.settings), std::invalid_argument);
  BOOST_CHECK_THROW((void)trackingQA(toy.rows, toy.truth, toy.tracks, toy.kinematics, toy.seeded, {}, nullptr, kQAMaxIterations + 1, toy.settings),
                    std::invalid_argument);
}

} // namespace o2::gtrk::study
