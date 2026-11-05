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

/// \file testTrkChainPulls.cxx
/// \brief Fit pulls of chains through the ALICE 3 barrel, curved VD included,
///        and the TRK's readout frames, and on to FT3's discs.
///
/// The chains are made as testDiscChainPulls.cxx makes its own: a straight
/// line from the origin picks the first VD sensor, then the helix crosses the
/// barrel layer by layer at 2 T -- or, forward, the VD, the ML layers it still
/// reaches and the +z discs: each step goes to the layer it crosses first --
/// scattering at every sensor but the last as a
/// thin scatterer carrying the material of the step outward of it, as the fit
/// models it. A VD crossing is bound to the sensor's tangent plane at the
/// crossing (BoundTrackState::surfacePhi), found by moving the plane to where
/// the track meets it until the two agree; a measurement there keeps its
/// unrolled coordinates, so the fit finds the same plane from it
/// (TrackingGeometry::planeOf). Each measurement is smeared by exactly the
/// covariance it is given, a pixel's pitch^2 / 12 -- not the hit's exact
/// crossing, whose pulls would be narrow however right the geometry
/// -- and the fitted state on the first measurement's plane is
/// held to truth[0].
///
/// Runs with `GTRK_ALICE3_GEOMETRY` naming an o2sim geometry with the TRK;
/// without it, it reports that it did not. The pulls are
/// asserted with scattering and without material, barrel and forward, to
/// pT 20 GeV (200 GeV of momentum at |eta| 3): the fit carries its covariance
/// as factors (CovarianceFactors), and holds what a float covariance could not.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TRK
#include <boost/test/unit_test.hpp>
#include "PullStatistics.h"

#include "BetweenBoundCheck.h"
#include "GridTestSupport.h"

#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Alice/TRK.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Tracking/LayerMaterial.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/TrackFitter.h"

#include "CommonUtils/ConfigurableParam.h"
#include "DetectorsBase/GeometryManager.h"
#include "UniformMaterial.h"

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <numbers>
#include <optional>
#include <random>
#include <span>
#include <string>
#include <unistd.h>
#include <vector>

namespace o2::gtrk
{

namespace
{

using detectors::TRKFT3;

const char* geometryFile() { return std::getenv("GTRK_ALICE3_GEOMETRY"); }

constexpr Scalar kBzKiloGauss = 20.F;
constexpr double kX2X0PerCm = 1e-3;
constexpr std::size_t kMinMeasurements = 8; ///< the VD's three and the ML's five
/// Forward: the VD's three and three discs. At |eta| 3 a track leaves the VD
/// at its end, passes the ML and falls through the inner holes of two discs.
constexpr std::size_t kMinForwardMeasurements = 6;

/// Where a state meets a curved sensor: the tangent plane at the angle the
/// track crosses it, found by moving the plane to where the propagated track
/// meets it. Each round halves the angle's error many times over (the plane
/// and the cylinder part by the sagitta), so a few rounds reach float.
std::optional<std::pair<BoundTrackState, Scalar>> toCylinder(const Propagator& propagator, const BoundTrackState& state, const TrackingSurface& from,
                                                             const TrackingSurface& cylinder, Scalar angle)
{
  BoundTrackState probe = state;
  for (int round = 0; round < 6; ++round) {
    probe = state;
    const auto plane = tangentPlane(cylinder, angle);
    const auto status = propagate(propagator, probe, from, plane, PropagationDirection::kForward).status;
    if (status != PropagationStatus::kOk && status != PropagationStatus::kOutOfBounds) {
      return std::nullopt;
    }
    Scalar global[3];
    surfaceToGlobal(plane, probe.loc0(), probe.loc1(), global);
    Scalar d[3];
    for (int i = 0; i < 3; ++i) {
      d[i] = global[i] - cylinder.origin[i];
    }
    const Scalar x = (d[0] * cylinder.axis0[0]) + (d[1] * cylinder.axis0[1]) + (d[2] * cylinder.axis0[2]);
    const Scalar y = (d[0] * cylinder.axis1[0]) + (d[1] * cylinder.axis1[1]) + (d[2] * cylinder.axis1[2]);
    angle = std::atan2(y, x);
  }
  return std::pair{probe, angle};
}

/// A sensor's centre in the global frame, and where its layer is crossed.
struct Candidate {
  SensorID id;
  double z;
  double phi;
};

struct BarrelChain {
  const Geometry<TRKFT3>* geometry;
  const TrackingGeometry<TRKFT3>* tracking;
  Propagator propagator;
  std::array<std::vector<Candidate>, Detector<TRKFT3>::nLayers> sensorsOfLayer;
  std::vector<BoundTrackState> truth;
  std::vector<Measurement> measurements;

  BarrelChain(const Geometry<TRKFT3>& geo, const TrackingGeometry<TRKFT3>& trk) : geometry(&geo), tracking(&trk), truth(1)
  {
    propagator.bzKiloGauss = kBzKiloGauss;
    for (SensorID id = 0; id < static_cast<SensorID>(TRKFT3::nSensors); ++id) {
      const auto& sensor = geo.sensor(id);
      sensorsOfLayer[sensor.layer].push_back({.id = id, .z = sensor.trans[2], .phi = std::atan2(sensor.trans[1], sensor.trans[0])});
    }
  }

  [[nodiscard]] static bool curved(LayerID layer) { return Geometry<TRKFT3>::layerSpec(layer).curved(); }

  /// The curved sensor of \a layer whose arc holds the global azimuth
  /// \a phi, with that azimuth as its local angle; nothing if none does.
  [[nodiscard]] std::optional<std::pair<SensorID, Scalar>> petalAt(LayerID layer, Scalar phi) const
  {
    for (const auto& candidate : sensorsOfLayer[layer]) {
      const auto& c = tracking->surface(candidate.id);
      // The direction at phi in the sensor's frame: against its local x and y.
      const Scalar x = (std::cos(phi) * c.axis0[0]) + (std::sin(phi) * c.axis0[1]);
      const Scalar y = (std::cos(phi) * c.axis1[0]) + (std::sin(phi) * c.axis1[1]);
      const Scalar angle = std::atan2(y, x);
      if (std::abs(angle) * c.radius <= c.bounds.halfX) {
        return std::pair{candidate.id, angle};
      }
    }
    return std::nullopt;
  }

  /// The state on the sensor of \a layer the track crosses next; nothing if
  /// it crosses none.
  [[nodiscard]] std::optional<BoundTrackState> cross(const BoundTrackState& state, LayerID layer) const
  {
    const auto from = tracking->planeOf(state);
    Scalar global[3];
    surfaceToGlobal(from, state.loc0(), state.loc1(), global);
    if (!Detector<TRKFT3>::isBarrelLayer(layer)) {
      return crossDisc(state, from, layer);
    }
    if (curved(layer)) {
      // Start the plane where the track points, then move it onto the crossing.
      const auto petal = petalAt(layer, state.phi());
      if (!petal) {
        return std::nullopt;
      }
      const auto& cylinder = tracking->surface(petal->first);
      const auto found = toCylinder(propagator, state, from, cylinder, petal->second);
      if (!found || std::abs(found->second) * cylinder.radius > cylinder.bounds.halfX) {
        return std::nullopt;
      }
      BoundTrackState result = state;
      if (propagate(propagator, result, from, tracking->plane(petal->first, found->second), PropagationDirection::kForward).status != PropagationStatus::kOk) {
        return std::nullopt;
      }
      return result;
    }
    // A flat layer: the chips near the line's crossing of the layer's radius,
    // the first one the helix reaches inside its bounds.
    const auto& envelope = PlacedEnvelopes<TRKFT3>::all()[layer];
    TrackingSurface ring{};
    ring.axis0[0] = 1.F;
    ring.axis1[1] = 1.F;
    ring.normal[2] = 1.F;
    ring.radius = static_cast<Scalar>(envelope.position);
    ring.bounds = {.halfX = 1e9F, .halfY = 1e9F};
    const auto atRing = toCylinder(propagator, state, from, ring, state.phi());
    if (!atRing) {
      return std::nullopt;
    }
    Scalar point[3];
    surfaceToGlobal(tangentPlane(ring, atRing->second), atRing->first.loc0(), atRing->first.loc1(), point);
    const auto pointPhi = static_cast<double>(std::atan2(point[1], point[0]));
    std::optional<BoundTrackState> best;
    Scalar bestPath = std::numeric_limits<Scalar>::infinity();
    for (const auto& candidate : sensorsOfLayer[layer]) {
      double dPhi = candidate.phi - pointPhi;
      math_utils::detail::bringToPMPi<double>(dPhi);
      if (std::abs(candidate.z - static_cast<double>(point[2])) > 2.5 || std::abs(dPhi) * envelope.position > 2.5) {
        continue;
      }
      BoundTrackState probe = state;
      const auto result = propagate(propagator, probe, from, tracking->surface(candidate.id), PropagationDirection::kForward);
      if (result.status == PropagationStatus::kOk && result.pathLength > 0.F && result.pathLength < bestPath) {
        bestPath = result.pathLength;
        best = probe;
      }
    }
    return best;
  }

  /// A disc: the chips near where the track meets the disc's plane, the first
  /// one the helix reaches inside its bounds.
  [[nodiscard]] std::optional<BoundTrackState> crossDisc(const BoundTrackState& state, const TrackingSurface& from, LayerID layer) const
  {
    TrackingSurface plane{};
    plane.axis0[0] = 1.F;
    plane.axis1[1] = 1.F;
    plane.normal[2] = 1.F;
    plane.origin[2] = static_cast<Scalar>(PlacedEnvelopes<TRKFT3>::all()[layer].position);
    plane.bounds = {.halfX = 1e9F, .halfY = 1e9F};
    BoundTrackState atPlane = state;
    if (propagate(propagator, atPlane, from, plane, PropagationDirection::kForward).status != PropagationStatus::kOk) {
      return std::nullopt;
    }
    Scalar point[3];
    surfaceToGlobal(plane, atPlane.loc0(), atPlane.loc1(), point);
    std::optional<BoundTrackState> best;
    Scalar bestPath = std::numeric_limits<Scalar>::infinity();
    for (const auto& candidate : sensorsOfLayer[layer]) {
      const auto& sensor = geometry->sensor(candidate.id);
      if (std::hypot(sensor.trans[0] - static_cast<double>(point[0]), sensor.trans[1] - static_cast<double>(point[1])) > 2.5) {
        continue;
      }
      BoundTrackState probe = state;
      const auto result = propagate(propagator, probe, from, tracking->surface(candidate.id), PropagationDirection::kForward);
      if (result.status == PropagationStatus::kOk && result.pathLength > 0.F && result.pathLength < bestPath) {
        bestPath = result.pathLength;
        best = probe;
      }
    }
    return best;
  }

  /// truth[0] on the VD sensor of layer 0 a straight line from (0, 0, z0)
  /// at (phi, theta) crosses; false if it crosses none.
  bool start(Scalar phi, Scalar theta, Scalar qOverP, Scalar z0)
  {
    const auto petal = petalAt(0, phi);
    if (!petal) {
      return false;
    }
    const auto plane = tracking->plane(petal->first, petal->second);
    const auto& cylinder = tracking->surface(petal->first);
    const Scalar t = cylinder.radius / std::sin(theta);
    const Scalar g[3]{t * std::sin(theta) * std::cos(phi), t * std::sin(theta) * std::sin(phi), z0 + (t * std::cos(theta))};
    Scalar l[3];
    surfaceToLocal(plane, g, l);
    if (!plane.bounds.isInside(l[0], l[1])) {
      return false;
    }
    truth[0] = BoundTrackState{};
    truth[0].par[kLoc0] = l[0];
    truth[0].par[kLoc1] = l[1];
    truth[0].par[kPhi] = phi;
    truth[0].par[kTheta] = theta;
    truth[0].par[kQOverP] = qOverP;
    truth[0].surface = petal->first;
    truth[0].surfacePhi = petal->second;
    return true;
  }

  /// Where \a state is, globally.
  [[nodiscard]] std::array<Scalar, 3> globalOf(const BoundTrackState& state) const
  {
    std::array<Scalar, 3> g{};
    surfaceToGlobal(tracking->planeOf(state), state.loc0(), state.loc1(), g.data());
    return g;
  }

  /// The layers a chain may cross after layer 0, in order: the barrel's, then
  /// with \a forward the +z discs.
  [[nodiscard]] static std::vector<LayerID> sequence(bool forward)
  {
    std::vector<LayerID> layers;
    for (LayerID layer = 1; std::cmp_less(layer, Detector<TRKFT3>::nBarrelLayers); ++layer) {
      layers.push_back(layer);
    }
    if (forward) {
      // FT3's order: six discs at -z, then their mirrors at +z.
      for (int disc = 6; disc < Detector<TRKFT3>::nDiscLayers; ++disc) {
        layers.push_back(static_cast<LayerID>(Detector<TRKFT3>::nBarrelLayers + disc));
      }
    }
    return layers;
  }

  /// The rest of the chain from truth[0], scattering if \a scatter: each step
  /// to the layer of \a layers, after the last one crossed, that the track
  /// crosses first; false if fewer than \a minMeasurements remain.
  bool generate(bool scatter, std::mt19937& rng, std::span<const LayerID> layers, std::size_t minMeasurements)
  {
    std::normal_distribution<double> gauss(0., 1.);
    truth.resize(1);
    std::size_t next = 0;
    while (next < layers.size()) {
      BoundTrackState outgoing = truth.back();
      const auto here = globalOf(outgoing);
      std::optional<BoundTrackState> unkicked;
      std::size_t chosen = layers.size();
      Scalar nearest = std::numeric_limits<Scalar>::infinity();
      for (std::size_t i = next; i < layers.size(); ++i) {
        const auto crossing = cross(outgoing, layers[i]);
        if (!crossing) {
          continue;
        }
        const auto there = globalOf(*crossing);
        const Scalar distance = std::hypot(there[0] - here[0], there[1] - here[1], there[2] - here[2]);
        if (distance < nearest) {
          nearest = distance;
          unkicked = crossing;
          chosen = i;
        }
      }
      if (!unkicked) {
        break;
      }
      const LayerID layer = layers[chosen];
      next = chosen + 1;
      Scalar theta0Half = 0.F;
      const auto kick = [&](BoundTrackState& state) {
        state.par[kTheta] += static_cast<Scalar>(gauss(rng)) * theta0Half;
        state.par[kPhi] += static_cast<Scalar>(gauss(rng)) * theta0Half / std::sin(state.theta());
      };
      if (scatter) {
        Scalar a[3]{};
        Scalar b[3]{};
        surfaceToGlobal(tracking->planeOf(outgoing), outgoing.loc0(), outgoing.loc1(), a);
        surfaceToGlobal(tracking->planeOf(*unkicked), unkicked->loc0(), unkicked->loc1(), b);
        const Scalar chord = std::hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
        theta0Half = highlandTheta0(0.5F * static_cast<Scalar>(kX2X0PerCm) * chord, betaMomentum(1.F / std::abs(outgoing.qOverP()), propagator.massGeV));
        kick(outgoing);
      }
      auto crossed = cross(outgoing, layer);
      if (!crossed) {
        break;
      }
      if (scatter) {
        kick(*crossed);
      }
      truth.push_back(*crossed);
    }
    if (truth.size() < minMeasurements) {
      return false;
    }
    measurements.clear();
    for (std::size_t i = 0; i < truth.size(); ++i) {
      const auto& spec = Geometry<TRKFT3>::layerSpec(tracking->surface(truth[i].surface).layer);
      const auto var = static_cast<Scalar>(spec.pitchX * spec.pitchX / 12.);
      measurements.push_back(Measurement{.loc0 = truth[i].loc0(),
                                         .loc1 = truth[i].loc1(),
                                         .timePS = test::toyTimeBox(test::kToyTrackTimeNS).start,
                                         .widthPS = test::toyTimeBox(test::kToyTrackTimeNS).width,
                                         .varLoc0 = var,
                                         .covLoc01 = 0.F,
                                         .varLoc1 = var,
                                         .clusterIndex = static_cast<std::uint32_t>(i),
                                         .surface = truth[i].surface});
    }
    return true;
  }
};

using namespace test::pulls;

PullPoint pullsAt(BarrelChain& chain, double ptGeV, double eta, bool withMaterial, int trials, std::uint32_t seed)
{
  const bool forward = eta > 1.5;
  const auto layers = BarrelChain::sequence(forward);
  const UniformMaterial material{.x2X0PerCm = withMaterial ? static_cast<Scalar>(kX2X0PerCm) : 0.F};
  const double theta = 2. * std::atan(std::exp(-eta));
  const double p = ptGeV / std::sin(theta);
  std::mt19937 rng(seed);
  std::normal_distribution<double> gauss(0., 1.);
  std::uniform_real_distribution<double> uniformPhi(-std::numbers::pi, std::numbers::pi);
  // Vertices along the beam, so that the chains cross the chips everywhere
  // and not all on the seams between them at z = 0.
  std::uniform_real_distribution<double> uniformZ(-5., 5.);
  return collectPulls(chain, ptGeV, eta, trials, [&](int trial) {
    const Scalar charge = (trial % 2) != 0 ? 1.F : -1.F;
    bool built = false;
    // Gaps between petals and chips lose some lines; draw again. At |eta| 3
    // most lines fall through the discs' inner holes, hence the many draws.
    for (int attempt = 0; attempt < 1000 && !built; ++attempt) {
      built = chain.start(static_cast<Scalar>(uniformPhi(rng)), static_cast<Scalar>(theta), static_cast<Scalar>(static_cast<double>(charge) / p),
                          static_cast<Scalar>(uniformZ(rng))) &&
              chain.generate(withMaterial, rng, layers, forward ? kMinForwardMeasurements : kMinMeasurements);
    }
    BOOST_REQUIRE_MESSAGE(built, "no chain in 1000 draws at pT " << ptGeV << " eta " << eta << (withMaterial ? "" : ", no material"));
    auto smeared = chain.measurements;
    for (auto& m : smeared) {
      m.loc0 += static_cast<Scalar>(gauss(rng)) * std::sqrt(m.varLoc0);
      m.loc1 += static_cast<Scalar>(gauss(rng)) * std::sqrt(m.varLoc1);
    }
    return fitMeasurements(chain.propagator, *chain.tracking, material, smeared.data(), static_cast<int>(smeared.size()));
  });
}

/// From where a 2 T track still crosses the ML (0.2 GeV curls at 33 cm), to
/// where it is nearly straight; at |eta| the OT still covers.
constexpr std::array<double, 5> kPtScan{0.2, 0.5, 1., 5., 20.};
constexpr std::array<double, 3> kEtas{0., 0.6, 1.1};
/// Forward: through the VD and the inner ML onto FT3's +z discs (their |eta|
/// 1.5-3.2). At high momentum there a float covariance cannot hold the track
/// (VD to VD at |eta| 3 moves loc1 by 71 cm per radian of theta); the fit
/// carries it as factors.
constexpr std::array<double, 5> kForwardPtScan{0.2, 0.5, 1., 5., 20.};
constexpr std::array<double, 3> kForwardEtas{2., 2.5, 3.};
/// The widths of n unit-Gaussian pulls scatter by 1 / sqrt(2 n) (plain) and
/// about 1.6 times that (robust, the median's): at 2000, 1.6 % and 2.6 %, so
/// 15 % is beyond 5 sigma of either over the scan's 600 widths.
constexpr int kTrials = 2000;
constexpr double kTolerance = 0.15;

struct Scans {
  std::vector<PullPoint> clean;
  std::vector<PullPoint> scattered;
  std::vector<PullPoint> forwardClean;
  std::vector<PullPoint> forwardScattered;
};

/// Null without `GTRK_ALICE3_GEOMETRY`.
const Scans* scans()
{
  static const std::unique_ptr<const Scans> result = []() -> std::unique_ptr<const Scans> {
    if (geometryFile() == nullptr) {
      return nullptr;
    }
    o2::base::GeometryManager::loadGeometry(geometryFile());
    const auto geometry = std::make_unique<const Geometry<TRKFT3>>(alice::placeChecked<TRKFT3>());
    const auto tracking = std::make_unique<const TrackingGeometry<TRKFT3>>(*geometry);
    BarrelChain chain(*geometry, *tracking);
    auto out = std::make_unique<Scans>();
    std::uint32_t seed = 20261007;
    std::uint32_t forwardSeed = 20261107; // apart from the barrel's 30 seeds
    for (const bool withMaterial : {false, true}) {
      auto& points = withMaterial ? out->scattered : out->clean;
      for (const double eta : kEtas) {
        for (const double pt : kPtScan) {
          points.push_back(pullsAt(chain, pt, eta, withMaterial, kTrials, seed++));
        }
      }
      auto& forward = withMaterial ? out->forwardScattered : out->forwardClean;
      for (const double eta : kForwardEtas) {
        for (const double pt : kForwardPtScan) {
          forward.push_back(pullsAt(chain, pt, eta, withMaterial, kTrials, forwardSeed++));
        }
      }
    }
    BOOST_TEST_MESSAGE(report("ALICE 3 barrel, no material", out->clean, false));
    BOOST_TEST_MESSAGE(report("ALICE 3 barrel, 1e-3 X0/cm", out->scattered, false));
    BOOST_TEST_MESSAGE(report("ALICE 3 forward, onto FT3, no material", out->forwardClean, false));
    BOOST_TEST_MESSAGE(report("ALICE 3 forward, onto FT3, 1e-3 X0/cm", out->forwardScattered, false));
    return out;
  }();
  return result.get();
}

void checkUnit(std::span<const PullPoint> points) { test::pulls::checkUnit(points, kTolerance); }

} // namespace

/// A tangent plane is right-handed, radial and touches its sensor where its
/// angle says, at loc0 = radius * angle.
/// Every one of the 12 VD sensors: the plane a measurement's unrolled
/// coordinates name contains that point of the sensor.
BOOST_AUTO_TEST_CASE(EveryPetalPlacesItsMeasurements)
{
  if (geometryFile() == nullptr) {
    BOOST_TEST_MESSAGE("GTRK_ALICE3_GEOMETRY not set: the TRK is not loaded");
    return;
  }
  o2::base::GeometryManager::loadGeometry(geometryFile());
  const auto geometry = std::make_unique<const Geometry<TRKFT3>>(alice::placeChecked<TRKFT3>());
  const auto tracking = std::make_unique<const TrackingGeometry<TRKFT3>>(*geometry);
  int petals = 0;
  for (SensorID id = 0; id < static_cast<SensorID>(TRKFT3::nSensors); ++id) {
    const auto& spec = geometry->specOf(id);
    if (!spec.curved()) {
      continue;
    }
    ++petals;
    for (const double u : {-0.99 * spec.halfX, -0.3 * spec.halfX, 0., 0.5 * spec.halfX, 0.99 * spec.halfX}) {
      const double v = 17.;
      double frame[3];
      spec.unrolledToFrame(u, v, 0., frame);
      double expected[3];
      geometry->sensor(id).toGlobal(frame, expected);
      const Measurement m{.loc0 = static_cast<Scalar>(u), .loc1 = static_cast<Scalar>(v), .surface = id};
      const auto plane = tracking->planeOf(m);
      Scalar global[3];
      surfaceToGlobal(plane, m.loc0, m.loc1, global);
      for (int i = 0; i < 3; ++i) {
        BOOST_TEST(std::abs(static_cast<double>(global[i]) - expected[i]) < 1e-5); // float, at 25 cm
      }
      // The plane's normal is radial there.
      const double r = std::hypot(expected[0], expected[1]);
      BOOST_TEST(std::abs((static_cast<double>(plane.normal[0]) * expected[0] / r) + (static_cast<double>(plane.normal[1]) * expected[1] / r) - 1.) < 1e-6);
    }
  }
  BOOST_TEST(petals == 12);
}

/// With scattering at the fit's own model, chains through the VD's tangent
/// planes and the ML/OT chips fit with unit pulls.
BOOST_AUTO_TEST_CASE(BarrelPullsWithScatteringAreUnit)
{
  const auto* result = scans();
  if (result == nullptr) {
    BOOST_TEST_MESSAGE("GTRK_ALICE3_GEOMETRY not set: the TRK's chains are not fitted");
    return;
  }
  checkUnit(result->scattered);
}

/// The same forward, through the VD and the inner ML onto FT3's discs.
BOOST_AUTO_TEST_CASE(DiscPullsWithScatteringAreUnit)
{
  const auto* result = scans();
  if (result == nullptr) {
    BOOST_TEST_MESSAGE("GTRK_ALICE3_GEOMETRY not set: the TRK's chains are not fitted");
    return;
  }
  checkUnit(result->forwardScattered);
}

/// No pair's kink between its layers exceeds its link's bound on any link with
/// one of the measured passives between its layers, at pT 0.1 GeV over a
/// 5 cm beam spot: the tracklet stage's test at the bound refuses only pairs
/// the pair's own kink would refuse too.
BOOST_AUTO_TEST_CASE(NoPairsKinkExceedsItsLinksBound)
{
  if (geometryFile() == nullptr) {
    BOOST_TEST_MESSAGE("GTRK_ALICE3_GEOMETRY not set: the TRK is not loaded");
    return;
  }
  o2::base::GeometryManager::loadGeometry(geometryFile());
  const auto geometry = std::make_unique<const Geometry<TRKFT3>>(alice::placeChecked<TRKFT3>());
  SearchWindowParams windows;
  windows.minPtGeV = 0.1;
  windows.bzKiloGauss = 20.;
  const auto scattering = scatteringTable<TRKFT3>(uniformLayerMaterial<TRKFT3>(1.e-3));
  const auto ratios = test::betweenOverBound<TRKFT3>(scattering, windows, VertexConstraint{.z = 0., .sigmaZ = 5., .sigmaTransverse = 0.005});
  BOOST_TEST(!ratios.empty());
  BOOST_TEST(test::segmentMaskMismatches<TRKFT3>(scattering) == 0);
  for (const auto& ratio : ratios) {
    BOOST_TEST_CONTEXT("link " << static_cast<int>(ratio.from) << " -> " << static_cast<int>(ratio.to))
    {
      BOOST_TEST(ratio.worst <= 1.);
      BOOST_TEST(ratio.maskMismatches == 0);
      BOOST_TEST(ratio.sourcePartMismatches == 0);
    }
  }
}

/// The material LUT the tracking reads against the geometry it was made from,
/// band by band as the layout splits the material: each passive's band
/// (`passiveMaterial`) and each layer's slab (`layerMaterial`, with the
/// passives its silicon envelope), sampled at the same azimuths and z on both.
/// The passives' measured x/X0 (TRKFT3.h) is reported beside them. Runs with
/// `GTRK_ALICE3_MATLUT` naming the LUT and `GTRK_ALICE3_GEOMETRY` the geometry.
BOOST_AUTO_TEST_CASE(LutHoldsTheGeometrysMaterial)
{
  const char* lutFile = std::getenv("GTRK_ALICE3_MATLUT");
  if (geometryFile() == nullptr || lutFile == nullptr) {
    BOOST_TEST_MESSAGE("GTRK_ALICE3_GEOMETRY or GTRK_ALICE3_MATLUT not set: the LUT is not held to the geometry");
    return;
  }
  o2::base::GeometryManager::loadGeometry(geometryFile());
  const auto geometry = std::make_unique<const Geometry<TRKFT3>>(alice::placeChecked<TRKFT3>());
  const std::unique_ptr<o2::base::MatLayerCylSet> lut{o2::base::MatLayerCylSet::loadFromFile(lutFile)};
  BOOST_REQUIRE(lut);
  // The same averages from the geometry, at random azimuths and z, with the
  // mean's standard error. Some gaps hold thin radial walls -- the petals'
  // side walls, the OT's quarter-barrel panels -- inside which a radial
  // segment runs its whole length, and structures that end along z: a regular
  // grid hits them by aliasing, and the mean is known to its statistical error.
  struct Estimate {
    double mean;
    double error;
  };
  std::mt19937 rng(20261007);
  std::uniform_real_distribution<double> uniformPhi(0., 2. * std::numbers::pi);
  const auto fromGeometry = [&](double low, double high, double zMin, double zMax) {
    constexpr int kSamples = 65536;
    std::uniform_real_distribution<double> uniformZ(zMin, zMax);
    double sum = 0.;
    double sum2 = 0.;
    {
      for (int sample = 0; sample < kSamples; ++sample) {
        const double phi = uniformPhi(rng);
        const double z = uniformZ(rng);
        const auto x =
          static_cast<double>(o2::base::GeometryManager::meanMaterialBudget(o2::math_utils::Point3D<double>{low * std::cos(phi), low * std::sin(phi), z},
                                                                            o2::math_utils::Point3D<double>{high * std::cos(phi), high * std::sin(phi), z})
                                .meanX2X0);
        sum += x;
        sum2 += x * x;
      }
    }
    const double n = kSamples;
    const double mean = sum / n;
    return Estimate{.mean = mean, .error = std::sqrt(std::max(0., (sum2 / n) - (mean * mean)) / n)};
  };
  // The LUT averages each cell over its own 10 x 10 samples; agreement within
  // 5 % or 3 standard errors of the geometry's mean, whichever is wider.
  const auto close = [](double fromLut, const Estimate& geo) { return std::abs(fromLut - geo.mean) <= std::max(0.05 * geo.mean, 3. * geo.error) + 1e-7; };
  const auto passives = passiveMaterial<TRKFT3>(*lut);
  std::string out = "passive: radius; x/X0 measured (TRKFT3.h), geometry +- its error, LUT\n";
  for (std::size_t i = 0; i < passives.size(); ++i) {
    const auto& spec = TRKFT3::passives[i];
    const auto geo = fromGeometry(spec.bandLow, spec.bandHigh, spec.low, spec.high);
    out += fmt::format("  {:2d}: {:8.4f} cm  {:.3e}  {:.3e} +- {:.1e}  {:.3e}\n", i, spec.position, spec.x2X0, geo.mean, geo.error, passives[i]);
    BOOST_TEST_CONTEXT("passive " << i << " geometry " << geo.mean << " +- " << geo.error << " LUT " << passives[i]) { BOOST_TEST(close(passives[i], geo)); }
  }
  // A disc's slab from the geometry as layerMaterial reads it: segments from
  // z = low to high along the ray from the nominal vertex through each radius
  // of the disc, scaled to normal incidence.
  const auto fromDisc = [&](double low, double high, double position, double rMin, double rMax) {
    constexpr int kSamples = 65536;
    std::uniform_real_distribution<double> uniformR(rMin, rMax);
    double sum = 0.;
    double sum2 = 0.;
    for (int sample = 0; sample < kSamples; ++sample) {
      const double phi = uniformPhi(rng);
      const double r = uniformR(rng);
      const double rLow = r * low / position;
      const double rHigh = r * high / position;
      const double length = std::hypot(rHigh - rLow, high - low);
      const auto x =
        static_cast<double>(o2::base::GeometryManager::meanMaterialBudget(o2::math_utils::Point3D<double>{rLow * std::cos(phi), rLow * std::sin(phi), low},
                                                                          o2::math_utils::Point3D<double>{rHigh * std::cos(phi), rHigh * std::sin(phi), high})
                              .meanX2X0) *
        std::abs(high - low) / length;
      sum += x;
      sum2 += x * x;
    }
    const double n = kSamples;
    const double mean = sum / n;
    return Estimate{.mean = mean, .error = std::sqrt(std::max(0., (sum2 / n) - (mean * mean)) / n)};
  };
  const auto layers = layerMaterial<TRKFT3>(*lut);
  out +=
    "layer: x/X0 of its slab (a barrel layer's silicon envelope, the outermost's mirrored about it; a disc's half-gaps to the discs\n"
    "beside it on its side, the innermost's mirrored), geometry +- its error, LUT\n";
  constexpr auto kBarrel = static_cast<std::size_t>(Detector<TRKFT3>::nBarrelLayers);
  for (std::size_t i = 0; i < layers.size(); ++i) {
    const auto& e = PlacedEnvelopes<TRKFT3>::all()[i];
    Estimate geo{};
    if (i < kBarrel) {
      const double low = TRKFT3::passives[i].bandHigh;
      const double high = i + 1 < kBarrel ? TRKFT3::passives[i + 1].bandLow : e.position + (e.position - low);
      geo = fromGeometry(low, high, e.zMin, e.zMax);
    } else {
      // The neighbours of the same side, as layerMaterial takes them.
      const auto sameSide = [&](std::size_t other) {
        return other >= kBarrel && other < layers.size() && PlacedEnvelopes<TRKFT3>::all()[other].position * e.position > 0.;
      };
      const bool hasBelow = sameSide(i - 1);
      const bool hasAbove = sameSide(i + 1);
      const double below = hasBelow ? PlacedEnvelopes<TRKFT3>::all()[i - 1].position : 0.;
      const double above = hasAbove ? PlacedEnvelopes<TRKFT3>::all()[i + 1].position : 0.;
      double low = hasBelow ? 0.5 * (below + e.position) : 0.;
      double high = hasAbove ? 0.5 * (above + e.position) : 0.;
      low = hasBelow ? low : e.position - (high - e.position);
      high = hasAbove ? high : e.position + (e.position - low);
      // The outermost discs' mirrored slabs end at the LUT's |z|, as the
      // LUT's reading does; beyond lies material behind the last disc.
      const auto zMax = static_cast<double>(lut->getZMax());
      geo = fromDisc(std::clamp(low, -zMax, zMax), std::clamp(high, -zMax, zMax), e.position, e.rMin, e.rMax);
    }
    out += fmt::format("  {:2d}: {:.3e} +- {:.1e}  {:.3e}\n", i, geo.mean, geo.error, layers[i]);
    BOOST_TEST_CONTEXT("layer " << i << " geometry " << geo.mean << " +- " << geo.error << " LUT " << layers[i]) { BOOST_TEST(close(layers[i], geo)); }
  }
  BOOST_TEST_MESSAGE(out);
}

/// Without material, barrel and forward: with nothing to scatter on, 2.9
/// and 5.8 um over 80 cm at 2 T fix the track far below what a float
/// covariance holds; the fit carries it as factors.
BOOST_AUTO_TEST_CASE(PullsWithoutMaterialAreUnit)
{
  const auto* result = scans();
  if (result == nullptr) {
    BOOST_TEST_MESSAGE("GTRK_ALICE3_GEOMETRY not set: the TRK's chains are not fitted");
    return;
  }
  checkUnit(result->clean);
  checkUnit(result->forwardClean);
}

/// The TRK's and FT3's readout frames come from TRKAlmiraParam, here loaded
/// from a digitiser's configuration, one per layer, the discs' after the
/// barrel's: an unset length is O2's 18 BC, and the first start is the bias
/// and the delay together.
BOOST_AUTO_TEST_CASE(TheTrksFramesAreTheDigitisersPerLayer)
{
  const auto path = (std::filesystem::temp_directory_path() / ("gtrk-test-almira-" + std::to_string(getpid()) + ".ini")).string();
  {
    std::ofstream out{path};
    out << "[TRKAlmiraParam]\n";
    // The VD's 22 BC, ML and OT's 81 BC, the last layer unset.
    constexpr std::array<int, 11> kLengths{22, 22, 22, 81, 81, 81, 81, 81, 81, 81, 0};
    for (std::size_t layer = 0; layer < kLengths.size(); ++layer) {
      out << "roFrameLengthInBCPerLayer[" << layer << "]=" << kLengths[layer] << '\n';
    }
    out << "roFrameBiasInBCPerLayer[4]=5\nroFrameDelayInBCPerLayer[4]=2\n";
    // FT3's first disc at ML and OT's 81 BC, the others unset.
    out << "roFrameLengthInBCPerLayer[11]=81\nroFrameBiasInBCPerLayer[11]=3\n";
  }
  o2::conf::ConfigurableParam::updateFromFile(path, "TRKAlmiraParam");
  std::filesystem::remove(path);
  const auto frames = alice::trkFrameTiming();
  BOOST_REQUIRE_EQUAL(frames.size(), 23U);
  BOOST_TEST((frames[0] == FrameTiming{.lengthBC = 22, .firstStartBC = 0}));
  BOOST_TEST((frames[3] == FrameTiming{.lengthBC = 81, .firstStartBC = 0}));
  BOOST_TEST((frames[4] == FrameTiming{.lengthBC = 81, .firstStartBC = 7}));
  BOOST_TEST((frames[10] == FrameTiming{.lengthBC = 18, .firstStartBC = 0}));
  BOOST_TEST((frames[11] == FrameTiming{.lengthBC = 81, .firstStartBC = 3}));
  BOOST_TEST((frames[22] == FrameTiming{.lengthBC = 18, .firstStartBC = 0}));
}

} // namespace o2::gtrk
