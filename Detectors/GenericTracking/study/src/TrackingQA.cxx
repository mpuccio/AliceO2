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

#include "GenericTracking/Study/TrackingQA.h"

#include "GenericTracking/Tracking/MathUtils.h"

#include "CommonUtils/TreeStreamRedirector.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <format>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace o2::gtrk::study
{

std::vector<HitCrossing> outwardLeg(std::vector<HitCrossing> crossings, const std::array<double, 3>& production)
{
  std::ranges::stable_sort(crossings, {}, &HitCrossing::time);
  double last = 0.;
  for (std::size_t i = 0; i < crossings.size(); ++i) {
    const auto& position = crossings[i].position;
    const double distance = std::hypot(position[0] - production[0], position[1] - production[1], position[2] - production[2]);
    if (distance < last) {
      crossings.resize(i);
      break;
    }
    last = distance;
  }
  return crossings;
}

std::uint8_t outwardLayers(std::span<const HitCrossing> leg)
{
  std::uint64_t layers = 0U;
  for (const auto& crossing : leg) {
    layers |= std::uint64_t{1} << crossing.layer;
  }
  return static_cast<std::uint8_t>(std::popcount(layers));
}

EtaRange etaReach(std::span<const LayerEnvelope> envelopes)
{
  double low = std::numeric_limits<double>::infinity();
  double high = -std::numeric_limits<double>::infinity();
  for (const auto& envelope : envelopes) {
    for (const double r : {envelope.rMin, envelope.rMax}) {
      for (const double z : {envelope.zMin, envelope.zMax}) {
        const double eta = std::asinh(z / r);
        low = std::min(low, eta);
        high = std::max(high, eta);
      }
    }
  }
  if (!(low < high)) {
    throw std::invalid_argument("QA: the envelopes span no eta");
  }
  return {.low = std::floor(low * 2.) / 2., .high = std::ceil(high * 2.) / 2.};
}

namespace
{

/// A helix's circle: its centre, its radius, and the sense of motion, +1
/// counter-clockwise; none for a neutral or no field.
struct Circle {
  double x = 0.;
  double y = 0.;
  double radius = 0.;
  double sense = 0.;
};

std::optional<Circle> circleOf(const TruthHelix& helix, double bzKiloGauss)
{
  const double pt = std::hypot(helix.momentum[0], helix.momentum[1]);
  const double charge = helix.charge;
  if (exactlyEqual(charge, 0.) || exactlyEqual(bzKiloGauss, 0.) || !(pt > 0.)) {
    return std::nullopt;
  }
  const double radius = curvatureRadius(pt, bzKiloGauss) / std::abs(charge);
  // F = q v x B: with q Bz > 0 the force points to the right of the
  // direction, (uy, -ux), and the motion is clockwise.
  const double side = std::copysign(1., charge * bzKiloGauss);
  const double ux = helix.momentum[0] / pt;
  const double uy = helix.momentum[1] / pt;
  return Circle{.x = helix.position[0] + (side * radius * uy), .y = helix.position[1] - (side * radius * ux), .radius = radius, .sense = -side};
}

/// Whether \a helix, \a bzKiloGauss and the point (\a x, \a y) are finite,
/// and the helix has some momentum.
bool followable(const TruthHelix& helix, double bzKiloGauss, double x, double y)
{
  const auto finite = [](const std::array<double, 3>& v) { return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]); };
  return finite(helix.position) && finite(helix.momentum) && std::isfinite(helix.charge) && std::isfinite(bzKiloGauss) && std::isfinite(x) && std::isfinite(y) && std::hypot(helix.momentum[0], helix.momentum[1], helix.momentum[2]) > 0.;
}

/// Whether \a helix moves transversely at all.
bool transverse(const TruthHelix& helix) { return std::hypot(helix.momentum[0], helix.momentum[1]) > 0.; }

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

} // namespace

HelixImpact helixImpact(const TruthHelix& helix, double bzKiloGauss, double x, double y)
{
  if (!followable(helix, bzKiloGauss, x, y) || !transverse(helix)) {
    return {.d0 = kNaN, .z = kNaN};
  }
  const double pt = std::hypot(helix.momentum[0], helix.momentum[1]);
  const double cotTheta = helix.momentum[2] / pt;
  const auto circle = circleOf(helix, bzKiloGauss);
  if (!circle) {
    // The straight line: the axis's side by the cross product of the
    // direction with it, the transverse path back to the closest point.
    const double ux = helix.momentum[0] / pt;
    const double uy = helix.momentum[1] / pt;
    const double dx = x - helix.position[0];
    const double dy = y - helix.position[1];
    const double along = (ux * dx) + (uy * dy);
    return {.d0 = (ux * dy) - (uy * dx), .z = helix.position[2] + (along * cotTheta)};
  }
  const double cx = x - circle->x;
  const double cy = y - circle->y;
  const double distance = std::hypot(cx, cy);
  // The closest point is on the centre's ray through the axis; the turn
  // from it to the production point in the sense of motion, in (-pi, pi].
  const double px = helix.position[0] - circle->x;
  const double py = helix.position[1] - circle->y;
  double turn = std::atan2((cx * py) - (cy * px), (cx * px) + (cy * py)) * circle->sense;
  if (turn <= -std::numbers::pi) {
    turn += 2. * std::numbers::pi; // the half turn's tie: the convention's +pi
  }
  return {.d0 = distance - circle->radius, .z = helix.position[2] - (circle->radius * turn * cotTheta)};
}

std::array<double, 3> helixDirectionAt(const TruthHelix& helix, double bzKiloGauss, double x, double y)
{
  if (!followable(helix, bzKiloGauss, x, y)) {
    return {kNaN, kNaN, kNaN};
  }
  const double p = std::hypot(helix.momentum[0], helix.momentum[1], helix.momentum[2]);
  const double pt = std::hypot(helix.momentum[0], helix.momentum[1]);
  const auto circle = circleOf(helix, bzKiloGauss);
  if (!circle) {
    return {helix.momentum[0] / p, helix.momentum[1] / p, helix.momentum[2] / p};
  }
  const double rx = x - circle->x;
  const double ry = y - circle->y;
  const double r = std::hypot(rx, ry);
  if (!(r > 0.)) {
    return {kNaN, kNaN, kNaN};
  }
  // Counter-clockwise, the tangent is the radius turned by +90 degrees.
  const double tx = -circle->sense * ry / r;
  const double ty = circle->sense * rx / r;
  return {tx * pt / p, ty * pt / p, helix.momentum[2] / p};
}

double helixMiss(const TruthHelix& helix, double bzKiloGauss, double x, double y)
{
  if (!followable(helix, bzKiloGauss, x, y)) {
    return kNaN;
  }
  if (!transverse(helix)) {
    return std::hypot(x - helix.position[0], y - helix.position[1]); // a line along z
  }
  const auto circle = circleOf(helix, bzKiloGauss);
  if (!circle) {
    return std::abs(helixImpact(helix, bzKiloGauss, x, y).d0);
  }
  return std::abs(std::hypot(x - circle->x, y - circle->y) - circle->radius);
}

TrackingQA trackingQA(std::span<const TruthRow> rows, const ParticleTable& truth, const TrackRows& tracks, std::span<const TrackKinematics> kinematics, const SeededBy& seeded, const ParticleExtrasTable& extras, const TrackRows* theirs, int iterations, const ComparisonSettings& settings)
{
  if (iterations < 1 || iterations > kQAMaxIterations) {
    throw std::invalid_argument(std::format("QA: {} iterations, 1 to {} supported", iterations, kQAMaxIterations));
  }
  if (kinematics.size() != tracks.size()) {
    throw std::invalid_argument(std::format("QA: {} tracks, {} kinematics", tracks.size(), kinematics.size()));
  }
  TrackingQA qa;

  // The tracks: each one's particle, a clone after the particle's first good
  // track in iteration and then track order.
  std::vector<std::size_t> order(tracks.size());
  for (std::size_t track = 0; track < order.size(); ++track) {
    order[track] = track;
    if (std::cmp_greater_equal(kinematics[track].iteration, iterations)) {
      throw std::invalid_argument(std::format("QA: a track of iteration {}, {} configured", kinematics[track].iteration, iterations));
    }
  }
  std::ranges::stable_sort(order, {}, [&](std::size_t track) { return kinematics[track].iteration; });
  std::unordered_map<std::uint64_t, std::uint64_t> foundBy;
  qa.tracks.resize(tracks.size());
  for (const auto track : order) {
    const auto& k = kinematics[track];
    const auto particle = goodParticle(rows, tracks[track]);
    auto& record = qa.tracks[track];
    record = {
      .pt = k.pt,
      .eta = k.eta,
      .phi = wrapPhi(static_cast<double>(k.phi)),
      .iteration = k.iteration,
      .fake = particle == kNoTrackParticle,
      .chi2 = k.chi2,
      .ndf = k.ndf,
    };
    if (!record.fake) {
      auto& mask = foundBy[particle];
      record.clone = mask != 0U;
      mask |= std::uint64_t{1} << k.iteration;
    }
  }

  // The other tracker's good particles.
  std::unordered_set<std::uint64_t> foundByTheirs;
  if (theirs != nullptr) {
    for (const auto& track : *theirs) {
      if (const auto particle = goodParticle(rows, track); particle != kNoTrackParticle) {
        foundByTheirs.insert(particle);
      }
    }
  }

  // The findable, by either rule.
  for (const auto& [id, layers] : particleLayers(rows, settings)) {
    const auto found = truth.find(id);
    if (found == truth.end()) {
      continue;
    }
    const auto more = extras.find(id);
    const ParticleExtras extra = more == extras.end() ? ParticleExtras{} : more->second;
    const bool findable = isFindable(found->second.pt, layers, settings);
    if (!findable && (static_cast<double>(found->second.pt) < settings.minPtGeV || std::cmp_less(extra.outwardLayers, settings.minOutwardLayers))) {
      continue;
    }
    const auto& particle = found->second;
    const auto entry = foundBy.find(id);
    const auto seeds = seeded.find(id);
    qa.particles.push_back({
      .pt = particle.pt,
      .eta = particle.eta,
      .phi = wrapPhi(static_cast<double>(particle.phi)),
      .origin = static_cast<std::uint8_t>(particle.origin),
      .foundBy = entry == foundBy.end() ? 0U : entry->second,
      .seededBy = seeds == seeded.end() ? 0U : seeds->second,
      .pdg = particle.pdg,
      .motherPdg = particle.motherPdg,
      .vz = particle.vz,
      .layers = layers.layers,
      .foundByO2 = foundByTheirs.contains(id),
      .findable = findable,
      .extras = extra,
    });
  }
  return qa;
}

void stream(o2::utils::TreeStreamRedirector& output, const TrackingQA& qa)
{
  for (const auto& p : qa.particles) {
    const auto& e = p.extras;
    output << "particles" << "pt=" << p.pt << "eta=" << p.eta << "phi=" << p.phi << "origin=" << p.origin << "foundBy=" << p.foundBy << "seededBy=" << p.seededBy << "pdg=" << p.pdg << "motherPdg=" << p.motherPdg << "vr=" << e.vr << "vz=" << p.vz << "d0=" << e.d0 << "z0=" << e.z0 << "openingAngle=" << e.openingAngle << "motherMiss=" << e.motherMiss << "layers=" << p.layers << "lostAt=" << e.lostAt << "pairLoss=" << e.pairLoss << "roadLoss=" << e.roadLoss << "fitLoss=" << e.fitLoss << "fitRoadLayers=" << e.fitRoadLayers << "roadLosses=" << e.roadLosses << "fitLosses=" << e.fitLosses << "raisedFit=" << e.raisedFit << "foundByO2=" << p.foundByO2 << "findable=" << p.findable << "outwardLayers=" << e.outwardLayers << "\n";
  }
  for (const auto& t : qa.tracks) {
    output << "tracks" << "pt=" << t.pt << "eta=" << t.eta << "phi=" << t.phi << "iteration=" << t.iteration << "fake=" << t.fake << "clone=" << t.clone << "chi2=" << t.chi2 << "ndf=" << t.ndf << "\n";
  }
}

void streamResolution(o2::utils::TreeStreamRedirector& output, std::span<const QAResolution> entries)
{
  for (const auto& r : entries) {
    output << "resolution" << "pt=" << r.pt << "eta=" << r.eta << "phi=" << r.phi << "charge=" << r.charge << "pdg=" << r.pdg << "origin=" << r.origin << "iteration=" << r.iteration << "kinks=" << r.kinks << "clusters=" << r.clusters << "layer=" << r.layer << "dLoc0=" << r.residual[0] << "dLoc1=" << r.residual[1] << "dPhi=" << r.residual[2] << "dTheta=" << r.residual[3] << "dQOverP=" << r.residual[4] << "sLoc0=" << r.sigma[0] << "sLoc1=" << r.sigma[1] << "sPhi=" << r.sigma[2] << "sTheta=" << r.sigma[3] << "sQOverP=" << r.sigma[4] << "fitPt=" << r.fitPt << "dPtRel=" << r.ptResidual << "sPt=" << r.ptSigma << "dD0=" << r.dca[0] << "dZ0=" << r.dca[1] << "sD0=" << r.dcaSigma[0] << "sZ0=" << r.dcaSigma[1] << "scatD0=" << r.dcaScattering[0] << "scatZ0=" << r.dcaScattering[1] << "productionOffset=" << r.productionOffset << "\n";
  }
}

void streamOracle(o2::utils::TreeStreamRedirector& output, std::span<const QAOracle> entries)
{
  for (const auto& o : entries) {
    output << "oracle" << "pdg=" << o.pdg << "charge=" << o.charge << "origin=" << o.origin << "eta=" << o.eta << "layer=" << o.layer << "kinks=" << o.kinks << "clusters=" << o.clusters << "productionOffset=" << o.productionOffset << "p=" << o.momentum << "pOuter=" << o.pOuter << "variant=" << o.variant << "status=" << o.status << "chi2=" << o.chi2 << "ndf=" << o.ndf << "dLoc0=" << o.residual[0] << "dLoc1=" << o.residual[1] << "dPhi=" << o.residual[2] << "dTheta=" << o.residual[3] << "dQOverP=" << o.residual[4] << "sLoc0=" << o.sigma[0] << "sLoc1=" << o.sigma[1] << "sPhi=" << o.sigma[2] << "sTheta=" << o.sigma[3] << "sQOverP=" << o.sigma[4] << "dPtRel=" << o.ptResidual << "\n";
  }
}

void streamMomenta(o2::utils::TreeStreamRedirector& output, std::span<const QAMomenta> entries)
{
  for (const auto& m : entries) {
    std::vector<float> pAtLayer = m.pAtLayer;
    output << "momenta" << "pdg=" << m.pdg << "origin=" << m.origin << "pt=" << m.pt << "eta=" << m.eta << "p0=" << m.p0 << "found=" << m.found << "lostAt=" << m.lostAt << "pairLoss=" << m.pairLoss << "layers=" << m.layers << "fakeWrong=" << m.fakeWrong << "fakeWrongLayer=" << m.fakeWrongLayer << "pAtLayer=" << pAtLayer << "\n";
  }
}

void streamRun(o2::utils::TreeStreamRedirector& output, int iterations, EtaRange eta, const ComparisonSettings& settings, std::size_t events) { output << "run" << "iterations=" << iterations << "etaLow=" << eta.low << "etaHigh=" << eta.high << "minLayers=" << settings.minLayers << "minPtGeV=" << settings.minPtGeV << "minOutwardLayers=" << settings.minOutwardLayers << "events=" << static_cast<std::uint64_t>(events) << "\n"; }

} // namespace o2::gtrk::study
