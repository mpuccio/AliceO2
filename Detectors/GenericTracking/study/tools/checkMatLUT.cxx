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

/// \file checkMatLUT.cxx
/// \brief Hold a material LUT to the geometry it was built from, along the
///        paths simulated particles take.
///
///     o2-gtrk-check-matlut --detector TRKFT3 --sim <dir> --lut gtrk_matlut_trk.root [--events 200]
///
/// A LUT holds each cell's material evenly, so a thin dense object -- a
/// layer's metal stack, a petal's side wall, a disc's support -- lands on
/// every track through its cell, crossing it or not. The radial total cannot
/// see that; the fit can. This tool measures what the fit
/// is given against what the particles crossed.
///
/// Each MC particle's crossings of the layers, in time order, collapsed to the
/// first crossing of each layer and cut where a layer comes back (a looper's
/// second leg), give the steps the fit takes: from the production vertex to
/// the first crossing, then from one crossing to the next, each from the
/// crossing's entry point to the next one's. The path between is the true
/// helix, without a field: the transverse circle tangent to the momentum at
/// the first point through the second, z linear in the arc. A step whose
/// circle turns the momentum more than --max-turn away from the true one at
/// its end is left out: more than half a turn, or a hard scatter. Along the
/// path, in --pieces straight pieces, the material from the loaded geometry
/// (`GeometryManager::meanMaterialBudget`, as the LUT builder fills its cells)
/// and from the LUT (`MatLayerCylSet::getMatBudget`, as the fit reads it).
///
/// The tables:
///
/// - **nodes**: per layer, the material the fit applies at the crossing --
///   the far half of the step in and the near half of the step out (its two
///   halves). The material of a layer's own band then counts in full, however
///   the band splits it between the two steps. The ratio of the means, and of
///   the particles' own ratios the median and the 5 % and 95 % quantiles: a
///   cell smeared over its neighbours widens these, a wrong total moves them.
/// - **steps**: per pair of layers, the material, the share of it in the
///   step's near half (where it sits), and the true change of the polar angle
///   (which the solenoid leaves alone) in Highland's theta0 of the geometry
///   and of the LUT: its core, 1.4826 MAD, near 1 for the geometry's says
///   the comparison holds; the LUT's is what the fit's angle pulls see.
///   Electrons are left out of this one (bremsstrahlung).
/// - **worst cells**: the nodes binned by layer and azimuth and by layer and
///   z (barrel) or r (disc), the cells with at least --min-cell nodes whose
///   ratio of the means is furthest from 1. A wall, a support or a stave
///   overlap shows up here without being looked for.
///
/// --output writes the steps as a tree for plots.

#include "GenericTracking/Study/Alice/Hits.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Study/Statistics.h"

#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/MatLayerCylSet.h"
#include "SimulationDataFormat/MCTrack.h"
#include "Steer/MCKinematicsReader.h"

#include <TFile.h>
#include <TGeoManager.h>
#include <TTree.h>

#include <boost/program_options.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <format>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace bpo = boost::program_options;

namespace
{

using o2::gtrk::LayerID;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

struct Options {
  std::string detector = "TRKFT3";
  std::string sim;
  std::string lut;
  std::string geometry;
  std::string output;
  int events = 200;
  double minPt = 0.1;
  bool secondaries = false;
  int pieces = 8;
  double maxTurn = 0.2;
  int minSteps = 200;
  int minCell = 30;
  int worst = 25;
  double phiBinDeg = 1.;
  double coordinateBin = 2.;
  int threads = -1;
};

/// A point on the path [cm] with the momentum there [GeV].
struct Point {
  std::array<double, 3> x{};
  std::array<double, 3> p{};
};

/// One step: from a crossing's entry (or the vertex) to the next one's.
struct Step {
  int event = 0;
  int track = 0;
  int pdg = 0;
  double pt = 0.;
  double eta = 0.;
  double mass = 0.;
  int from = -1; ///< layer, -1 the vertex
  int to = -1;
  Point a;
  Point b;
  double turn = 0.;            ///< the circle's direction at b against the true one [rad]
  std::array<double, 2> geo{}; ///< x/X0, near and far half
  std::array<double, 2> lut{};
  bool lastOfLeg = false;
};

/// The step's path, \a pieces + 1 points: the transverse circle tangent to a's
/// momentum through b, z linear in the arc. Also the circle's direction at b.
std::vector<std::array<double, 3>> path(const Step& step, int pieces, double& directionAtB)
{
  const auto& [ax, ap] = step.a;
  const auto& bx = step.b.x;
  const double tNorm = std::hypot(ap[0], ap[1]);
  const std::array<double, 2> t{ap[0] / tNorm, ap[1] / tNorm};
  const std::array<double, 2> n{-t[1], t[0]};
  const std::array<double, 2> d{bx[0] - ax[0], bx[1] - ax[1]};
  const double dd = (d[0] * d[0]) + (d[1] * d[1]);
  const double dn = (d[0] * n[0]) + (d[1] * n[1]);
  std::vector<std::array<double, 3>> points;
  points.reserve(static_cast<std::size_t>(pieces) + 1U);
  if (std::abs(dn) < 1e-9 * dd) {
    for (int i = 0; i <= pieces; ++i) {
      const double f = static_cast<double>(i) / pieces;
      points.push_back({ax[0] + (f * d[0]), ax[1] + (f * d[1]), ax[2] + (f * (bx[2] - ax[2]))});
    }
    directionAtB = std::atan2(t[1], t[0]);
    return points;
  }
  const double radius = dd / (2. * dn); // signed, along n
  const std::array<double, 2> centre{ax[0] + (radius * n[0]), ax[1] + (radius * n[1])};
  const double start = std::atan2(ax[1] - centre[1], ax[0] - centre[0]);
  double sweep = std::atan2(bx[1] - centre[1], bx[0] - centre[0]) - start;
  // the direction of travel: counter-clockwise for a positive radius
  if (radius > 0.) {
    while (sweep < 0.) {
      sweep += 2. * std::numbers::pi;
    }
  } else {
    while (sweep > 0.) {
      sweep -= 2. * std::numbers::pi;
    }
  }
  const double r = std::abs(radius);
  for (int i = 0; i <= pieces; ++i) {
    const double f = static_cast<double>(i) / pieces;
    const double angle = start + (f * sweep);
    points.push_back({centre[0] + (r * std::cos(angle)), centre[1] + (r * std::sin(angle)), ax[2] + (f * (bx[2] - ax[2]))});
  }
  directionAtB = std::atan2(t[1], t[0]) + sweep;
  return points;
}

double wrapped(double angle) { return std::remainder(angle, 2. * std::numbers::pi); }

/// The step's material along its path, near and far half, from the geometry
/// (on \a navigator) and from \a lut.
void measure(Step& step, int pieces, const o2::base::MatLayerCylSet& lut, TGeoNavigator* navigator)
{
  double direction = 0.;
  const auto points = path(step, pieces, direction);
  step.turn = std::abs(wrapped(direction - std::atan2(step.b.p[1], step.b.p[0])));
  for (int i = 0; i < pieces; ++i) {
    const auto& p = points[static_cast<std::size_t>(i)];
    const auto& q = points[static_cast<std::size_t>(i) + 1U];
    const std::size_t half = 2 * i < pieces ? 0U : 1U;
    const auto geo = o2::base::GeometryManager::meanMaterialBudget(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]), static_cast<float>(q[0]), static_cast<float>(q[1]), static_cast<float>(q[2]), navigator);
    if (geo.length > 0.F) {
      step.geo[half] += static_cast<double>(geo.meanX2X0);
    }
    const auto fromLut = lut.getMatBudget(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]), static_cast<float>(q[0]), static_cast<float>(q[1]), static_cast<float>(q[2]));
    step.lut[half] += static_cast<double>(fromLut.meanX2X0);
  }
}

/// Highland's theta0, the tracking's formula (`Tracking/Material.h`,
/// `highlandTheta0`) in double; x/X0 <= 0 gives 0.
double highland(double x2X0, double betaP) { return x2X0 > 0. ? 0.0136 / betaP * std::sqrt(x2X0) * (1. + (0.038 * std::log(x2X0))) : 0.; }

/// beta p of momentum \a p and \a mass [GeV], the tracking's `betaMomentum` in double.
double betaMomentum(const std::array<double, 3>& momentum, double mass)
{
  const double p = std::sqrt((momentum[0] * momentum[0]) + (momentum[1] * momentum[1]) + (momentum[2] * momentum[2]));
  return p * p / std::sqrt((p * p) + (mass * mass));
}

double polar(const std::array<double, 3>& p) { return std::atan2(std::hypot(p[0], p[1]), p[2]); }

double quantile(std::vector<double> v, double q)
{
  if (v.empty()) {
    return kNaN;
  }
  return o2::gtrk::study::quantile(std::span{v}, q);
}

/// 1.4826 MAD.
double core(std::vector<double> v)
{
  if (v.size() < 20U) {
    return kNaN;
  }
  const double median = quantile(v, 0.5);
  for (auto& x : v) {
    x = std::abs(x - median);
  }
  return 1.4826 * quantile(v, 0.5);
}

double beyond(const std::vector<double>& v, double width, double k)
{
  if (v.empty() || !std::isfinite(width)) {
    return kNaN;
  }
  const auto n = std::ranges::count_if(v, [&](double x) { return std::abs(x) > k * width; });
  return 100. * static_cast<double>(n) / static_cast<double>(v.size());
}

template <o2::gtrk::DetectorLayout Det>
std::string layerName(int layer)
{
  if (layer < 0) {
    return "vertex";
  }
  return std::format("{}{}", o2::gtrk::Detector<Det>::isBarrelLayer(static_cast<LayerID>(layer)) ? "B" : "D", layer);
}

/// The steps of every particle above --min-pt in the first --events events;
/// \a events, the events read.
template <o2::gtrk::DetectorLayout Det>
std::vector<Step> collect(const Options& options, const o2::gtrk::Geometry<Det>& geometry, int& events)
{
  constexpr auto kChips = o2::gtrk::study::chipsOf<Det>();
  const std::string hitName = [] {
    if constexpr (kChips == o2::gtrk::study::Chips::kITS) {
      return "ITS";
    } else if constexpr (kChips == o2::gtrk::study::Chips::kMFT) {
      return "MFT";
    } else {
      return "TRK";
    }
  }();
  const std::string discPath = kChips == o2::gtrk::study::Chips::kTRK ? options.sim + "/sgn_HitsFT3.root" : std::string{};
  o2::gtrk::study::HitFile hits(options.sim + "/sgn_Hits" + hitName + ".root", kChips, discPath);
  // The kinematics through O2's reader, each event released once its steps are taken.
  if (!std::filesystem::exists(options.sim + "/sgn_Kine.root")) {
    throw std::runtime_error("no " + options.sim + "/sgn_Kine.root");
  }
  o2::steer::MCKinematicsReader kine(options.sim + "/sgn", o2::steer::MCKinematicsReader::Mode::kMCKine);
  std::vector<Step> steps;
  events = std::min({options.events, hits.nEvents(), static_cast<int>(kine.getNEvents(0))});
  for (int event = 0; event < events; ++event) {
    const auto& tracks = kine.getTracks(0, event);
    std::map<int, std::vector<const o2::gtrk::study::HitSegment*>> byTrack;
    for (const auto& segment : hits.event(event)) {
      byTrack[segment.track].push_back(&segment);
    }
    for (auto& [id, segments] : byTrack) {
      if (id < 0 || std::cmp_greater_equal(id, tracks.size())) {
        continue;
      }
      const auto& mc = tracks[static_cast<std::size_t>(id)];
      if ((!options.secondaries && !mc.isPrimary()) || mc.GetPt() < options.minPt) {
        continue;
      }
      std::ranges::sort(segments, {}, &o2::gtrk::study::HitSegment::time);
      Point previous{.x = {mc.Vx(), mc.Vy(), mc.Vz()}, .p = {mc.Px(), mc.Py(), mc.Pz()}};
      int previousLayer = -1;
      std::vector<bool> seen(static_cast<std::size_t>(o2::gtrk::Detector<Det>::nLayers), false);
      const auto first = steps.size();
      for (const auto* segment : segments) {
        const int layer = geometry.sensor(segment->chip).layer;
        if (layer == previousLayer) {
          continue; // a second chip of the same layer: an overlap
        }
        if (seen[static_cast<std::size_t>(layer)]) {
          break; // the layer again: a looper's next leg
        }
        seen[static_cast<std::size_t>(layer)] = true;
        const Point entry{.x = segment->start, .p = segment->momentum};
        steps.push_back({
          .event = event,
          .track = id,
          .pdg = mc.GetPdgCode(),
          .pt = mc.GetPt(),
          .eta = mc.GetEta(),
          .mass = mc.GetMass(),
          .from = previousLayer,
          .to = layer,
          .a = previous,
          .b = entry,
        });
        previous = entry;
        previousLayer = layer;
      }
      if (steps.size() > first) {
        steps.back().lastOfLeg = true;
      }
    }
    kine.releaseTracksForSourceAndEvent(0, event);
  }
  return steps;
}

struct Sums {
  std::vector<double> geo;
  std::vector<double> lut;
  double geoSum = 0.;
  double lutSum = 0.;
  void add(double g, double l)
  {
    geo.push_back(g);
    lut.push_back(l);
    geoSum += g;
    lutSum += l;
  }
  [[nodiscard]] std::size_t size() const { return geo.size(); }
  [[nodiscard]] double ratio() const { return geoSum > 0. ? lutSum / geoSum : kNaN; }
  [[nodiscard]] std::vector<double> ratios() const
  {
    std::vector<double> r;
    r.reserve(geo.size());
    for (std::size_t i = 0; i < geo.size(); ++i) {
      if (geo[i] > 1e-6) {
        r.push_back(lut[i] / geo[i]);
      }
    }
    return r;
  }
};

template <o2::gtrk::DetectorLayout Det>
void report(const Options& options, const std::vector<Step>& steps)
{
  // nodes: the far half in and the near half out
  std::map<int, Sums> nodes;
  std::map<std::tuple<int, int, bool>, Sums> cells;
  for (std::size_t i = 0; i + 1 < steps.size(); ++i) {
    const auto& in = steps[i];
    const auto& out = steps[i + 1];
    if (in.lastOfLeg || in.event != out.event || in.track != out.track || in.turn > options.maxTurn || out.turn > options.maxTurn) {
      continue;
    }
    const double g = in.geo[1] + out.geo[0];
    const double l = in.lut[1] + out.lut[0];
    nodes[in.to].add(g, l);
    const auto& x = in.b.x;
    double phi = std::atan2(x[1], x[0]) * 180. / std::numbers::pi;
    phi = phi < 0. ? phi + 360. : phi;
    const bool barrel = o2::gtrk::Detector<Det>::isBarrelLayer(static_cast<LayerID>(in.to));
    const double coordinate = barrel ? x[2] : std::hypot(x[0], x[1]);
    cells[{in.to, static_cast<int>(std::floor(phi / options.phiBinDeg)), true}].add(g, l);
    cells[{in.to, static_cast<int>(std::floor(coordinate / options.coordinateBin)), false}].add(g, l);
  }
  std::cout << std::format("\n== nodes: the material the fit applies at each layer's crossing, x/X0 [%]\n");
  std::cout << std::format("{:>6} {:>8} | {:>8} {:>8} {:>7} | LUT/geometry per node: {:>6} {:>6} {:>6} | {:>9}\n", "layer", "n", "geometry", "LUT", "ratio", "5 %", "median", "95 %", "off 25 %");
  for (const auto& [layer, sums] : nodes) {
    if (std::cmp_less(sums.size(), options.minSteps)) {
      continue;
    }
    const auto r = sums.ratios();
    const auto off = std::ranges::count_if(r, [](double x) { return x < 0.8 || x > 1.25; });
    std::cout << std::format("{:>6} {:>8} | {:8.4f} {:8.4f} {:7.3f} | {:>22} {:6.3f} {:6.3f} {:6.3f} | {:8.1f}%\n", layerName<Det>(layer), sums.size(), 100. * sums.geoSum / static_cast<double>(sums.size()), 100. * sums.lutSum / static_cast<double>(sums.size()), sums.ratio(), "", quantile(r, 0.05), quantile(r, 0.5), quantile(r, 0.95), r.empty() ? kNaN : 100. * static_cast<double>(off) / static_cast<double>(r.size()));
  }

  struct StepStats {
    Sums material;
    double geoNear = 0.;
    double lutNear = 0.;
    std::vector<double> pullGeo;
    std::vector<double> pullLut;
  };
  std::map<std::pair<int, int>, StepStats> pairs;
  std::size_t turned = 0U;
  for (const auto& step : steps) {
    if (step.turn > options.maxTurn) {
      ++turned;
      continue;
    }
    auto& s = pairs[{step.from, step.to}];
    s.material.add(step.geo[0] + step.geo[1], step.lut[0] + step.lut[1]);
    s.geoNear += step.geo[0];
    s.lutNear += step.lut[0];
    if (std::abs(step.pdg) == 11) {
      continue;
    }
    const double betaP = betaMomentum(step.a.p, step.mass);
    const double dTheta = polar(step.b.p) - polar(step.a.p);
    // the fit's two halves
    const double geoTheta = std::hypot(highland(step.geo[0], betaP), highland(step.geo[1], betaP));
    const double lutTheta = std::hypot(highland(step.lut[0], betaP), highland(step.lut[1], betaP));
    if (geoTheta > 0. && lutTheta > 0.) {
      s.pullGeo.push_back(dTheta / geoTheta);
      s.pullLut.push_back(dTheta / lutTheta);
    }
  }
  std::cout << std::format("\n== steps ({} of {} left out: the circle turned more than {} rad from the truth), x/X0 [%]\n", turned, steps.size(), options.maxTurn);
  std::cout << std::format("{:>13} {:>8} | {:>8} {:>8} {:>7} | near half: {:>5} {:>5} | dtheta/theta0 core, >3: geometry {:>5} {:>5}  LUT {:>5} {:>5}\n", "step", "n", "geometry", "LUT", "ratio", "geo", "LUT", "", "", "", "");
  for (const auto& [key, s] : pairs) {
    if (std::cmp_less(s.material.size(), options.minSteps)) {
      continue;
    }
    const double cg = core(s.pullGeo);
    const double cl = core(s.pullLut);
    const auto n = static_cast<double>(s.material.size());
    std::cout << std::format("{:>13} {:>8} | {:8.4f} {:8.4f} {:7.3f} | {:>11} {:5.2f} {:5.2f} | {:>32} {:5.3f} {:5.2f}  {:>3} {:5.3f} {:5.2f}\n", layerName<Det>(key.first) + "->" + layerName<Det>(key.second), s.material.size(), 100. * s.material.geoSum / n, 100. * s.material.lutSum / n, s.material.ratio(), "", s.material.geoSum > 0. ? s.geoNear / s.material.geoSum : kNaN, s.material.lutSum > 0. ? s.lutNear / s.material.lutSum : kNaN, "", cg, beyond(s.pullGeo, cg, 3.), "", cl, beyond(s.pullLut, cl, 3.));
  }

  std::vector<std::pair<double, std::tuple<int, int, bool>>> ranked;
  for (const auto& [key, sums] : cells) {
    if (std::cmp_less(sums.size(), options.minCell) || !(sums.geoSum > 0.)) {
      continue;
    }
    ranked.emplace_back(std::abs(std::log(sums.ratio())), key);
  }
  std::ranges::sort(ranked, [](const auto& a, const auto& b) { return a.first > b.first; });
  std::cout << std::format("\n== worst cells: nodes by layer and azimuth ({} deg) or z / r ({} cm), at least {} nodes, x/X0 [%]\n", options.phiBinDeg, options.coordinateBin, options.minCell);
  std::cout << std::format("{:>6} {:>24} {:>6} | {:>8} {:>8} {:>7} | LUT/geometry per node {:>6} {:>6} {:>6}\n", "layer", "cell", "n", "geometry", "LUT", "ratio", "5 %", "median", "95 %");
  for (std::size_t i = 0; i < ranked.size() && std::cmp_less(i, options.worst); ++i) {
    const auto& [layer, bin, azimuth] = ranked[i].second;
    const auto& sums = cells.at(ranked[i].second);
    const bool barrel = o2::gtrk::Detector<Det>::isBarrelLayer(static_cast<LayerID>(layer));
    const double width = azimuth ? options.phiBinDeg : options.coordinateBin;
    const char* axis = barrel ? "z" : "r";
    if (azimuth) {
      axis = "phi";
    }
    const std::string cell = std::format("{} {:.1f} .. {:.1f} {}", axis, bin * width, (bin + 1) * width, azimuth ? "deg" : "cm");
    const auto r = sums.ratios();
    const auto n = static_cast<double>(sums.size());
    std::cout << std::format("{:>6} {:>24} {:>6} | {:8.4f} {:8.4f} {:7.3f} | {:>21} {:6.3f} {:6.3f} {:6.3f}\n", layerName<Det>(layer), cell, sums.size(), 100. * sums.geoSum / n, 100. * sums.lutSum / n, sums.ratio(), "", quantile(r, 0.05), quantile(r, 0.5), quantile(r, 0.95));
  }
}

void write(const std::string& path, const std::vector<Step>& steps)
{
  const TFile file(path.c_str(), "RECREATE");
  TTree tree("steps", "o2-gtrk-check-matlut steps");
  Step s;
  float ax = 0.F;
  float ay = 0.F;
  float az = 0.F;
  float bx = 0.F;
  float by = 0.F;
  float bz = 0.F;
  float pt = 0.F;
  float eta = 0.F;
  float turn = 0.F;
  float dTheta = 0.F;
  float betaP = 0.F;
  std::array<float, 2> geo{};
  std::array<float, 2> lut{};
  tree.Branch("event", &s.event);
  tree.Branch("track", &s.track);
  tree.Branch("pdg", &s.pdg);
  tree.Branch("from", &s.from);
  tree.Branch("to", &s.to);
  tree.Branch("pt", &pt);
  tree.Branch("eta", &eta);
  tree.Branch("ax", &ax);
  tree.Branch("ay", &ay);
  tree.Branch("az", &az);
  tree.Branch("bx", &bx);
  tree.Branch("by", &by);
  tree.Branch("bz", &bz);
  tree.Branch("turn", &turn);
  tree.Branch("dTheta", &dTheta);
  tree.Branch("betaP", &betaP);
  tree.Branch("geo", geo.data(), "geo[2]/F");
  tree.Branch("lut", lut.data(), "lut[2]/F");
  for (const auto& step : steps) {
    s = step;
    ax = static_cast<float>(step.a.x[0]);
    ay = static_cast<float>(step.a.x[1]);
    az = static_cast<float>(step.a.x[2]);
    bx = static_cast<float>(step.b.x[0]);
    by = static_cast<float>(step.b.x[1]);
    bz = static_cast<float>(step.b.x[2]);
    pt = static_cast<float>(step.pt);
    eta = static_cast<float>(step.eta);
    turn = static_cast<float>(step.turn);
    dTheta = static_cast<float>(polar(step.b.p) - polar(step.a.p));
    betaP = static_cast<float>(betaMomentum(step.a.p, step.mass));
    geo = {static_cast<float>(step.geo[0]), static_cast<float>(step.geo[1])};
    lut = {static_cast<float>(step.lut[0]), static_cast<float>(step.lut[1])};
    tree.Fill();
  }
  tree.Write();
}

template <o2::gtrk::DetectorLayout Det>
int runFor(const Options& options)
{
  // O2 exports the ideal geometry before applying alignment, then transports
  // with the aligned geometry (O2MCApplicationBase::MisalignGeometry).
  auto geometryPath = options.geometry;
  if (geometryPath.empty()) {
    const auto aligned = options.sim + "/sgn_geometry-aligned.root";
    geometryPath = std::filesystem::exists(aligned) ? aligned : options.sim + "/sgn_geometry.root";
  }
  o2::base::GeometryManager::loadGeometry(geometryPath);
  const auto geometry = std::make_unique<const o2::gtrk::Geometry<Det>>(o2::gtrk::alice::placeChecked<Det>());
  const std::unique_ptr<o2::base::MatLayerCylSet> lut{o2::base::MatLayerCylSet::loadFromFile(options.lut)};
  if (!lut) {
    throw std::runtime_error("cannot load the LUT " + options.lut);
  }
  int events = 0;
  auto steps = collect<Det>(options, *geometry, events);
  std::cout << std::format("{} steps of the particles above {} GeV{} in {} events\n", steps.size(), options.minPt, options.secondaries ? "" : " (primaries)", events);

  int threads = options.threads < 0 ? static_cast<int>(std::thread::hardware_concurrency()) : options.threads;
  threads = std::max(threads, 1);
  gGeoManager->SetMaxThreads(threads);
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> workers;
  workers.reserve(static_cast<std::size_t>(threads));
  for (int thread = 0; thread < threads; ++thread) {
    workers.emplace_back([&] {
      TGeoNavigator* navigator = gGeoManager->AddNavigator();
      for (std::size_t i = next++; i < steps.size(); i = next++) {
        measure(steps[i], options.pieces, *lut, navigator);
      }
      gGeoManager->RemoveNavigator(navigator);
    });
  }
  for (auto& worker : workers) {
    worker.join();
  }
  report<Det>(options, steps);
  if (!options.output.empty()) {
    write(options.output, steps);
  }
  return 0;
}

/// The body of main, so that clang-tidy can see that nothing escapes it.
int run(int argc, char** argv)
{
  Options options;
  bpo::options_description description("Hold a material LUT to the geometry along simulated particles' paths");
  description.add_options()                                                                                                                            //
    ("help,h", "this message")                                                                                                                         //
    ("detector,d", bpo::value(&options.detector)->default_value(options.detector), "ITS, MFT or TRKFT3")                                               //
    ("sim", bpo::value(&options.sim), "the simulation's directory: sgn_Hits*.root, sgn_Kine.root, sgn_geometry[-aligned].root")                        //
    ("lut", bpo::value(&options.lut), "the material LUT")                                                                                              //
    ("geometry", bpo::value(&options.geometry), "override geometry; default <sim>/sgn_geometry-aligned.root if present, otherwise sgn_geometry.root")  //
    ("output,o", bpo::value(&options.output), "write the steps to this ROOT file")                                                                     //
    ("events,n", bpo::value(&options.events)->default_value(options.events), "the first MC events")                                                    //
    ("min-pt", bpo::value(&options.minPt)->default_value(options.minPt), "the particles' pT floor [GeV]")                                              //
    ("secondaries", bpo::bool_switch(&options.secondaries), "every particle, not only primaries")                                                      //
    ("pieces", bpo::value(&options.pieces)->default_value(options.pieces), "straight pieces per step (even)")                                          //
    ("max-turn", bpo::value(&options.maxTurn)->default_value(options.maxTurn), "leave out steps whose circle misses the true direction by more [rad]") //
    ("min-steps", bpo::value(&options.minSteps)->default_value(options.minSteps), "the fewest steps or nodes a table row shows")                       //
    ("min-cell", bpo::value(&options.minCell)->default_value(options.minCell), "the fewest nodes a worst cell holds")                                  //
    ("worst", bpo::value(&options.worst)->default_value(options.worst), "how many worst cells")                                                        //
    ("phi-bin", bpo::value(&options.phiBinDeg)->default_value(options.phiBinDeg), "the worst cells' azimuth bin [deg]")                                //
    ("coordinate-bin", bpo::value(&options.coordinateBin)->default_value(options.coordinateBin), "the worst cells' z (barrel) or r (disc) bin [cm]")   //
    ("threads,j", bpo::value(&options.threads)->default_value(options.threads), "threads; negative: all");
  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(argc, argv, description), vm);
  bpo::notify(vm);
  if (vm.contains("help")) {
    std::cout << description << '\n';
    return 0;
  }
  const char* missing = options.sim.empty() ? "--sim is required" : options.lut.empty() ? "--lut is required"
                                                                                        : nullptr;
  if (missing == nullptr && (options.pieces < 2 || options.pieces % 2 != 0)) {
    missing = "--pieces must be even and at least 2";
  }
  if (missing != nullptr) {
    std::cerr << missing << "\n\n"
              << description << '\n';
    return 1;
  }
  if (options.detector == "TRKFT3") {
    return runFor<o2::gtrk::detectors::TRKFT3>(options);
  }
  if (options.detector == "ITS") {
    return runFor<o2::gtrk::detectors::ITS>(options);
  }
  if (options.detector == "MFT") {
    return runFor<o2::gtrk::detectors::MFT>(options);
  }
  std::cerr << "unknown detector " << options.detector << '\n';
  return 1;
}

} // namespace

int main(int argc, char** argv)
{
  try {
    return run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
