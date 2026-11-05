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

/// \file GeneratedHelixFixture.h
/// \brief Shared generated-helix support for tracking-cut tests.

#ifndef GTRK_TEST_GENERATEDHELIXFIXTURE_H_
#define GTRK_TEST_GENERATEDHELIXFIXTURE_H_

#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/TripletCut.h"

#include "GridTestSupport.h"
#include "UniformMaterial.h"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <random>
#include <stdexcept>

namespace o2::gtrk
{

template <DetectorLayout Det>
struct GeneratedHelixFixture {
  static constexpr double kMaterialX2X0 = 0.02;

  /// A position on the helix and its direction there.
  struct State {
    double x = 0.;
    double y = 0.;
    double z = 0.;
    double phi = 0.;   ///< transverse direction
    double slope = 0.; ///< tan(lambda)
  };

  struct Errors {
    double varZ = 0.;
    double varR = 0.;
    double varRPhi = 0.;
  };

  static constexpr Errors kPixelErrors{.varZ = 2.25e-6, .varR = 1.e-6, .varRPhi = 2.25e-6};
  static constexpr Errors kLongZErrors{.varZ = 1.e-2, .varR = 1.e-6, .varRPhi = 4.e-6};
  static constexpr Errors kDiscErrors{.varZ = 1.e-6, .varR = 2.25e-6, .varRPhi = 2.25e-6};

  /// \a kinkSigma is the iteration's kink allowance, 1 for none.
  [[nodiscard]] static IterationParameters<Det> parametersAt(double minPtGeV, double kinkSigma = 1.)
  {
    IterationSettings settings;
    settings.windows.minPtGeV = minPtGeV;
    settings.kinkSigma = kinkSigma;
    return iterationParameters<Det>(settings, VertexConstraint{}, scatteringTable<Det>(uniformLayerMaterial<Det>(kMaterialX2X0)));
  }

  [[nodiscard]] static State propagate(const State& start, double curvature, double arc)
  {
    return {.x = start.x + ((std::sin(start.phi + (curvature * arc)) - std::sin(start.phi)) / curvature),
            .y = start.y + ((std::cos(start.phi) - std::cos(start.phi + (curvature * arc))) / curvature),
            .z = start.z + (start.slope * arc),
            .phi = start.phi + (curvature * arc),
            .slope = start.slope};
  }

  /// To the layer: a barrel's radius by bisection on the outgoing half-turn,
  /// or a disc's z directly.
  [[nodiscard]] static State toLayer(const State& start, double curvature, LayerID layer)
  {
    if (!Detector<Det>::isBarrelLayer(layer)) {
      return propagate(start, curvature, (discZ<Det>(layer) - start.z) / start.slope);
    }
    const double radius = layerDistance<Det>(layer);
    double low = 0.;
    double high = std::numbers::pi_v<double> / curvature;
    for (int i = 0; i < 100; ++i) {
      const double middle = 0.5 * (low + high);
      const State at = propagate(start, curvature, middle);
      (std::hypot(at.x, at.y) < radius ? low : high) = middle;
    }
    return propagate(start, curvature, 0.5 * (low + high));
  }

  /// Scatter at the layer at the true momentum: Highland with
  /// p = pT sqrt(1 + slope^2) and the layer's path.
  [[nodiscard]] static State scatter(const State& at, LayerID layer, double ptGeV, const IterationParameters<Det>& parameters, std::mt19937_64& generator)
  {
    std::normal_distribution<double> normal;
    const double cosine = 1. / std::sqrt(1. + (at.slope * at.slope));
    const double sine = std::abs(at.slope) * cosine;
    const double path = Detector<Det>::isBarrelLayer(layer) ? 1. / cosine : 1. / sine;
    const auto k = static_cast<double>(parameters.scattering.layers[layer].k);
    const double theta0 = std::sqrt(k * path) * cosine / ptGeV;
    State kinked = at;
    kinked.slope = std::tan(std::atan(at.slope) + (normal(generator) * theta0));
    kinked.phi += normal(generator) * theta0 / cosine;
    return kinked;
  }

  [[nodiscard]] static TimedPoint measure(const State& state, const Errors& errors, double timeNS, std::mt19937_64& generator)
  {
    std::normal_distribution<double> normal;
    const double r = std::hypot(state.x, state.y);
    const double radial = normal(generator) * std::sqrt(errors.varR);
    const double tangential = normal(generator) * std::sqrt(errors.varRPhi);
    const double x = state.x + (((radial * state.x) - (tangential * state.y)) / r);
    const double y = state.y + (((radial * state.y) + (tangential * state.x)) / r);
    const double z = state.z + (normal(generator) * std::sqrt(errors.varZ));
    return {.point = {.x = static_cast<Scalar>(x),
                      .y = static_cast<Scalar>(y),
                      .z = static_cast<Scalar>(z),
                      .r = static_cast<Scalar>(std::hypot(x, y)),
                      .varZ = static_cast<Scalar>(errors.varZ),
                      .varR = static_cast<Scalar>(errors.varR),
                      .varRPhi = static_cast<Scalar>(errors.varRPhi)},
            .time = test::toyTimeBox(timeNS)};
  }
};

} // namespace o2::gtrk

#endif // GTRK_TEST_GENERATEDHELIXFIXTURE_H_
