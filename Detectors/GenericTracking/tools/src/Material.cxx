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

#include "GenericTracking/Tools/Material.h"

#include "DetectorsBase/MatLayerCylSet.h"
#include "TRKBase/GeometryTGeo.h"

#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoNavigator.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <format>
#include <functional>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace o2::gtrk::tools
{

int fillThreads(int threads) { return std::max(threads < 0 ? o2::base::MatLayerCylSet::getNThreadsFromEnv() : threads, 1); }

void forEachOnNavigators(std::size_t count, int threads, const std::function<void(std::size_t, TGeoNavigator*)>& work)
{
  threads = std::max(threads, 1);
  if (threads == 1) {
    TGeoNavigator* navigator = gGeoManager->GetCurrentNavigator();
    for (std::size_t i = 0; i < count; ++i) {
      work(i, navigator);
    }
    return;
  }
  gGeoManager->SetMaxThreads(threads);
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> workers;
  workers.reserve(static_cast<std::size_t>(threads));
  for (int thread = 0; thread < threads; ++thread) {
    workers.emplace_back([&] {
      TGeoNavigator* navigator = gGeoManager->AddNavigator();
      for (std::size_t i = next++; i < count; i = next++) {
        work(i, navigator);
      }
      gGeoManager->RemoveNavigator(navigator);
    });
  }
  for (auto& worker : workers) {
    worker.join();
  }
}

SideWalls trkSideWalls()
{
  if (gGeoManager == nullptr || gGeoManager->GetTopVolume() == nullptr) {
    throw std::runtime_error("trkSideWalls: no geometry loaded");
  }
  SideWalls walls;
  TGeoIterator next(gGeoManager->GetTopVolume());
  while (const TGeoNode* node = next()) {
    if (std::string_view{node->GetVolume()->GetName()} != "VD_SideWall") {
      continue;
    }
    const auto* box = dynamic_cast<const TGeoBBox*>(node->GetVolume()->GetShape());
    if (box == nullptr) {
      throw std::runtime_error("trkSideWalls: a VD_SideWall is not a box");
    }
    const TGeoMatrix* matrix = next.GetCurrentMatrix();
    const auto radius = [&](double x) {
      const std::array<double, 3> local{x, 0., 0.};
      std::array<double, 3> global{};
      matrix->LocalToMaster(local.data(), global.data());
      return std::hypot(global[0], global[1]);
    };
    const double inner = std::min(radius(-box->GetDX()), radius(box->GetDX()));
    const double outer = std::max(radius(-box->GetDX()), radius(box->GetDX()));
    const double thickness = 2. * box->GetDY();
    walls.thickness = walls.thickness > 0. ? std::min(walls.thickness, thickness) : thickness;
    walls.rMin = walls.rMax > 0. ? std::min(walls.rMin, inner) : inner;
    walls.rMax = std::max(walls.rMax, outer);
  }
  return walls;
}

std::vector<double> beamLineCylinderFaces(double rMax, double zMax, double dz, int threads)
{
  if (gGeoManager == nullptr || gGeoManager->GetTopVolume() == nullptr) {
    throw std::runtime_error("beamLineCylinderFaces: no geometry loaded");
  }
  constexpr int kAzimuths = 256;
  constexpr double kResolution = 1e-4; // cm
  constexpr double kGasDensity = 0.01; // g/cm^3
  const double golden = std::numbers::pi * (3. - std::sqrt(5.));
  const auto nz = static_cast<int>(std::floor(2. * zMax / dz)) + 1;
  const auto facesAt = [&](double z, TGeoNavigator* navigator) {
    std::map<long, int> seen; // per radius, the azimuths meeting it
    for (int ia = 0; ia < kAzimuths; ++ia) {
      const double phi = ia * golden;
      const std::array<double, 3> direction{std::cos(phi), std::sin(phi), 0.};
      std::array<double, 3> start{0., 0., z};
      navigator->InitTrack(start.data(), direction.data());
      std::set<long> mine;
      double r = 0.;
      bool solid = false;
      while (r < rMax && !navigator->IsOutside()) {
        const bool inSolid = navigator->GetCurrentVolume()->GetMaterial()->GetDensity() >= kGasDensity;
        if (inSolid != solid) {
          mine.insert(std::lround(r / kResolution));
          solid = inSolid;
        }
        navigator->FindNextBoundaryAndStep(rMax - r);
        const double step = navigator->GetStep();
        if (step <= 0.) {
          break;
        }
        r += step;
      }
      for (const long face : mine) {
        ++seen[face];
      }
    }
    std::vector<long> found;
    for (const auto& [face, count] : seen) {
      if (2 * count > kAzimuths && face > 0) {
        found.push_back(face);
      }
    }
    return found;
  };
  std::vector<std::vector<long>> byZ(static_cast<std::size_t>(nz));
  forEachOnNavigators(byZ.size(), threads, [&](std::size_t iz, TGeoNavigator* navigator) {
    byZ[iz] = facesAt(-zMax + (static_cast<double>(iz) * dz), navigator);
  });
  std::set<long> faces;
  for (const auto& found : byZ) {
    faces.insert(found.begin(), found.end());
  }
  std::vector<double> result;
  result.reserve(faces.size());
  for (const long face : faces) {
    result.push_back(static_cast<double>(face) * kResolution);
  }
  return result;
}

std::vector<ChipStack> trkStackThickness(std::span<const Sensor> sensors, std::span<const SensorSpec> specs, int nBarrelLayers)
{
  if (gGeoManager == nullptr) {
    throw std::runtime_error("trkStackThickness: no geometry loaded");
  }
  const std::string sensorPattern = o2::trk::GeometryTGeo::getTRKSensorPattern();
  const std::string stackPattern = o2::trk::GeometryTGeo::getTRKMetalStackPattern();
  const auto named = [](const std::string& pattern) { return std::string{gGeoManager->GetCurrentVolume()->GetName()}.find(pattern) != std::string::npos; };
  std::vector<ChipStack> thickness(static_cast<std::size_t>(nBarrelLayers));
  for (const auto& sensor : sensors) {
    if (std::cmp_greater_equal(sensor.layer, nBarrelLayers)) {
      continue;
    }
    // The centre of the silicon, and the normal away from the beam: a flat
    // sensor's local z, a curved one's radius from its axis.
    const auto& spec = specs[sensor.layer];
    double local[3];
    spec.unrolledToFrame(0., 0., 0.5 * (spec.surfaceMin() + spec.surfaceMax()), local);
    double centre[3];
    sensor.toGlobal(local, centre);
    double normal[3];
    if (spec.curved()) {
      const double dx = centre[0] - sensor.trans[0];
      const double dy = centre[1] - sensor.trans[1];
      const double norm = std::hypot(dx, dy);
      normal[0] = dx / norm;
      normal[1] = dy / norm;
      normal[2] = 0.;
    } else {
      normal[0] = sensor.rot[2];
      normal[1] = sensor.rot[5];
      normal[2] = sensor.rot[8];
      if ((normal[0] * centre[0]) + (normal[1] * centre[1]) < 0.) {
        for (double& n : normal) {
          n = -n;
        }
      }
    }
    // Out of the sensor along the normal, then through the stack's volumes.
    const auto stackAlong = [&](double sign) {
      const double direction[3] = {sign * normal[0], sign * normal[1], sign * normal[2]};
      gGeoManager->InitTrack(centre, direction);
      if (!named(sensorPattern)) {
        throw std::runtime_error(std::format("trkStackThickness: the centre of a sensor of layer {} is in {}, not a sensor", sensor.layer,
                                             gGeoManager->GetCurrentVolume()->GetName()));
      }
      gGeoManager->FindNextBoundaryAndStep();
      double stack = 0.;
      while (named(stackPattern)) {
        gGeoManager->FindNextBoundaryAndStep();
        stack += gGeoManager->GetStep();
      }
      return stack;
    };
    auto& layer = thickness[sensor.layer];
    layer.inner = std::max(layer.inner, stackAlong(-1.));
    layer.outer = std::max(layer.outer, stackAlong(1.));
  }
  return thickness;
}

std::vector<LayerEnvelope> withStacks(std::span<const LayerEnvelope> envelopes, std::span<const Sensor> sensors, std::span<const SensorSpec> specs,
                                      std::span<const ChipStack> stack)
{
  // A flat sensor's silicon is centred on its frame: a face moves by its
  // stack with the half-thickness grown by all of it, and each face is
  // measured on its own. A curved one's starts at its radius, which the inner
  // stack moves in.
  const auto thicker = [&](bool outer) {
    std::vector<SensorSpec> result(specs.begin(), specs.end());
    for (std::size_t layer = 0; layer < stack.size(); ++layer) {
      auto& spec = result[layer];
      if (spec.curved()) {
        spec.radius -= stack[layer].inner;
        spec.halfThickness += 0.5 * (stack[layer].inner + stack[layer].outer);
      } else {
        spec.halfThickness += outer ? stack[layer].outer : stack[layer].inner;
      }
    }
    return measureEnvelopes(sensors, result, static_cast<int>(stack.size()), kEnvelopeResolution);
  };
  const auto inward = thicker(false);
  const auto outward = thicker(true);
  std::vector<LayerEnvelope> result(envelopes.begin(), envelopes.end());
  for (std::size_t layer = 0; layer < stack.size(); ++layer) {
    result[layer].rMin = inward[layer].rMin;
    result[layer].rMax = outward[layer].rMax;
  }
  return result;
}

} // namespace o2::gtrk::tools
