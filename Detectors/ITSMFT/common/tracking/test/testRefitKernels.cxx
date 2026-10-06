// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE ITSMFT RefitKernels
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "ITSMFTTrackingGPU/TimeFrameGPU.h"

using namespace o2::itsmft::tracking;

BOOST_AUTO_TEST_CASE(RefitParityAndFailureRecovery)
{
  const int device = gpu::currentDevice();
  gpu::TimeFrameGPU frame{device};
  for (int layout : {0, 1, 2}) {
    const auto kind = layout == 1 ? SurfaceKind::Disk : SurfaceKind::Cylinder;
    for (float bz : {0.f, 5.f, -5.f}) {
      for (bool repeat : {false, true}) {
        bounded_vector<SurfaceDescriptor> surfaces(7);
        RefitJob job{};
        job.nSlots = 7;
        SurfaceTrackState state{};
        state.kind = kind;
        state.absCharge = 1;
        state.referenceCoordinate = kind == SurfaceKind::Cylinder ? 2.f : 10.f;
        const auto parametersArray = kind == SurfaceKind::Cylinder ? std::array<float, 5>{0.1f, 1.f, 0.03f, 0.5f, 1.f} : std::array<float, 5>{2.f, 0.1f, 0.03f, 2.f, 1.f};
        for (int i = 0; i < 5; ++i) {
          state.parameters[i] = parametersArray[i];
          state.covariance[packedCovarianceIndex(i, i)] = 0.01f;
        }
        LayerMask mask{};
        for (int layer = 0; layer < 7; ++layer) {
          const auto measurementKind = layout == 2 && layer % 2 ? SurfaceKind::Disk : kind;
          BOOST_REQUIRE(Propagator::convertKind(state, measurementKind, bz));
          surfaces[layer].kind = measurementKind;
          // Small but nonzero material exercises each leg's direction.
          surfaces[layer].material = {1.e-4f, 1.e-3f};
          if (layer) {
            BOOST_REQUIRE(Propagator::propagateToReference(state, state.referenceCoordinate + 1.f, bz));
          }
          if (layer == 3) { // a genuine hole between attached surfaces
            continue;
          }
          job.slots[layer] = {{{state.referenceCoordinate, state.parameters[0], state.parameters[1], state.alpha},
                                {1.e-4f, 0.f, 1.e-4f}}, LayerId{static_cast<uint16_t>(layer)}, true};
          const float x = measurementKind == SurfaceKind::Cylinder ? state.referenceCoordinate : state.parameters[0];
          const float y = measurementKind == SurfaceKind::Cylinder ? state.parameters[0] : state.parameters[1];
          job.points[job.nPoints++] = {x, y, 1.e-4f, 0.f, 1.e-4f};
          job.seed.setCluster(layer, 0);
          mask.set(layer);
        }
        job.seed.state() = state;
        job.seed.setHitLayerMask(mask);
        bounded_vector<RefitJob> jobs(129, job);
        RefitParameters parameters{bz, 100.f, 100.f, true, repeat};
        const SurfaceCatalogView catalog{surfaces.data(), static_cast<uint32_t>(surfaces.size())};
        BOOST_REQUIRE(evaluateRefit(job, catalog, parameters).accepted);
        jobs[1].minPt = 100.f;
        jobs[2].nSlots = 0;
        jobs[3].slots[0].surface = LayerId{31};
        jobs[4].slots[0].measurement.covariance.uu = -1.f;
        bounded_vector<RefitResult> output;
        gpu::refitTracks(frame, {jobs, surfaces, parameters}, output);
        BOOST_REQUIRE_EQUAL(output.size(), jobs.size());
        for (size_t i = 0; i < jobs.size(); ++i) {
          BOOST_CHECK(refitsEquivalent(evaluateRefit(jobs[i], catalog, parameters), output[i]));
        }
        BOOST_CHECK(!output[1].accepted);
        BOOST_CHECK(!output[2].accepted);
        BOOST_CHECK(!output[3].accepted);
        BOOST_CHECK(!output[4].accepted);
        BoundedMemoryResource tiny{1};
        bounded_vector<RefitResult> limited{&tiny};
        BOOST_CHECK_THROW(gpu::refitTracks(frame, {jobs, surfaces, parameters}, limited), BoundedMemoryResource::MemoryLimitExceeded);
        jobs[0].nPoints = MaxLayoutSurfaces + 1;
        BOOST_CHECK_THROW(gpu::refitTracks(frame, {jobs, surfaces, parameters}, output), std::out_of_range);
        jobs[0] = job;
        gpu::refitTracks(frame, {jobs, surfaces, parameters}, output);
        BOOST_REQUIRE(output[0].accepted);
        gpu::refitTracks(frame, {{}, surfaces, parameters}, output);
        BOOST_CHECK(output.empty());

        RoadJob road{};
        road.initialize = true;
        road.endTarget = 1;
        road.initial.cell = Triplet{4, 0, 0, 0, 1, 2, {100, 20}};
        road.initial.cell.setLevel(2);
        for (int hit = 0; hit < 3; ++hit) {
          const auto& measurement = job.slots[hit + 4].measurement;
          road.initial.measurements[hit] = measurement;
          road.initial.surfaces[hit] = surfaces[hit + 4];
          const auto measurementKind = surfaces[hit + 4].kind;
          auto& position = road.initial.globals[hit].position;
          position.x = measurementKind == SurfaceKind::Cylinder ? measurement.frame.q : measurement.frame.u;
          position.y = measurementKind == SurfaceKind::Cylinder ? measurement.frame.u : measurement.frame.v;
          position.z = measurementKind == SurfaceKind::Cylinder ? measurement.frame.v : measurement.frame.q;
        }
        RoadTarget target{};
        target.cell = Triplet{LayerMask{2, 4, 5}, 0, 0, 0, 0, 1, {100, 20}};
        target.available = true;
        target.surface = surfaces[2];
        target.measurement = job.slots[2].measurement;
        target.cellId = 2;
        target.cellPathId = 0;
        bounded_vector<RoadJob> roads(129, road);
        for (size_t i = 0; i < roads.size(); ++i) {
          roads[i].source = i;
        }
        bounded_vector<RoadTarget> targets(1, target);
        bounded_vector<RoadResult> roadOutput;
        gpu::computeRoads(frame, {roads, targets, bz, 100.f}, roadOutput);
        bounded_vector<RoadSeedEmission> reference;
        forEachRoad(road.seed, &road.initial, targets.data(), road.firstTarget, road.endTarget, bz, 100.f, [&](const auto& emission) { reference.push_back(emission); });
        BOOST_REQUIRE_EQUAL(reference.size(), 1u);
        BOOST_REQUIRE_EQUAL(roadOutput.size(), roads.size());
        for (size_t i = 0; i < roadOutput.size(); ++i) {
          BOOST_CHECK_EQUAL(roadOutput[i].job, i);
          BOOST_CHECK(roadSeedsEquivalent(reference[0], roadOutput[i].emission));
        }
      }
    }
  }
}
