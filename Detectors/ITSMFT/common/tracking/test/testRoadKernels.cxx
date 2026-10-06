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

#define BOOST_TEST_MODULE CommonTrackingRoadKernels
#include <boost/test/unit_test.hpp>
#include "ITSMFTTrackingGPU/TimeFrameGPU.h"

using namespace o2::itsmft::tracking;

BOOST_AUTO_TEST_CASE(RoadBranchingPropagationAndRecovery)
{
  const auto device = gpu::currentDevice();
  gpu::TimeFrameGPU frame{device};
  for (auto kind : {SurfaceKind::Cylinder, SurfaceKind::Disk}) {
    for (auto targetKind : {SurfaceKind::Cylinder, SurfaceKind::Disk}) {
      for (float bz : {0.f, 5.f, -5.f}) {
        SurfaceTrackState state{};
        state.kind = kind;
        state.absCharge = 1;
        state.referenceCoordinate = kind == SurfaceKind::Cylinder ? 4.f : 8.f;
        state.parameters[0] = kind == SurfaceKind::Cylinder ? 0.1f : 4.f;
        state.parameters[1] = kind == SurfaceKind::Cylinder ? 8.f : 0.1f;
        state.parameters[2] = 0.01f;
        // Real batch_015 material regression: float log rounding differed on
        // CUDA and changed q/pT by one bit, surviving the final CPU refit.
        state.parameters[3] = 0x1.10cf4ep+0f;
        state.parameters[4] = -0x1.f8d8a4p+0f;
        for (int i = 0; i < 5; ++i) {
          state.covariance[packedCovarianceIndex(i, i)] = 0.01f;
        }
        Triplet cell{1, 0, 0, 0, 1, 2, {100, 20}};
        cell.setLevel(2);
        TrackSeed seed{cell, state, 0.f};
        RoadTarget target{};
        target.cell = Triplet{0, 0, 0, 0, 0, 1, {110, 20}};
        target.cellId = 0;
        target.cellPathId = 0;
        target.available = true;
        target.surface.kind = targetKind;
        target.surface.material = {0x1.47ae14p-7f, 0x1.bea4eap-3f};
        auto prediction = state;
        BOOST_REQUIRE(Propagator::convertKind(prediction, targetKind, bz));
        const float coordinate = prediction.referenceCoordinate - 0.1f;
        BOOST_REQUIRE(Propagator::propagateToReference(prediction, coordinate, bz));
        target.measurement = {{coordinate, prediction.parameters[0], prediction.parameters[1], prediction.alpha}, {1.e-4f, 0.f, 1.e-4f}};
        RoadSeedEmission initialReference;
        BOOST_REQUIRE(extendRoad(seed, target, bz, 100.f, initialReference));
        bounded_vector<RoadJob> jobs(257, RoadJob{seed, 0, 3});
        // Sparse sources: results must report job.source, not the job position.
        for (size_t i = 0; i < jobs.size(); ++i) {
          jobs[i].source = 3 * i;
        }
        bounded_vector<RoadTarget> targets(3, target);
        targets[1].available = false;
        targets[2].cellId = 2;
        const auto launch = [&](const gpu::RoadInput& input, auto& output) { gpu::computeRoads(frame, input, output); };
        gpu::RoadInput input{jobs, targets, bz, 100.f};
        bounded_vector<RoadResult> output;
        launch(input, output);
        BOOST_REQUIRE_EQUAL(output.size(), 514u);
        for (size_t i = 0; i < output.size(); ++i) {
          RoadSeedEmission reference;
          BOOST_REQUIRE(extendRoad(seed, targets[i % 2 == 0 ? 0 : 2], bz, input.maxChi2, reference));
          BOOST_CHECK_EQUAL(output[i].job, 3 * (i / 2));
          BOOST_CHECK(roadSeedsEquivalent(reference, output[i].emission));
          BOOST_CHECK_EQUAL(output[i].emission.seed.getHitLayerMask().count(), 4);
          BOOST_CHECK_EQUAL(output[i].emission.seed.getLevel(), 1);
        }
        // Only every other predecessor has a job: the device must read each
        // continuing state from the predecessor buffer at job.source.
        bounded_vector<RoadJob> nextJobs;
        for (size_t i = 0; i < output.size(); i += 2) {
          RoadJob next{output[i].emission.seed, 0, 1};
          next.source = i;
          nextJobs.push_back(next);
        }
        auto nextTarget = target;
        nextTarget.cell.setLevel(0);
        nextTarget.cell.setSecondTrackletIndex(0);
        bounded_vector<RoadTarget> nextTargets{nextTarget};
        gpu::RoadInput nextInput{nextJobs, nextTargets, bz, 100.f, true, output.size()};
        const auto previousUpload = frame.statistics().uploadedBytes;
        bounded_vector<RoadResult> nextOutput;
        launch(nextInput, nextOutput);
        BOOST_REQUIRE_EQUAL(nextOutput.size(), nextJobs.size());
        for (size_t i = 0; i < nextOutput.size(); ++i) {
          RoadSeedEmission reference;
          BOOST_REQUIRE(extendRoad(nextJobs[i].seed, nextTarget, bz, 100.f, reference));
          BOOST_CHECK_EQUAL(nextOutput[i].job, nextJobs[i].source);
          BOOST_CHECK(roadSeedsEquivalent(reference, nextOutput[i].emission));
        }
        // Continuing states were not uploaded: descriptors plus one target
        // are smaller than the seeds alone.
        BOOST_CHECK_LT(frame.statistics().uploadedBytes - previousUpload, nextJobs.size() * sizeof(TrackSeed));
        nextInput.nSources = nextOutput.size() + 1;
        BOOST_CHECK_THROW(launch(nextInput, nextOutput), std::invalid_argument);
        targets[0].cell.getTimeStamp() = {120, 20};
        targets[2].cell.setLevel(2);
        launch(input, output);
        BOOST_CHECK(output.empty());
        targets[0] = target;
        targets[2] = target;
        targets[2].cell.setSecondTrackletIndex(99);
        launch(input, output);
        BOOST_CHECK_EQUAL(output.size(), 257u);
        jobs[0].endTarget = 4;
        BOOST_CHECK_THROW(launch(input, output), std::out_of_range);
        jobs[0].endTarget = 3;
        BoundedMemoryResource tiny{1};
        bounded_vector<RoadResult> limited{&tiny};
        BOOST_CHECK_THROW(launch(input, limited), BoundedMemoryResource::MemoryLimitExceeded);
        BOOST_CHECK_NO_THROW(launch(input, output));
        input.maxChi2 = -1.f;
        launch(input, output);
        BOOST_CHECK(output.empty());
        input.jobs = {};
        launch(input, output);
        BOOST_CHECK(output.empty());
      }
    }
  }
}
