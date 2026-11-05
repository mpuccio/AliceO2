# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

# The toy layouts, the oracles and the fixtures live beside the tests.
function(gtrk_tracking_test name source)
  o2_add_test(${name} SOURCES ${source} COMPONENT_NAME generic-tracking
              PUBLIC_LINK_LIBRARIES O2::GTrkTracking ${ARGN}
              LABELS generic-tracking TARGETVARNAME target)
  if(target)
    # the study's header-only helpers (digests, the point-only stage functions):
    # these tests build without the upgrades too
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../../study/include)
  endif()
endfunction()

gtrk_tracking_test(math-utils testMathUtils.cxx)
gtrk_tracking_test(placed-layouts testPlacedLayouts.cxx)
gtrk_tracking_test(tracking-surface testTrackingSurface.cxx)
gtrk_tracking_test(timeframe testTimeFrame.cxx)
gtrk_tracking_test(timeframe-scratch testTimeFrameScratch.cxx)
gtrk_tracking_test(stage-output testStageOutput.cxx)
gtrk_tracking_test(device-arena testDeviceArena.cxx)
gtrk_tracking_test(profile testProfile.cxx)
gtrk_tracking_test(tracklet-stage testTrackletStage.cxx)
gtrk_tracking_test(tracklet-sort testTrackletSort.cxx)
gtrk_tracking_test(triplet-fit testTripletFit.cxx)
gtrk_tracking_test(seed-cuts testSeedCuts.cxx)
gtrk_tracking_test(triplet-stage testTripletStage.cxx)
gtrk_tracking_test(tracklet-paths testTrackletPaths.cxx)
gtrk_tracking_test(triplet-end-index testTripletEndIndex.cxx)
gtrk_tracking_test(host-algorithms testHostAlgorithms.cxx)
gtrk_tracking_test(neighbour-stage testNeighbourStage.cxx)
gtrk_tracking_test(roads testRoads.cxx)
gtrk_tracking_test(road-fit testRoadFit.cxx)
gtrk_tracking_test(selection testSelection.cxx)
gtrk_tracking_test(iterations testIterations.cxx)
gtrk_tracking_test(tracker testTracker.cxx)
gtrk_tracking_test(vertex-stage testVertexStage.cxx)
gtrk_tracking_test(vertexing-pass testVertexingPass.cxx)
gtrk_tracking_test(measurement-grid testMeasurementGrid.cxx)
gtrk_tracking_test(layer-links testLayerLinks.cxx)
gtrk_tracking_test(iteration testIteration.cxx)
gtrk_tracking_test(measurement-update testMeasurementUpdate.cxx)
gtrk_tracking_test(propagator testPropagator.cxx)
gtrk_tracking_test(propagation-precision testPropagationPrecision.cxx)
gtrk_tracking_test(covariance-transport testCovarianceTransport.cxx)
gtrk_tracking_test(track-fitter testTrackFitter.cxx O2::ReconstructionDataFormats)
gtrk_tracking_test(passive-scattering testPassiveScattering.cxx)
gtrk_tracking_test(displaced-pair testDisplacedPair.cxx)
gtrk_tracking_test(follow testFollow.cxx)

add_subdirectory(device)
