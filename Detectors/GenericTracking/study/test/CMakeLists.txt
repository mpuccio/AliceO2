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

# The study's tests; the toy layouts and the fixtures are the tracking's.
function(gtrk_study_test name source)
  o2_add_test(${name} SOURCES ${source} COMPONENT_NAME generic-tracking
              PUBLIC_LINK_LIBRARIES O2::GTrkStudy ${ARGN}
              LABELS generic-tracking TARGETVARNAME target)
  if(target)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../tracking/test)
  endif()
endfunction()

gtrk_study_test(mc-truth testMCTruth.cxx)
gtrk_study_test(track-comparison testTrackComparison.cxx)
gtrk_study_test(stage-funnel testStageFunnel.cxx)
gtrk_study_test(tracking-qa testTrackingQA.cxx)
gtrk_study_test(tracking-truth testTrackingTruth.cxx)
gtrk_study_test(vertex-qa testVertexQA.cxx)
gtrk_study_test(study-profile testProfile.cxx)
gtrk_study_test(stage-comparison testStageComparison.cxx)
