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

o2_add_library(GTrkTOF
               SOURCES src/Matcher.cxx src/Clusters.cxx src/Workflow.cxx
               PUBLIC_LINK_LIBRARIES O2::GTrkWorkflow O2::DataFormatsIOTOF
               PRIVATE_LINK_LIBRARIES O2::IOTOFBase O2::DPLUtils ROOT::RIO ROOT::Tree
               TARGETVARNAME targetName)
target_include_directories(${targetName} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(DIRECTORY include/GenericTracking DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
o2_target_root_dictionary(GTrkTOF HEADERS include/GenericTracking/TOF/Match.h LINKDEF src/LinkDef.h)

o2_add_test(tof-matching SOURCES test/testTOFMatching.cxx COMPONENT_NAME generic-tracking
            PUBLIC_LINK_LIBRARIES O2::GTrkTOF O2::GTrkTracking ROOT::RIO
            LABELS generic-tracking TARGETVARNAME testTarget)
if(testTarget)
  target_include_directories(${testTarget} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../tracking/test)
endif()

# Frozen tracking inputs and transported TOF hits, for an initial diagnostic.
# GTrkStudy is declared later by GenericTracking's parent CMakeLists.
o2_add_executable(tof-hit-study
                  SOURCES study/tofHitStudy.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkTOF O2::GTrkStudy O2::ITSMFTSimulation
                                        Boost::program_options ROOT::RIO ROOT::Tree)

o2_add_executable(prepare-tof-clusters
                  SOURCES study/prepareTOFClusters.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::IOTOFReconstruction O2::IOTOFSimulation O2::DetectorsBase
                                        Boost::program_options ROOT::RIO ROOT::Tree nlohmann_json::nlohmann_json)
