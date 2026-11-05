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

# What measures or checks the tracking: the bench and its analyses, the replay
# format of the dumps, the study workflow's sink and its QA, the MC vertex
# reader, the checks of the material LUT, the clusters and the hits. Built
# with the upgrades (the bench, the hits and the replay format name the TRK).
# The library, the sink and the benches instantiate the header tracker
# (GTRK_TRACKER_OPTIMIZATION).
function(gtrk_study_o3 target)
  target_compile_options(${target} PRIVATE ${GTRK_TRACKER_OPTIMIZATION})
endfunction()

# What the study's sources and public headers include, outside GenericTracking.
set(studyLibraries O2::GTrkTools O2::GTrkWorkflow O2::GTrkTracking O2::Framework O2::CommonUtils O2::CommonConstants
                   O2::DetectorsBase O2::SimulationDataFormat O2::ReconstructionDataFormats)
o2_add_library(GTrkStudy
               SOURCES src/Hits.cxx
                       src/Tracks.cxx
                       src/MCVertexParam.cxx
                       src/MCVertexReaderSpec.cxx
                       src/ProfileReport.cxx
                       src/TrackComparison.cxx
                       src/TrackPulls.cxx
                       src/TrackingInput.cxx
                       src/TrackingQA.cxx
               PUBLIC_LINK_LIBRARIES ${studyLibraries}
               PRIVATE_LINK_LIBRARIES O2::ITSMFTSimulation O2::DataFormatsTRKFT3 O2::DataFormatsITS O2::DataFormatsMFT O2::Steer O2::DetectorsRaw
                                      ROOT::Core ROOT::RIO ROOT::Tree
               TARGETVARNAME targetName)
gtrk_study_o3(${targetName})

# o2_add_library installs the headers of include/<library name> only.
target_include_directories(${targetName} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(DIRECTORY include/GenericTracking DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# GTRKMCVertex, a ConfigurableParam: a dictionary for --configKeyValues.
o2_target_root_dictionary(GTrkStudy
                          HEADERS include/GenericTracking/Study/MCVertexParam.h
                          LINKDEF src/GTrkStudyLinkDef.h)

# The sink, one unit per detector: each instantiates the tracker and its QA.
foreach(detector ITS MFT TRKFT3)
  add_library(gtrk_sink_${detector} OBJECT src/MeasurementSinkSpec.cxx)
  set_property(TARGET gtrk_sink_${detector} PROPERTY POSITION_INDEPENDENT_CODE ON)
  target_compile_definitions(gtrk_sink_${detector} PRIVATE GTRK_SINK_DETECTOR=${detector})
  target_include_directories(gtrk_sink_${detector} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  target_link_libraries(gtrk_sink_${detector} PRIVATE ${studyLibraries} O2::Steer)
  gtrk_study_o3(gtrk_sink_${detector})
  target_sources(${targetName} PRIVATE $<TARGET_OBJECTS:gtrk_sink_${detector}>)
endforeach()

# The study workflows, one per detector (tools/StudyWorkflow.h).
o2_add_executable(its-study-workflow
                  SOURCES tools/gtrk-its-study-workflow.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkStudy O2::DetectorsRaw O2::ITSMFTWorkflow)
o2_add_executable(mft-study-workflow
                  SOURCES tools/gtrk-mft-study-workflow.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkStudy O2::DetectorsRaw O2::ITSMFTWorkflow)
o2_add_executable(trk-study-workflow
                  SOURCES tools/gtrk-trk-study-workflow.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkStudy O2::DetectorsRaw O2::TRKWorkflow)

# One unit per detector's run (TrackBench.h), the options and the dispatch in
# trackBench.cxx: the detectors compile in parallel.
o2_add_executable(track-bench
                  SOURCES tools/trackBench.cxx tools/trackBenchITS.cxx tools/trackBenchMFT.cxx tools/trackBenchTRKFT3.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkStudy O2::Steer O2::CCDB Boost::program_options
                  TARGETVARNAME benchTarget)
gtrk_study_o3(${benchTarget})

o2_add_executable(fit-bench
                  SOURCES tools/fitBench.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkStudy Boost::program_options
                  TARGETVARNAME fitBenchTarget)
gtrk_study_o3(${fitBenchTarget})

o2_add_executable(check-matlut
                  SOURCES tools/checkMatLUT.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkStudy O2::DetectorsBase O2::Steer ROOT::Geom ROOT::RIO ROOT::Tree Boost::program_options)

o2_add_executable(cluster-closure
                  SOURCES tools/clusterClosure.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkStudy O2::ITSBase O2::MFTBase O2::DataFormatsITSMFT O2::TRKBase O2::FT3Base
                                        O2::TRKReconstruction O2::DataFormatsTRKFT3 O2::DetectorsRaw O2::CommonUtils ROOT::RIO ROOT::Tree
                                        Boost::program_options)

o2_add_executable(check-hits
                  SOURCES tools/checkHits.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkStudy Boost::program_options)

add_subdirectory(test)
