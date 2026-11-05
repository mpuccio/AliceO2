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

# The ALICE adapters (the chips as our sensors, O2's clusters as ours, the
# fields and readout frames), the cluster record, the parameters, the DPL
# devices that make the tracking's measurements, each detector's in its own
# unit, ALICE 3's (TRK.cxx) with the upgrades, and the tracking device
# (TrackerDevice.h), one unit per detector, each instantiating the tracker.
set(sources src/Clusters.cxx src/Placement.cxx src/RecoParam.cxx src/ITS.cxx src/MFT.cxx src/TrackerITS.cxx src/TrackerMFT.cxx src/AdaptorITS.cxx
            src/AdaptorMFT.cxx)
set(trackerUnits src/TrackerITS.cxx src/TrackerMFT.cxx)
set(publicLibraries O2::GTrkTracking O2::Framework O2::CommonUtils O2::DataFormatsCalibration O2::DetectorsBase
                    O2::DataFormatsITSMFT O2::ITSBase O2::MFTBase O2::DetectorsRaw ROOT::Core)
set(privateLibraries O2::ITSMFTBase O2::ITSMFTTracking O2::DataFormatsITS O2::DataFormatsMFT O2::SimulationDataFormat O2::CommonConstants
                     O2::MathUtils O2::ReconstructionDataFormats ROOT::Geom ROOT::RIO ROOT::Tree)
if(ENABLE_UPGRADES)
  list(APPEND sources src/TRK.cxx src/TrackerTRKFT3.cxx)
  list(APPEND trackerUnits src/TrackerTRKFT3.cxx)
  list(APPEND publicLibraries O2::DataFormatsTRKFT3 O2::TRKReconstruction O2::TRKBase)
  list(APPEND privateLibraries O2::FT3Base O2::Field)
endif()

o2_add_library(GTrkWorkflow
               SOURCES ${sources}
               PUBLIC_LINK_LIBRARIES ${publicLibraries}
               PRIVATE_LINK_LIBRARIES ${privateLibraries}
               TARGETVARNAME targetName)

# The tracker's units instantiate the header tracker (GTRK_TRACKER_OPTIMIZATION).
set_source_files_properties(${trackerUnits} PROPERTIES COMPILE_OPTIONS "${GTRK_TRACKER_OPTIMIZATION}")

# o2_add_library installs the headers of include/<library name> only.
target_include_directories(${targetName} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(DIRECTORY include/GenericTracking DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# GTRKReco, a ConfigurableParam, which needs a dictionary for
# --configKeyValues to reach its fields by name, and the cluster-shape band,
# a CCDB object.
o2_target_root_dictionary(GTrkWorkflow
                          HEADERS include/GenericTracking/Workflow/RecoParam.h
                                  include/GenericTracking/Workflow/ClusterShapeBand.h
                          LINKDEF src/GTrkWorkflowLinkDef.h)

# The tracking workflows, one per detector: O2's clusters in, GenericTracking's
# tracks and vertices out; the ITS's and the MFT's also as O2's formats and files.
o2_add_executable(its-reco-workflow
                  SOURCES src/gtrk-its-reco-workflow.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkWorkflow O2::DetectorsRaw O2::ITSMFTWorkflow O2::ITSWorkflow O2::GlobalTrackingWorkflowWriters)
o2_add_executable(mft-reco-workflow
                  SOURCES src/gtrk-mft-reco-workflow.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkWorkflow O2::DetectorsRaw O2::ITSMFTWorkflow O2::ITSMFTCAWriter)
if(ENABLE_UPGRADES)
  o2_add_executable(trk-reco-workflow
                    SOURCES src/gtrk-trk-reco-workflow.cxx
                    COMPONENT_NAME gtrk
                    PUBLIC_LINK_LIBRARIES O2::GTrkWorkflow O2::DetectorsRaw O2::TRKWorkflow)
endif()

add_subdirectory(test)
