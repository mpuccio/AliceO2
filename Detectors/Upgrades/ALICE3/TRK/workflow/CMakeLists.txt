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

set(trkGpuTargets "")
if(CUDA_ENABLED)
  list(APPEND trkGpuTargets O2::ITStrackingCUDA)
elseif(HIP_ENABLED)
  list(APPEND trkGpuTargets O2::ITStrackingHIP)
endif()

o2_add_library(TRKWorkflow
               TARGETVARNAME targetName
               SOURCES src/DigitReaderSpec.cxx
                       src/DigitWriterSpec.cxx
                       src/ClustererSpec.cxx
                       src/ClusterWriterSpec.cxx
                       src/TrackerSpec.cxx
                       src/TrackWriterSpec.cxx
                       src/RecoWorkflow.cxx
               PUBLIC_LINK_LIBRARIES O2::Framework
                                     O2::GPUWorkflow
                                     O2::SimConfig
                                     O2::DataFormatsITSMFT
                                     O2::DataFormatsTRK
                                     O2::SimulationDataFormat
                                     O2::DPLUtils
                                     O2::TRKBase
                                     O2::TRKSimulation
                                     O2::TRKReconstruction
                                     ${trkGpuTargets}
                                     nlohmann_json::nlohmann_json)

if(trkGpuTargets)
  target_compile_definitions(${targetName} PRIVATE TRK_HAS_GPU_TRACKING)
endif()

if(CUDA_ENABLED)
  target_compile_definitions(${targetName} PRIVATE TRK_HAS_CUDA_TRACKING)
elseif(HIP_ENABLED)
  target_compile_definitions(${targetName} PRIVATE TRK_HAS_HIP_TRACKING)
endif()

o2_add_executable(reco-workflow
                  SOURCES src/trk-reco-workflow.cxx
                  COMPONENT_NAME alice3-trk
                  PUBLIC_LINK_LIBRARIES O2::Framework
                                        O2::TRKWorkflow
                                        O2::TRKSimulation
                                        O2::TRKReconstruction
                                        O2::ITStracking)
