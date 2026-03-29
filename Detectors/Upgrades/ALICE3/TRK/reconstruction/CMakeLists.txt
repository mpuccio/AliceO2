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

if(Acts_FOUND)
  set(actsTarget Acts::Core)
endif()

set(trkGpuSources "")
set(trkGpuTargets "")
set(trkGpuPrivateTargets "")
if(CUDA_ENABLED)
  find_package(CUDAToolkit REQUIRED)
  list(APPEND trkGpuSources src/TimeFrameGPU.cxx src/GPUExternalAllocator.cxx)
  list(APPEND trkGpuTargets O2::ITStrackingCUDA)
  list(APPEND trkGpuPrivateTargets CUDA::cudart)
elseif(HIP_ENABLED)
  list(APPEND trkGpuSources src/TimeFrameGPU.cxx src/GPUExternalAllocator.cxx)
  list(APPEND trkGpuTargets O2::ITStrackingHIP)
  list(APPEND trkGpuPrivateTargets hip::host)
endif()

o2_add_library(TRKReconstruction
               TARGETVARNAME targetName
               SOURCES src/TimeFrame.cxx
                       ${trkGpuSources}
                       src/Clusterer.cxx
                       $<$<BOOL:${Acts_FOUND}>:src/ClustererACTS.cxx>
                       $<$<BOOL:${Acts_FOUND}>:src/TrackerACTS.cxx>
               PUBLIC_LINK_LIBRARIES
                       O2::ITStracking
                       O2::GPUCommon
                       Microsoft.GSL::GSL
                       O2::CommonConstants
                       O2::DataFormatsITSMFT
                       O2::DataFormatsTRK
                       O2::SimulationDataFormat
                       O2::ITSBase
                       O2::ITSReconstruction
                       O2::ITSMFTReconstruction
                       O2::DataFormatsITS
                       O2::TRKSimulation
                       nlohmann_json::nlohmann_json
                       ${trkGpuTargets}
                       ${actsTarget}
               PRIVATE_LINK_LIBRARIES
                       O2::Steer
                       TBB::tbb
                       ${trkGpuPrivateTargets})

if(trkGpuTargets)
  target_compile_definitions(${targetName} PUBLIC TRK_HAS_GPU_TRACKING)
endif()

if(CUDA_ENABLED)
  target_include_directories(${targetName} PRIVATE ${CUDAToolkit_INCLUDE_DIRS})
endif()

if(CUDA_ENABLED)
  target_compile_definitions(${targetName} PUBLIC TRK_HAS_CUDA_TRACKING)
elseif(HIP_ENABLED)
  target_compile_definitions(${targetName} PUBLIC TRK_HAS_HIP_TRACKING)
endif()

if(Acts_FOUND)
  target_compile_definitions(${targetName} PUBLIC O2_WITH_ACTS)
endif()

set(dictHeaders include/TRKReconstruction/TimeFrame.h
                include/TRKReconstruction/Clusterer.h)

if(Acts_FOUND)
  list(APPEND dictHeaders include/TRKReconstruction/ClustererACTS.h
                          include/TRKReconstruction/TrackerACTS.h)
endif()

o2_target_root_dictionary(TRKReconstruction
                          HEADERS ${dictHeaders}
                          LINKDEF src/TRKReconstructionLinkDef.h)
