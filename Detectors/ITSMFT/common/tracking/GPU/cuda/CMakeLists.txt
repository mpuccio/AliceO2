# Copyright 2019-2026 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

# Triplet Jacobians are sensitive to approximate transcendental/division math.
# Override the global fast-math defaults only in this library's directory.
foreach(config "" _DEBUG _RELEASE _RELWITHDEBINFO _MINSIZEREL)
  string(REPLACE "--use_fast_math" "" CMAKE_CUDA_FLAGS${config} "${CMAKE_CUDA_FLAGS${config}}")
  string(REPLACE "-use_fast_math" "" CMAKE_CUDA_FLAGS${config} "${CMAKE_CUDA_FLAGS${config}}")
  string(REPLACE "--ftz=true" "" CMAKE_CUDA_FLAGS${config} "${CMAKE_CUDA_FLAGS${config}}")
endforeach()
o2_add_library(ITSMFTTrackingCUDA
               SOURCES TrackerTraitsGPU.cxx
                       TrackingKernels.cu
               PUBLIC_INCLUDE_DIRECTORIES ../
               PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
               TARGETVARNAME targetName)
set_target_gpu_arch("CUDA" ${targetName})
# Match the host's unfused selection arithmetic at cut boundaries.
target_compile_options(${targetName} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false;--ftz=false;--prec-div=true;--prec-sqrt=true>)
