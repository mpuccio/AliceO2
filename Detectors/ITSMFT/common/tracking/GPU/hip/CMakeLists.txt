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

# Keep one implementation, as in the legacy ITS tracker.
o2_add_hipified_library(ITSMFTTrackingHIP
                       SOURCES ../cuda/TrackerTraitsGPU.cxx
                               ../cuda/TrackingKernels.cu
                       PUBLIC_INCLUDE_DIRECTORIES ../
                       PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking hip::host
                       TARGETVARNAME targetName)
set_target_gpu_arch("HIP" ${targetName})
target_compile_options(${targetName} PRIVATE -fno-fast-math -ffp-contract=off)
# Selects the backend namespace in the (not hipified) headers, for users too.
target_compile_definitions(${targetName} PUBLIC ITSMFT_TRACKING_HIP)
