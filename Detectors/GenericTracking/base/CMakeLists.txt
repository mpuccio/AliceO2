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

# The base: what includes no tracking code -- the geometry and time model
# (Base/, the detectors' layouts in Detectors/), the bounded allocator, the SoA
# columns, the host algorithms, the lanes, the portable math and the capacity
# estimator. Header code but for the three sources; host only.
o2_add_library(GTrkBase
               SOURCES src/BoundedAllocator.cxx
                       src/CapacityEstimator.cxx
                       src/LayerEnvelopes.cxx
               PUBLIC_LINK_LIBRARIES O2::GPUCommon
               TARGETVARNAME targetName)

# o2_add_library installs the headers of include/<library name> only.
target_include_directories(${targetName} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_include_directories(${targetName} SYSTEM PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/External>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/GenericTracking/External>)
install(DIRECTORY include/GenericTracking DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY External/memlayout DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/GenericTracking/External)

add_subdirectory(test)
