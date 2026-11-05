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

# The tracking: the track model, the stages, the tracker and the vertexer,
# header code compiled into whichever target includes it. The compiled part is
# host only; TBB stays private behind Parallel.h.
o2_add_library(GTrkTracking
               SOURCES src/Parallel.cxx
                       src/Profile.cxx
                       src/WorkCounters.cxx
                       src/GPU/DeviceArena.cxx
                       src/GPU/DeviceContext.cxx
                       src/GPU/GPUChainGenericTracking.cxx
                       src/GPU/RuntimeFactory.cxx
                       src/GPU/DeviceStream.cxx
                       src/GPU/MaterialPort.cxx
               PUBLIC_LINK_LIBRARIES O2::GTrkBase O2::GPUCommon O2::CommonConstants O2::MathUtils O2::DetectorsBase
               PRIVATE_LINK_LIBRARIES TBB::tbb O2::GPUTracking
               TARGETVARNAME targetName)

# The work counters (Tracking/WorkCounters.h) count only where this is on: a
# tree of its own for the studies, since a count touches the hot loops. Public,
# so that every unit that instantiates the tracking counts alike.
option(GTRK_WORK_COUNTERS "Count the work inside the generic tracking's stages" OFF)
if(GTRK_WORK_COUNTERS)
  target_compile_definitions(${targetName} PUBLIC GTRK_WORK_COUNTERS)
endif()

# o2_add_library installs the headers of include/<library name> only.
target_include_directories(${targetName} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(DIRECTORY include/GenericTracking DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_subdirectory(GPU)
add_subdirectory(test)
