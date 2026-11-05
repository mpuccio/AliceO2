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

# What makes an input the tracking reads: the layer link tables, the layouts'
# layer tables, the material LUT, and the sensor-frame and material helpers
# they need. The ALICE 3 parts name the TRK, so all of it is built with the
# upgrades, the ITS's and the MFT's tables too, until a build without the
# upgrades splits them off.
o2_add_library(GTrkTools
               SOURCES src/Material.cxx
               PUBLIC_LINK_LIBRARIES O2::GTrkBase
               PRIVATE_LINK_LIBRARIES O2::TRKBase O2::DetectorsBase ROOT::Geom
               TARGETVARNAME targetName)

# o2_add_library installs the headers of include/<library name> only.
target_include_directories(${targetName} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(DIRECTORY include/GenericTracking DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

o2_add_executable(layer-links
                  SOURCES src/layerLinks.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkTools O2::GTrkWorkflow O2::DetectorsBase Boost::program_options)

o2_add_executable(import-geometry
                  SOURCES src/importGeometry.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkTools O2::GTrkWorkflow O2::TRKBase O2::FT3Base O2::MathUtils Boost::program_options)

o2_add_executable(build-matlut
                  SOURCES src/buildMatLUT.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkTools O2::GTrkWorkflow O2::DetectorsBase ROOT::Geom Boost::program_options)
