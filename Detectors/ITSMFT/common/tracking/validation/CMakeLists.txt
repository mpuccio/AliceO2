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

# Gate 4 B1 Slice 1: detached nominal-geometry validator infrastructure.
#
# Deliberately outside O2::ITSMFTTracking's link graph: this directory does
# not link that library, and ../CMakeLists.txt's O2::ITSMFTTracking target
# does not link anything from here either -- the one-way, zero-coupling
# relationship the Gate 4 B1/B2 design requires. The existing
# DetectorSurfaceCatalogAggregation.{h,cxx} (../src) is compiled a second
# time, directly, into ITSMFTNominalGeometryValidatorCore below; it is not
# moved or duplicated -- exactly one implementation exists, currently owned
# by ../CMakeLists.txt's O2::ITSMFTTracking SOURCES list, until a future B2
# slice physically moves it here once static-plan activation and replay
# parity are confirmed (see the accepted B1/B2 aggregation-ownership delta).
o2_add_library(ITSMFTNominalGeometryValidatorCore
               TARGETVARNAME coreTarget
               SOURCES src/NominalGeometryReport.cxx
                       src/ValidatorArgs.cxx
                       src/GeometryTGeoSurfaceSource.cxx
                       ../src/DetectorSurfaceCatalogAggregation.cxx
               PUBLIC_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_LIST_DIR}/src
                                         ${CMAKE_CURRENT_LIST_DIR}/../src
               PUBLIC_LINK_LIBRARIES O2::DetectorsBase
                                     O2::ITSBase
                                     O2::MFTBase
                                     O2::ITSMFTBase
                                     O2::CommonUtils
                                     O2::MathUtils)

o2_add_executable(nominal-geometry-validator
                  SOURCES src/nominal-geometry-validator.cxx
                  COMPONENT_NAME itsmft
                  NO_INSTALL
                  PUBLIC_LINK_LIBRARIES O2::ITSMFTNominalGeometryValidatorCore)

add_subdirectory(test)
