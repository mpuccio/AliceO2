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

o2_add_test(boundedmemoryresource
            SOURCES testBoundedMemoryResource.cxx
            COMPONENT_NAME its-tracking
            LABELS "its;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITStracking)

o2_add_test(roflookuptables
            SOURCES testROFLookupTables.cxx
            COMPONENT_NAME its-tracking
            LABELS "its;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITStracking)

o2_add_test(trackingtopology
            SOURCES testTrackingTopology.cxx
            COMPONENT_NAME its-tracking
            LABELS "its;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITStracking)

o2_add_test(commontrackingparameters
            SOURCES testCommonTrackingParameters.cxx
            COMPONENT_NAME its-tracking
            LABELS "its;tracking;itsmft;gate3"
            PUBLIC_LINK_LIBRARIES O2::ITSCommonTrackingAdapter)

o2_add_executable(invalid-common-tracking-mode
                  SOURCES invalidCommonTrackingMode.cxx
                  COMPONENT_NAME its-tracking
                  NO_INSTALL
                  PUBLIC_LINK_LIBRARIES O2::ITSCommonTrackingAdapter)
add_dependencies(O2test-its-tracking-commontrackingparameters
                 O2exe-its-tracking-invalid-common-tracking-mode)

# Link-only Gate 3 onboarding proof. It initializes the translated ITS
# parameters, seven-cylinder layout/topology, Tracker<7>, and the
# CylinderCylinder policy boundary without processing or publishing an event.
o2_add_executable(its-common-tracker-onboarding-compile
                  SOURCES compileCommonTrackerOnboarding.cxx
                  COMPONENT_NAME its-tracking
                  NO_INSTALL
                  PUBLIC_LINK_LIBRARIES O2::ITSCommonTrackingAdapter
                                        O2::MFTTracking)
