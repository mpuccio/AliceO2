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

# Gate 3 workflow-onboarding Slice 2: isolated ITS common-CA opt-in workflow
# library. Deliberately its own library/executable, not an addition to
# ITSWorkflow -- the frozen legacy o2-its-reco-workflow's link graph (built
# from ITSWorkflow, which links O2::ITSTrackingInterface/O2::ITSReconstruction)
# has no overlap with this target, which links O2::ITSMFTTracking instead.

o2_add_library(ITSCAWorkflow
               TARGETVARNAME targetName
               SOURCES src/ConfigPreflight.cxx
                       src/CATrackerSpec.cxx
                       src/TrackWriterSpec.cxx
               PUBLIC_LINK_LIBRARIES O2::Framework
                                     O2::SimulationDataFormat
                                     O2::DataFormatsITS
                                     O2::DataFormatsITSMFT
                                     O2::ITSBase
                                     O2::ITSMFTTracking
                                     # ITSMFTTracking/TimeFrame.h transitively
                                     # includes MFTCATrack.h, which needs
                                     # MFTTracking/Constants.h; O2::ITSMFTTracking
                                     # links O2::MFTTracking only PRIVATE-ly
                                     # (Detectors/ITSMFT/common/tracking/CMakeLists.txt),
                                     # so an ITS-only consumer must still add it
                                     # directly to see that public header's own
                                     # include path (same as O2::MFTWorkflow does).
                                     O2::MFTTracking
                                     O2::CCDB)

o2_add_executable(ca-tracker-workflow
                  SOURCES src/its-ca-tracker-workflow.cxx
                  COMPONENT_NAME its
                  PUBLIC_LINK_LIBRARIES O2::ITSCAWorkflow)

o2_add_test(its-ca-config-preflight
            COMPONENT_NAME its
            LABELS "its;workflow;itsmft"
            SOURCES test/testITSCAConfigPreflight.cxx
            PUBLIC_LINK_LIBRARIES O2::ITSCAWorkflow)

o2_add_test(its-ca-trackits-conversion
            COMPONENT_NAME its
            LABELS "its;workflow;itsmft"
            SOURCES test/testITSCATrackITSConversion.cxx
            PUBLIC_LINK_LIBRARIES O2::ITSCAWorkflow)

o2_add_test(its-ca-tracker-dpl-contract
            COMPONENT_NAME its
            LABELS "its;workflow;itsmft"
            SOURCES test/testITSCATrackerDPLContract.cxx
            PUBLIC_LINK_LIBRARIES O2::ITSCAWorkflow)
