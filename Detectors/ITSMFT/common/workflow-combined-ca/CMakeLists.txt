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

# Gate 4 C4: isolated opt-in combined ITS+MFT common-CA DPL workflow library.
# Deliberately its own library/executable, not an addition to ITSCAWorkflow
# or MFTWorkflow, and does not link either -- same isolation rationale as
# ITSCAWorkflow's own CMakeLists.txt comment (D009, AgentCoordination.md):
# no link-graph overlap with either single-detector opt-in workflow or
# either frozen legacy reco workflow. Its own ITS/MFT CommonTrack writer
# specs (src/ITSCommonTrackWriterSpec.cxx/src/MFTCommonTrackWriterSpec.cxx)
# are therefore small local reimplementations of ITSCAWorkflow/
# MFTWorkflow's own TrackWriterSpec.cxx -- same DPLUtils::MakeRootTreeWriterSpec
# utility, same default output filenames/tree/branch contract, so downstream
# consumers of o2trac_its_ca.root/mfttracks.root are unaffected -- not a
# link dependency on either forbidden library.
o2_add_library(ITSMFTCombinedCAWorkflow
               TARGETVARNAME targetName
               SOURCES src/ConfigPreflight.cxx
                       src/CombinedCATrackerSpec.cxx
                       src/ITSCommonTrackWriterSpec.cxx
                       src/MFTCommonTrackWriterSpec.cxx
               PUBLIC_LINK_LIBRARIES O2::Framework
                                     O2::SimulationDataFormat
                                     O2::DataFormatsITS
                                     O2::DataFormatsITSMFT
                                     O2::DataFormatsMFT
                                     O2::ITSBase
                                     O2::MFTBase
                                     O2::ITSMFTTracking
                                     O2::MFTTracking
                                     O2::CCDB)

o2_add_executable(combined-ca-tracker-workflow
                  SOURCES src/itsmft-combined-ca-tracker-workflow.cxx
                  COMPONENT_NAME itsmft
                  PUBLIC_LINK_LIBRARIES O2::ITSMFTCombinedCAWorkflow)

o2_add_test(itsmft-combined-ca-config-preflight
            COMPONENT_NAME itsmft
            LABELS "itsmft;workflow"
            SOURCES test/testITSMFTCombinedCAConfigPreflight.cxx
            PUBLIC_LINK_LIBRARIES O2::ITSMFTCombinedCAWorkflow)

o2_add_test(itsmft-combined-ca-tracker-dpl-contract
            COMPONENT_NAME itsmft
            LABELS "itsmft;workflow"
            SOURCES test/testITSMFTCombinedCATrackerDPLContract.cxx
            PUBLIC_LINK_LIBRARIES O2::ITSMFTCombinedCAWorkflow)
