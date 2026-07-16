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

o2_add_test(layermask
            SOURCES testLayerMask.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(trackingtopology
            SOURCES testTrackingTopology.cxx

            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(surfacelayout
            SOURCES testSurfaceLayout.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(transitionpolicystate
            SOURCES testTransitionPolicyState.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(transitionpolicydispatch
            SOURCES testTransitionPolicyDispatch.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(surfacemeasurement
            SOURCES testSurfaceMeasurement.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(surfacemeasurement-adapters
            SOURCES testSurfaceMeasurementAdapters.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(surfacetiming
            SOURCES testSurfaceTiming.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(multisourceloading
            SOURCES testMultiSourceLoading.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(timeframenormalizedsource
            SOURCES testTimeFrameNormalizedSource.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

# The Gate 0 extractor needs generated tracker/MC dictionaries and external
# fixture files at execution time. Keep its runnable replay in the documented
# baseline protocol, while compiling the macro here so clean O2 configurations
# validate its includes and types.
o2_add_test_root_macro(gate0-baseline/extract_metrics.C
                       COMPILE_ONLY
                       LABELS "itsmft;tracking;gate0"
                       PUBLIC_LINK_LIBRARIES
                         O2::DataFormatsITS
                         O2::DataFormatsITSMFT
                         O2::DataFormatsMFT
                         O2::DetectorsBase
                         O2::ITSBase
                         O2::SimulationDataFormat)
