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

o2_add_test(mclabelaccumulator
            SOURCES testMCLabelAccumulator.cxx
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

o2_add_test(detectorlayoutbuilder
            SOURCES testDetectorLayoutBuilder.cxx
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

o2_add_test(indextableconfigbinding
            SOURCES testIndexTableConfigBinding.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking O2::DetectorsCommonDataFormats)

o2_add_test(transitionpolicyoperations
            SOURCES testTransitionPolicyOperations.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(detectortraits-track-representation
            SOURCES testDetectorTraitsTrackRepresentation.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(surfacemeasurement
            SOURCES testSurfaceMeasurement.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(surfacespec
            SOURCES testSurfaceSpec.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(surfacekinematicstate
            SOURCES testSurfaceKinematicState.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(forward-surface-state-operations
            SOURCES testForwardSurfaceStateOperations.cxx
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

o2_add_test(timeframe-load-failure
            SOURCES testTimeFrameLoadFailure.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(rof-timing-uniformity
            SOURCES testROFTimingUniformity.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::ITStracking)

o2_add_test(tracking-interface-load-failure-contract
            SOURCES testTrackingInterfaceLoadFailureContract.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking O2::DetectorsBase O2::Field)

o2_add_test(timeframe-detector-layouts
            SOURCES testTimeFrameDetectorLayouts.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(timeframe-lifecycle
            SOURCES testTimeFrameLifecycle.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(catracker-failure-contract
            SOURCES testCATrackerFailureContract.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking O2::Field)

o2_add_test(indextableactivation
            SOURCES testIndexTableActivation.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(computelayercells-orchestration
            SOURCES testComputeLayerCellsOrchestration.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(computelayertracklets-orchestration
            SOURCES testComputeLayerTrackletsOrchestration.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(detector-surface-catalog-aggregation
            SOURCES testDetectorSurfaceCatalogAggregation.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(geometry-surface-catalog-providers
            SOURCES testGeometrySurfaceCatalogProviders.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_executable(geometry-surface-catalog-validation
                  SOURCES validateGeometrySurfaceCatalogs.cxx
                  COMPONENT_NAME itsmft-tracking
                  NO_INSTALL
                  PUBLIC_LINK_LIBRARIES
                    O2::DetectorsBase
                    O2::ITSBase
                    O2::ITSMFTTracking
                    O2::MFTBase)

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

o2_add_test_root_macro(gate0-baseline/extract_metrics_common_ca.C
                       COMPILE_ONLY
                       LABELS "itsmft;tracking;gate0"
                       PUBLIC_LINK_LIBRARIES
                         O2::DataFormatsITSMFT
                         O2::DataFormatsMFT
                         O2::SimulationDataFormat)

o2_add_test_root_macro(validateGeometrySurfaceCatalogs.C
                       COMPILE_ONLY
                       LABELS "itsmft;tracking;geometry"
                       PUBLIC_LINK_LIBRARIES
                         O2::DetectorsBase
                         O2::ITSBase
                         O2::ITSMFTTracking
                         O2::MFTBase)
