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

o2_add_test(detectorlayoutview
            SOURCES testDetectorLayoutView.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(transitionpolicystate
            SOURCES testTransitionPolicyState.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(mft-normalized-refit
            SOURCES testMFTNormalizedRefit.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(transitionpolicydispatch
            SOURCES testTransitionPolicyDispatch.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(material-correction-mode-preflight
            SOURCES testMaterialCorrectionModePreflight.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

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

o2_add_test(transitionpolicyoperations-native
            SOURCES testTransitionPolicyOperationsNative.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(detectortraits-track-representation
            SOURCES testDetectorTraitsTrackRepresentation.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(cell-representation
            SOURCES testCellRepresentation.cxx
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

# Gate 4 B1 Slice C3a: ITSSurfaceSpec/MFTSurfaceSpec candidate table tests.
o2_add_test(its-mft-surfacespec-projection
            SOURCES testITSMFTSurfaceSpecProjection.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(surfacekinematicstate
            SOURCES testSurfaceKinematicState.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(commontrack
            SOURCES testCommonTrack.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(forward-surface-state-operations
            SOURCES testForwardSurfaceStateOperations.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(barrel-surface-state-operations
            SOURCES testBarrelSurfaceStateOperations.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(surfacelinearizationreference
            SOURCES testSurfaceLinearizationReference.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(barrel-linearization-operations
            SOURCES testBarrelLinearizationOperations.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(forward-linearization-operations
            SOURCES testForwardLinearizationOperations.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(refithit
            SOURCES testRefitHit.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(drive-refit-leg
            SOURCES testDriveRefitLeg.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(refit-leg-assembly
            SOURCES testRefitLegAssembly.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

o2_add_test(native-cylindercylinder-refit-driver
            SOURCES testNativeCylinderCylinderRefitDriver.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking O2::ITStracking)

o2_add_test(material-physics
            SOURCES testMaterialPhysics.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(family-material-operations
            SOURCES testFamilyMaterialOperations.cxx
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

o2_add_test(timeframe-covariance-lifecycle
            SOURCES testTimeFrameCovarianceLifecycle.cxx
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

o2_add_test(combined-surface-catalog-builder
            SOURCES testCombinedSurfaceCatalogBuilder.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(its-common-ca-tracking-mode-configuration
            SOURCES testITSCommonCATrackingModeConfiguration.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
                                  FairLogger::FairLogger)

o2_add_test(its-common-ca-configkeyvalues-preflight
            SOURCES testITSCommonCAConfigKeyValuesPreflight.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTrackingParams)

o2_add_test(static-diamond-rof-compatibility
            SOURCES testStaticDiamondROFCompatibility.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITStracking)

o2_add_test(its-common-ca-nthreads
            SOURCES testITSCommonCANThreads.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
                                  O2::MFTTracking
                                  FairLogger::FairLogger)

o2_add_test(tracking-mode-ownership
            SOURCES testTrackingModeOwnership.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
                                  O2::MFTTracking
                                  O2::DetectorsBase
                                  O2::Field
                                  FairLogger::FairLogger)

o2_add_test(beam-position-ownership
            SOURCES testBeamPositionOwnership.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
                                  O2::MFTTracking
                                  O2::DataFormatsCalibration)

o2_add_test(ioutils-syserror-ownership
            SOURCES testIOUtilsSysErrorOwnership.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
                                  O2::MFTTracking)

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

# Gate 3 Slice 3 ITS common-CA characterization extractor (see
# gate3-slice3-its-ca-validation/README.md). Sibling to the gate0-baseline
# macros above: compiled here so a clean configure validates its includes
# and types, while its runnable replay/extraction stays in that directory's
# own documented protocol (external fixture/replay artifacts, not present
# at configure time).
o2_add_test_root_macro(gate3-slice3-its-ca-validation/extract_metrics_its_common_ca.C
                       COMPILE_ONLY
                       LABELS "itsmft;tracking;gate3"
                       PUBLIC_LINK_LIBRARIES
                         O2::DataFormatsITS
                         O2::DataFormatsITSMFT
                         O2::DetectorsBase
                         O2::ITSBase
                         O2::SimulationDataFormat)

# Gate 3 Slice 3b legacy track-content-hash extractor (see
# gate3-slice3b-its-ca-thread-validation/README.md). Same treatment as the
# gate3-slice3 macro above: compiled here so a clean configure validates its
# includes and types, while its runnable replay/extraction stays in that
# directory's own documented protocol (external replay artifacts, not
# present at configure time).
o2_add_test_root_macro(gate3-slice3b-its-ca-thread-validation/extract_track_content_hash_its_legacy.C
                       COMPILE_ONLY
                       LABELS "itsmft;tracking;gate3"
                       PUBLIC_LINK_LIBRARIES
                         O2::DataFormatsITS)

# Publication/CMake correction handoff: cluster-attachment identity
# checkers (see cluster-attachment-validation/ -- proving every emitted
# flattened cluster reference has the correct hit layer, external index,
# published cluster size, and associated ROF; siblings to the
# gate3-slice3-its-ca-validation and gate0-baseline extractors above, not a
# replacement for either). Compiled here so a clean configure validates
# their includes and types; their runnable replay/check protocol stays
# external (fixture/replay artifacts not present at configure time).
o2_add_test_root_macro(cluster-attachment-validation/check_cluster_attachment_its_common_ca.C
                       COMPILE_ONLY
                       LABELS "itsmft;tracking;gate3"
                       PUBLIC_LINK_LIBRARIES
                         O2::CCDB
                         O2::DataFormatsITS
                         O2::DataFormatsITSMFT
                         O2::DetectorsBase
                         O2::ITSBase)

o2_add_test_root_macro(cluster-attachment-validation/check_cluster_attachment_mft_common_ca.C
                       COMPILE_ONLY
                       LABELS "itsmft;tracking;gate3"
                       PUBLIC_LINK_LIBRARIES
                         O2::CCDB
                         O2::DataFormatsITSMFT
                         O2::DataFormatsMFT
                         O2::DetectorsBase
                         O2::MFTBase)
