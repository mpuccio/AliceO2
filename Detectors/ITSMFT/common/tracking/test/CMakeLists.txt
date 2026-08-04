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

o2_add_test(surfacemask
            SOURCES testSurfaceMask.cxx
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

o2_add_test(transitionpolicytag-containment
            SOURCES testTransitionPolicyTagContainment.cxx
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

o2_add_test(no-legacy-fitting-dependency
            SOURCES testNoLegacyFittingDependency.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(propagator
            SOURCES testPropagator.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(native-cylindercylinder-refit-driver
            SOURCES testNativeCylinderCylinderRefitDriver.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking O2::ITStracking)

o2_add_test(its-native-vs-legacy-refit-characterization-harness
            SOURCES testITSNativeVsLegacyRefitCharacterizationHarness.cxx
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

o2_add_test(covariance-sanitization
            SOURCES testCovarianceSanitization.cxx
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

# M2b: MultiSourceTimeFrameLoader's participant-count-generic atomic loading
# transaction (loadEvent()) and its loadITSAndMFT() compatibility wrapper.
o2_add_test(multisource-timeframe-loader
            SOURCES testMultiSourceTimeFrameLoader.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
                                  O2::MFTTracking)

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

o2_add_test(gate4-b31-ownership-contract
            SOURCES testGate4B31OwnershipContract.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking O2::DetectorsBase)

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

o2_add_test(traversaloperationbinding-containment
            SOURCES testTraversalOperationBindingContainment.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

# Gate 4 B2 Slice 1: replacement combined-layout test over the static,
# compile-time-concatenated ITS+MFT catalog (StaticDetectorCatalogs.h). Gate
# 4 B2 Slice 3 deleted CombinedSurfaceCatalogBuilder, its test, the runtime
# catalog providers, and moved DetectorSurfaceCatalogAggregation into the
# detached validator target (see ../validation/test/CMakeLists.txt) -- this
# is now the sole combined-catalog coverage.
o2_add_test(combined-static-surface-catalog-topology
            SOURCES testCombinedStaticSurfaceCatalogTopology.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(detector-traversal-binding
            SOURCES testDetectorTraversalBinding.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

# Gate 4 C2 Slice 1: TrackerTraits<MFTNLayers> hot-loop orchestration bound to
# a DetectorTraversalBinding built against the combined 17-surface ITS+MFT
# catalog (global MFT SurfaceIds 7..16, ClusterSourceId{1}).
o2_add_test(detector-traversal-binding-orchestration
            SOURCES testDetectorTraversalBindingOrchestration.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

# M6b (doc/design/0002-m6-generic-workspace-migration.md Sec 9): the additive,
# detector-neutral successor to DetectorTraversalBinding. Zero production
# callers -- these are the only two things that build against
# detail/SurfacePlanBinding.h today.
o2_add_test(surface-plan-binding
            SOURCES testSurfacePlanBinding.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(surface-plan-binding-no-detector-dependency
            SOURCES testSurfacePlanBindingNoDetectorDependency.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

# M6c (doc/design/0002-m6-generic-workspace-migration.md Sec 4, 9): the
# additive, detector-neutral successor to the generic portions of
# LegacyTrackerScratch<NLayers>. Zero production callers.
o2_add_test(surface-tracking-scratch
            SOURCES testSurfaceTrackingScratch.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            # M6d: SurfaceTrackingScratch.h now pulls in MFTCATrack.h/
            # IndexTableUtils.h's own MFTTracking/Constants.h dependency
            # directly (the file-doc's own narrow MFTNLayers exception) --
            # O2::ITSMFTTracking links O2::MFTTracking only PRIVATE_LINK_LIBRARIES,
            # so consumers of this now-MFT-aware public header need it too.
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking O2::MFTTracking)

# M6c: the additive, GPU-portable, non-templated successor to
# TrackSeedTpl<NLayers>. Zero production callers.
o2_add_test(track-seed
            SOURCES testTrackSeed.cxx
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

# GenericTrackingEngineMigration.md M1: the detector-neutral TrackingEngine/
# TrackingParticipant contract (ADR 0007).
o2_add_test(tracking-engine-contract
            SOURCES testTrackingEngineContract.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

o2_add_test(tracking-engine-dependency-boundary
            SOURCES testTrackingEngineDependencyBoundary.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking)

# M2c: the one explicit ITS/MFT-named application-layer participant set
# (ITSMFTLegacyParticipantSet) the combined DPL task composes construction/
# configuration through.
o2_add_test(itsmft-legacy-participant-set
            SOURCES testITSMFTLegacyParticipantSet.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
                                  O2::MFTTracking)

# M3: the combined ITS+MFT load/track/publish composition, ported from the
# deleted host-only combined coordinator's own test suite.
o2_add_test(combined-tracking-composition
            SOURCES testCombinedTrackingComposition.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
                                  O2::MFTTracking
                                  O2::DetectorsBase
                                  O2::Field)

# M6d (doc/design/0002-m6-generic-workspace-migration.md Sec 9): focused
# coverage for the MFT common-CA participant's migration onto
# SurfaceTrackingScratch/SurfacePlanBinding, exercised through the real
# production ITSMFTLegacyParticipantSet wiring.
o2_add_test(m6d-mft-migration
            SOURCES testM6dMFTMigration.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
                                  O2::MFTTracking
                                  O2::DetectorsBase
                                  O2::Field)

# M6e1 (doc/design/0002-m6-generic-workspace-migration.md): the standalone
# MFT common-CA interface (ITSMFTTrackingInterfaceMFT, o2-mft-ca-tracker-workflow)
# migrated onto SurfaceTrackingScratch/SurfacePlanBinding.
o2_add_test(m6e1-standalone-mft-migration
            SOURCES testM6e1StandaloneMFTMigration.cxx
            COMPONENT_NAME itsmft-tracking
            LABELS "itsmft;tracking"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTTracking
                                  O2::MFTTracking
                                  O2::DetectorsBase
                                  O2::Field)

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

# Gate 4 B4.2b persisted-output comparator. Its runnable A/B protocol uses
# external fixture files; compile it here to validate the ROOT and data-format
# dependencies in normal configurations.
o2_add_test_root_macro(gate4-b42b-validation/compare_common_ca_output.C
                       COMPILE_ONLY
                       LABELS "itsmft;tracking;gate4"
                       PUBLIC_LINK_LIBRARIES
                         O2::DataFormatsITSMFT
                         O2::SimulationDataFormat)

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
