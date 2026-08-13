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

# Lightweight param registration so MFTCATrackerParam keys in --configKeyValues
# are known to workflows that do not link O2::ITSMFTTracking (e.g. ctf-reader).
# Workflow-specific --configKeyValues policy checks live at their workflow
# boundary; this library only owns shared tracking parameter registration.
o2_add_library(ITSMFTTrackingParams
               SOURCES src/TrackingConfigParam.cxx
               PUBLIC_LINK_LIBRARIES O2::CommonUtils
                                     O2::DetectorsCommonDataFormats)

o2_target_root_dictionary(ITSMFTTrackingParams
                          HEADERS include/ITSMFTTracking/TrackingConfigParam.h
                          LINKDEF src/ITSMFTTrackingLinkDef.h)

# Shared CA tracking infrastructure for ITS and MFT. Currently wired from the
# Standalone and combined CA workflows compose Tracker, Loader, and TimeFrame
# at their workflow edges; frozen legacy ITS tracking remains separate.
# Dependency is one-way (ITSMFTTracking -> ITStracking) to avoid a circular link.
o2_add_library(ITSMFTTracking
               TARGETVARNAME targetName
               SOURCES src/IOUtils.cxx
                       src/Propagator.cxx
                       src/Configuration.cxx
                       src/TimeFrame.cxx
                       src/SurfaceTrackingScratch.cxx
                       src/TrackerTraits.cxx
                       src/CandidateFinding.cxx
                       src/TripletFitting.cxx
                       src/BarrelSurfaceStateOperations.cxx
                       src/ForwardSurfaceStateOperations.cxx
                       src/SurfaceLinearizationReference.cxx
                       src/MaterialPhysics.cxx
                       src/FamilyMaterialOperations.cxx
                       src/IndexTableConfiguration.cxx
                       src/SurfaceGraphBuilder.cxx
                       src/Tracker.cxx
               PUBLIC_LINK_LIBRARIES
                       O2::ITSMFTTrackingParams
                       O2::ITStracking
                       O2::GPUCommon
                       O2::CommonConstants
                       O2::DetectorsCommonDataFormats
                       O2::DataFormatsITSMFT
                       O2::ITSMFTBase
                       O2::CommonUtils
                       O2::DetectorsBase
                       O2::MathUtils
                       Microsoft.GSL::GSL
                       O2::SimulationDataFormat
                       O2::ReconstructionDataFormats
                       O2::DataFormatsCalibration
                       O2::DataFormatsMFT
                       TBB::tbb
               PRIVATE_LINK_LIBRARIES
                       O2::Framework
                       O2::ITSBase
                       O2::MFTBase
                       O2::MFTTracking
                       O2::DataFormatsITS)

# Gate 4 B1 Slice 1: standalone nominal-geometry validator, deliberately not
# an O2::ITSMFTTracking dependency and not a dependency of it -- see
# validation/CMakeLists.txt.
add_subdirectory(validation)

add_subdirectory(test)
