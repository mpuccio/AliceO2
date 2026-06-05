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

o2_add_library(ITSMFTTracking
               TARGETVARNAME targetName
               SOURCES src/IOUtils.cxx
                       src/TrackingConfigParam.cxx
                       src/Configuration.cxx
                       src/Cluster.cxx
                       src/ClusterLines.cxx
                       src/TimeFrame.cxx
               PUBLIC_LINK_LIBRARIES
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
               PRIVATE_LINK_LIBRARIES
                       O2::Framework
                       O2::ITSBase
                       O2::MFTBase
                       O2::DataFormatsITS)

o2_target_root_dictionary(ITSMFTTracking
                          HEADERS include/ITSMFTTracking/IOUtils.h
                                  include/ITSMFTTracking/TrackingConfigParam.h
                                  include/ITSMFTTracking/Configuration.h
                                  include/ITSMFTTracking/IndexTableUtils.h
                                  include/ITSMFTTracking/TimeFrame.h
                          LINKDEF src/ITSMFTTrackingLinkDef.h)
