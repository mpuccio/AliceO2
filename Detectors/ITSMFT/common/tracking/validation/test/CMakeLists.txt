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

o2_add_test(nominal-geometry-validator-args
            SOURCES testValidatorArgs.cxx
            COMPONENT_NAME itsmft
            LABELS "itsmft;validation"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTNominalGeometryValidatorCore)

o2_add_test(nominal-geometry-report
            SOURCES testNominalGeometryReport.cxx
            COMPONENT_NAME itsmft
            LABELS "itsmft;validation"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTNominalGeometryValidatorCore)

# RapidJSON is test-only here (genuine JSON parsing of captured stdout, not
# just substring inspection); neither ITSMFTNominalGeometryValidatorCore nor
# the o2-itsmft-nominal-geometry-validator executable link it.
o2_add_test(nominal-geometry-json-stdout-contract
            SOURCES testJsonStdoutContract.cxx
            COMPONENT_NAME itsmft
            LABELS "itsmft;validation"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTNominalGeometryValidatorCore
                                  RapidJSON::RapidJSON)

# Gate 4 acceptance-cleanup C1: bit-exact float round-trip proof
# (source float32 -> JSON token -> RapidJSON-parsed float32).
o2_add_test(nominal-geometry-lossless-float-round-trip
            SOURCES testLosslessFloatRoundTrip.cxx
            COMPONENT_NAME itsmft
            LABELS "itsmft;validation"
            PUBLIC_LINK_LIBRARIES O2::ITSMFTNominalGeometryValidatorCore
                                  RapidJSON::RapidJSON)
