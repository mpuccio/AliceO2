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

o2_add_dpl_workflow(table-maker
                    SOURCES tableMaker.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisDataModel O2::DetectorsBase O2::AnalysisCore O2::PWGDQCore
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(table-maker-inpbpb
                    SOURCES tableMaker_PbPb.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisDataModel O2::DetectorsBase O2::AnalysisCore O2::PWGDQCore
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(table-reader
                    SOURCES tableReader.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisDataModel O2::DetectorsBase O2::AnalysisCore O2::PWGDQCore
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(dilepton-ee
                    SOURCES dileptonEE.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisDataModel O2::DetectorsBase O2::AnalysisCore O2::PWGDQCore
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(dilepton-mumu
                    SOURCES dileptonMuMu.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisDataModel O2::DetectorsBase O2::AnalysisCore O2::PWGDQCore
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(dq-filter-pp
                    SOURCES filterPP.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisDataModel O2::DetectorsBase O2::AnalysisCore O2::PWGDQCore
                    COMPONENT_NAME Analysis)

o2_add_dpl_workflow(v0-selector
                    SOURCES v0selector.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::AnalysisDataModel O2::DetectorsBase O2::AnalysisCore O2::PWGDQCore  O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)
