// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file gtrk-trk-study-workflow.cxx
/// \brief ALICE 3's TRK and FT3 clusters from O2's readers or upstream, with the study's sink (StudyWorkflow.h).

#include "StudyWorkflow.h"

#include "CommonUtils/ConfigurableParam.h"
#include "Framework/CallbacksPolicy.h"

using namespace o2::framework;

// DPL looks both of these up by symbol before runDataProcessing.h is included,
// so they have to keep external linkage.
// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::CallbacksPolicy>& policies) { o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies); }

// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions) { o2::gtrk::study::addStudyOptions(workflowOptions); }

#include "GenericTracking/Detectors/TRKFT3.h"

#include "TRKWorkflow/ClusterReaderSpec.h"

#include <string>

#include "Framework/runDataProcessing.h" // IWYU pragma: keep

WorkflowSpec defineDataProcessing(const ConfigContext& context)
{
  o2::conf::ConfigurableParam::updateFromString(context.options().get<std::string>("configKeyValues"));

  const bool useMC = !context.options().get<bool>("disable-mc");
  auto specs = o2::gtrk::study::studyChain<o2::gtrk::detectors::TRKFT3>(context, {o2::trk::getTRKClusterReaderSpec(useMC), o2::trk::getFT3ClusterReaderSpec(useMC)});

  const o2::raw::HBFUtilsInitializer hbfIni(context, specs);

  return specs;
}
