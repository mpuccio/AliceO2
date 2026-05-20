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

/// @file   CATrackerSpec.h

#ifndef O2_MFT_CATRACKERSPEC_H_
#define O2_MFT_CATRACKERSPEC_H_

#include <memory>
#include <utility>

#include <gsl/span>

#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITStracking/ROFLookupTables.h"

namespace o2::mft
{

class CATrackerDPL : public o2::framework::Task
{
 public:
  static constexpr int NCALayers = 5;
  using ROFOverlapTable = o2::its::ROFOverlapTable<NCALayers>;

  CATrackerDPL(std::shared_ptr<o2::base::GRPGeomRequest> gr, bool useMC) : mGGCCDBRequest(std::move(gr)), mUseMC(useMC) {}
  ~CATrackerDPL() override = default;

  void init(framework::InitContext& ic) final;
  void run(framework::ProcessingContext& pc) final;
  void finaliseCCDB(framework::ConcreteDataMatcher& matcher, void* obj) final;

 private:
  void updateTimeDependentParams(framework::ProcessingContext& pc);
  void initialiseROFTable(gsl::span<const o2::itsmft::ROFRecord> rofs);

  bool mUseMC = false;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  const o2::itsmft::TopologyDictionary* mDict = nullptr;
  ROFOverlapTable mROFTable;
  ROFOverlapTable::View mROFTableView;
};

/// create a processor spec that reads all MFT tracker inputs and provides the
/// future insertion point for wiring them into a generalized ITS TimeFrame.
o2::framework::DataProcessorSpec getCATrackerSpec(bool useMC, bool useGeom, bool useIRFrames);

} // namespace o2::mft

#endif // O2_MFT_CATRACKERSPEC_H_
