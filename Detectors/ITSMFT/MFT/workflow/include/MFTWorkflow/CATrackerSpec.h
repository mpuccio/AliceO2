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

#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITSMFTTracking/Configuration.h"
#include "ITSMFTTracking/TrackingInterface.h"

namespace o2::mft
{

/// MFT CA tracker DPL task. Delegates reconstruction to ITSMFTTrackingInterfaceMFT.
class CATrackerDPL : public o2::framework::Task
{
 public:
  CATrackerDPL(std::shared_ptr<o2::base::GRPGeomRequest> gr,
               bool useMC,
               o2::itsmft::TrackingMode::Type trMode)
    : mGGCCDBRequest(std::move(gr)), mUseMC(useMC), mTracking(useMC, trMode, false)
  {
  }
  ~CATrackerDPL() override = default;

  void init(framework::InitContext& ic) final;
  void run(framework::ProcessingContext& pc) final;
  void finaliseCCDB(framework::ConcreteDataMatcher& matcher, void* obj) final;

 private:
  void updateTimeDependentParams(framework::ProcessingContext& pc);

  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  bool mUseMC = false;
  bool mTrackingInitialised = false;
  o2::itsmft::tracking::ITSMFTTrackingInterfaceMFT mTracking;
};

o2::framework::DataProcessorSpec getCATrackerSpec(bool useMC, bool useGeom, bool useIRFrames,
                                                  o2::itsmft::TrackingMode::Type trMode);

} // namespace o2::mft

#endif // O2_MFT_CATRACKERSPEC_H_
