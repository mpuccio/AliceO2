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

/// \file GPURecoTracking.h
/// \brief A detector's tracking inside O2's GPU reconstruction workflow
///        (GPU/Workflow's GPURecoWorkflowSpec): the tracking device's
///        timeframe, after the TPC's, on the same GPU, in the pool's free
///        memory, borrowed per timeframe (DeviceContext). Its inputs are the
///        measurement producer's message, its products the tracking device's
///        (TrackerSpec.h); the converter, the producer and the adaptor stay
///        devices of their own.

#ifndef GTRK_WORKFLOW_GPURECOTRACKING_H_
#define GTRK_WORKFLOW_GPURECOTRACKING_H_

#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InitContext.h"
#include "Framework/InputSpec.h"
#include "Framework/OutputSpec.h"
#include "Framework/ProcessingContext.h"

#include <memory>
#include <string_view>
#include <vector>

namespace o2::gpu
{
class GPUReconstruction;
} // namespace o2::gpu

namespace o2::gtrk::workflow
{

/// The calls GPURecoWorkflowSpec makes, in its own: the inputs and outputs
/// when the workflow is defined; init after its reconstruction is
/// initialised, on the same thread; finaliseCCDB before GRPGeomHelper's;
/// run after the TPC's tracking and before the pool is cleared, on the
/// thread that ran it; stop and endOfStream before the reconstruction goes.
class GPURecoTracking
{
 public:
  GPURecoTracking() = default;
  GPURecoTracking(const GPURecoTracking&) = delete;
  GPURecoTracking(GPURecoTracking&&) = delete;
  GPURecoTracking& operator=(const GPURecoTracking&) = delete;
  GPURecoTracking& operator=(GPURecoTracking&&) = delete;
  virtual ~GPURecoTracking() = default;

  /// The detector's inputs, appended to \a inputs: the measurements and the
  /// detector's own conditions (the run-wide ones are the reconstruction's).
  virtual void appendInputs(std::vector<o2::framework::InputSpec>& inputs) = 0;
  /// The detector's products, appended to \a outputs.
  virtual void appendOutputs(std::vector<o2::framework::OutputSpec>& outputs) = 0;
  /// On \a reconstruction's GPU, initialised.
  virtual void init(o2::framework::InitContext& ic, o2::gpu::GPUReconstruction& reconstruction) = 0;
  /// Whether \a matcher was one of the detector's conditions, then taken.
  virtual bool finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) = 0;
  /// The timeframe; \a reconstructed false if the TPC's failed, then its
  /// products are published empty and nothing is borrowed.
  virtual void run(o2::framework::ProcessingContext& pc, bool reconstructed) = 0;
  /// The GPU released.
  virtual void stop() = 0;
  virtual void endOfStream() = 0;
};

/// The tracking of \a detector, "its" or "mft"; std::invalid_argument for
/// another.
[[nodiscard]] std::unique_ptr<GPURecoTracking> makeGPURecoTracking(std::string_view detector, bool useMC = false);

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_GPURECOTRACKING_H_
