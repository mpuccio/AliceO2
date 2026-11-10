// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef GTRK_WORKFLOW_ITSROFSELECTION_H_
#define GTRK_WORKFLOW_ITSROFSELECTION_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "Framework/InputSpec.h"

#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace o2::framework
{
class ProcessingContext;
}

namespace o2::gtrk::workflow
{

/// ITS's native fastMultConfig selection: cluster/random selection before
/// vertexing, vertex-contributor selection before tracking. It consumes raw
/// cluster ROFs so that expanded measurement time margins cannot change a
/// cluster's ROF identity. Both host and device receive the selected rows.
class ITSROFSelection
{
 public:
  static void addInputs(std::vector<o2::framework::InputSpec>& inputs);
  ITSROFSelection(o2::framework::ProcessingContext& pc, std::span<const FrameTiming> timing, int orbitsPerTF);
  ~ITSROFSelection();
  ITSROFSelection(const ITSROFSelection&) = delete;
  ITSROFSelection& operator=(const ITSROFSelection&) = delete;

  /// Vertices remain in the output, as in native ITS; the cut disables their
  /// overlapping ROFs for tracking. Returns whether a vertex cut is active.
  bool selectWithVertices(std::span<const VertexSeed> seeds, std::span<const std::uint32_t> contributors);
  [[nodiscard]] std::vector<std::uint8_t> clockMask() const;

  struct Measurements {
    MeasurementStore columns;
    std::vector<std::uint32_t> offsets;
    std::vector<std::uint32_t> originalRows; ///< selected row -> original input measurement row
  };

  /// Preserve bucket ordering and original measurement identities. No copy
  /// when every row passes. originalRows composes mappings after vertex cuts.
  [[nodiscard]] std::optional<Measurements> filter(MeasurementView source, std::span<const std::uint32_t> offsets, std::span<const std::uint32_t> originalRows = {}) const;

 private:
  struct State;
  std::unique_ptr<State> mState;
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_ITSROFSELECTION_H_
