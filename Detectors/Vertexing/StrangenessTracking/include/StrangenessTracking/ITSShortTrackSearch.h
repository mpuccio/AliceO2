// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// Distributed under the GNU GPL v3. CERN privileges and immunities are preserved.

#ifndef O2_STRANGENESS_TRACKING_ITSSHORTTRACKSEARCH_H
#define O2_STRANGENESS_TRACKING_ITSSHORTTRACKSEARCH_H

#include <array>
#include <cstddef>
#include <set>

namespace o2::strangeness_tracking::shorttrack
{
// Search hypotheses, in GeV/c, not cuts on a track's fitted momentum.
// Keep this schedule independent of the configured final-fit momentum floor.
inline constexpr std::array<float, 3> progressiveSearchMinPt()
{
  return {.150f, .100f, .050f};
}

// Final fitting depends on the PV, original hits and fixed final cuts, never
// on which search pass found a triplet. Remember failures as well as successes.
// This unbounded set avoids repeated work; it never reserves a cluster.
class VisitedShortTriplets
{
 public:
  bool markTested(int pv, int l0, int l1, int l2)
  {
    return mKeys.insert({pv, l0, l1, l2}).second;
  }
  std::size_t size() const { return mKeys.size(); }

 private:
  std::set<std::array<int, 4>> mKeys;
};
} // namespace o2::strangeness_tracking::shorttrack
#endif
