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

#ifndef ALICEO2_ITSMFT_TRACKING_MCLABELACCUMULATOR_H_
#define ALICEO2_ITSMFT_TRACKING_MCLABELACCUMULATOR_H_

#ifndef GPUCA_GPUCODE

#include <cstddef>
#include <vector>

#include <gsl/gsl>

#include "SimulationDataFormat/MCCompLabel.h"

namespace o2::itsmft::tracking
{

// Host-side post-tracking MC-label reduction. Each call represents one attached
// cluster, including an empty label span. Identity is MCCompLabel::operator==,
// so the first encountered object remains the representative of that identity.
class MCLabelAccumulator
{
 public:
  void addCluster(gsl::span<const MCCompLabel> labels)
  {
    ++mAttachedClusters;
    for (const auto& label : labels) {
      bool foundCandidate = false;
      for (auto& candidate : mCandidates) {
        if (label == candidate.representative) {
          // Repeated identities in this span see the same cluster ordinal and
          // therefore contribute only once without per-cluster storage.
          if (candidate.lastSeenCluster != mAttachedClusters) {
            ++candidate.count;
            candidate.lastSeenCluster = mAttachedClusters;
          }
          foundCandidate = true;
          break;
        }
      }
      if (!foundCandidate) {
        mCandidates.push_back({label, 1, mAttachedClusters});
      }
    }
  }

  MCCompLabel finalize() const
  {
    MCCompLabel winner;
    if (mAttachedClusters == 0 || mCandidates.empty()) {
      winner.setFakeFlag();
      return winner;
    }

    const Candidate* best = &mCandidates.front();
    for (const auto& candidate : mCandidates) {
      if (candidate.count > best->count) {
        best = &candidate;
      }
    }

    winner = best->representative;
    if (best->count != mAttachedClusters) {
      winner.setFakeFlag();
    }
    return winner;
  }

  size_t getAttachedClusterCount() const { return mAttachedClusters; }

 private:
  struct Candidate {
    MCCompLabel representative;
    size_t count;
    size_t lastSeenCluster;
  };

  size_t mAttachedClusters = 0;
  std::vector<Candidate> mCandidates;
};

} // namespace o2::itsmft::tracking

#endif // GPUCA_GPUCODE
#endif
