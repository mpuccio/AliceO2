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

    // A cluster may contain the same identity several times. Keep a local
    // identity set so it contributes at most once to the global counts.
    std::vector<MCCompLabel> clusterIdentities;
    clusterIdentities.reserve(labels.size());
    for (const auto& label : labels) {
      bool seenInCluster = false;
      for (const auto& identity : clusterIdentities) {
        if (label == identity) {
          seenInCluster = true;
          break;
        }
      }
      if (seenInCluster) {
        continue;
      }
      clusterIdentities.push_back(label);

      bool foundCandidate = false;
      for (auto& candidate : mCandidates) {
        if (label == candidate.representative) {
          ++candidate.count;
          foundCandidate = true;
          break;
        }
      }
      if (!foundCandidate) {
        mCandidates.push_back({label, 1});
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
  };

  size_t mAttachedClusters = 0;
  std::vector<Candidate> mCandidates;
};

} // namespace o2::itsmft::tracking

#endif
