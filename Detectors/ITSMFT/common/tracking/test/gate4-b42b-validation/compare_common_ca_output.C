// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Field-level comparison for the B4.2b legacy/CommonTrack output A/B gate.
// It deliberately compares the persisted TTree leaves and flattened vector
// products, rather than ROOT file bytes (which contain run-specific metadata).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "DataFormatsITSMFT/ROFRecord.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "TBranch.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TObjArray.h"
#include "TTree.h"

namespace
{
template <typename T>
bool compareVectorProduct(TTree& left, TTree& right, const char* name)
{
  left.ResetBranchAddresses();
  right.ResetBranchAddresses();
  std::vector<T>* leftValue = nullptr;
  std::vector<T>* rightValue = nullptr;
  if (left.SetBranchAddress(name, &leftValue) < 0 || right.SetBranchAddress(name, &rightValue) < 0) {
    std::cerr << "missing vector product " << name << '\n';
    return false;
  }
  for (Long64_t entry = 0; entry < left.GetEntries(); ++entry) {
    left.GetBranch(name)->GetEntry(entry);
    right.GetBranch(name)->GetEntry(entry);
    if (leftValue == nullptr || rightValue == nullptr || *leftValue != *rightValue) {
      std::cerr << "mismatch in " << name << " at entry " << entry << '\n';
      return false;
    }
  }
  return true;
}

bool compareROFs(TTree& left, TTree& right, const char* name)
{
  left.ResetBranchAddresses();
  right.ResetBranchAddresses();
  auto* leftBranch = left.GetBranch(name);
  auto* rightBranch = right.GetBranch(name);
  if (leftBranch == nullptr || rightBranch == nullptr) {
    std::cerr << "missing ROF product " << name << '\n';
    return false;
  }
  std::vector<o2::itsmft::ROFRecord>* leftValue = nullptr;
  std::vector<o2::itsmft::ROFRecord>* rightValue = nullptr;
  if (left.SetBranchAddress(name, &leftValue) < 0 || right.SetBranchAddress(name, &rightValue) < 0) {
    std::cerr << "cannot read ROF product " << name << '\n';
    return false;
  }
  for (Long64_t entry = 0; entry < left.GetEntries(); ++entry) {
    leftBranch->GetEntry(entry);
    rightBranch->GetEntry(entry);
    if (leftValue == nullptr || rightValue == nullptr || leftValue->size() != rightValue->size()) {
      std::cerr << "mismatch in " << name << " at entry " << entry << '\n';
      return false;
    }
    for (size_t index = 0; index < leftValue->size(); ++index) {
      const auto& legacy = (*leftValue)[index];
      const auto& common = (*rightValue)[index];
      if (legacy.getBCData().bc != common.getBCData().bc || legacy.getBCData().orbit != common.getBCData().orbit ||
          legacy.getROFrame() != common.getROFrame() || legacy.getFlags() != common.getFlags() ||
          legacy.getFirstEntry() != common.getFirstEntry() || legacy.getNEntries() != common.getNEntries()) {
        std::cerr << "mismatch in " << name << " at entry " << entry << " ROF " << index << '\n';
        return false;
      }
    }
  }
  return true;
}

bool compareLabels(TTree& left, TTree& right, const char* name)
{
  left.ResetBranchAddresses();
  right.ResetBranchAddresses();
  auto* leftBranch = left.GetBranch(name);
  auto* rightBranch = right.GetBranch(name);
  if (leftBranch == nullptr || rightBranch == nullptr) {
    if (leftBranch == rightBranch) {
      return true; // Optional MC is absent from both products.
    }
    std::cerr << "MC branch presence differs for " << name << '\n';
    return false;
  }
  std::vector<o2::MCCompLabel>* leftValue = nullptr;
  std::vector<o2::MCCompLabel>* rightValue = nullptr;
  if (left.SetBranchAddress(name, &leftValue) < 0 || right.SetBranchAddress(name, &rightValue) < 0) {
    std::cerr << "cannot read MC product " << name << '\n';
    return false;
  }
  for (Long64_t entry = 0; entry < left.GetEntries(); ++entry) {
    leftBranch->GetEntry(entry);
    rightBranch->GetEntry(entry);
    if (leftValue == nullptr || rightValue == nullptr || leftValue->size() != rightValue->size()) {
      std::cerr << "mismatch in " << name << " at entry " << entry << '\n';
      return false;
    }
    for (size_t index = 0; index < leftValue->size(); ++index) {
      if ((*leftValue)[index].getRawValue() != (*rightValue)[index].getRawValue()) {
        std::cerr << "mismatch in " << name << " at entry " << entry << " label " << index << '\n';
        return false;
      }
    }
  }
  return true;
}

bool isMFTForwardStateOrCovarianceLeaf(const std::string& name)
{
  return name.find("MFTTrack") != std::string::npos &&
         (name.find("mParameters") != std::string::npos || name.find("mCovariances") != std::string::npos ||
          name.find(".mZ") != std::string::npos || name.find("mTrackChi2") != std::string::npos);
}

struct ProjectionStats {
  double maxAbsolute{};
  double maxRelative{};
  size_t values{};
};

bool compareLeaves(TTree& left, TTree& right, bool allowMFTFloatProjection, ProjectionStats& projection)
{
  if (left.GetEntries() != right.GetEntries()) {
    std::cerr << "tree entry count differs\n";
    return false;
  }
  auto* leftLeaves = left.GetListOfLeaves();
  auto* rightLeaves = right.GetListOfLeaves();
  if (leftLeaves->GetEntries() != rightLeaves->GetEntries()) {
    std::cerr << "leaf count differs\n";
    return false;
  }
  for (int leafIndex = 0; leafIndex < leftLeaves->GetEntries(); ++leafIndex) {
    auto* leftLeaf = static_cast<TLeaf*>(leftLeaves->At(leafIndex));
    auto* rightLeaf = right.GetLeaf(leftLeaf->GetName());
    if (rightLeaf == nullptr) {
      std::cerr << "missing leaf " << leftLeaf->GetName() << '\n';
      return false;
    }
    for (Long64_t entry = 0; entry < left.GetEntries(); ++entry) {
      left.GetEntry(entry);
      right.GetEntry(entry);
      if (leftLeaf->GetNdata() != rightLeaf->GetNdata()) {
        std::cerr << "array length mismatch in " << leftLeaf->GetName() << " at entry " << entry << '\n';
        return false;
      }
      for (int value = 0; value < leftLeaf->GetNdata(); ++value) {
        const double legacy = leftLeaf->GetValue(value);
        const double common = rightLeaf->GetValue(value);
        const bool projected = allowMFTFloatProjection && isMFTForwardStateOrCovarianceLeaf(leftLeaf->GetName());
        const double expected = projected ? static_cast<double>(static_cast<float>(legacy)) : legacy;
        if (projected) {
          const double absolute = std::abs(common - legacy);
          projection.maxAbsolute = std::max(projection.maxAbsolute, absolute);
          projection.maxRelative = std::max(projection.maxRelative, absolute / std::max(1., std::abs(legacy)));
          ++projection.values;
        }
        if (common != expected) {
          std::cerr << "mismatch in " << leftLeaf->GetName() << " at entry " << entry << " element " << value << ": legacy=" << legacy
                    << ", expected=" << expected << ", common=" << common << '\n';
          return false;
        }
      }
    }
  }
  return true;
}
} // namespace

int compare_common_ca_output(const char* leftFile, const char* rightFile, bool mft = false)
{
  TFile leftFileHandle{leftFile};
  TFile rightFileHandle{rightFile};
  auto* left = dynamic_cast<TTree*>(leftFileHandle.Get("o2sim"));
  auto* right = dynamic_cast<TTree*>(rightFileHandle.Get("o2sim"));
  if (left == nullptr || right == nullptr) {
    std::cerr << "missing o2sim tree\n";
    return 1;
  }
  ProjectionStats projection;
  bool equal = compareLeaves(*left, *right, mft, projection);
  equal &= compareVectorProduct<int>(*left, *right, mft ? "MFTTrackClusIdx" : "ITSTrackClusIdx");
  if (mft) {
    equal &= compareVectorProduct<unsigned short>(*left, *right, "MFTTrackSeedPattern");
    equal &= compareROFs(*left, *right, "MFTTracksROF");
    equal &= compareLabels(*left, *right, "MFTTrackMCTruth");
  } else {
    equal &= compareROFs(*left, *right, "ITSTracksROF");
    equal &= compareLabels(*left, *right, "ITSTrackMCTruth");
  }
  if (mft) {
    std::cout << "MFT forward float projection: " << projection.values << " values, max abs " << projection.maxAbsolute << ", max rel "
              << projection.maxRelative << '\n';
  }
  std::cout << (equal ? "Common-CA output products match field-by-field\n" : "Common-CA output products differ\n");
  return equal ? 0 : 1;
}
