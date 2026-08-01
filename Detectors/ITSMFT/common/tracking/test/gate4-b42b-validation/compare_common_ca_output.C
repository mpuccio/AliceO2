// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

// Field-level comparison for the B4.2b legacy/CommonTrack output A/B gate.
// It deliberately compares the persisted TTree leaves and flattened vector
// products, rather than ROOT file bytes (which contain run-specific metadata).

#include <iostream>
#include <string>
#include <vector>

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
  std::vector<T>* leftValue = nullptr;
  std::vector<T>* rightValue = nullptr;
  if (left.SetBranchAddress(name, &leftValue) < 0 || right.SetBranchAddress(name, &rightValue) < 0) {
    std::cerr << "missing vector product " << name << '\n';
    return false;
  }
  for (Long64_t entry = 0; entry < left.GetEntries(); ++entry) {
    left.GetEntry(entry);
    right.GetEntry(entry);
    if (leftValue == nullptr || rightValue == nullptr || *leftValue != *rightValue) {
      std::cerr << "mismatch in " << name << " at entry " << entry << '\n';
      return false;
    }
  }
  return true;
}

bool compareLeaves(TTree& left, TTree& right)
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
        if (leftLeaf->GetValue(value) != rightLeaf->GetValue(value)) {
          std::cerr << "mismatch in " << leftLeaf->GetName() << " at entry " << entry << " element " << value << '\n';
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
  bool equal = compareLeaves(*left, *right);
  equal &= compareVectorProduct<int>(*left, *right, mft ? "MFTTrackClusIdx" : "ITSTrackClusIdx");
  if (mft) {
    equal &= compareVectorProduct<unsigned short>(*left, *right, "MFTTrackSeedPattern");
  }
  std::cout << (equal ? "Common-CA output products match field-by-field\n" : "Common-CA output products differ\n");
  return equal ? 0 : 1;
}
