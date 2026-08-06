// Analysis-only persisted-output comparison for the P1 native-Propagator
// physics-sign-off campaign. It deliberately observes written products only:
// no tracker code, configuration, or reconstruction decision is changed.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "SimulationDataFormat/MCCompLabel.h"
#endif

namespace
{
struct TrackRecord {
  int index{};
  int rof{-1};
  std::vector<int> refs;
  bool labelSet{};
  bool labelFake{};
  uint64_t label{};
  std::vector<double> state;
  std::vector<double> covarianceDiagonal;
};

struct RunningSummary {
  int n{};
  double sum{};
  double maximum{};

  void add(double value)
  {
    ++n;
    sum += value;
    maximum = std::max(maximum, value);
  }
};

[[noreturn]] void fail(const std::string& message)
{
  std::cerr << "[compare_native_frozen_tracks] " << message << '\n';
  std::exit(1);
}

int findROF(const std::vector<o2::itsmft::ROFRecord>& rofs, int track)
{
  for (int index = 0; index < static_cast<int>(rofs.size()); ++index) {
    const auto& rof = rofs[index];
    if (track >= rof.getFirstEntry() && track < rof.getFirstEntry() + rof.getNEntries()) {
      return index;
    }
  }
  return -1;
}

int overlap(const std::vector<int>& left, const std::vector<int>& right)
{
  int result = 0;
  size_t leftIndex = 0;
  size_t rightIndex = 0;
  while (leftIndex < left.size() && rightIndex < right.size()) {
    if (left[leftIndex] == right[rightIndex]) {
      ++result;
      ++leftIndex;
      ++rightIndex;
    } else if (left[leftIndex] < right[rightIndex]) {
      ++leftIndex;
    } else {
      ++rightIndex;
    }
  }
  return result;
}

void writeSummary(std::ostream& output, const RunningSummary& summary)
{
  output << "{\"n\": " << summary.n << ", \"meanAbs\": "
         << (summary.n == 0 ? 0. : summary.sum / summary.n) << ", \"maxAbs\": " << summary.maximum << '}';
}

std::vector<TrackRecord> loadITS(const std::string& path)
{
  TFile file(path.c_str());
  if (file.IsZombie()) {
    fail("cannot open ITS file " + path);
  }
  auto* tree = dynamic_cast<TTree*>(file.Get("o2sim"));
  if (tree == nullptr) {
    fail("missing o2sim tree in " + path);
  }
  std::vector<o2::its::TrackITS>* tracks = nullptr;
  std::vector<int>* refs = nullptr;
  std::vector<o2::itsmft::ROFRecord>* rofs = nullptr;
  std::vector<o2::MCCompLabel>* labels = nullptr;
  if (tree->SetBranchAddress("ITSTrack", &tracks) < 0 || tree->SetBranchAddress("ITSTrackClusIdx", &refs) < 0 ||
      tree->SetBranchAddress("ITSTracksROF", &rofs) < 0 || tree->SetBranchAddress("ITSTrackMCTruth", &labels) < 0 ||
      tree->GetEntry(0) <= 0 || tracks == nullptr || refs == nullptr || rofs == nullptr || labels == nullptr) {
    fail("cannot read required ITS products from " + path);
  }
  if (tracks->size() != labels->size()) {
    fail("ITS track/label cardinality differs in " + path);
  }
  std::vector<TrackRecord> result;
  result.reserve(tracks->size());
  for (int index = 0; index < static_cast<int>(tracks->size()); ++index) {
    const auto& track = (*tracks)[index];
    TrackRecord record;
    record.index = index;
    record.rof = findROF(*rofs, index);
    const int first = track.getFirstClusterEntry();
    const int size = track.getNumberOfClusters();
    if (first < 0 || size < 0 || first + size > static_cast<int>(refs->size())) {
      fail("invalid ITS cluster-reference range in " + path);
    }
    record.refs.assign(refs->begin() + first, refs->begin() + first + size);
    std::sort(record.refs.begin(), record.refs.end());
    const auto& label = (*labels)[index];
    record.labelSet = label.isSet();
    record.labelFake = record.labelSet && label.isFake();
    record.label = label.getRawValue();
    record.state = {track.getX(), track.getAlpha(), track.getY(), track.getZ(), track.getSnp(), track.getTgl(), track.getQ2Pt(), track.getChi2()};
    for (int diagonal = 0; diagonal < 5; ++diagonal) {
      record.covarianceDiagonal.push_back(track.getCovarElem(diagonal, diagonal));
    }
    result.push_back(std::move(record));
  }
  return result;
}

std::vector<TrackRecord> loadMFT(const std::string& path)
{
  TFile file(path.c_str());
  if (file.IsZombie()) {
    fail("cannot open MFT file " + path);
  }
  auto* tree = dynamic_cast<TTree*>(file.Get("o2sim"));
  if (tree == nullptr) {
    fail("missing o2sim tree in " + path);
  }
  std::vector<o2::mft::TrackMFT>* tracks = nullptr;
  std::vector<int>* refs = nullptr;
  std::vector<o2::itsmft::ROFRecord>* rofs = nullptr;
  std::vector<o2::MCCompLabel>* labels = nullptr;
  if (tree->SetBranchAddress("MFTTrack", &tracks) < 0 || tree->SetBranchAddress("MFTTrackClusIdx", &refs) < 0 ||
      tree->SetBranchAddress("MFTTracksROF", &rofs) < 0 || tree->SetBranchAddress("MFTTrackMCTruth", &labels) < 0 ||
      tree->GetEntry(0) <= 0 || tracks == nullptr || refs == nullptr || rofs == nullptr || labels == nullptr) {
    fail("cannot read required MFT products from " + path);
  }
  if (tracks->size() != labels->size()) {
    fail("MFT track/label cardinality differs in " + path);
  }
  std::vector<TrackRecord> result;
  result.reserve(tracks->size());
  for (int index = 0; index < static_cast<int>(tracks->size()); ++index) {
    const auto& track = (*tracks)[index];
    TrackRecord record;
    record.index = index;
    record.rof = findROF(*rofs, index);
    const int first = track.getExternalClusterIndexOffset();
    const int size = track.getNumberOfPoints();
    if (first < 0 || size < 0 || first + size > static_cast<int>(refs->size())) {
      fail("invalid MFT cluster-reference range in " + path);
    }
    record.refs.assign(refs->begin() + first, refs->begin() + first + size);
    std::sort(record.refs.begin(), record.refs.end());
    const auto& label = (*labels)[index];
    record.labelSet = label.isSet();
    record.labelFake = record.labelSet && label.isFake();
    record.label = label.getRawValue();
    record.state = {track.getX(), track.getY(), track.getZ(), track.getPhi(), track.getTanl(), track.getInvQPt(), track.getTrackChi2()};
    for (int diagonal = 0; diagonal < 5; ++diagonal) {
      record.covarianceDiagonal.push_back(track.getCovariances()(diagonal, diagonal));
    }
    result.push_back(std::move(record));
  }
  return result;
}

void compare(const std::vector<TrackRecord>& native, const std::vector<TrackRecord>& frozen, int nominalSurfaces, std::ostream& output)
{
  struct Candidate {
    int nativeIndex{};
    int frozenIndex{};
    int shared{};
    bool sameLabel{};
    double jaccard{};
  };
  std::vector<Candidate> candidates;
  for (int nativeIndex = 0; nativeIndex < static_cast<int>(native.size()); ++nativeIndex) {
    for (int frozenIndex = 0; frozenIndex < static_cast<int>(frozen.size()); ++frozenIndex) {
      if (native[nativeIndex].rof < 0 || native[nativeIndex].rof != frozen[frozenIndex].rof) {
        continue;
      }
      const int shared = overlap(native[nativeIndex].refs, frozen[frozenIndex].refs);
      if (shared == 0) {
        continue;
      }
      const int unionSize = static_cast<int>(native[nativeIndex].refs.size() + frozen[frozenIndex].refs.size()) - shared;
      candidates.push_back({nativeIndex, frozenIndex, shared,
                            native[nativeIndex].labelSet && frozen[frozenIndex].labelSet && native[nativeIndex].label == frozen[frozenIndex].label,
                            static_cast<double>(shared) / unionSize});
    }
  }
  std::sort(candidates.begin(), candidates.end(), [](const Candidate& left, const Candidate& right) {
    if (left.sameLabel != right.sameLabel) {
      return left.sameLabel > right.sameLabel;
    }
    if (left.shared != right.shared) {
      return left.shared > right.shared;
    }
    if (left.jaccard != right.jaccard) {
      return left.jaccard > right.jaccard;
    }
    return std::tie(left.nativeIndex, left.frozenIndex) < std::tie(right.nativeIndex, right.frozenIndex);
  });
  std::vector<bool> nativeUsed(native.size());
  std::vector<bool> frozenUsed(frozen.size());
  RunningSummary overlapSummary;
  RunningSummary jaccardSummary;
  RunningSummary nativeOnlyHoles;
  RunningSummary frozenOnlyHoles;
  std::vector<RunningSummary> state(native.empty() ? 0 : native.front().state.size());
  std::vector<RunningSummary> covariance(native.empty() ? 0 : native.front().covarianceDiagonal.size());
  int matched = 0;
  int sameLabel = 0;
  int exactReferences = 0;
  for (const auto& candidate : candidates) {
    if (nativeUsed[candidate.nativeIndex] || frozenUsed[candidate.frozenIndex]) {
      continue;
    }
    nativeUsed[candidate.nativeIndex] = true;
    frozenUsed[candidate.frozenIndex] = true;
    ++matched;
    sameLabel += candidate.sameLabel;
    exactReferences += native[candidate.nativeIndex].refs == frozen[candidate.frozenIndex].refs;
    overlapSummary.add(candidate.shared);
    jaccardSummary.add(candidate.jaccard);
    for (size_t index = 0; index < state.size(); ++index) {
      state[index].add(std::abs(native[candidate.nativeIndex].state[index] - frozen[candidate.frozenIndex].state[index]));
    }
    for (size_t index = 0; index < covariance.size(); ++index) {
      covariance[index].add(std::abs(native[candidate.nativeIndex].covarianceDiagonal[index] - frozen[candidate.frozenIndex].covarianceDiagonal[index]));
    }
  }
  for (int index = 0; index < static_cast<int>(native.size()); ++index) {
    if (!nativeUsed[index]) {
      nativeOnlyHoles.add(nominalSurfaces - native[index].refs.size());
    }
  }
  for (int index = 0; index < static_cast<int>(frozen.size()); ++index) {
    if (!frozenUsed[index]) {
      frozenOnlyHoles.add(nominalSurfaces - frozen[index].refs.size());
    }
  }
  output << "{\n  \"nativeTracks\": " << native.size() << ",\n  \"frozenTracks\": " << frozen.size()
         << ",\n  \"association\": {\"method\": \"same output ROF; greedy maximum tuple (same MC raw label, shared cluster references, Jaccard overlap, native index, frozen index); pairs require at least one shared reference\", \"matched\": "
         << matched << ", \"nativeOnly\": " << native.size() - matched << ", \"frozenOnly\": " << frozen.size() - matched
         << ", \"sameMCLabel\": " << sameLabel << ", \"exactReferenceSet\": " << exactReferences << ", \"sharedReferences\": ";
  writeSummary(output, overlapSummary);
  output << ", \"jaccard\": ";
  writeSummary(output, jaccardSummary);
  output << "},\n  \"matchedAbsoluteDelta\": {\"state\": [";
  for (size_t index = 0; index < state.size(); ++index) {
    if (index != 0) {
      output << ',';
    }
    writeSummary(output, state[index]);
  }
  output << "], \"covarianceDiagonal\": [";
  for (size_t index = 0; index < covariance.size(); ++index) {
    if (index != 0) {
      output << ',';
    }
    writeSummary(output, covariance[index]);
  }
  output << "]},\n  \"unmatchedHoles\": {\"nativeOnly\": ";
  writeSummary(output, nativeOnlyHoles);
  output << ", \"frozenOnly\": ";
  writeSummary(output, frozenOnlyHoles);
  output << "},\n  \"limitations\": \"Persisted writer products expose accepted tracks, references, MC labels, ROFs, states and covariance, but not CA-stage populations, rejected refit reasons, per-step material, or pre-publication timestamps. Those require separately authorized diagnostic instrumentation.\"\n}\n";
}
} // namespace

void compare_native_frozen_tracks(std::string detector, std::string nativeFile, std::string frozenFile, std::string outputFile)
{
  std::ofstream output(outputFile);
  if (!output) {
    fail("cannot open output " + outputFile);
  }
  if (detector == "ITS") {
    compare(loadITS(nativeFile), loadITS(frozenFile), 7, output);
  } else if (detector == "MFT") {
    compare(loadMFT(nativeFile), loadMFT(frozenFile), 10, output);
  } else {
    fail("detector must be ITS or MFT");
  }
  std::cout << "[compare_native_frozen_tracks] wrote " << outputFile << '\n';
}
