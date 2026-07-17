// Gate 2 common-CA MFT baseline: test-only metric extraction from an
// o2-mft-ca-tracker-workflow replay output directory (replay_tracking_common_ca.sh).
//
// This is a SIBLING to extract_metrics.C, not a replacement or a reuse of
// its ITS+MFT-combined entry point: extract_metrics.C's extract_metrics()
// unconditionally reads o2trac_its.root from its replayDir, but the common-CA
// replay directory only ever contains an MFT common-CA mfttracks.root (no
// ITS component runs in this bounded replay), so calling it directly here
// would fail on a missing ITS file. The counting/denominator logic below is
// otherwise an unmodified copy of extract_metrics.C's extractMFT(): the
// o2-mft-ca-tracker-workflow output is written by the same TrackWriterSpec
// (identical "mfttracks.root" / "o2sim" tree / "MFTTrack" and
// "MFTTrackMCTruth" branch names, same o2::mft::TrackMFT type, same
// getNumberOfPoints()/getTrackChi2() accessors) as the legacy MFT tracker,
// so the semantics genuinely match and no new extraction logic was needed --
// only a smaller driver that does not also require an ITS replay output.
//
// Every ROOT file/tree/branch/entry access below is explicitly validated
// (open succeeded, tree/branch exists, GetEntry succeeded, resulting pointer
// non-null) before use, and the output JSON is only written after the full
// extraction succeeds -- a failure prints a message to stderr and exits
// non-zero instead of dereferencing a null pointer or leaving a partial
// output file on disk.
//
// Usage:
//   root -l -b -q 'extract_metrics_common_ca.C("<fixtureDir>", "<out.json>", "<replayDir>")'
//
// <fixtureDir> (from generate_fixture.sh) must contain mftclusters.root,
// o2sim_Kine.root. <replayDir> (from replay_tracking_common_ca.sh) must
// contain mfttracks.root.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TMD5.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#endif

namespace
{
struct Summary {
  int n = 0;
  double mean = 0, median = 0, min = 0, max = 0;
};

Summary summarize(std::vector<double> v)
{
  Summary s;
  s.n = (int)v.size();
  if (v.empty()) {
    return s;
  }
  std::sort(v.begin(), v.end());
  s.min = v.front();
  s.max = v.back();
  double sum = 0;
  for (auto x : v) {
    sum += x;
  }
  s.mean = sum / v.size();
  s.median = v[v.size() / 2];
  return s;
}

void writeSummary(std::ostream& out, const char* name, const Summary& s, bool last = false)
{
  out << "    \"" << name << "\": {\"n\": " << s.n << ", \"mean\": " << s.mean
      << ", \"median\": " << s.median << ", \"min\": " << s.min << ", \"max\": " << s.max
      << "}" << (last ? "\n" : ",\n");
}

// Every validation failure below goes through this single choke point: print
// a specific message to stderr and exit non-zero immediately. Never returns,
// so callers do not need to guard against a "failed but kept going" state --
// there is no path from a validation failure back into extraction logic that
// would dereference an unchecked pointer or reach the output-writing stage.
[[noreturn]] void failExtract(const std::string& msg)
{
  std::cerr << "[extract_metrics_common_ca] " << msg << std::endl;
  std::exit(1);
}

void requireOpen(TFile& f, const std::string& path)
{
  if (f.IsZombie() || !f.IsOpen()) {
    failExtract("failed to open ROOT file: " + path);
  }
}

TTree* requireTree(TFile& f, const char* name, const std::string& fileDesc)
{
  auto* t = dynamic_cast<TTree*>(f.Get(name));
  if (!t) {
    failExtract("missing tree '" + std::string(name) + "' in " + fileDesc);
  }
  return t;
}

template <typename T>
void requireBranch(TTree* t, const char* name, T*& ptr, const std::string& treeDesc)
{
  if (!t->GetBranch(name)) {
    failExtract("missing branch '" + std::string(name) + "' in " + treeDesc);
  }
  if (t->SetBranchAddress(name, &ptr) < 0) {
    failExtract("SetBranchAddress failed for branch '" + std::string(name) + "' in " + treeDesc);
  }
}

void requireEntry(TTree* t, Long64_t entry, const std::string& treeDesc)
{
  if (entry >= t->GetEntries()) {
    failExtract(treeDesc + " has no entry " + std::to_string(entry) +
                " (nEntries=" + std::to_string(t->GetEntries()) + ")");
  }
  if (t->GetEntry(entry) <= 0) {
    failExtract("GetEntry(" + std::to_string(entry) + ") failed for " + treeDesc);
  }
}

template <typename T>
void requireNonNull(T* ptr, const std::string& what)
{
  if (!ptr) {
    failExtract("null branch object after read: " + what);
  }
}
} // namespace

// Copy of extract_metrics.C's extractMFT(), with explicit validation added
// at every file/tree/branch/entry access, except for the JSON key
// ("mftCommonCA" instead of "mft") and the denominatorDefinition string, so
// common-CA and legacy MFT results are never confused when both JSON files
// are read side by side.
void extractMFTCommonCA(const std::string& fixtureDir, const std::string& replayDir, std::ostream& out)
{
  using namespace o2::mft;
  const int minClustersMFT = 4; // matches MFTTrackingParam.MinTrackPointsCA default

  const std::string clusPath = fixtureDir + "/mftclusters.root";
  TFile fClus(clusPath.c_str());
  requireOpen(fClus, clusPath);
  auto clusTree = requireTree(fClus, "o2sim", clusPath);
  std::vector<o2::itsmft::CompClusterExt>* clusArr = nullptr;
  std::vector<o2::itsmft::ROFRecord>* clusROF = nullptr;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArr = nullptr;
  requireBranch(clusTree, "MFTClusterComp", clusArr, clusPath + ":o2sim");
  requireBranch(clusTree, "MFTClustersROF", clusROF, clusPath + ":o2sim");
  requireBranch(clusTree, "MFTClusterMCTruth", clusLabArr, clusPath + ":o2sim");
  requireEntry(clusTree, 0, clusPath + ":o2sim");
  requireNonNull(clusArr, clusPath + ":o2sim:MFTClusterComp");
  requireNonNull(clusROF, clusPath + ":o2sim:MFTClustersROF");
  requireNonNull(clusLabArr, clusPath + ":o2sim:MFTClusterMCTruth");

  int nClusters = (int)clusArr->size();
  int nROF = (int)clusROF->size();

  const std::string tracPath = replayDir + "/mfttracks.root";
  TFile fTrac(tracPath.c_str());
  requireOpen(fTrac, tracPath);
  auto trTree = requireTree(fTrac, "o2sim", tracPath);
  std::vector<o2::mft::TrackMFT>* trkArr = nullptr;
  std::vector<o2::MCCompLabel>* trkLabArr = nullptr;
  requireBranch(trTree, "MFTTrack", trkArr, tracPath + ":o2sim");
  requireBranch(trTree, "MFTTrackMCTruth", trkLabArr, tracPath + ":o2sim");
  requireEntry(trTree, 0, tracPath + ":o2sim");
  requireNonNull(trkArr, tracPath + ":o2sim:MFTTrack");
  requireNonNull(trkLabArr, tracPath + ":o2sim:MFTTrackMCTruth");

  int nTracks = (int)trkArr->size();
  std::vector<double> clustersPerTrack, chi2;
  // Ordered-content hash over (nPoints, chi2, x, y, z, phi, tanl, invQPt) per
  // track, independent of the mfttracks.root file's own bytes: a ROOT TFile
  // embeds a per-write UUID/timestamp in its metadata, so hashing the raw
  // file bytes is NOT repeatability evidence across separately-produced runs
  // even when this tuple is identical. This hash covers exactly these eight
  // fields per track, in track order -- it is not a hash of every branch
  // written by the track writer (e.g. MFTTrackROF, MFTTrackClusIdx,
  // MFTTrackSeedPattern are not included).
  TMD5 contentMD5;
  char buf[256];
  for (auto& t : *trkArr) {
    clustersPerTrack.push_back(t.getNumberOfPoints());
    chi2.push_back(t.getTrackChi2());
    int len = snprintf(buf, sizeof(buf), "%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g;",
                        t.getNumberOfPoints(), t.getTrackChi2(), t.getX(), t.getY(), t.getZ(),
                        t.getPhi(), t.getTanl(), t.getInvQPt());
    contentMD5.Update((UChar_t*)buf, len);
  }
  contentMD5.Final();
  std::string trackContentHash = contentMD5.AsString();

  const std::string kinePath = fixtureDir + "/o2sim_Kine.root";
  TFile fKine(kinePath.c_str());
  requireOpen(fKine, kinePath);
  auto kineTree = requireTree(fKine, "o2sim", kinePath);
  std::vector<o2::MCTrack>* mcArr = nullptr;
  requireBranch(kineTree, "MCTrack", mcArr, kinePath + ":o2sim");
  int nev = (int)kineTree->GetEntries();
  if (nev <= 0) {
    failExtract(kinePath + ":o2sim has no entries (nev=" + std::to_string(nev) + ")");
  }

  struct PInfo {
    bool isPrimary = false;
    int nClusters = 0;
    int isReco = 0;
    int isFake = 0;
  };
  std::vector<std::vector<PInfo>> info(nev);
  for (int n = 0; n < nev; n++) {
    requireEntry(kineTree, n, kinePath + ":o2sim");
    requireNonNull(mcArr, kinePath + ":o2sim:MCTrack@entry" + std::to_string(n));
    info[n].resize(mcArr->size());
    for (size_t i = 0; i < mcArr->size(); ++i) {
      info[n][i].isPrimary = mcArr->at(i).isPrimary();
    }
  }

  // Simpler denominator than ITS: raw count of correct cluster labels, no
  // per-disk uniqueness requirement (no MFT GeometryTGeo lookup here).
  for (unsigned int iC = 0; iC < clusArr->size(); ++iC) {
    auto labs = clusLabArr->getLabels(iC);
    if (labs.empty()) {
      continue;
    }
    auto lab = labs[0];
    if (!lab.isValid() || lab.getSourceID() != 0 || !lab.isCorrect()) {
      continue;
    }
    int trackID, evID, srcID;
    bool fake;
    lab.get(trackID, evID, srcID, fake);
    if (evID < 0 || evID >= nev || trackID < 0 || trackID >= (int)info[evID].size()) {
      continue;
    }
    info[evID][trackID].nClusters++;
  }

  for (unsigned int iT = 0; iT < trkLabArr->size(); ++iT) {
    auto lab = trkLabArr->at(iT);
    if (!lab.isSet()) {
      continue;
    }
    int trackID, evID, srcID;
    bool fake;
    lab.get(trackID, evID, srcID, fake);
    if (evID < 0 || evID >= nev || trackID < 0 || trackID >= (int)info[evID].size()) {
      continue;
    }
    info[evID][trackID].isReco += !fake;
    info[evID][trackID].isFake += fake;
  }

  int reconstructable = 0, matched = 0, cloneTracks = 0, fakeTracks = 0;
  for (auto& ev : info) {
    for (auto& p : ev) {
      if (p.nClusters < minClustersMFT || !p.isPrimary) {
        continue;
      }
      reconstructable++;
      if (p.isReco) {
        matched++;
        if (p.isReco > 1) {
          cloneTracks += (p.isReco - 1);
        }
      }
    }
  }
  for (unsigned int iT = 0; iT < trkLabArr->size(); ++iT) {
    if (trkLabArr->at(iT).isSet() && trkLabArr->at(iT).isFake()) {
      fakeTracks++;
    }
  }

  double efficiency = reconstructable > 0 ? double(matched) / reconstructable : 0.0;
  double fakeRate = nTracks > 0 ? double(fakeTracks) / nTracks : 0.0;
  double cloneRate = reconstructable > 0 ? double(cloneTracks) / reconstructable : 0.0;

  out << "  \"mftCommonCA\": {\n";
  out << "    \"inputClusters\": " << nClusters << ",\n";
  out << "    \"inputROFs\": " << nROF << ",\n";
  out << "    \"outputTracks\": " << nTracks << ",\n";
  out << "    \"trackContentHash\": \"" << trackContentHash << "\",\n";
  out << "    \"trackContentHashDefinition\": \"MD5 over ordered per-track (nPoints,chi2,x,y,z,phi,tanl,invQPt), %.9g-formatted; covers exactly this tuple per track in track order, NOT every branch written by the track writer (MFTTrackROF/MFTTrackClusIdx/MFTTrackSeedPattern are excluded) and NOT a hash of mfttracks.root's bytes, which vary run-to-run due to ROOT TFile UUID/timestamp metadata even for identical content\",\n";
  writeSummary(out, "clustersPerTrack", summarize(clustersPerTrack));
  writeSummary(out, "chi2", summarize(chi2));
  out << "    \"mcReconstructable\": " << reconstructable << ",\n";
  out << "    \"matched\": " << matched << ",\n";
  out << "    \"efficiency\": " << efficiency << ",\n";
  out << "    \"fakeTracks\": " << fakeTracks << ",\n";
  out << "    \"fakeRate\": " << fakeRate << ",\n";
  out << "    \"cloneTracks\": " << cloneTracks << ",\n";
  out << "    \"cloneRate\": " << cloneRate << ",\n";
  out << "    \"denominatorDefinition\": \"primary MCTrack with >=4 correct MFT cluster labels (no per-disk uniqueness check); identical convention to extract_metrics.C's extractMFT, applied to the o2-mft-ca-tracker-workflow (common o2::itsmft::tracking core) output instead of the legacy o2-mft-reco-workflow output\"\n";
  out << "  }";
}

void extract_metrics_common_ca(std::string fixtureDir, std::string outFile, std::string replayDir = "")
{
  if (replayDir.empty()) {
    replayDir = fixtureDir;
  }
  // Extraction is built into an in-memory buffer first and outFile is only
  // touched after it fully succeeds, so a validation failure partway through
  // (via failExtract's exit(1)) never leaves a truncated/partial JSON file
  // on disk -- there is simply no output file at all on failure.
  std::ostringstream body;
  extractMFTCommonCA(fixtureDir, replayDir, body);

  std::ofstream out(outFile);
  if (!out) {
    failExtract("failed to open output file for writing: " + outFile);
  }
  out << "{\n" << body.str() << "\n}\n";
  out.close();
  if (out.fail()) {
    failExtract("failed while writing output file: " + outFile);
  }
  std::cout << "[extract_metrics_common_ca] wrote " << outFile << std::endl;
}
