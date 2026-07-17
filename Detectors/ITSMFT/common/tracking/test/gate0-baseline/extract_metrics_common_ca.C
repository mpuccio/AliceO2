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
#include <fstream>
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

void writeSummary(std::ofstream& out, const char* name, const Summary& s, bool last = false)
{
  out << "    \"" << name << "\": {\"n\": " << s.n << ", \"mean\": " << s.mean
      << ", \"median\": " << s.median << ", \"min\": " << s.min << ", \"max\": " << s.max
      << "}" << (last ? "\n" : ",\n");
}
} // namespace

// Unmodified copy of extract_metrics.C's extractMFT(), except for the JSON
// key ("mftCommonCA" instead of "mft") and the denominatorDefinition string,
// so common-CA and legacy MFT results are never confused when both JSON
// files are read side by side.
void extractMFTCommonCA(const std::string& fixtureDir, const std::string& replayDir, std::ofstream& out)
{
  using namespace o2::mft;
  const int minClustersMFT = 4; // matches MFTTrackingParam.MinTrackPointsCA default

  TFile fClus((fixtureDir + "/mftclusters.root").c_str());
  auto clusTree = (TTree*)fClus.Get("o2sim");
  std::vector<o2::itsmft::CompClusterExt>* clusArr = nullptr;
  std::vector<o2::itsmft::ROFRecord>* clusROF = nullptr;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArr = nullptr;
  clusTree->SetBranchAddress("MFTClusterComp", &clusArr);
  clusTree->SetBranchAddress("MFTClustersROF", &clusROF);
  clusTree->SetBranchAddress("MFTClusterMCTruth", &clusLabArr);
  clusTree->GetEntry(0);

  int nClusters = (int)clusArr->size();
  int nROF = (int)clusROF->size();

  TFile fTrac((replayDir + "/mfttracks.root").c_str());
  auto trTree = (TTree*)fTrac.Get("o2sim");
  std::vector<o2::mft::TrackMFT>* trkArr = nullptr;
  std::vector<o2::MCCompLabel>* trkLabArr = nullptr;
  trTree->SetBranchAddress("MFTTrack", &trkArr);
  trTree->SetBranchAddress("MFTTrackMCTruth", &trkLabArr);
  trTree->GetEntry(0);

  int nTracks = (int)trkArr->size();
  std::vector<double> clustersPerTrack, chi2;
  // Ordered-content hash over (nPoints, chi2, x, y, z, phi, tanl, invQPt) per
  // track, independent of the mfttracks.root file's own bytes: a ROOT TFile
  // embeds a per-write UUID/timestamp in its metadata, so hashing the raw
  // file bytes is NOT repeatability evidence across separately-produced runs
  // even when the tracking output is bit-identical. This hash is.
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

  TFile fKine((fixtureDir + "/o2sim_Kine.root").c_str());
  auto kineTree = (TTree*)fKine.Get("o2sim");
  std::vector<o2::MCTrack>* mcArr = nullptr;
  kineTree->SetBranchAddress("MCTrack", &mcArr);
  int nev = (int)kineTree->GetEntries();

  struct PInfo {
    bool isPrimary = false;
    int nClusters = 0;
    int isReco = 0;
    int isFake = 0;
  };
  std::vector<std::vector<PInfo>> info(nev);
  for (int n = 0; n < nev; n++) {
    kineTree->GetEntry(n);
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
  out << "    \"trackContentHashDefinition\": \"MD5 over ordered per-track (nPoints,chi2,x,y,z,phi,tanl,invQPt), %.9g-formatted; NOT a hash of mfttracks.root's bytes, which vary run-to-run due to ROOT TFile UUID/timestamp metadata even for identical content\",\n";
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
  std::ofstream out(outFile);
  out << "{\n";
  extractMFTCommonCA(fixtureDir, replayDir, out);
  out << "\n}\n";
  out.close();
  std::cout << "[extract_metrics_common_ca] wrote " << outFile << std::endl;
}
