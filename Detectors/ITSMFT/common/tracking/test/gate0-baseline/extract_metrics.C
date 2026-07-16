// Gate 0 baseline: test-only metric extraction from a tracker-replay output
// directory. Not a production macro; reused where practical from the
// counting/denominator conventions already established in
// Detectors/ITSMFT/ITS/macros/test/CheckTracksCA.C (see extractITS below),
// which prints plots interactively rather than machine-readable numbers, so
// this file only re-derives the counting logic and writes JSON instead.
//
// Usage:
//   root -l -b -q 'extract_metrics.C("<fixtureDir>", "<out.json>", "<replayDir>")'
//
// <fixtureDir> (from generate_fixture.sh) must contain o2clus_its.root,
// mftclusters.root, o2sim_Kine.root, o2sim_geometry.root /
// o2sim_geometry-aligned.root. <replayDir> (from replay_tracking.sh, may
// equal fixtureDir) must contain o2trac_its.root, mfttracks.root.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DetectorsBase/GeometryManager.h"
#include "ITSBase/GeometryTGeo.h"
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

void extractITS(const std::string& fixtureDir, const std::string& replayDir, std::ofstream& out)
{
  using namespace o2::its;

  auto gman = GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::TransformType::T2L);

  TFile fClus((fixtureDir + "/o2clus_its.root").c_str());
  auto clusTree = (TTree*)fClus.Get("o2sim");
  std::vector<o2::itsmft::CompClusterExt>* clusArr = nullptr;
  std::vector<o2::itsmft::ROFRecord>* clusROF = nullptr;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterComp", &clusArr);
  clusTree->SetBranchAddress("ITSClustersROF", &clusROF);
  clusTree->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);
  clusTree->GetEntry(0);

  int nClusters = (int)clusArr->size();
  int nROF = (int)clusROF->size();

  TFile fTrac((replayDir + "/o2trac_its.root").c_str());
  auto trTree = (TTree*)fTrac.Get("o2sim");
  std::vector<o2::its::TrackITS>* trkArr = nullptr;
  std::vector<o2::MCCompLabel>* trkLabArr = nullptr;
  trTree->SetBranchAddress("ITSTrack", &trkArr);
  trTree->SetBranchAddress("ITSTrackMCTruth", &trkLabArr);
  trTree->GetEntry(0);

  int nTracks = (int)trkArr->size();
  std::vector<double> clustersPerTrack, chi2;
  for (auto& t : *trkArr) {
    clustersPerTrack.push_back(t.getNumberOfClusters());
    chi2.push_back(t.getChi2());
  }

  TFile fKine((fixtureDir + "/o2sim_Kine.root").c_str());
  auto kineTree = (TTree*)fKine.Get("o2sim");
  std::vector<o2::MCTrack>* mcArr = nullptr;
  kineTree->SetBranchAddress("MCTrack", &mcArr);
  int nev = (int)kineTree->GetEntries();

  struct PInfo {
    bool isPrimary = false;
    unsigned char clusters = 0;
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

  // Denominator/association convention matches CheckTracksCA.C: take the
  // first label of each cluster, require it valid/correct/source-0.
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
    auto layer = gman->getLayer((*clusArr)[iC].getSensorID());
    info[evID][trackID].clusters |= (1 << layer);
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
      // Reconstructable: primary MC particle with a correct cluster label on
      // all 7 ITS layers (bitmask == 0x7f). Same cut as CheckTracksCA.C.
      if ((p.clusters & 0x7f) != 0x7f || !p.isPrimary) {
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

  out << "  \"its\": {\n";
  out << "    \"inputClusters\": " << nClusters << ",\n";
  out << "    \"inputROFs\": " << nROF << ",\n";
  out << "    \"outputTracks\": " << nTracks << ",\n";
  writeSummary(out, "clustersPerTrack", summarize(clustersPerTrack));
  writeSummary(out, "chi2", summarize(chi2));
  out << "    \"mcReconstructable\": " << reconstructable << ",\n";
  out << "    \"matched\": " << matched << ",\n";
  out << "    \"efficiency\": " << efficiency << ",\n";
  out << "    \"fakeTracks\": " << fakeTracks << ",\n";
  out << "    \"fakeRate\": " << fakeRate << ",\n";
  out << "    \"cloneTracks\": " << cloneTracks << ",\n";
  out << "    \"cloneRate\": " << cloneRate << ",\n";
  out << "    \"denominatorDefinition\": \"primary MCTrack with a correct cluster label on all 7 ITS layers (bitmask==0x7f); matches CheckTracksCA.C\"\n";
  out << "  }";
}

void extractMFT(const std::string& fixtureDir, const std::string& replayDir, std::ofstream& out)
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
  for (auto& t : *trkArr) {
    clustersPerTrack.push_back(t.getNumberOfPoints());
    chi2.push_back(t.getTrackChi2());
  }

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

  out << "  \"mft\": {\n";
  out << "    \"inputClusters\": " << nClusters << ",\n";
  out << "    \"inputROFs\": " << nROF << ",\n";
  out << "    \"outputTracks\": " << nTracks << ",\n";
  writeSummary(out, "clustersPerTrack", summarize(clustersPerTrack));
  writeSummary(out, "chi2", summarize(chi2));
  out << "    \"mcReconstructable\": " << reconstructable << ",\n";
  out << "    \"matched\": " << matched << ",\n";
  out << "    \"efficiency\": " << efficiency << ",\n";
  out << "    \"fakeTracks\": " << fakeTracks << ",\n";
  out << "    \"fakeRate\": " << fakeRate << ",\n";
  out << "    \"cloneTracks\": " << cloneTracks << ",\n";
  out << "    \"cloneRate\": " << cloneRate << ",\n";
  out << "    \"denominatorDefinition\": \"primary MCTrack with >=4 correct MFT cluster labels (no per-disk uniqueness check; simpler than the ITS 7-layer bitmask)\"\n";
  out << "  }";
}

void extract_metrics(std::string fixtureDir, std::string outFile, std::string replayDir = "")
{
  if (replayDir.empty()) {
    replayDir = fixtureDir;
  }
  o2::base::GeometryManager::loadGeometry(fixtureDir + "/o2sim"); // prefix, not a filename: appends _geometry-aligned.root

  std::ofstream out(outFile);
  out << "{\n";
  extractITS(fixtureDir, replayDir, out);
  out << ",\n";
  extractMFT(fixtureDir, replayDir, out);
  out << "\n}\n";
  out.close();
  std::cout << "[extract_metrics] wrote " << outFile << std::endl;
}
