// Gate 3 Slice 3 common-CA ITS characterization: test-only metric
// extraction from an o2-its-ca-tracker-workflow replay output directory
// (replay_tracking_its_common_ca.sh).
//
// This is a SIBLING to gate0-baseline/extract_metrics.C and
// gate0-baseline/extract_metrics_common_ca.C (the MFT common-CA
// characterization extractor), not a replacement or a modification of
// either -- no file under gate0-baseline is read or written by this macro.
// The counting/denominator convention below (primary MCTrack with a
// correct cluster label on all 7 ITS layers, bitmask==0x7f) is an
// unmodified copy of extract_metrics.C's extractITS() convention, applied
// to the o2-its-ca-tracker-workflow (common o2::itsmft::tracking core)
// output (o2trac_its_ca.root) instead of the legacy o2-its-reco-workflow
// output (o2trac_its.root). Per task instruction #5, this is
// characterization, not a bitwise-parity gate: no physics tolerance is
// invented here, and no claim is made that these numbers match the legacy
// leg's.
//
// Every ROOT file/tree/branch/entry access below is explicitly validated
// (open succeeded, tree/branch exists, GetEntry succeeded, resulting
// pointer non-null) before use, and the output JSON is only written after
// the full extraction succeeds -- a failure prints a message to stderr and
// exits non-zero instead of dereferencing a null pointer or leaving a
// partial output file on disk.
//
// Usage:
//   root -l -b -q 'extract_metrics_its_common_ca.C("<fixtureDir>", "<out.json>", "<replayDir>")'
//
// <fixtureDir> (from gate0-baseline/generate_fixture.sh) must contain
// o2clus_its.root, o2sim_Kine.root, o2sim_geometry.root /
// o2sim_geometry-aligned.root. <replayDir> (from
// replay_tracking_its_common_ca.sh) must contain o2trac_its_ca.root.

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

#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
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

void writeSummary(std::ostream& out, const char* name, const Summary& s, bool last = false)
{
  out << "    \"" << name << "\": {\"n\": " << s.n << ", \"mean\": " << s.mean
      << ", \"median\": " << s.median << ", \"min\": " << s.min << ", \"max\": " << s.max
      << "}" << (last ? "\n" : ",\n");
}

// Single choke point for every validation failure below: print a specific
// message to stderr and exit non-zero immediately. Never returns, so
// callers do not need to guard against a "failed but kept going" state --
// there is no path from a validation failure back into extraction logic
// that would dereference an unchecked pointer or reach the output-writing
// stage.
[[noreturn]] void failExtract(const std::string& msg)
{
  std::cerr << "[extract_metrics_its_common_ca] " << msg << std::endl;
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

void extractITSCommonCA(const std::string& fixtureDir, const std::string& replayDir, std::ostream& out, bool mcExpected)
{
  using namespace o2::its;

  auto gman = GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::TransformType::T2L);

  const std::string clusPath = fixtureDir + "/o2clus_its.root";
  TFile fClus(clusPath.c_str());
  requireOpen(fClus, clusPath);
  auto clusTree = requireTree(fClus, "o2sim", clusPath);
  std::vector<o2::itsmft::CompClusterExt>* clusArr = nullptr;
  std::vector<o2::itsmft::ROFRecord>* clusROF = nullptr;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArr = nullptr;
  requireBranch(clusTree, "ITSClusterComp", clusArr, clusPath + ":o2sim");
  requireBranch(clusTree, "ITSClustersROF", clusROF, clusPath + ":o2sim");
  requireBranch(clusTree, "ITSClusterMCTruth", clusLabArr, clusPath + ":o2sim");
  requireEntry(clusTree, 0, clusPath + ":o2sim");
  requireNonNull(clusArr, clusPath + ":o2sim:ITSClusterComp");
  requireNonNull(clusROF, clusPath + ":o2sim:ITSClustersROF");
  requireNonNull(clusLabArr, clusPath + ":o2sim:ITSClusterMCTruth");

  int nClusters = (int)clusArr->size();
  int nROF = (int)clusROF->size();

  const std::string tracPath = replayDir + "/o2trac_its_ca.root";
  TFile fTrac(tracPath.c_str());
  requireOpen(fTrac, tracPath);
  auto trTree = requireTree(fTrac, "o2sim", tracPath);
  std::vector<o2::its::TrackITS>* trkArr = nullptr;
  requireBranch(trTree, "ITSTrack", trkArr, tracPath + ":o2sim");
  requireEntry(trTree, 0, tracPath + ":o2sim");
  requireNonNull(trkArr, tracPath + ":o2sim:ITSTrack");

  // MC labels: the writer only creates the ITSTrackMCTruth branch when MC
  // was enabled at replay time (see TrackWriterSpec.cxx's useMC-gated
  // BranchDefinition); absence of the branch when MC was NOT requested is
  // expected, not a failure. When MC WAS requested, the branch's presence,
  // non-null read, and count-vs-track-count alignment are all recorded
  // explicitly rather than assumed.
  std::vector<o2::MCCompLabel>* trkLabArr = nullptr;
  bool mcLabelsAvailable = false;
  if (trTree->GetBranch("ITSTrackMCTruth")) {
    if (trTree->SetBranchAddress("ITSTrackMCTruth", &trkLabArr) < 0) {
      failExtract("SetBranchAddress failed for branch 'ITSTrackMCTruth' in " + tracPath + ":o2sim");
    }
    requireEntry(trTree, 0, tracPath + ":o2sim (re-read for ITSTrackMCTruth)");
    mcLabelsAvailable = (trkLabArr != nullptr);
  } else if (mcExpected) {
    failExtract("MC was requested for this replay but 'ITSTrackMCTruth' branch is missing in " + tracPath + ":o2sim");
  }

  int nTracks = (int)trkArr->size();
  std::vector<double> clustersPerTrack, chi2;
  // Ordered-content hash over (nClusters, chi2, x, alpha, y, z, snp, tgl,
  // q2pt) per track, independent of o2trac_its_ca.root's own bytes: a ROOT
  // TFile embeds a per-write UUID/timestamp in its metadata, so hashing the
  // raw file bytes is NOT repeatability evidence across separately-produced
  // runs even when this tuple is identical. These nine fields are exactly
  // o2::its::TrackITS's own documented state: getNumberOfClusters(),
  // getChi2(), and the seven o2::track::TrackParametrization fields
  // (getX/getAlpha/getY/getZ/getSnp/getTgl/getQ2Pt) that
  // TrackParametrization::hash() itself hashes over -- so this is not an
  // ad hoc field choice but the framework's own notion of a track's
  // parametrized state, plus the two ITS-specific summary fields
  // (nClusters, chi2) extract_metrics.C's extractITS() already reports.
  // This is not a hash of every branch written by the track writer
  // (ITSTrackClusIdx, ITSTracksROF are not included).
  TMD5 contentMD5;
  char buf[320];
  for (auto& t : *trkArr) {
    clustersPerTrack.push_back(t.getNumberOfClusters());
    chi2.push_back(t.getChi2());
    int len = snprintf(buf, sizeof(buf), "%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g;",
                        t.getNumberOfClusters(), t.getChi2(), t.getX(), t.getAlpha(), t.getY(), t.getZ(),
                        t.getSnp(), t.getTgl(), t.getQ2Pt());
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
    unsigned char clusters = 0;
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

  // Denominator/association convention matches extract_metrics.C's
  // extractITS() (itself matching CheckTracksCA.C): take the first label of
  // each cluster, require it valid/correct/source-0.
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

  if (mcLabelsAvailable) {
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
  }

  int reconstructable = 0, matched = 0, cloneTracks = 0, fakeTracks = 0;
  for (auto& ev : info) {
    for (auto& p : ev) {
      // Reconstructable: primary MC particle with a correct cluster label on
      // all 7 ITS layers (bitmask == 0x7f). Same cut as extractITS()/CheckTracksCA.C.
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
  if (mcLabelsAvailable) {
    for (unsigned int iT = 0; iT < trkLabArr->size(); ++iT) {
      if (trkLabArr->at(iT).isSet() && trkLabArr->at(iT).isFake()) {
        fakeTracks++;
      }
    }
  }

  double efficiency = reconstructable > 0 ? double(matched) / reconstructable : 0.0;
  double fakeRate = nTracks > 0 ? double(fakeTracks) / nTracks : 0.0;
  double cloneRate = reconstructable > 0 ? double(cloneTracks) / reconstructable : 0.0;
  bool mcLabelsAligned = mcLabelsAvailable && (int)trkLabArr->size() == nTracks;

  out << "  \"itsCommonCA\": {\n";
  out << "    \"inputClusters\": " << nClusters << ",\n";
  out << "    \"inputROFs\": " << nROF << ",\n";
  out << "    \"outputTracks\": " << nTracks << ",\n";
  out << "    \"trackContentHash\": \"" << trackContentHash << "\",\n";
  out << "    \"trackContentHashDefinition\": \"MD5 over ordered per-track (nClusters,chi2,x,alpha,y,z,snp,tgl,q2pt), %.9g-formatted; covers exactly o2::its::TrackITS::getNumberOfClusters()/getChi2() plus the seven o2::track::TrackParametrization fields TrackParametrization::hash() itself hashes over, in track order. NOT every branch written by the track writer (ITSTrackClusIdx/ITSTracksROF excluded) and NOT a hash of o2trac_its_ca.root's bytes, which vary run-to-run due to ROOT TFile UUID/timestamp metadata even for identical content\",\n";
  writeSummary(out, "clustersPerTrack", summarize(clustersPerTrack));
  writeSummary(out, "chi2", summarize(chi2));
  out << "    \"mcLabelsAvailable\": " << (mcLabelsAvailable ? "true" : "false") << ",\n";
  out << "    \"mcLabelCount\": " << (mcLabelsAvailable ? (long long)trkLabArr->size() : 0LL) << ",\n";
  out << "    \"mcLabelsAlignedWithTracks\": " << (mcLabelsAligned ? "true" : "false") << ",\n";
  out << "    \"mcReconstructable\": " << reconstructable << ",\n";
  out << "    \"matched\": " << matched << ",\n";
  out << "    \"efficiency\": " << efficiency << ",\n";
  out << "    \"fakeTracks\": " << fakeTracks << ",\n";
  out << "    \"fakeRate\": " << fakeRate << ",\n";
  out << "    \"cloneTracks\": " << cloneTracks << ",\n";
  out << "    \"cloneRate\": " << cloneRate << ",\n";
  out << "    \"denominatorDefinition\": \"primary MCTrack with a correct cluster label on all 7 ITS layers (bitmask==0x7f); identical convention to extract_metrics.C's extractITS(), applied to the o2-its-ca-tracker-workflow (common o2::itsmft::tracking core) output instead of the legacy o2-its-reco-workflow output. efficiency/fakeRate/cloneRate/matched/fakeTracks/cloneTracks are all 0 when mcLabelsAvailable is false (no MC label pass was run), which is distinct from a genuine zero-efficiency result -- check mcLabelsAvailable before interpreting these fields\"\n";
  out << "  }";
}

void extract_metrics_its_common_ca(std::string fixtureDir, std::string outFile, std::string replayDir = "", bool mcExpected = true)
{
  if (replayDir.empty()) {
    replayDir = fixtureDir;
  }
  o2::base::GeometryManager::loadGeometry(fixtureDir + "/o2sim"); // prefix, not a filename: appends _geometry-aligned.root

  // Extraction is built into an in-memory buffer first and outFile is only
  // touched after it fully succeeds, so a validation failure partway
  // through (via failExtract's exit(1)) never leaves a truncated/partial
  // JSON file on disk -- there is simply no output file at all on failure.
  std::ostringstream body;
  extractITSCommonCA(fixtureDir, replayDir, body, mcExpected);

  std::ofstream out(outFile);
  if (!out) {
    failExtract("failed to open output file for writing: " + outFile);
  }
  out << "{\n" << body.str() << "\n}\n";
  out.close();
  if (out.fail()) {
    failExtract("failed while writing output file: " + outFile);
  }
  std::cout << "[extract_metrics_its_common_ca] wrote " << outFile << std::endl;
}
