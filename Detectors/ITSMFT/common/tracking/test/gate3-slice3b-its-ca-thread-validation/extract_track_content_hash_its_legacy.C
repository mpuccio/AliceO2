// Gate 3 Slice 3b: semantic track-content digest for the legacy ITS labelled
// reference (replay_tracking_its_legacy_sync.sh's o2trac_its.root).
//
// gate0-baseline/extract_metrics.C's extractITS() (reused unmodified via
// gate3-slice3-its-ca-validation/extract_metrics_its_legacy_diamond.sh for
// this leg's track-count/efficiency/fake-rate/etc. metrics) does not compute
// a track-content hash. This macro adds ONLY that digest, using the exact
// same definition as
// gate3-slice3-its-ca-validation/extract_metrics_its_common_ca.C's
// trackContentHash (MD5 over ordered per-track
// (nClusters,chi2,x,alpha,y,z,snp,tgl,q2pt), %.9g-formatted) -- so the two
// legs' hashes are computed identically, even though (per task scope) no
// equivalence or parity claim is made between their values: the legacy leg
// exercises the real per-event Sync vertexer, the common-CA leg a static
// diamond, and different tracker implementations are involved. This is
// provided purely so a "semantic digest" exists for the legacy leg per task
// instruction #5, and so the SAME digest definition can be used to prove
// (or refute) determinism of repeated legacy replays if ever re-run.
//
// Neither gate0-baseline nor gate3-slice3-its-ca-validation is read, written,
// or modified by this macro -- it opens replayDir's own o2trac_its.root
// directly. No fixture/geometry access is needed (the digest depends only on
// the track branch itself), unlike extractITS()/extractITSCommonCA() which
// also need geometry+MC for efficiency/fake-rate/denominator computation.
//
// Usage:
//   root -l -b -q 'extract_track_content_hash_its_legacy.C("<replayDir>", "<out.json>")'

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TMD5.h>
#include <TTree.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "DataFormatsITS/TrackITS.h"
#endif

namespace
{
[[noreturn]] void failExtract(const std::string& msg)
{
  std::cerr << "[extract_track_content_hash_its_legacy] " << msg << std::endl;
  std::exit(1);
}
} // namespace

void extract_track_content_hash_its_legacy(std::string replayDir, std::string outFile)
{
  const std::string tracPath = replayDir + "/o2trac_its.root";
  TFile fTrac(tracPath.c_str());
  if (fTrac.IsZombie() || !fTrac.IsOpen()) {
    failExtract("failed to open ROOT file: " + tracPath);
  }
  auto* trTree = dynamic_cast<TTree*>(fTrac.Get("o2sim"));
  if (!trTree) {
    failExtract("missing tree 'o2sim' in " + tracPath);
  }
  std::vector<o2::its::TrackITS>* trkArr = nullptr;
  if (!trTree->GetBranch("ITSTrack")) {
    failExtract("missing branch 'ITSTrack' in " + tracPath + ":o2sim");
  }
  if (trTree->SetBranchAddress("ITSTrack", &trkArr) < 0) {
    failExtract("SetBranchAddress failed for branch 'ITSTrack' in " + tracPath + ":o2sim");
  }
  if (trTree->GetEntries() < 1) {
    failExtract(tracPath + ":o2sim has no entries");
  }
  if (trTree->GetEntry(0) <= 0) {
    failExtract("GetEntry(0) failed for " + tracPath + ":o2sim");
  }
  if (!trkArr) {
    failExtract("null branch object after read: " + tracPath + ":o2sim:ITSTrack");
  }

  // Identical field list/formatting/order to extract_metrics_its_common_ca.C's
  // trackContentHash -- see that file's header comment for why these nine
  // fields (not the full set of writer branches) were chosen.
  TMD5 contentMD5;
  char buf[320];
  int nTracks = (int)trkArr->size();
  for (auto& t : *trkArr) {
    int len = snprintf(buf, sizeof(buf), "%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g;",
                        t.getNumberOfClusters(), t.getChi2(), t.getX(), t.getAlpha(), t.getY(), t.getZ(),
                        t.getSnp(), t.getTgl(), t.getQ2Pt());
    contentMD5.Update((UChar_t*)buf, len);
  }
  contentMD5.Final();
  std::string trackContentHash = contentMD5.AsString();

  std::ofstream out(outFile);
  if (!out) {
    failExtract("failed to open output file for writing: " + outFile);
  }
  out << "{\n"
      << "  \"outputTracks\": " << nTracks << ",\n"
      << "  \"trackContentHash\": \"" << trackContentHash << "\",\n"
      << "  \"trackContentHashDefinition\": \"identical definition to gate3-slice3-its-ca-validation/extract_metrics_its_common_ca.C's trackContentHash, applied to this leg's o2trac_its.root instead\"\n"
      << "}\n";
  out.close();
  if (out.fail()) {
    failExtract("failed while writing output file: " + outFile);
  }
  std::cout << "[extract_track_content_hash_its_legacy] wrote " << outFile << std::endl;
}
