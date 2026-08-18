// Publication/CMake correction handoff: focused output checker proving the
// ITS common-CA workflow's identity contract for every emitted flattened
// cluster reference (o2-its-ca-tracker-workflow's o2trac_its_ca.root, from
// replay_tracking_its_common_ca.sh). This is a SIBLING to
// extract_metrics_its_common_ca.C (gate3-slice3-its-ca-validation/), not a
// replacement or modification of it: no file under gate3-slice3-its-ca-
// validation or gate0-baseline is read or written by this macro, and no
// historical metric recorded there is altered.
//
// extract_metrics_its_common_ca.C's trackContentHash intentionally never
// covered cluster identity/size/ROF content (see its own
// trackContentHashDefinition string) -- it hashes only
// (nClusters,chi2,x,alpha,y,z,snp,tgl,q2pt), the o2::track::TrackParametrization
// state. This macro exists specifically to validate the content that hash
// excludes: for every (track, layer) cluster attachment, independently
// derived ground truth for
//   - the correct hit layer (from the referenced cluster's own sensorID via
//     GeometryTGeo::getLayer(), cross-checked against
//     TrackITS::hasHitOnLayer());
//   - the external cluster index itself, bounds-checked against the input
//     cluster array;
//   - the published cluster size (TrackITS::getClusterSize(layer)),
//     compared against an independently recomputed pixel count from the raw
//     cluster pattern + CCDB TopologyDictionary, mirroring exactly the
//     dict->isGroup(pattID) ? ClusterPattern(pattIt).getNPixels()
//     : dict->getNpixels(pattID) convention TimeFrame::loadNormalizedSource()
//     itself uses (ITSMFTTracking/TimeFrame.cxx) -- never re-deriving size
//     from the track/TimeFrame's own bookkeeping, which is exactly what the
//     local-vs-external cluster-size index-domain bug corrupted;
//   - the associated ROF, cross-checked between the cluster's input ROF
//     (ITSClustersROF) and the track's output ROF (ITSTrackROF), which are
//     the same vector by construction (CATrackerSpec.cxx's fillITSOutputs:
//     trackROFs.assign(inputROFs.begin(), inputROFs.end())), so index i in
//     one must equal index i in the other for a given track/cluster pair.
//
// Because the ground-truth size is recomputed independently from the raw
// cluster file + dictionary (never from TrackITS::getClusterSize() itself),
// this check would have failed loudly under the old external-index lookup bug: for every
// track cluster NOT on layer 0, the published size would have come from an
// unrelated slot of layer 0's own (much smaller) per-layer vector instead
// of the correct layer's, so sizeMismatch (and, depending on layer 0's
// vector length that run, an out-of-bounds vector access) would fire for
// nearly every reference.
//
// Usage:
//   root -l -b -q 'check_cluster_attachment_its_common_ca.C("<fixtureDir>", "<replayDir>", "<out.json>", <timestamp>)'
//
// <fixtureDir> must contain o2clus_its.root, o2sim_geometry-aligned.root.
// <replayDir> must contain o2trac_its_ca.root. <timestamp> must be the
// ACTUAL resolved CCDB condition timestamp the replay itself used for
// ITS/Calib/ClusterDictionary -- read it from the replay's own log (grep
// for "ccdb reads .*ClusterDictionary.*for <timestamp>"), NOT the
// --condition-not-after CLI value passed to the workflow: DPL resolves the
// real per-run condition timestamp internally (from GRP/run metadata), and
// --condition-not-after is only an upper-bound ceiling on it. Passing the
// wrong timestamp here silently fetches a different (but still valid-
// looking) TopologyDictionary object and produces spurious sizeMismatch
// output that has nothing to do with the publication code under test.
//
// Exits non-zero (after writing the full JSON report) if ANY reference
// fails any of the four checks above -- this is a gate, not just a report.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TMD5.h>
#include <TTree.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsBase/GeometryManager.h"
#include "ITSBase/GeometryTGeo.h"
#endif

namespace
{
[[noreturn]] void failCheck(const std::string& msg)
{
  std::cerr << "[check_cluster_attachment_its_common_ca] " << msg << std::endl;
  std::exit(1);
}

void requireOpen(TFile& f, const std::string& path)
{
  if (f.IsZombie() || !f.IsOpen()) {
    failCheck("failed to open ROOT file: " + path);
  }
}

TTree* requireTree(TFile& f, const char* name, const std::string& fileDesc)
{
  auto* t = dynamic_cast<TTree*>(f.Get(name));
  if (!t) {
    failCheck("missing tree '" + std::string(name) + "' in " + fileDesc);
  }
  return t;
}

template <typename T>
void requireBranch(TTree* t, const char* name, T*& ptr, const std::string& treeDesc)
{
  if (!t->GetBranch(name)) {
    failCheck("missing branch '" + std::string(name) + "' in " + treeDesc);
  }
  if (t->SetBranchAddress(name, &ptr) < 0) {
    failCheck("SetBranchAddress failed for branch '" + std::string(name) + "' in " + treeDesc);
  }
}

void requireEntry(TTree* t, Long64_t entry, const std::string& treeDesc)
{
  if (entry >= t->GetEntries()) {
    failCheck(treeDesc + " has no entry " + std::to_string(entry));
  }
  if (t->GetEntry(entry) <= 0) {
    failCheck("GetEntry(" + std::to_string(entry) + ") failed for " + treeDesc);
  }
}

// Find the ROF index whose [firstEntry, firstEntry+nEntries) range contains
// `pos`, or -1 if none does.
int findROF(const std::vector<o2::itsmft::ROFRecord>& rofs, int pos)
{
  for (int i = 0; i < (int)rofs.size(); ++i) {
    const auto first = rofs[i].getFirstEntry();
    const auto n = rofs[i].getNEntries();
    if (pos >= first && pos < first + n) {
      return i;
    }
  }
  return -1;
}
} // namespace

void check_cluster_attachment_its_common_ca(std::string fixtureDir, std::string replayDir, std::string outFile, long timestamp)
{
  using namespace o2::its;

  o2::base::GeometryManager::loadGeometry(fixtureDir + "/o2sim");
  auto gman = GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::TransformType::T2L);

  auto& mgr = o2::ccdb::BasicCCDBManager::instance();
  mgr.setURL("http://alice-ccdb.cern.ch");
  mgr.setTimestamp(timestamp);
  const auto* dict = mgr.get<o2::itsmft::TopologyDictionary>("ITS/Calib/ClusterDictionary");
  if (!dict) {
    failCheck("failed to fetch ITS/Calib/ClusterDictionary from CCDB at timestamp " + std::to_string(timestamp));
  }

  const std::string clusPath = fixtureDir + "/o2clus_its.root";
  TFile fClus(clusPath.c_str());
  requireOpen(fClus, clusPath);
  auto clusTree = requireTree(fClus, "o2sim", clusPath);
  std::vector<o2::itsmft::CompClusterExt>* clusArr = nullptr;
  std::vector<o2::itsmft::ROFRecord>* clusROF = nullptr;
  std::vector<unsigned char>* patternsPtr = nullptr;
  requireBranch(clusTree, "ITSClusterComp", clusArr, clusPath + ":o2sim");
  requireBranch(clusTree, "ITSClustersROF", clusROF, clusPath + ":o2sim");
  requireBranch(clusTree, "ITSClusterPatt", patternsPtr, clusPath + ":o2sim");
  requireEntry(clusTree, 0, clusPath + ":o2sim");
  if (!clusArr || !clusROF || !patternsPtr) {
    failCheck("null branch object after read in " + clusPath + ":o2sim");
  }
  const int nClusters = (int)clusArr->size();

  // Ground-truth pixel count per raw cluster row, decoded exactly as
  // TimeFrame::loadNormalizedSource() computes m.shape.nPixels (identical
  // dict->isGroup(...) branch to IOUtils.h's extractClusterDataBounded()):
  // patterns are consumed sequentially from the shared byte stream only for
  // rows that need it (InvalidPatternID or a grouped pattern ID), so this
  // must iterate every cluster row in file order to keep the iterator in
  // sync, even though only some rows are ever attached to a track.
  //
  // Clamped to 15, not 255: o2::its::TrackITS::setClusterSize() packs the
  // published size into a 4-bit field per layer (DataFormatsITS/TrackITS.h),
  // so any real cluster wider than 15 pixels is clamped by production
  // before publication -- the ground truth here must apply the identical
  // clamp or a handful of genuinely->15px clusters would show as spurious
  // mismatches.
  std::vector<uint8_t> groundTruthSize(nClusters, 0);
  {
    auto pattIt = patternsPtr->cbegin();
    for (int i = 0; i < nClusters; ++i) {
      const auto& c = (*clusArr)[i];
      const auto pattID = c.getPatternID();
      unsigned int nPixels = 0;
      if (pattID != o2::itsmft::CompCluster::InvalidPatternID && !dict->isGroup(pattID)) {
        nPixels = dict->getNpixels(pattID);
      } else {
        o2::itsmft::ClusterPattern patt(pattIt);
        nPixels = (unsigned int)patt.getNPixels();
      }
      groundTruthSize[i] = (uint8_t)std::clamp(nPixels, 0u, 15u);
    }
  }

  const std::string tracPath = replayDir + "/o2trac_its_ca.root";
  TFile fTrac(tracPath.c_str());
  requireOpen(fTrac, tracPath);
  auto trTree = requireTree(fTrac, "o2sim", tracPath);
  std::vector<o2::its::TrackITS>* trkArr = nullptr;
  std::vector<int>* clsIdxArr = nullptr;
  std::vector<o2::itsmft::ROFRecord>* trkROF = nullptr;
  requireBranch(trTree, "ITSTrack", trkArr, tracPath + ":o2sim");
  requireBranch(trTree, "ITSTrackClusIdx", clsIdxArr, tracPath + ":o2sim");
  requireBranch(trTree, "ITSTracksROF", trkROF, tracPath + ":o2sim");
  requireEntry(trTree, 0, tracPath + ":o2sim");
  if (!trkArr || !clsIdxArr || !trkROF) {
    failCheck("null branch object after read in " + tracPath + ":o2sim");
  }

  long long nRefs = 0, nBoundsFail = 0, nLayerMismatch = 0, nPatternBitMissing = 0,
            nSizeZero = 0, nSizeMismatch = 0, nRofMismatch = 0, nRofUnresolved = 0;
  TMD5 digestMD5;
  char buf[128];

  for (int iTrk = 0; iTrk < (int)trkArr->size(); ++iTrk) {
    auto& t = (*trkArr)[iTrk];
    const int outputRofIdx = findROF(*trkROF, iTrk);
    if (outputRofIdx < 0) {
      failCheck("track " + std::to_string(iTrk) + " is not covered by any ITSTracksROF entry");
    }
    const int first = t.getFirstClusterEntry();
    const int n = t.getNumberOfClusters();
    for (int k = 0; k < n; ++k) {
      const int pos = first + k;
      if (pos < 0 || pos >= (int)clsIdxArr->size()) {
        failCheck("track " + std::to_string(iTrk) + " cluster range index " + std::to_string(pos) +
                   " is out of bounds of ITSTrackClusIdx (size " + std::to_string(clsIdxArr->size()) + ")");
      }
      const int extIdx = (*clsIdxArr)[pos];
      ++nRefs;

      const bool boundsOK = extIdx >= 0 && extIdx < nClusters;
      if (!boundsOK) {
        ++nBoundsFail;
        std::cerr << "[check_cluster_attachment_its_common_ca] track " << iTrk
                  << " pos " << pos << ": external index " << extIdx
                  << " out of bounds [0," << nClusters << ")" << std::endl;
        continue; // cannot look up ground truth for an out-of-bounds index
      }

      const int layer = gman->getLayer((*clusArr)[extIdx].getSensorID());
      const bool layerOK = layer >= 0 && layer < 7;
      const bool patternBitOK = layerOK && t.hasHitOnLayer(layer);
      if (!layerOK) {
        ++nLayerMismatch;
      }
      if (layerOK && !patternBitOK) {
        ++nPatternBitMissing;
      }

      const int publishedSize = layerOK ? t.getClusterSize(layer) : 0;
      const uint8_t expectedSize = groundTruthSize[extIdx];
      if (publishedSize == 0) {
        ++nSizeZero;
      }
      if (!layerOK || publishedSize != (int)expectedSize) {
        ++nSizeMismatch;
        std::cerr << "[check_cluster_attachment_its_common_ca] track " << iTrk
                  << " pos " << pos << " extIdx " << extIdx << " layer " << layer
                  << ": published size " << publishedSize << " != expected " << (int)expectedSize
                  << std::endl;
      }

      const int inputRofIdx = findROF(*clusROF, extIdx);
      if (inputRofIdx < 0) {
        ++nRofUnresolved;
      } else if (inputRofIdx != outputRofIdx) {
        ++nRofMismatch;
        std::cerr << "[check_cluster_attachment_its_common_ca] track " << iTrk
                  << " (output ROF " << outputRofIdx << ") pos " << pos << " extIdx " << extIdx
                  << ": input ROF " << inputRofIdx << " != output ROF " << outputRofIdx << std::endl;
      }

      const int len = snprintf(buf, sizeof(buf), "%d,%d,%d,%d,%d;", iTrk, layer, extIdx, publishedSize, outputRofIdx);
      digestMD5.Update((UChar_t*)buf, len);
    }
  }
  digestMD5.Final();
  const std::string expandedDigest = digestMD5.AsString();

  const bool allOK = nBoundsFail == 0 && nLayerMismatch == 0 && nPatternBitMissing == 0 &&
                      nSizeZero == 0 && nSizeMismatch == 0 && nRofMismatch == 0 && nRofUnresolved == 0;

  std::ofstream out(outFile);
  if (!out) {
    failCheck("failed to open output file for writing: " + outFile);
  }
  out << "{\n"
      << "  \"itsCommonCAClusterAttachment\": {\n"
      << "    \"tracks\": " << trkArr->size() << ",\n"
      << "    \"clusterReferences\": " << nRefs << ",\n"
      << "    \"expandedDigest\": \"" << expandedDigest << "\",\n"
      << "    \"expandedDigestDefinition\": \"MD5 over ordered (trackIdx,layer,externalClusterIndex,publishedClusterSize,outputRofIdx) per emitted cluster reference, in (track,position) order -- covers exactly the identity/size/ROF content extract_metrics_its_common_ca.C's trackContentHash excludes\",\n"
      << "    \"boundsFail\": " << nBoundsFail << ",\n"
      << "    \"layerMismatch\": " << nLayerMismatch << ",\n"
      << "    \"patternBitMissing\": " << nPatternBitMissing << ",\n"
      << "    \"sizeZero\": " << nSizeZero << ",\n"
      << "    \"sizeMismatch\": " << nSizeMismatch << ",\n"
      << "    \"rofMismatch\": " << nRofMismatch << ",\n"
      << "    \"rofUnresolved\": " << nRofUnresolved << ",\n"
      << "    \"allOK\": " << (allOK ? "true" : "false") << "\n"
      << "  }\n"
      << "}\n";
  out.close();
  if (out.fail()) {
    failCheck("failed while writing output file: " + outFile);
  }
  std::cout << "[check_cluster_attachment_its_common_ca] wrote " << outFile
            << " (" << nRefs << " cluster references, allOK=" << (allOK ? "true" : "false") << ")" << std::endl;
  if (!allOK) {
    std::exit(1);
  }
}
