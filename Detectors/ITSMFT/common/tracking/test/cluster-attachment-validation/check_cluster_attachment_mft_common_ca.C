// Publication/CMake correction handoff: focused output checker proving the
// MFT common-CA workflow's identity contract for every emitted flattened
// cluster reference (o2-mft-ca-tracker-workflow's mfttracks.root, from
// gate0-baseline/replay_tracking_common_ca.sh). Sibling to
// check_cluster_attachment_its_common_ca.C (same directory) and to
// gate0-baseline/extract_metrics_common_ca.C -- no file under gate0-baseline
// is read or written by this macro, and no historical metric recorded there
// is altered.
//
// extract_metrics_common_ca.C's trackContentHash covers only
// (nPoints,chi2,x,y,z,phi,tanl,invQPt) (see its own trackContentHashDefinition
// string) -- it never covers cluster identity/size/ROF content. This macro
// validates exactly that excluded content: for every (track, cluster
// reference), independently derived ground truth for
//   - the correct hit layer: o2::mft::TrackMFT carries no per-track
//     pattern/hit-layer bitmask once published (unlike o2::its::TrackITS),
//     so the layer is derived solely from the referenced cluster's own
//     sensorID via MFT GeometryTGeo::getLayer() -- and self-consistency is
//     checked by requiring each track never attaches two references to the
//     same layer;
//   - the external cluster index itself, bounds-checked against the input
//     cluster array;
//   - the published cluster size, decoded from TrackMFT::getClusterSizes()'s
//     packed 6-bits-per-layer field at the geometry-derived layer, compared
//     against an independently recomputed pixel count from the raw cluster
//     pattern + CCDB TopologyDictionary (identical dict->isGroup(pattID)
//     convention to the ITS checker and to
//     ITSMFTTracking/TimeFrame.cxx's loadNormalizedSource()) -- never
//     re-deriving size from the track/TimeFrame's own bookkeeping, which is
//     exactly what the local-vs-external cluster-size index-domain bug
//     corrupted (MFTFwdTrackHelpers.cxx's refitTrackFwdImpl and
//     MFT/workflow/CATrackerSpec.cxx's fillMFTOutputs);
//   - the associated ROF, cross-checked between the cluster's input ROF
//     (MFTClustersROF) and the track's output ROF (MFTTracksROF), which are
//     the same vector by construction (CATrackerSpec.cxx's fillMFTOutputs:
//     trackROFs.assign(inputROFs.begin(), inputROFs.end())).
//
// Usage:
//   root -l -b -q 'check_cluster_attachment_mft_common_ca.C("<fixtureDir>", "<replayDir>", "<out.json>", <timestamp>)'
//
// <fixtureDir> must contain mftclusters.root, o2sim_geometry-aligned.root.
// <replayDir> must contain mfttracks.root. <timestamp> must be the ACTUAL
// resolved CCDB condition timestamp the replay itself used for
// MFT/Calib/ClusterDictionary -- read it from the replay's own log (grep
// for "ccdb reads .*ClusterDictionary.*for <timestamp>"), NOT the
// --condition-not-after CLI value: see check_cluster_attachment_its_common_ca.C's
// header comment for why those two values differ.
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
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DetectorsBase/GeometryManager.h"
#include "MFTBase/GeometryTGeo.h"
#endif

namespace
{
constexpr int MFTNLayers = 10;

[[noreturn]] void failCheck(const std::string& msg)
{
  std::cerr << "[check_cluster_attachment_mft_common_ca] " << msg << std::endl;
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

// TrackMFT has no public getClusterSize(layer) accessor (unlike TrackITS) --
// only the raw packed getClusterSizes() and setClusterSize(). Decode the
// same 6-bits-per-layer packing setClusterSize() itself uses.
int decodeClusterSize(uint64_t packed, int layer)
{
  return (int)((packed >> (layer * 6)) & 0x3fULL);
}
} // namespace

void check_cluster_attachment_mft_common_ca(std::string fixtureDir, std::string replayDir, std::string outFile, long timestamp)
{
  o2::base::GeometryManager::loadGeometry(fixtureDir + "/o2sim");
  auto gman = o2::mft::GeometryTGeo::Instance();

  auto& mgr = o2::ccdb::BasicCCDBManager::instance();
  mgr.setURL("http://alice-ccdb.cern.ch");
  mgr.setTimestamp(timestamp);
  const auto* dict = mgr.get<o2::itsmft::TopologyDictionary>("MFT/Calib/ClusterDictionary");
  if (!dict) {
    failCheck("failed to fetch MFT/Calib/ClusterDictionary from CCDB at timestamp " + std::to_string(timestamp));
  }

  const std::string clusPath = fixtureDir + "/mftclusters.root";
  TFile fClus(clusPath.c_str());
  requireOpen(fClus, clusPath);
  auto clusTree = requireTree(fClus, "o2sim", clusPath);
  std::vector<o2::itsmft::CompClusterExt>* clusArr = nullptr;
  std::vector<o2::itsmft::ROFRecord>* clusROF = nullptr;
  std::vector<unsigned char>* patternsPtr = nullptr;
  requireBranch(clusTree, "MFTClusterComp", clusArr, clusPath + ":o2sim");
  requireBranch(clusTree, "MFTClustersROF", clusROF, clusPath + ":o2sim");
  requireBranch(clusTree, "MFTClusterPatt", patternsPtr, clusPath + ":o2sim");
  requireEntry(clusTree, 0, clusPath + ":o2sim");
  if (!clusArr || !clusROF || !patternsPtr) {
    failCheck("null branch object after read in " + clusPath + ":o2sim");
  }
  const int nClusters = (int)clusArr->size();

  // Clamped to 63, not 255: o2::mft::TrackMFT::setClusterSize() (and
  // MFTCATrack::setClusterSize()) pack the published size into a 6-bit
  // field per layer, so ground truth must apply the identical clamp.
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
      groundTruthSize[i] = (uint8_t)std::clamp(nPixels, 0u, 63u);
    }
  }

  const std::string tracPath = replayDir + "/mfttracks.root";
  TFile fTrac(tracPath.c_str());
  requireOpen(fTrac, tracPath);
  auto trTree = requireTree(fTrac, "o2sim", tracPath);
  std::vector<o2::mft::TrackMFT>* trkArr = nullptr;
  std::vector<int>* clsIdxArr = nullptr;
  std::vector<o2::itsmft::ROFRecord>* trkROF = nullptr;
  requireBranch(trTree, "MFTTrack", trkArr, tracPath + ":o2sim");
  requireBranch(trTree, "MFTTrackClusIdx", clsIdxArr, tracPath + ":o2sim");
  requireBranch(trTree, "MFTTracksROF", trkROF, tracPath + ":o2sim");
  requireEntry(trTree, 0, tracPath + ":o2sim");
  if (!trkArr || !clsIdxArr || !trkROF) {
    failCheck("null branch object after read in " + tracPath + ":o2sim");
  }

  long long nRefs = 0, nBoundsFail = 0, nLayerMismatch = 0, nDuplicateLayer = 0,
            nSizeZero = 0, nSizeMismatch = 0, nRofMismatch = 0, nRofUnresolved = 0;
  TMD5 digestMD5;
  char buf[128];

  for (int iTrk = 0; iTrk < (int)trkArr->size(); ++iTrk) {
    auto& t = (*trkArr)[iTrk];
    const int outputRofIdx = findROF(*trkROF, iTrk);
    if (outputRofIdx < 0) {
      failCheck("track " + std::to_string(iTrk) + " is not covered by any MFTTracksROF entry");
    }
    const int first = t.getExternalClusterIndexOffset();
    const int n = t.getNumberOfPoints();
    const uint64_t packedSizes = t.getClusterSizes();
    std::set<int> seenLayers;
    for (int k = 0; k < n; ++k) {
      const int pos = first + k;
      if (pos < 0 || pos >= (int)clsIdxArr->size()) {
        failCheck("track " + std::to_string(iTrk) + " cluster range index " + std::to_string(pos) +
                   " is out of bounds of MFTTrackClusIdx (size " + std::to_string(clsIdxArr->size()) + ")");
      }
      const int extIdx = (*clsIdxArr)[pos];
      ++nRefs;

      const bool boundsOK = extIdx >= 0 && extIdx < nClusters;
      if (!boundsOK) {
        ++nBoundsFail;
        std::cerr << "[check_cluster_attachment_mft_common_ca] track " << iTrk
                  << " pos " << pos << ": external index " << extIdx
                  << " out of bounds [0," << nClusters << ")" << std::endl;
        continue;
      }

      const int layer = gman->getLayer((*clusArr)[extIdx].getSensorID());
      const bool layerOK = layer >= 0 && layer < MFTNLayers;
      if (!layerOK) {
        ++nLayerMismatch;
      } else if (!seenLayers.insert(layer).second) {
        ++nDuplicateLayer;
        std::cerr << "[check_cluster_attachment_mft_common_ca] track " << iTrk
                  << ": layer " << layer << " attached more than once" << std::endl;
      }

      const int publishedSize = layerOK ? decodeClusterSize(packedSizes, layer) : 0;
      const uint8_t expectedSize = groundTruthSize[extIdx];
      if (publishedSize == 0) {
        ++nSizeZero;
      }
      if (!layerOK || publishedSize != (int)expectedSize) {
        ++nSizeMismatch;
        std::cerr << "[check_cluster_attachment_mft_common_ca] track " << iTrk
                  << " pos " << pos << " extIdx " << extIdx << " layer " << layer
                  << ": published size " << publishedSize << " != expected " << (int)expectedSize
                  << std::endl;
      }

      const int inputRofIdx = findROF(*clusROF, extIdx);
      if (inputRofIdx < 0) {
        ++nRofUnresolved;
      } else if (inputRofIdx != outputRofIdx) {
        ++nRofMismatch;
        std::cerr << "[check_cluster_attachment_mft_common_ca] track " << iTrk
                  << " (output ROF " << outputRofIdx << ") pos " << pos << " extIdx " << extIdx
                  << ": input ROF " << inputRofIdx << " != output ROF " << outputRofIdx << std::endl;
      }

      const int len = snprintf(buf, sizeof(buf), "%d,%d,%d,%d,%d;", iTrk, layer, extIdx, publishedSize, outputRofIdx);
      digestMD5.Update((UChar_t*)buf, len);
    }
  }
  digestMD5.Final();
  const std::string expandedDigest = digestMD5.AsString();

  const bool allOK = nBoundsFail == 0 && nLayerMismatch == 0 && nDuplicateLayer == 0 &&
                      nSizeZero == 0 && nSizeMismatch == 0 && nRofMismatch == 0 && nRofUnresolved == 0;

  std::ofstream out(outFile);
  if (!out) {
    failCheck("failed to open output file for writing: " + outFile);
  }
  out << "{\n"
      << "  \"mftCommonCAClusterAttachment\": {\n"
      << "    \"tracks\": " << trkArr->size() << ",\n"
      << "    \"clusterReferences\": " << nRefs << ",\n"
      << "    \"expandedDigest\": \"" << expandedDigest << "\",\n"
      << "    \"expandedDigestDefinition\": \"MD5 over ordered (trackIdx,layer,externalClusterIndex,publishedClusterSize,outputRofIdx) per emitted cluster reference, in (track,position) order -- covers exactly the identity/size/ROF content extract_metrics_common_ca.C's trackContentHash excludes\",\n"
      << "    \"boundsFail\": " << nBoundsFail << ",\n"
      << "    \"layerMismatch\": " << nLayerMismatch << ",\n"
      << "    \"duplicateLayer\": " << nDuplicateLayer << ",\n"
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
  std::cout << "[check_cluster_attachment_mft_common_ca] wrote " << outFile
            << " (" << nRefs << " cluster references, allOK=" << (allOK ? "true" : "false") << ")" << std::endl;
  if (!allOK) {
    std::exit(1);
  }
}
