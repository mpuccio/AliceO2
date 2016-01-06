//-*- Mode: C++ -*-
// **************************************************************************
// This file is property of and copyright by the ALICE ITSU Project         *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Author: Ruben Shahoyan                                           *
//                                                                          *
// Adapted to ITSU: Maximiliano Puccio <maximiliano.puccio@cern.ch>         *
//                  for the ITS Upgrade project                             *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

#include "CATrackingStation.h"
#include <cmath>

using namespace AliceO2::ITS::CA;

TrackingStation::TrackingStation()
:mVIDOffset(0)
,mNClusters(0)
,mZMin(0)
,mZMax(0)
,mDZInv(-1)
,mDPhiInv(-1)
,mNZBins(20)
,mNPhiBins(20)
,mQueryZBmin(-1)
,mQueryZBmax(-1)
,mQueryPhiBmin(-1)
,mQueryPhiBmax(-1)
,mBins(0)
,mOccBins(0)
,mNOccBins(0)
,mNFoundClusters(0)
,mFoundClusterIterator(0)
,mFoundBinIterator(0)
,mIndex()
,mFoundBins(0)
,mSortedClInfo(0)
,mDetectors() {
  // def. c-tor
}

TrackingStation::TrackingStation(int nzbins,int nphibins)
:mVIDOffset()
,mNClusters(0)
,mZMin(0.f)
,mZMax(0.f)
,mDZInv(-1)
,mDPhiInv(-1)
,mNZBins(nzbins)
,mNPhiBins(nphibins)
,mQueryZBmin(-1)
,mQueryZBmax(-1)
,mQueryPhiBmin(-1)
,mQueryPhiBmax(-1)
,mBins(0)
,mOccBins(0)
,mNOccBins(0)
,mNFoundClusters(0)
,mFoundClusterIterator(0)
,mFoundBinIterator(0)
,mIndex()
,mFoundBins(0)
,mSortedClInfo(0)
,mDetectors() {
  // c-tor
  Init();
}

TrackingStation::~TrackingStation() {
  // d-tor
  delete[] mBins;
  delete[] mOccBins;
}

void TrackingStation::Init() {

  mClusters = *lr->GetClustersAddress();
  mZMin = lr->GetZMin();
  mZMax = lr->GetZMax();
  if (mNZBins < 1)   mNZBins = 2;
  if (mNPhiBins < 1) mNPhiBins = 1;
  mDZInv   = mNZBins / (mZMax - mZMin);
  mDPhiInv = mNPhiBins / TMath::TwoPi();
  //
  mBins = new ClBinInfo_t[mNZBins * mNPhiBins];
  mOccBins = new int[mNZBins * mNPhiBins];
  mNClusters = mClusters->GetEntriesFast();
  if(mNClusters == 0) return;
  mSortedClInfo.reserve(mClusters->GetEntriesFast());
  mVIDOffset = ((AliITSUClusterPix*)mClusters->UncheckedAt(0))->GetVolumeId();
  // prepare detectors info
  int detID = -1;
  mIndex.resize(lr->GetNSensors(),-1);
  mDetectors.reserve(lr->GetNSensors());
  for (int iCl = 0; iCl < mClusters->GetEntriesFast(); ++iCl) { //Fill this layer with detectors
    AliITSUClusterPix* c = (AliITSUClusterPix*)mClusters->UncheckedAt(iCl);
    if (detID == c->GetVolumeId()) {
      continue;
    }
    detID = c->GetVolumeId();
    ITSDetInfo_t det;
    det.index = iCl;
    //
    TGeoHMatrix m;
    geo->GetOrigMatrix(detID,m);
    //
    mIndex[detID - mVIDOffset] = mDetectors.size();
    AliITSURecoSens* sens = lr->GetSensorFromID(detID);
    const TGeoHMatrix *tm = geo->GetMatrixT2L(detID);
    m.Multiply(tm);
    double txyz[3] = {0.,0.,0.}, xyz[3] = {0.,0.,0.};
    m.LocalToMaster(txyz,xyz);
    det.xTF = sens->GetXTF(); // TMath::Sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);

    det.phiTF = sens->GetPhiTF();//TMath::ATan2(xyz[1],xyz[0]);
    det.sinTF = TMath::Sin(det.phiTF);
    det.cosTF = TMath::Cos(det.phiTF);
    //
    // compute the real radius (with misalignment)
    TGeoHMatrix mmisal(*(geo->GetMatrix(detID)));
    mmisal.Multiply(tm);
    xyz[0] = 0.;
    xyz[1] = 0.;
    xyz[2] = 0.;
    mmisal.LocalToMaster(txyz,xyz);
    det.xTFmisal = TMath::Sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
    mDetectors.push_back(det);
  } // end loop on detectors
}

void TrackingStation::SortClusters(const float vtx[3]) {
  // sort clusters and build fast lookup table
  //
  ClearSortedInfo();
  mSortedClInfo.reserve(mNClusters);
  //
  ClsInfo_t cl;
  for (int icl = mNClusters;icl--;) {
    AliITSUClusterPix* cluster = (AliITSUClusterPix*)mClusters->UncheckedAt(icl);
    cluster->GetGlobalXYZ( (float*)&cl );
    //
    if (vtx) { // phi and r will be computed wrt vertex
      cl.x -= vtx[0];
      cl.y -= vtx[1];
    }
    //
    cl.r = TMath::Sqrt(cl.x*cl.x + cl.y*cl.y);
    cl.phi = TMath::ATan2(cl.y,cl.x);
    BringTo02Pi(cl.phi);
    cl.index = icl;
    cl.zphibin = GetBinIndex(GetZBin(cl.z),GetPhiBin(cl.phi));
    cl.detid = cluster->GetVolumeId() - mVIDOffset;
    //
    mSortedClInfo.push_back(cl);
    //
  }
  sort(mSortedClInfo.begin(), fSortedClInfo.end()); // sort in phi, z
  //
  // fill cells in phi,z
  int currBin = -1;
  for (int icl = 0;icl < mNClusters; ++icl) {
    ClsInfo_t &t = mSortedClInfo[icl];
    if (t.zphibin > currBin) { // register new occupied bin
      currBin = t.zphibin;
      mBins[currBin].first = icl;
      mBins[currBin].index = mNOccBins;
      mOccBins[mNOccBins++] = currBin;
    }
    mBins[currBin].ncl++;
  }
}

void TrackingStation::Clear(Option_t *) {
  // clear cluster info
  ClearSortedInfo();
  mIndex.clear();
  mNClusters = 0;
  if (mClusters) fClusters = 0x0;
  //
}

void TrackingStation::ClearSortedInfo() {
  // clear cluster info
  mSortedClInfo.clear();
  memset(mBins,0,mNZBins * mNPhiBins * sizeof(ClBinInfo_t));
  memset(mOccBins,0,mNZBins * mNPhiBins * sizeof(int));
  mNOccBins = 0;
}

/*void TrackingStation::Print(Option_t *opt) const {
// dump cluster bins info
TString opts = opt;
opts.ToLower();
printf("Stored %d clusters in %d occupied bins\n",mNClusters,mNOccBins);
//
if (opts.Contains("c")) {
printf("\nCluster info\n");
for (int i = 0; i < mNClusters;i++) {
const ClsInfo_t &t = mSortedClInfo[i];
printf("#%5d Bin(phi/z):%03d/%03d Z:%+8.3f Phi:%+6.3f R:%7.3f Ind:%d ",
i,t.zphibin/mNZBins,t.zphibin%fNZBins,t.z,t.phi,t.r,t.index);
if (opts.Contains("l")) { // mc labels
AliITSUClusterPix* rp = (AliITSUClusterPix*)mClusters->UncheckedAt(t.index);
for (int l = 0;l < 3; l++) if (rp->GetLabel(l) >= 0) printf("| %d ",rp->GetLabel(l));
}
printf("\n");
}
}
//
if (opts.Contains("b")) {
printf("\nBins info (occupied only)\n");
for (int i=0;i<mNOccBins;i++) {
printf("%4d %5d(phi/z: %03d/%03d) -> %3d cl from %d\n",i,mOccBins[i],fOccBins[i]/mNZBins,fOccBins[i]%fNZBins,
mBins[mOccBins[i]].ncl,fBins[fOccBins[i]].first);
}
}
//
}*/

int TrackingStation::SelectClusters(float zmin,float zmax,float phimin,float phimax) {
  // prepare occupied bins in the requested region
  if (!mNOccBins) return 0;
  if (zmax < mZMin || zmin > mZMax || zmin > zmax) return 0;
  mFoundBins.clear();
  mQueryZBmin = GetZBin(zmin);
  if (mQueryZBmin < 0) fQueryZBmin = 0;
  mQueryZBmax = GetZBin(zmax);
  if (mQueryZBmax >= mNZBins) fQueryZBmax = fNZBins - 1;
  BringTo02Pi(phimin);
  BringTo02Pi(phimax);
  mQueryPhiBmin = GetPhiBin(phimin);
  mQueryPhiBmax = GetPhiBin(phimax);
  int dbz = 0;
  mNFoundClusters = 0;
  int nbcheck = mQueryPhiBmax - fQueryPhiBmin + 1; //TODO:(MP) check if a circular buffer is feasible
  if (nbcheck > 0) { // no wrapping around 0-2pi, fast case
    for (int ip = mQueryPhiBmin;ip <= fQueryPhiBmax;ip++) {
      int binID = GetBinIndex(mQueryZBmin,ip);
      if ( !(dbz = (mQueryZBmax-fQueryZBmin)) ) { // just one Z bin in the query range
        ClBinInfo_t& binInfo = mBins[binID];
        if (!binInfo.ncl) continue;
        mNFoundClusters += binInfo.ncl;
        mFoundBins.push_back(binID);
        continue;
      }
      int binMax = binID + dbz;
      for ( ; binID <= binMax; binID++) {
        ClBinInfo_t& binInfo = mBins[binID];
        if (!binInfo.ncl) continue;
        mNFoundClusters += binInfo.ncl;
        mFoundBins.push_back(binID);
      }
    }
  } else {  // wrapping
    nbcheck += mNPhiBins;
    for (int ip0 = 0;ip0 <= nbcheck;ip0++) {
      int ip = mQueryPhiBmin + ip0;
      if (ip >= mNPhiBins) ip -= fNPhiBins;
      int binID = GetBinIndex(mQueryZBmin,ip);
      if ( !(dbz = (mQueryZBmax - fQueryZBmin)) ) { // just one Z bin in the query range
        ClBinInfo_t& binInfo = mBins[binID];
        if (!binInfo.ncl) continue;
        mNFoundClusters += binInfo.ncl;
        mFoundBins.push_back(binID);
        continue;
      }
      int binMax = binID + dbz;
      for (;binID <= binMax;binID++) {
        ClBinInfo_t& binInfo = mBins[binID];
        if (!binInfo.ncl) continue;
        mNFoundClusters += binInfo.ncl;
        mFoundBins.push_back(binID);
      }
    }
  }
  mFoundClusterIterator = mFoundBinIterator = 0;
  return mNFoundClusters;
}

int TrackingStation::GetNextClusterInfoID() {
  if (mFoundBinIterator < 0) return 0;
  int currBin = mFoundBins[mFoundBinIterator];
  if (mFoundClusterIterator < mBins[currBin].ncl) { // same bin
    return mBins[currBin].first + mFoundClusterIterator++;
  }
  if (++mFoundBinIterator < int(mFoundBins.size())) {  // need to change bin
    currBin = mFoundBins[mFoundBinIterator];
    mFoundClusterIterator = 1;
    return mBins[currBin].first;
  }
  mFoundBinIterator = -1;
  return -1;
}

void TrackingStation::ResetFoundIterator() {
  // prepare for a new loop over found clusters
  if (mNFoundClusters)  mFoundClusterIterator = mFoundBinIterator = 0;
}

