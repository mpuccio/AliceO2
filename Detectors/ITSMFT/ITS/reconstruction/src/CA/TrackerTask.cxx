// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CA/TrackerTask.cxx
/// \brief Implementation of the ITS Cellular Automaton tracker task
/// \author maximiliano.puccio@cern.ch

#include "ITSReconstruction/CA/TrackerTask.h"
#include "ITSReconstruction/Cluster.h"

#include "FairLogger.h"      // for LOG
#include "FairRootManager.h" // for FairRootManager
#include "TClonesArray.h"    // for TClonesArray

ClassImp(o2::ITS::CA::TrackerTask)

using namespace o2::ITS::CA;

//_____________________________________________________________________
TrackerTask::TrackerTask() :
  FairTask{"ITS::CA::TrackerTask"},
  mGeometry{},
  mEvent{0,&mGeometry},
  mTracker{mEvent},
  mClustersArray{nullptr},
  mTracksArray{nullptr}
  {}

//_____________________________________________________________________
TrackerTask::~TrackerTask()
{
  delete mTracksArray;
}

//_____________________________________________________________________
/// \brief Init function
/// Inititializes the tracker and connects input and output container
InitStatus TrackerTask::Init()
{
  FairRootManager* mgr = FairRootManager::Instance();
  if (!mgr) {
    LOG(ERROR) << "Could not instantiate FairRootManager. Exiting ..." << FairLogger::endl;
    return kERROR;
  }

  mClustersArray = dynamic_cast<const TClonesArray*>(mgr->GetObject("ITSCluster"));
  if (!mClustersArray) {
    LOG(ERROR) << "ITS clusters not registered in the FairRootManager. Exiting ..." << FairLogger::endl;
    return kERROR;
  }

  // Register output container
  mTracksArray = new TClonesArray("o2::ITS::CookedTrack");
  mgr->Register("ITSTrack", "ITS", mTracksArray, kTRUE);

  mGeometry.Build(kTRUE);
  ITS::Cluster::setGeom(&mGeometry);

  return kSUCCESS;
}

//_____________________________________________________________________
void TrackerTask::Exec(Option_t* option)
{
  mTracksArray->Clear();
  LOG(DEBUG) << "Running "<< this->GetName() <<" on new event" << FairLogger::endl;

  mTracker.clustersToTracks();
}
