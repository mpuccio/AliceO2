// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TrackerTask.h
/// \brief Definition of the ITS Cellular Automaton tracker task
/// \author maximiliano.puccio@cern.ch

#ifndef O2_ITSMFT_RECONSTRUCTION_CA_TRACKERTASK_H_
#define O2_ITSMFT_RECONSTRUCTION_CA_TRACKERTASK_H_

#include "FairTask.h" 

#include "ITSBase/GeometryTGeo.h"
#include "ITSReconstruction/CA/Event.h"
#include "ITSReconstruction/CA/Tracker.h"

class TClonesArray;

namespace o2
{
namespace ITS
{
namespace CA
{
class TrackerTask : public FairTask
{
 public:
  TrackerTask();
  ~TrackerTask() override;

  InitStatus Init() override;
  void Exec(Option_t* option) override;

 private:
  GeometryTGeo mGeometry; ///< ITS geometry
  Event        mEvent;    ///< CA tracker event
  Tracker      mTracker;  ///< Track finder

  const TClonesArray* mClustersArray;   ///< Array of clusters
  TClonesArray*       mTracksArray;     ///< Array of tracks

  ClassDefOverride(TrackerTask, 1)
};
}
}
}
#endif /* O2_ITSMFT_RECONSTRUCTION_CA_TRACKERTASK_H_ */
