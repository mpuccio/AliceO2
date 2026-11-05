// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file testMCTruth.cxx
/// \brief The sink's MC label rules.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK MCTruth
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Study/MCTruth.h"

#include <vector>

namespace o2::gtrk::study
{

BOOST_AUTO_TEST_CASE(AParticleIsAPrimaryADecaysOrTheMaterials)
{
  using study::ParticleOrigin;
  // The generator's: no mother, produced in the primary interaction.
  const o2::MCTrack generated;
  BOOST_TEST(generated.isPrimary());
  BOOST_TEST((particleOrigin(generated, true) == ParticleOrigin::kPhysicalPrimary));
  // Feed-down the generator decayed.
  BOOST_TEST((particleOrigin(generated, false) == ParticleOrigin::kDecay));

  o2::MCTrack decayed;
  decayed.SetMotherTrackId(0);
  decayed.setProcess(kPDecay);
  BOOST_TEST((particleOrigin(decayed, false) == ParticleOrigin::kDecay));
  // A heavy-flavour daughter of transport is a physical primary.
  BOOST_TEST((particleOrigin(decayed, true) == ParticleOrigin::kPhysicalPrimary));

  for (const auto process : {kPHadronic, kPPair, kPCompton, kPDeltaRay}) {
    o2::MCTrack secondary;
    secondary.SetMotherTrackId(0);
    secondary.setProcess(process);
    BOOST_TEST((particleOrigin(secondary, false) == ParticleOrigin::kMaterial));
  }
}

BOOST_AUTO_TEST_CASE(AnOrphanOfTransportIsNoPhysicalPrimary)
{
  // Its mother is not in the event: O2's test would read through a null one.
  o2::MCTrack orphan;
  orphan.SetMotherTrackId(5);
  orphan.setProcess(kPHadronic);
  BOOST_TEST((particleOrigin(orphan, std::vector<o2::MCTrack>{orphan}) == study::ParticleOrigin::kMaterial));
}

} // namespace o2::gtrk::study
