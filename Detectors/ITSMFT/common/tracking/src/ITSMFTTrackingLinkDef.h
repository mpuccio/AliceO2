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

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class o2::itsmft::TrackerParamConfig < o2::detectors::DetID::MFT> + ;
#pragma link C++ class o2::conf::ConfigurableParamHelper < o2::itsmft::TrackerParamConfig < o2::detectors::DetID::MFT>> + ;

// workflow-onboarding Slice 2 fix: ITSCommonCATrackerParam (added in Slice 1)
// was never added here, so it had no ROOT dictionary and was invisible to
// o2::conf::ConfigurableParam::updateFromString()/setValue() at runtime --
// undetected until a test actually exercised the string-keyed path (Slice 1's
// own tests only ever read ITSCommonCATrackerParam::Instance() directly in
// C++). The opt-in ITS common-CA workflow's --configKeyValues contract
// requires this to work.
#pragma link C++ class o2::itsmft::ITSCommonCATrackerParam + ;
#pragma link C++ class o2::conf::ConfigurableParamHelper < o2::itsmft::ITSCommonCATrackerParam> + ;

#endif
