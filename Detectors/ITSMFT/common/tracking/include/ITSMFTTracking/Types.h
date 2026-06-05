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
///
/// \file Types.h
/// \brief Shared data-format types for ITSMFT CA tracking
///

#ifndef ALICEO2_ITSMFT_TRACKING_TYPES_H_
#define ALICEO2_ITSMFT_TRACKING_TYPES_H_

#include <cstdint>

#include "CommonDataFormat/TimeStamp.h"
#include "DataFormatsITS/TimeEstBC.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITS/Vertex.h"

namespace o2::itsmft::tracking
{

// BC time types used by ROFLookupTables (same layout as DataFormatsITS/TimeEstBC.h)
using TimeStampType = uint32_t;
using TimeStampErrorType = uint16_t;
using TimeStamp = o2::dataformats::TimeStampWithError<float, float>;

using TimeEstBC = o2::its::TimeEstBC;
using Vertex = o2::its::Vertex;
using VertexLabel = o2::its::VertexLabel;
using TrackITS = o2::its::TrackITS;
using TrackITSExt = o2::its::TrackITSExt;

} // namespace o2::itsmft::tracking

#endif /* ALICEO2_ITSMFT_TRACKING_TYPES_H_ */
