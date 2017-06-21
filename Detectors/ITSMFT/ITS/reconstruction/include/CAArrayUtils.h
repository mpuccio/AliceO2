// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CAArrayUtils.h
/// \brief
/// \author iacopo.colonnelli@cern.ch
/// \author maximiliano.puccio@cern.ch

#ifndef TRACKINGITSU_INCLUDE_CAARRAYUTILS_H_
#define TRACKINGITSU_INCLUDE_CAARRAYUTILS_H_

#include <array>
#include <cstddef>
#include <utility>

namespace CAArrayUtils {
template<typename T, std::size_t ...Is, typename Initializer>
constexpr std::array<T, sizeof...(Is)> fillArray(Initializer, std::index_sequence<Is...>);
template<typename T, std::size_t N, typename Initializer>
constexpr std::array<T, N> fillArray(Initializer);
}

template<typename T, std::size_t ...Is, typename Initializer>
constexpr std::array<T, sizeof...(Is)> CAArrayUtils::fillArray(Initializer initializer, std::index_sequence<Is...>)
{
  return std::array<T, sizeof...(Is)> { { initializer(Is)... } };
}

template<typename T, std::size_t N, typename Initializer>
constexpr std::array<T, N> CAArrayUtils::fillArray(Initializer initializer)
{
  return CAArrayUtils::fillArray<T>(initializer, std::make_index_sequence<N> { });
}

#endif /* TRACKINGITSU_INCLUDE_CAARRAYUTILS_H_ */
