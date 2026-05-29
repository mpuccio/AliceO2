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

/// \file TrackFwd.h
/// \brief Base forward track model, params only, w/o covariance
///
/// \author Philippe Pillot, Subatech; adapted by Rafael Pezzi, UFRGS

#ifndef ALICEO2_BASE_TRACKFWD
#define ALICEO2_BASE_TRACKFWD

#include <array>
#include "GPUCommonRtypes.h"
#ifndef GPUCA_GPUCODE_DEVICE
#include <ostream>
#endif
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/TrackUtils.h"
#include "ReconstructionDataFormats/TrackParametrizationInterface.h"
#include "MathUtils/Primitive2D.h"

namespace o2::track
{

template <typename value_T>
class TrackParametrization; // fwd declaration for conversion method

template <typename value_T>
class TrackParametrizationWithError; // fwd declaration for conversion method

class TrackParFwd : public TrackParametrizationData<float, 5>, public TrackParFwdInterface<TrackParFwd, float>
{ // Forward track parameterization, kinematics only.
 public:
  using value_t = float;
  using base_t = TrackParametrizationData<value_t, 5>;
  using fwd_interface_t = TrackParFwdInterface<TrackParFwd, value_t>;
  using params_t = std::array<value_t, 5>;

  TrackParFwd() = default;
  ~TrackParFwd() = default;

  TrackParFwd(const TrackParFwd& tp) = default;
  TrackParFwd& operator=(const TrackParFwd& tp) = default;
  TrackParFwd(TrackParFwd&&) = delete;
  TrackParFwd& operator=(TrackParFwd&&) = delete;

  template <typename T>
  void toBarrelTrackPar(TrackParametrization<T>& t) const;

  // Track parameter propagation
  void propagateParamToZlinear(double zEnd);
  void propagateParamToZquadratic(double zEnd, double zField);
  void propagateParamToZhelix(double zEnd, double zField);
  void getCircleParams(float bz, o2::math_utils::CircleXY<float>& c, float& sna, float& csa) const;

 protected:
  using base_t::mP;
  using base_t::mX;

  /// Track parameters ordered as follow:      <pre>
  /// X       = X coordinate   (cm)
  /// Y       = Y coordinate   (cm)
  /// PHI     = azimutal angle
  /// TANL    = tangent of \lambda (dip angle)
  /// INVQPT    = Inverse transverse momentum (GeV/c ** -1) times charge (assumed forward motion)  </pre>

  // v2: storage entirely relocated into TrackParametrizationData<float, 5> base
  //     (was SMatrix5 mParameters + Double_t mZ + Double_t mTrackChi2). Old persistent
  //     data requires a TBufferIO read rule (SMatrix5 → mP[5], mZ → mX) to migrate.
  ClassDefNV(TrackParFwd, 2);
};

class TrackParCovFwd : public TrackParFwd, public TrackCovarianceData<TrackParFwd::value_t, 15>, public TrackParCovFwdInterface<TrackParCovFwd, TrackParFwd::value_t>
{ // Forward track+error parameterization
 public:
  using cov_base_t = TrackCovarianceData<value_t, 15>;
  using cov_interface_t = TrackParCovFwdInterface<TrackParCovFwd, value_t>;
  using covMat_t = std::array<value_t, 15>;
  using typename cov_interface_t::cov_view_t;
  using TrackParFwd::TrackParFwd; // inherit base constructors

  TrackParCovFwd() = default;
  ~TrackParCovFwd() = default;
  TrackParCovFwd& operator=(const TrackParCovFwd& tpf) = default;
  TrackParCovFwd(value_t z, const params_t& parameters, const covMat_t& covariances, value_t chi2);

  template <typename T>
  void toBarrelTrackParCov(TrackParametrizationWithError<T>& t) const;

  // Propagate parameters and covariances matrix
  void propagateToZlinear(double zEnd);
  void propagateToZquadratic(double zEnd, double zField);
  void propagateToZhelix(double zEnd, double zField);
  void propagateToZ(double zEnd, double zField); // Parameters: helix; errors: quadratic
  void propagateToDCAhelix(double zField, const std::array<double, 3>& p, std::array<double, 3>& dca);

  // Add Multiple Coulomb Scattering effects
  void addMCSEffect(double x2X0);

  // Kalman filter/fitting
  bool update(const std::array<float, 2>& p, const std::array<float, 2>& cov);

  // Propagate fwd track to vertex including MCS effects
  bool propagateToVtxhelixWithMCS(double z, const std::array<float, 2>& p, const std::array<float, 2>& cov, double field, double x_over_X0);
  bool propagateToVtxlinearWithMCS(double z, const std::array<float, 2>& p, const std::array<float, 2>& cov, double x_over_X0);

  bool getCovXYZPxPyPzGlo(std::array<float, 21>& cv) const;

 protected:
  using cov_base_t::mC;

 private:
  /// Covariance matrix of track parameters, ordered as follows:    <pre>
  ///  <X,X>         <Y,X>           <PHI,X>       <TANL,X>        <INVQPT,X>
  ///  <X,Y>         <Y,Y>           <PHI,Y>       <TANL,Y>        <INVQPT,Y>
  /// <X,PHI>       <Y,PHI>         <PHI,PHI>     <TANL,PHI>      <INVQPT,PHI>
  /// <X,TANL>      <Y,TANL>       <PHI,TANL>     <TANL,TANL>     <INVQPT,TANL>
  /// <X,INVQPT>   <Y,INVQPT>     <PHI,INVQPT>   <TANL,INVQPT>   <INVQPT,INVQPT>  </pre>
  // v2: cov storage relocated into TrackCovarianceData<float, 15> base (was SMatrix55Sym
  //     mCovariances); chi2 also moved to that base. Old persistent data requires a
  //     TBufferIO read rule (SMatrix55Sym → mC[15] packed) to migrate.
  ClassDefNV(TrackParCovFwd, 2);
};

static_assert(sizeof(TrackParFwd) == sizeof(TrackParFwd::base_t));
static_assert(sizeof(TrackParCovFwd) == sizeof(TrackParFwd) + sizeof(TrackParCovFwd::cov_base_t));

#ifndef GPUCA_GPUCODE_DEVICE
inline std::ostream& operator<<(std::ostream& os, const TrackParCovFwd::cov_view_t& cov)
{
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      os << cov(i, j);
      if (i != 4 || j != 4) {
        os << ' ';
      }
    }
  }
  return os;
}
#endif

} // namespace o2::track

#endif
