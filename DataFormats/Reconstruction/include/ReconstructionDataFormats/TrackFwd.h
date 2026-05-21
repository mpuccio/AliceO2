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

  /// return Z coordinate (cm)
  value_t getZ() const { return mX; }
  /// set Z coordinate (cm)
  void setZ(value_t z) { mX = z; }
  value_t getX() const { return mP[0]; }
  void setX(value_t x) { mP[0] = x; }

  value_t getY() const { return mP[1]; }
  void setY(value_t y) { mP[1] = y; }

  void setPhi(value_t phi) { mP[2] = phi; }
  value_t getPhi() const { return mP[2]; }

  value_t getSnp() const
  {
    return o2::math_utils::sin(mP[2]);
  }

  value_t getCsp2() const
  {
    auto snp = o2::math_utils::sin(mP[2]);
    value_t csp = o2::math_utils::sqrt((value_t(1) - snp) * (value_t(1) + snp));
    return csp * csp;
  }

  void setTanl(value_t tanl) { mP[3] = tanl; }
  value_t getTanl() const { return mP[3]; }

  value_t getTgl() const { return mP[3]; } // for the sake of helixhelper

  void setInvQPt(value_t invqpt) { mP[4] = invqpt; }
  value_t getInvQPt() const { return mP[4]; } // return Inverse charged pt
  value_t getPt() const { return o2::math_utils::abs(value_t(1) / mP[4]); }
  value_t getInvPt() const { return o2::math_utils::abs(mP[4]); }
  value_t getPx() const { return o2::math_utils::cos(getPhi()) * getPt(); } // return px
  value_t getPy() const { return o2::math_utils::sin(getPhi()) * getPt(); } // return py
  value_t getPz() const { return getTanl() * getPt(); }                     // return pz
  value_t getP() const { return getPt() * o2::math_utils::sqrt(value_t(1) + getTanl() * getTanl()); }
  value_t getInverseMomentum() const { return value_t(1) / getP(); }

  value_t getTheta() const { return value_t(0.5) * o2::math_utils::pi() - o2::math_utils::atan(getTanl()); }
  value_t getEta() const { return -o2::math_utils::log(o2::math_utils::tan(getTheta() / value_t(2))); } // return total momentum

  value_t getCurvature(value_t b) const
  {
    auto invqpt = getInvQPt();
    return o2::constants::math::B2C * b * invqpt;
  }

  /// return the charge (assumed forward motion)
  int getCharge() const { return mP[4] >= 0.f ? 1 : -1; }
  /// set the charge (assumed forward motion)
  void setCharge(int charge)
  {
    if (charge * mP[4] < 0.) {
      mP[4] *= -1.;
    }
  }

  /// return track parameters
  const value_t* getParameters() const { return mP; }
  const value_t* getParams() const { return mP; }
  params_t getParametersArray() const
  {
    params_t params;
    for (int i = 0; i < 5; ++i) {
      params[i] = mP[i];
    }
    return params;
  }
  /// set track parameters
  void setParameters(const params_t& parameters) { setParameters(parameters.data()); }
  void setParameters(const value_t* parameters)
  {
    for (int i = 0; i < 5; ++i) {
      mP[i] = parameters[i];
    }
  }
  /// add track parameters
  void addParameters(const params_t& parameters) { addParameters(parameters.data()); }
  void addParameters(const value_t* parameters)
  {
    for (int i = 0; i < 5; ++i) {
      mP[i] += parameters[i];
    }
  }

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

  ClassDefNV(TrackParFwd, 1);
};

class TrackParCovFwd : public TrackParFwd, public TrackCovarianceData<TrackParFwd::value_t, 15>, public TrackParCovFwdInterface<TrackParCovFwd, TrackParFwd::value_t>
{ // Forward track+error parameterization
 public:
  using cov_base_t = TrackCovarianceData<value_t, 15>;
  using cov_interface_t = TrackParCovFwdInterface<TrackParCovFwd, value_t>;
  using covMat_t = std::array<value_t, 15>;
  struct cov_view_t {
    const value_t* data = nullptr;
    value_t operator()(int i, int j) const { return data[TrackParCovFwd::covIndex(i, j)]; }
  };
  using TrackParFwd::TrackParFwd; // inherit base constructors

  TrackParCovFwd() = default;
  ~TrackParCovFwd() = default;
  TrackParCovFwd& operator=(const TrackParCovFwd& tpf) = default;
  TrackParCovFwd(value_t z, const params_t& parameters, const covMat_t& covariances, value_t chi2);

  template <typename T>
  void toBarrelTrackParCov(TrackParametrizationWithError<T>& t) const;

  cov_view_t getCovariances() const { return {mC.data()}; }
  const covMat_t& getCov() const { return mC; }
  value_t getCovarElem(int i, int j) const { return mC[covIndex(i, j)]; }
  void setCovariances(const covMat_t& covariances) { mC = covariances; }
  void setCovariances(const value_t* covariances);
  void deleteCovariances() { mC = {}; }

  value_t getSigma2X() const { return mC[covIndex(0, 0)]; }
  value_t getSigma2Y() const { return mC[covIndex(1, 1)]; }
  value_t getSigmaXY() const { return mC[covIndex(0, 1)]; }
  value_t getSigma2Phi() const { return mC[covIndex(2, 2)]; }
  value_t getSigma2Tanl() const { return mC[covIndex(3, 3)]; }
  value_t getSigma2InvQPt() const { return mC[covIndex(4, 4)]; }

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

  static constexpr int covIndex(int i, int j)
  {
    return i >= j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i;
  }

 protected:
  using cov_base_t::mC;

 private:
  /// Covariance matrix of track parameters, ordered as follows:    <pre>
  ///  <X,X>         <Y,X>           <PHI,X>       <TANL,X>        <INVQPT,X>
  ///  <X,Y>         <Y,Y>           <PHI,Y>       <TANL,Y>        <INVQPT,Y>
  /// <X,PHI>       <Y,PHI>         <PHI,PHI>     <TANL,PHI>      <INVQPT,PHI>
  /// <X,TANL>      <Y,TANL>       <PHI,TANL>     <TANL,TANL>     <INVQPT,TANL>
  /// <X,INVQPT>   <Y,INVQPT>     <PHI,INVQPT>   <TANL,INVQPT>   <INVQPT,INVQPT>  </pre>
  ClassDefNV(TrackParCovFwd, 1);
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
