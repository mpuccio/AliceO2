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

/// @file TrackParametrizationInterface.h
/// @brief Storage-less CRTP interfaces for track parameterizations.

#ifndef INCLUDE_RECONSTRUCTIONDATAFORMATS_TRACKPARAMETRIZATIONINTERFACE_H_
#define INCLUDE_RECONSTRUCTIONDATAFORMATS_TRACKPARAMETRIZATIONINTERFACE_H_

#include "ReconstructionDataFormats/TrackParametrizationData.h"
#include "CommonConstants/MathConstants.h"
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/PID.h"

#include <array>
#include <cstdint>

namespace o2::track
{

constexpr float MaxPT = 100000.;       // cap pT to avoid NaNs in derived kinematics
constexpr float MinPTInv = 1. / MaxPT; // floor on |q/pT|

/// CRTP mixin: stateless accessors over TrackParametrizationData.
///
/// Note: `derived_T` is used only to locate the right `TrackParametrizationData` base
/// via static_cast. When a class C inherits from B which mixes in this interface (e.g.
/// `TrackParCovFwd : TrackParFwd : TrackParFwdInterface<TrackParFwd,...>`), the
/// interface methods see `derived_T == B`, not the actual most-derived `C`. Do NOT add
/// methods to this mixin that need the actual leaf type (e.g. returning `derived_T` by
/// value, or calling methods only defined on `C`) — they would slice / miss overrides.
template <typename derived_T, typename value_T, int nParams>
class TrackParametrizationInterface
{
 public:
  using value_t = value_T;
  using params_t = std::array<value_t, nParams>;

  GPUd() const value_t* getParams() const { return data().mP; }
  GPUd() value_t getParam(int i) const { return data().mP[i]; }
  GPUd() value_t getS() const { return data().mX; }
  GPUd() void setS(value_t s) { data().mX = s; }
  GPUd() void setParam(value_t v, int i) { data().mP[i] = v; }
  GPUd() void setParams(const params_t& params) { setParams(params.data()); }
  GPUd() void setParams(const value_t* params)
  {
    for (int i = 0; i < nParams; ++i) {
      data().mP[i] = params[i];
    }
  }

  GPUd() value_t getAlpha() const { return data().mAlpha; }
  GPUd() void setAlpha(value_t alpha)
  {
    data().mAlpha = alpha;
    math_utils::detail::bringToPMPi<value_t>(data().mAlpha);
  }
  GPUd() int getAbsCharge() const { return data().mAbsCharge; }
  GPUd() void setAbsCharge(int q) { data().mAbsCharge = static_cast<char>(q < 0 ? -q : q); }
  GPUd() PID getPID() const { return data().mPID; }
  GPUd() void setPID(const PID pid, bool passCharge = false)
  {
    data().mPID = pid;
    if (passCharge) {
      setAbsCharge(pid.getCharge());
    }
  }
  GPUhd() uint16_t getUserField() const { return data().mUserField; }
  GPUhd() void setUserField(uint16_t v) { data().mUserField = v; }

  // Frame-agnostic kinematics, all in terms of common storage:
  //   mP[nParams-2] = tan(lambda) (same for barrel/fwd)
  //   mP[nParams-1] = signed inverse pT (same for barrel/fwd)
  //   mAbsCharge    = absolute charge (>1 only for multi-charged particles)
  GPUd() value_t getPtInv() const
  {
    value_t ptInv = o2::math_utils::abs(data().mP[nParams - 1]);
    if (ptInv < MinPTInv) {
      ptInv = MinPTInv;
    }
    return (data().mAbsCharge > 1) ? ptInv / data().mAbsCharge : ptInv;
  }
  GPUd() value_t getPt() const { return value_t(1) / getPtInv(); }
  GPUd() value_t getQ2P2() const
  {
    value_t q2pt2 = data().mP[nParams - 1] * data().mP[nParams - 1];
    if (q2pt2 < MinPTInv * MinPTInv) {
      q2pt2 = MinPTInv * MinPTInv;
    }
    const value_t tgl = data().mP[nParams - 2];
    return q2pt2 / (value_t(1) + tgl * tgl);
  }
  GPUd() value_t getP2Inv() const
  {
    const value_t p2 = getPtInv();
    const value_t tgl = data().mP[nParams - 2];
    return p2 * p2 / (value_t(1) + tgl * tgl);
  }
  GPUd() value_t getP2() const { return value_t(1) / getP2Inv(); }
  GPUd() value_t getPInv() const
  {
    const value_t tgl = data().mP[nParams - 2];
    return getPtInv() / o2::math_utils::sqrt(value_t(1) + tgl * tgl);
  }
  GPUd() value_t getP() const { return value_t(1) / getPInv(); }
  GPUd() value_t getInverseMomentum() const { return getPInv(); }
  GPUd() value_t getTheta() const
  {
    const value_t tgl = data().mP[nParams - 2];
    return value_t(0.5) * o2::math_utils::pi() - o2::math_utils::atan(tgl);
  }
  GPUd() value_t getEta() const
  {
    return -o2::math_utils::log(o2::math_utils::tan(getTheta() / value_t(2)));
  }
  GPUd() value_t getE2() const { return getP2() + getPID().getMass2(); }
  GPUd() value_t getE() const { return o2::math_utils::sqrt(getE2()); }

 protected:
  GPUdDefault() TrackParametrizationInterface() = default;
  GPUdDefault() TrackParametrizationInterface(const TrackParametrizationInterface&) = default;
  GPUdDefault() TrackParametrizationInterface(TrackParametrizationInterface&&) = default;
  GPUhdDefault() TrackParametrizationInterface& operator=(const TrackParametrizationInterface&) = default;
  GPUhdDefault() TrackParametrizationInterface& operator=(TrackParametrizationInterface&&) = default;
  GPUdDefault() ~TrackParametrizationInterface() = default;

  GPUd() TrackParametrizationData<value_t, nParams>& data()
  {
    return static_cast<TrackParametrizationData<value_t, nParams>&>(derived());
  }
  GPUd() const TrackParametrizationData<value_t, nParams>& data() const
  {
    return static_cast<const TrackParametrizationData<value_t, nParams>&>(derived());
  }

 private:
  GPUd() derived_T& derived() { return static_cast<derived_T&>(*this); }
  GPUd() const derived_T& derived() const { return static_cast<const derived_T&>(*this); }
};

/// CRTP mixin: stateless accessors over TrackCovarianceData. Same `derived_T` caveat as
/// TrackParametrizationInterface above.
template <typename derived_T, typename value_T, int nCov>
class TrackCovarianceInterface
{
 public:
  using value_t = value_T;
  using covMat_t = std::array<value_t, nCov>;

  struct cov_view_t {
    const value_t* data = nullptr;
    GPUd() value_t operator()(int i, int j) const { return data[covIndex(i, j)]; }
  };

  GPUd() cov_view_t getCovariances() const { return {covData().mC.data()}; }
  GPUd() const covMat_t& getCov() const { return covData().mC; }
  GPUd() value_t getCovarElem(int i, int j) const { return covData().mC[covIndex(i, j)]; }
  GPUd() void setCov(const covMat_t& covariances) { covData().mC = covariances; }
  GPUd() void setCov(const value_t* covariances)
  {
    for (int i = 0; i < nCov; ++i) {
      covData().mC[i] = covariances[i];
    }
  }
  GPUd() void deleteCovariances() { covData().mC = {}; }

  GPUd() value_t getTrackChi2() const { return covData().mTrackChi2; }
  GPUd() void setTrackChi2(value_t chi2) { covData().mTrackChi2 = chi2; }

  GPUdi() static constexpr int covIndex(int i, int j)
  {
    return i >= j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i;
  }

 protected:
  GPUdDefault() TrackCovarianceInterface() = default;
  GPUdDefault() TrackCovarianceInterface(const TrackCovarianceInterface&) = default;
  GPUdDefault() TrackCovarianceInterface(TrackCovarianceInterface&&) = default;
  GPUhdDefault() TrackCovarianceInterface& operator=(const TrackCovarianceInterface&) = default;
  GPUhdDefault() TrackCovarianceInterface& operator=(TrackCovarianceInterface&&) = default;
  GPUdDefault() ~TrackCovarianceInterface() = default;

  GPUd() TrackCovarianceData<value_t, nCov>& covData()
  {
    return static_cast<TrackCovarianceData<value_t, nCov>&>(derived());
  }
  GPUd() const TrackCovarianceData<value_t, nCov>& covData() const
  {
    return static_cast<const TrackCovarianceData<value_t, nCov>&>(derived());
  }

 private:
  GPUd() derived_T& derived() { return static_cast<derived_T&>(*this); }
  GPUd() const derived_T& derived() const { return static_cast<const derived_T&>(*this); }
};

template <typename derived_T, typename value_T>
class TrackParBarrelInterface : public TrackParametrizationInterface<derived_T, value_T, 5>
{
 public:
  using base_interface_t = TrackParametrizationInterface<derived_T, value_T, 5>;
  using typename base_interface_t::params_t;
  using typename base_interface_t::value_t;

  GPUd() value_t getX() const { return this->data().mX; }
  GPUd() void setX(value_t x) { this->data().mX = x; }
  GPUd() value_t getY() const { return this->data().mP[0]; }
  GPUd() void setY(value_t y) { this->data().mP[0] = y; }
  GPUd() value_t getZ() const { return this->data().mP[1]; }
  GPUd() void setZ(value_t z) { this->data().mP[1] = z; }
  GPUd() value_t getSnp() const { return this->data().mP[2]; }
  GPUd() void setSnp(value_t snp) { this->data().mP[2] = snp; }
  GPUd() value_t getTgl() const { return this->data().mP[3]; }
  GPUd() void setTgl(value_t tgl) { this->data().mP[3] = tgl; }
  GPUhd() value_t getQ2Pt() const { return this->data().mP[4]; }
  GPUd() void setQ2Pt(value_t q2pt) { this->data().mP[4] = q2pt; }
  GPUd() value_t getCharge2Pt() const { return this->data().mAbsCharge ? this->data().mP[4] : value_t(0); }

  GPUd() value_t getCsp2() const
  {
    const value_t snp = getSnp();
    const value_t csp2 = (value_t(1) - snp) * (value_t(1) + snp);
    return csp2 > o2::constants::math::Almost0 ? csp2 : value_t(o2::constants::math::Almost0);
  }
  GPUd() value_t getCsp() const { return o2::math_utils::sqrt(getCsp2()); }

  GPUd() value_t getCurvature(value_t b) const
  {
    return this->data().mAbsCharge ? this->data().mP[4] * b * o2::constants::math::B2C : value_t(0);
  }
};

template <typename derived_T, typename value_T>
class TrackParFwdInterface : public TrackParametrizationInterface<derived_T, value_T, 5>
{
 public:
  using base_interface_t = TrackParametrizationInterface<derived_T, value_T, 5>;
  using typename base_interface_t::params_t;
  using typename base_interface_t::value_t;

  GPUd() value_t getZ() const { return this->data().mX; }
  GPUd() void setZ(value_t z) { this->data().mX = z; }
  GPUd() value_t getX() const { return this->data().mP[0]; }
  GPUd() void setX(value_t x) { this->data().mP[0] = x; }
  GPUd() value_t getY() const { return this->data().mP[1]; }
  GPUd() void setY(value_t y) { this->data().mP[1] = y; }
  GPUd() value_t getPhi() const { return this->data().mP[2]; }
  GPUd() void setPhi(value_t phi) { this->data().mP[2] = phi; }
  GPUd() value_t getSnp() const { return o2::math_utils::sin(this->data().mP[2]); }
  GPUd() value_t getCsp2() const
  {
    auto snp = getSnp();
    value_t csp = o2::math_utils::sqrt((value_t(1) - snp) * (value_t(1) + snp));
    return csp * csp;
  }
  GPUd() value_t getTanl() const { return this->data().mP[3]; }
  GPUd() void setTanl(value_t tanl) { this->data().mP[3] = tanl; }
  GPUd() value_t getTgl() const { return this->data().mP[3]; }
  GPUd() value_t getInvQPt() const { return this->data().mP[4]; }
  GPUd() void setInvQPt(value_t invqpt) { this->data().mP[4] = invqpt; }
  GPUd() value_t getInvPt() const { return this->getPtInv(); } // fwd-side alias
  GPUd() value_t getPx() const { return o2::math_utils::cos(getPhi()) * this->getPt(); }
  GPUd() value_t getPy() const { return o2::math_utils::sin(getPhi()) * this->getPt(); }
  GPUd() value_t getPz() const { return getTanl() * this->getPt(); }
  GPUd() value_t getCurvature(value_t b) const { return o2::constants::math::B2C * b * getInvQPt(); }
  GPUd() int getCharge() const { return this->data().mP[4] >= 0.f ? 1 : -1; }
  GPUd() void setCharge(int charge)
  {
    if (charge * this->data().mP[4] < 0.) {
      this->data().mP[4] *= -1.;
    }
  }
  GPUd() const value_t* getParameters() const { return this->data().mP; }
  GPUd() params_t getParametersArray() const
  {
    params_t params;
    for (int i = 0; i < 5; ++i) {
      params[i] = this->data().mP[i];
    }
    return params;
  }
  GPUd() void setParameters(const params_t& parameters) { setParameters(parameters.data()); }
  GPUd() void setParameters(const value_t* parameters) { this->setParams(parameters); }
  GPUd() void addParameters(const params_t& parameters) { addParameters(parameters.data()); }
  GPUd() void addParameters(const value_t* parameters)
  {
    for (int i = 0; i < 5; ++i) {
      this->data().mP[i] += parameters[i];
    }
  }
};

template <typename derived_T, typename value_T>
class TrackParCovFwdInterface : public TrackCovarianceInterface<derived_T, value_T, 15>
{
 public:
  using base_interface_t = TrackCovarianceInterface<derived_T, value_T, 15>;
  using typename base_interface_t::covMat_t;
  using typename base_interface_t::cov_view_t;
  using typename base_interface_t::value_t;
  using base_interface_t::covIndex;
  using base_interface_t::getCov;
  using base_interface_t::getCovariances;
  using base_interface_t::getCovarElem;

  GPUd() void setCovariances(const covMat_t& covariances) { this->setCov(covariances); }
  GPUd() void setCovariances(const value_t* covariances) { this->setCov(covariances); }
  GPUd() value_t getSigma2X() const { return getCovarElem(0, 0); }
  GPUd() value_t getSigma2Y() const { return getCovarElem(1, 1); }
  GPUd() value_t getSigmaXY() const { return getCovarElem(0, 1); }
  GPUd() value_t getSigma2Phi() const { return getCovarElem(2, 2); }
  GPUd() value_t getSigma2Tanl() const { return getCovarElem(3, 3); }
  GPUd() value_t getSigma2InvQPt() const { return getCovarElem(4, 4); }
};

} // namespace o2::track

#endif
