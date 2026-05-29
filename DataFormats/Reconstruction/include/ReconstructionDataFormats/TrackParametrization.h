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

/// @file   TrackParametrization.h
/// @author ruben.shahoyan@cern.ch, michael.lettrich@cern.ch
/// @since  Oct 1, 2020
/// @brief

/*
  24/09/2020: Added new data member for abs. charge. This is needed for uniform treatment of tracks with non-standard
  charge: 0 (for V0s) and e.g. 2 for hypernuclei.
  In the aliroot AliExternalTrackParam this was treated by derived classes using virtual methods, which we don't use in O2.
  The meaning of mP[kQ2Pt] remains exactly the same, except for q=0 case: in this case the mP[kQ2Pt] is just an alias to
  1/pT, regardless of its sign, and the getCurvature() will 0 (because mAbsCharge is 0).
  The methods returning lab momentum or its combination account for eventual q>1.
 */

#ifndef INCLUDE_RECONSTRUCTIONDATAFORMATS_TRACKPARAMETRIZATION_H_
#define INCLUDE_RECONSTRUCTIONDATAFORMATS_TRACKPARAMETRIZATION_H_

#include "GPUCommonDef.h"
#include "GPUCommonRtypes.h"
#include "GPUCommonMath.h"
#include "GPUROOTCartesianFwd.h"

#include <array>

#ifndef GPUCA_GPUCODE_DEVICE
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <iosfwd>
#include <type_traits>
#endif

#ifndef GPUCA_ALIGPUCODE // Used only by functions that are hidden on the GPU
#include "ReconstructionDataFormats/BaseCluster.h"
#include <string>
#endif

#include "CommonConstants/MathConstants.h"
#include "MathUtils/Utils.h"
#include "MathUtils/Primitive2D.h"
#include "ReconstructionDataFormats/PID.h"

#include "ReconstructionDataFormats/TrackParametrizationData.h"
#include "ReconstructionDataFormats/TrackParametrizationInterface.h"
#include "ReconstructionDataFormats/TrackUtils.h"

namespace o2
{
template <typename T>
class BaseCluster;

namespace dataformats
{
class VertexBase;
class DCA;
} // namespace dataformats

namespace track
{

class TrackParFwd; // fwd declaration for conversion method

// aliases for track elements
enum ParLabels : int { kY,
                       kZ,
                       kSnp,
                       kTgl,
                       kQ2Pt };
enum CovLabels : int {
  kSigY2,
  kSigZY,
  kSigZ2,
  kSigSnpY,
  kSigSnpZ,
  kSigSnp2,
  kSigTglY,
  kSigTglZ,
  kSigTglSnp,
  kSigTgl2,
  kSigQ2PtY,
  kSigQ2PtZ,
  kSigQ2PtSnp,
  kSigQ2PtTgl,
  kSigQ2Pt2
};

enum DirType : int { DirInward = -1,
                     DirAuto = 0,
                     DirOutward = 1 };

constexpr int kNParams = 5, kCovMatSize = 15, kLabCovMatSize = 21;

constexpr float kCY2max = 100 * 100, // SigmaY<=100cm
  kCZ2max = 100 * 100,               // SigmaZ<=100cm
  kCSnp2max = 1 * 1,                 // SigmaSin<=1
  kCTgl2max = 1 * 1,                 // SigmaTan<=1
  kC1Pt2max = 100 * 100,             // Sigma1/Pt<=100 1/GeV
  kMostProbablePt = 0.6f,            // Most Probable Pt (GeV), for running with Bz=0
  kCalcdEdxAuto = -999.f;            // value indicating request for dedx calculation

// access to covariance matrix by row and column
GPUconstexpr() int CovarMap[kNParams][kNParams] = {{0, 1, 3, 6, 10},
                                                   {1, 2, 4, 7, 11},
                                                   {3, 4, 5, 8, 12},
                                                   {6, 7, 8, 9, 13},
                                                   {10, 11, 12, 13, 14}};

// access to covariance matrix diagonal elements
GPUconstexpr() int DiagMap[kNParams] = {0, 2, 5, 9, 14};

constexpr float HugeF = o2::constants::math::VeryBig;
constexpr float ELoss2EKinThreshInv = 1. / 0.025; // do not allow E.Loss correction step with dE/Ekin above the inverse of this value
constexpr int MaxELossIter = 50;                  // max number of iteration for the ELoss to account for BB dependence on beta*gamma
constexpr float DefaultDCA = 999.f;               // default DCA value
constexpr float DefaultDCACov = 999.f;            // default DCA cov value

// uncomment this to enable correction for BB dependence on beta*gamma via BB derivative
// #define _BB_NONCONST_CORR_

template <typename value_T = float>
class TrackParametrization : public TrackParametrizationData<value_T, kNParams>,
                              public TrackParBarrelInterface<TrackParametrization<value_T>, value_T>
{ // track parameterization, kinematics only.

 public:
  using base_t = TrackParametrizationData<value_T, kNParams>;
  using base_interface_t = TrackParBarrelInterface<TrackParametrization<value_T>, value_T>;
  using value_t = value_T;
  using dim2_t = std::array<value_t, 2>;
  using dim3_t = std::array<value_t, 3>;
  using params_t = std::array<value_t, kNParams>;

  // Expose interface accessors in the derived scope (the dependent base hides
  // them from unqualified lookup otherwise).
  using base_interface_t::getAlpha;
  using base_interface_t::getAbsCharge;
  using base_interface_t::getCharge2Pt;
  using base_interface_t::getCsp;
  using base_interface_t::getCsp2;
  using base_interface_t::getCurvature;
  using base_interface_t::getE;
  using base_interface_t::getE2;
  using base_interface_t::getEta;
  using base_interface_t::getInverseMomentum;
  using base_interface_t::getP;
  using base_interface_t::getP2;
  using base_interface_t::getP2Inv;
  using base_interface_t::getPInv;
  using base_interface_t::getPt;
  using base_interface_t::getPtInv;
  using base_interface_t::getQ2P2;
  using base_interface_t::getTheta;
  using base_interface_t::getParam;
  using base_interface_t::getParams;
  using base_interface_t::getPID;
  using base_interface_t::getQ2Pt;
  using base_interface_t::getSnp;
  using base_interface_t::getTgl;
  using base_interface_t::getUserField;
  using base_interface_t::getX;
  using base_interface_t::getY;
  using base_interface_t::getZ;
  using base_interface_t::setAbsCharge;
  using base_interface_t::setAlpha;
  using base_interface_t::setParam;
  using base_interface_t::setParams;
  using base_interface_t::setPID;
  using base_interface_t::setQ2Pt;
  using base_interface_t::setSnp;
  using base_interface_t::setTgl;
  using base_interface_t::setUserField;
  using base_interface_t::setX;
  using base_interface_t::setY;
  using base_interface_t::setZ;

  struct yzerr_t { // 2 measurement with error
    dim2_t yz;
    dim3_t yzerr;
  };

#ifndef GPUCA_GPUCODE_DEVICE
  static_assert(std::is_floating_point_v<value_t>);
#endif

  GPUdDefault() TrackParametrization() = default;
  GPUd() TrackParametrization(value_t x, value_t alpha, const params_t& par, int charge = 1, const PID pid = PID::Pion);
  GPUd() TrackParametrization(const dim3_t& xyz, const dim3_t& pxpypz, int charge, bool sectorAlpha = true, const PID pid = PID::Pion);
  GPUdDefault() TrackParametrization(const TrackParametrization&) = default;
  GPUdDefault() TrackParametrization(TrackParametrization&&) = default;
  GPUhdDefault() TrackParametrization& operator=(const TrackParametrization& src) = default;
  GPUhdDefault() TrackParametrization& operator=(TrackParametrization&& src) = default;
  GPUdDefault() ~TrackParametrization() = default;

  GPUd() void set(value_t x, value_t alpha, const params_t& par, int charge = 1, const PID pid = PID::Pion);
  GPUd() void set(value_t x, value_t alpha, const value_t* par, int charge = 1, const PID pid = PID::Pion);
  GPUd() value_t getR2() const;
  GPUd() value_t getR() const;

  // derived getters
  GPUd() bool getXatLabR(value_t r, value_t& x, value_t bz, DirType dir = DirAuto) const;
  GPUd() void getCircleParamsLoc(value_t bz, o2::math_utils::CircleXY<value_t>& circle) const;
  GPUd() void getCircleParams(value_t bz, o2::math_utils::CircleXY<value_t>& circle, value_t& sna, value_t& csa) const;
  GPUd() void getLineParams(o2::math_utils::IntervalXY<value_t>& line, value_t& sna, value_t& csa) const;
  GPUd() int getCharge() const;
  GPUd() int getSign() const;
  GPUd() value_t getPhi() const;
  GPUd() value_t getPhiPos() const;

  GPUdi() static value_t getdEdxBB(value_t betagamma) { return BetheBlochSolid(betagamma); }
  GPUdi() static value_t getdEdxBBOpt(value_t betagamma) { return BetheBlochSolidOpt(betagamma); }
  GPUdi() static value_t getBetheBlochSolidDerivativeApprox(value_T dedx, value_T bg) { return BetheBlochSolidDerivative(dedx, bg); }

  GPUd() math_utils::Point3D<value_t> getXYZGlo() const;
  GPUd() void getXYZGlo(dim3_t& xyz) const;
  GPUd() bool getPxPyPzGlo(dim3_t& pxyz) const;
  GPUd() bool getPosDirGlo(std::array<value_t, 9>& posdirp) const;

  // methods for track params estimate at other point
  GPUd() bool getYZAt(value_t xk, value_t b, value_t& y, value_t& z) const;
  GPUd() value_t getZAt(value_t xk, value_t b) const;
  GPUd() value_t getYAt(value_t xk, value_t b) const;
  GPUd() value_t getSnpAt(value_t xk, value_t b) const;
  GPUd() value_t getSnpAt(value_t alpha, value_t xk, value_t b) const;
  GPUd() value_t getPhiAt(value_t xk, value_t b) const;
  GPUd() value_t getPhiPosAt(value_t xk, value_t b) const;
  GPUd() value_t getDCAYtoMV(value_t b, value_t xmv = 0.f, value_t ymv = 0.f, value_t zmv = 0.f) const;
  GPUd() value_t getDCAZtoMV(value_t b, value_t xmv = 0.f, value_t ymv = 0.f, value_t zmv = 0.f) const;
  GPUd() math_utils::Point3D<value_t> getXYZGloAt(value_t xk, value_t b, bool& ok) const;

  // parameters manipulation
  GPUd() bool correctForELoss(value_t xrho, bool anglecorr = false);
  GPUd() bool rotateParam(value_t alpha);
  GPUd() bool rotateParam(value_t& alpha, value_t& ca, value_t& sa);
  GPUd() bool propagateParamTo(value_t xk, value_t b);
  GPUd() bool propagateParamTo(value_t xk, const dim3_t& b);
  GPUd() void invertParam();
  GPUd() bool propagateParamToDCA(const math_utils::Point3D<value_t>& vtx, value_t b, dim2_t* dca = nullptr, value_t maxD = 999.f);
  // aliases
  GPUd() bool rotate(value_t alpha) { return rotateParam(alpha); }
  GPUd() bool propagateTo(value_t xk, value_t b) { return propagateParamTo(xk, b); }
  GPUd() bool propagateTo(value_t xk, const dim3_t& b) { return propagateParamTo(xk, b); }
  GPUd() void invert() { invertParam(); }
  GPUd() bool propagateToDCA(const math_utils::Point3D<value_t>& vtx, value_t b, dim2_t* dca = nullptr, value_t maxD = 999.f) { return propagateParamToDCA(vtx, b, dca, maxD); }

  GPUd() bool isValid() const;
  GPUd() void invalidate();

  GPUd() void printParam() const;
  GPUd() void printParamHexadecimal();
#ifndef GPUCA_ALIGPUCODE
  void toFwdTrackPar(TrackParFwd& t) const;
  std::string asString() const;
  std::string asStringHexadecimal();
  size_t hash() const { return hash(getX(), getAlpha(), getY(), getZ(), getSnp(), getTgl(), getQ2Pt()); }
  static size_t hash(float x, float alp, float y, float z, float snp, float tgl, float q2pt);
#endif

  GPUd() void updateParam(value_t delta, int i);
  GPUd() void updateParams(const params_t& delta);
  GPUd() void updateParams(const value_t* delta);

  GPUd() yzerr_t getVertexInTrackFrame(const o2::dataformats::VertexBase& vtx) const;

 protected:
  using base_t::mP;
  using base_t::mX;
  using base_t::mAlpha;
  using base_t::mAbsCharge;
  using base_t::mPID;
  using base_t::mUserField;

 private:
  static constexpr value_t InvalidX = -99999.f;

  // v4: mX/mP/mAlpha/mAbsCharge/mPID/mUserField relocated into TrackParametrizationData<value_T, kNParams> base.
  ClassDefNV(TrackParametrization, 4);
};

//____________________________________________________________
template <typename value_T>
GPUdi() TrackParametrization<value_T>::TrackParametrization(value_t x, value_t alpha, const params_t& par, int charge, const PID pid)
  : base_t{x, alpha, charge, pid}
{
  for (int i = 0; i < kNParams; i++) {
    mP[i] = par[i];
  }
}

//____________________________________________________________
template <typename value_T>
GPUdi() void TrackParametrization<value_T>::set(value_t x, value_t alpha, const params_t& par, int charge, const PID pid)
{
  set(x, alpha, par.data(), charge, pid);
}

//____________________________________________________________
template <typename value_T>
GPUdi() void TrackParametrization<value_T>::set(value_t x, value_t alpha, const value_t* par, int charge, const PID pid)
{
  mX = x;
  mAlpha = alpha;
  math_utils::detail::bringToPMPi<value_t>(mAlpha);
  mAbsCharge = char(gpu::CAMath::Abs(charge));
  for (int i = 0; i < kNParams; i++) {
    mP[i] = par[i];
  }
  mPID = pid;
}

//____________________________________________________________
template <typename value_T>
GPUdi() auto TrackParametrization<value_T>::getR2() const -> value_t
{
  return mX * mX + mP[kY] * mP[kY];
}

//____________________________________________________________
template <typename value_T>
GPUdi() auto TrackParametrization<value_T>::getR() const -> value_t
{
  return gpu::CAMath::Sqrt(getR2());
}

//_______________________________________________________
template <typename value_T>
GPUdi() void TrackParametrization<value_T>::getCircleParamsLoc(value_t bz, o2::math_utils::CircleXY<value_t>& c) const
{
  // get circle params in track local frame, for straight line just set to local coordinates
  c.rC = getCurvature(bz);
  // treat as straight track if sagitta between the vertex and middle of TPC is below 0.01 cm
  constexpr value_t MinSagitta = 0.01f, TPCMidR = 160.f, MinCurv = 8 * MinSagitta / (TPCMidR * TPCMidR);
  if (gpu::CAMath::Abs(c.rC) > MinCurv) {
    c.rC = 1.f / getCurvature(bz);
    value_t sn = getSnp(), cs = gpu::CAMath::Sqrt((1.f - sn) * (1.f + sn));
    c.xC = getX() - sn * c.rC; // center in tracking
    c.yC = getY() + cs * c.rC; // frame. Note: r is signed!!!
    c.rC = gpu::CAMath::Abs(c.rC);
  } else {
    c.rC = 0.f; // signal straight line
    c.xC = getX();
    c.yC = getY();
  }
}

//_______________________________________________________
template <typename value_T>
GPUdi() void TrackParametrization<value_T>::getCircleParams(value_t bz, o2::math_utils::CircleXY<value_t>& c, value_t& sna, value_t& csa) const
{
  // get circle params in loc and lab frame, for straight line just set to global coordinates
  getCircleParamsLoc(bz, c);
  o2::math_utils::detail::sincos(getAlpha(), sna, csa);
  o2::math_utils::detail::rotateZ<value_t>(c.xC, c.yC, c.xC, c.yC, sna, csa); // center in global frame
}

//_______________________________________________________
template <typename value_T>
GPUdi() void TrackParametrization<value_T>::getLineParams(o2::math_utils::IntervalXY<value_t>& ln, value_t& sna, value_t& csa) const
{
  // get line parameterization as { x = x0 + xSlp*t, y = y0 + ySlp*t }
  o2::math_utils::detail::sincos(getAlpha(), sna, csa);
  o2::math_utils::detail::rotateZ<value_t>(getX(), getY(), ln.getX0(), ln.getY0(), sna, csa); // reference point in global frame
  value_t snp = getSnp(), csp = gpu::CAMath::Sqrt((1.f - snp) * (1.f + snp));
  ln.setDX(csp * csa - snp * sna);
  ln.setDY(snp * csa + csp * sna);
}

//____________________________________________________________
template <typename value_T>
GPUdi() int TrackParametrization<value_T>::getCharge() const
{
  return getSign() > 0 ? mAbsCharge : -mAbsCharge;
}

//____________________________________________________________
template <typename value_T>
GPUdi() int TrackParametrization<value_T>::getSign() const
{
  return mAbsCharge ? (mP[kQ2Pt] > 0.f ? 1 : -1) : 0;
}

//_______________________________________________________
template <typename value_T>
GPUdi() auto TrackParametrization<value_T>::getPhi() const -> value_t
{
  // track pt direction phi (in 0:2pi range)
  value_t phi = gpu::CAMath::ASin(getSnp()) + getAlpha();
  math_utils::detail::bringTo02Pi<value_t>(phi);
  return phi;
}

//_______________________________________________________
template <typename value_T>
GPUdi() auto TrackParametrization<value_T>::getPhiPos() const -> value_t
{
  // angle of track position (in -pi:pi range)
  value_t phi = gpu::CAMath::ATan2(getY(), getX()) + getAlpha();
  math_utils::detail::bringTo02Pi<value_t>(phi);
  return phi;
}

//_______________________________________________________
template <typename value_T>
GPUdi() auto TrackParametrization<value_T>::getXYZGlo() const -> math_utils::Point3D<value_t>
{
#ifndef GPUCA_ALIGPUCODE
  return math_utils::Rotation2D<value_t>(getAlpha())(math_utils::Point3D<value_t>(getX(), getY(), getZ()));
#else // mockup on GPU without ROOT
  float sina, cosa;
  gpu::CAMath::SinCos(getAlpha(), sina, cosa);
  return math_utils::Point3D<value_t>(cosa * getX() - sina * getY(), cosa * getY() + sina * getX(), getZ());
#endif
}

//_______________________________________________________
template <typename value_T>
GPUdi() void TrackParametrization<value_T>::getXYZGlo(dim3_t& xyz) const
{
  // track coordinates in lab frame
  xyz[0] = getX();
  xyz[1] = getY();
  xyz[2] = getZ();
  math_utils::detail::rotateZ<value_t>(xyz, getAlpha());
}

//_______________________________________________________
template <typename value_T>
GPUdi() auto TrackParametrization<value_T>::getXYZGloAt(value_t xk, value_t b, bool& ok) const -> math_utils::Point3D<value_t>
{
  //----------------------------------------------------------------
  // estimate global X,Y,Z in global frame at given X
  //----------------------------------------------------------------
  value_t y = 0.f, z = 0.f;
  ok = getYZAt(xk, b, y, z);
  if (ok) {
#ifndef GPUCA_ALIGPUCODE
    return math_utils::Rotation2D<value_t>(getAlpha())(math_utils::Point3D<value_t>(xk, y, z));
#else // mockup on GPU without ROOT
    float sina, cosa;
    gpu::CAMath::SinCos(getAlpha(), sina, cosa);
    return math_utils::Point3D<value_t>(cosa * xk - sina * y, cosa * y + sina * xk, z);
#endif
  } else {
    return math_utils::Point3D<value_t>();
  }
}

//____________________________________________________________
template <typename value_T>
GPUdi() bool TrackParametrization<value_T>::isValid() const
{
  return mX != InvalidX;
}

//____________________________________________________________
template <typename value_T>
GPUdi() void TrackParametrization<value_T>::invalidate()
{
  mX = InvalidX;
}

//____________________________________________________________
template <typename value_T>
GPUdi() void TrackParametrization<value_T>::updateParam(value_t delta, int i)
{
  mP[i] += delta;
}

//____________________________________________________________
template <typename value_T>
GPUdi() void TrackParametrization<value_T>::updateParams(const params_t& delta)
{
  updateParams(delta.data());
}

//____________________________________________________________
template <typename value_T>
GPUdi() void TrackParametrization<value_T>::updateParams(const value_t* delta)
{
  for (int i = kNParams; i--;) {
    mP[i] += delta[i];
  }
  // make sure that snp is in the valid range
  if (mP[kSnp] > constants::math::Almost1) {
    mP[kSnp] = constants::math::Almost1;
  } else if (mP[kSnp] < -constants::math::Almost1) {
    mP[kSnp] = -constants::math::Almost1;
  }
}

#ifndef GPUCA_ALIGPUCODE
template <typename value_T>
size_t TrackParametrization<value_T>::hash(float x, float alp, float y, float z, float snp, float tgl, float q2pt)
{
  size_t h = std::hash<float>{}(o2::math_utils::detail::truncateFloatFraction(x, 0xFFFFFFF0));
  h ^= std::hash<float>{}(o2::math_utils::detail::truncateFloatFraction(alp, 0xFFFFFFF0)) << 1;
  h ^= std::hash<float>{}(o2::math_utils::detail::truncateFloatFraction(y, 0xFFFFFFF0)) << 1;
  h ^= std::hash<float>{}(o2::math_utils::detail::truncateFloatFraction(z, 0xFFFFFFF0)) << 1;
  h ^= std::hash<float>{}(o2::math_utils::detail::truncateFloatFraction(snp, 0xFFFFFF00)) << 1;
  h ^= std::hash<float>{}(o2::math_utils::detail::truncateFloatFraction(tgl, 0xFFFFFF00)) << 1;
  h ^= std::hash<float>{}(o2::math_utils::detail::truncateFloatFraction(q2pt, 0xFFFFFC00)) << 1;
  return h;
}
#endif

} // namespace track
} // namespace o2

#endif /* INCLUDE_RECONSTRUCTIONDATAFORMATS_TRACKPARAMETRIZATION_H_ */
