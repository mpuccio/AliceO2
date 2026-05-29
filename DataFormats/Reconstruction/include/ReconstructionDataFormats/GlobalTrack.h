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

/// \file GlobalTrack.h
/// \brief Tagged-union track holding either a barrel or a forward parameterization.
///
/// A GlobalTrack holds exactly one of {TrackPar, TrackParFwd} at any time, plus
/// a Frame tag identifying which. The storage is a real C++ union (active-member
/// switching via placement-new) so that asBarrel() / asFwd() return genuine
/// references to the currently-active leaf type without aliasing UB.
///
/// GlobalTrack is *not* persisted directly: callers should stream the underlying
/// barrel or fwd track on disk and rebuild a GlobalTrack on read via the
/// converting constructors or assignment operators ("adoption").

#ifndef ALICEO2_GLOBALTRACK_H
#define ALICEO2_GLOBALTRACK_H

#include <cassert>
#include <cstdint>
#include <new>
#include <utility>

#include "GPUCommonDef.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackFwd.h"

namespace o2::track
{

enum class Frame : uint8_t {
  Barrel = 0,
  Fwd = 1,
};

/// Param-only global track.
class GlobalTrack
{
 public:
  using value_t = float;

  GPUd() GlobalTrack() noexcept : mFrame(Frame::Barrel) { ::new (&mRep.barrel) TrackPar(); }
  GPUd() GlobalTrack(const TrackPar& src) noexcept : mFrame(Frame::Barrel) { ::new (&mRep.barrel) TrackPar(src); }
  GPUd() GlobalTrack(const TrackParFwd& src) noexcept : mFrame(Frame::Fwd) { ::new (&mRep.fwd) TrackParFwd(src); }

  GPUd() GlobalTrack(const GlobalTrack& o) noexcept : mFrame(o.mFrame)
  {
    if (isBarrel()) {
      ::new (&mRep.barrel) TrackPar(o.mRep.barrel);
    } else {
      ::new (&mRep.fwd) TrackParFwd(o.mRep.fwd);
    }
  }
  GPUd() GlobalTrack(GlobalTrack&& o) noexcept : mFrame(o.mFrame)
  {
    if (isBarrel()) {
      ::new (&mRep.barrel) TrackPar(std::move(o.mRep.barrel));
    } else {
      ::new (&mRep.fwd) TrackParFwd(std::move(o.mRep.fwd));
    }
  }
  GPUhd() GlobalTrack& operator=(const GlobalTrack& o) noexcept
  {
    if (this != &o) {
      destroyActive();
      mFrame = o.mFrame;
      if (isBarrel()) {
        ::new (&mRep.barrel) TrackPar(o.mRep.barrel);
      } else {
        ::new (&mRep.fwd) TrackParFwd(o.mRep.fwd);
      }
    }
    return *this;
  }
  GPUhd() GlobalTrack& operator=(GlobalTrack&& o) noexcept
  {
    destroyActive();
    mFrame = o.mFrame;
    if (isBarrel()) {
      ::new (&mRep.barrel) TrackPar(std::move(o.mRep.barrel));
    } else {
      ::new (&mRep.fwd) TrackParFwd(std::move(o.mRep.fwd));
    }
    return *this;
  }

  /// Adoption: replace the current rep with a barrel/fwd track read from disk
  /// (or constructed elsewhere). Flips the frame tag if necessary.
  GPUhd() GlobalTrack& operator=(const TrackPar& src) noexcept
  {
    if (isBarrel()) {
      mRep.barrel = src;
    } else {
      mRep.fwd.~TrackParFwd();
      ::new (&mRep.barrel) TrackPar(src);
      mFrame = Frame::Barrel;
    }
    return *this;
  }
  GPUhd() GlobalTrack& operator=(const TrackParFwd& src) noexcept
  {
    if (isFwd()) {
      mRep.fwd = src;
    } else {
      mRep.barrel.~TrackPar();
      ::new (&mRep.fwd) TrackParFwd(src);
      mFrame = Frame::Fwd;
    }
    return *this;
  }

  GPUd() ~GlobalTrack() noexcept { destroyActive(); }

  GPUd() Frame getFrame() const noexcept { return mFrame; }
  GPUd() bool isBarrel() const noexcept { return mFrame == Frame::Barrel; }
  GPUd() bool isFwd() const noexcept { return mFrame == Frame::Fwd; }

  /// Return a reference to the active barrel rep. With autoConvert=true the
  /// current fwd rep (if any) is converted in place first; otherwise a frame
  /// mismatch trips an assertion.
  TrackPar& asBarrel(bool autoConvert = false)
  {
    if (!isBarrel()) {
      assert(autoConvert && "GlobalTrack::asBarrel on fwd-frame without autoConvert");
      convertToBarrel();
    }
    return mRep.barrel;
  }
  GPUd() const TrackPar& asBarrel() const
  {
    assert(isBarrel() && "GlobalTrack::asBarrel const on fwd-frame");
    return mRep.barrel;
  }
  TrackParFwd& asFwd(bool autoConvert = false)
  {
    if (!isFwd()) {
      assert(autoConvert && "GlobalTrack::asFwd on barrel-frame without autoConvert");
      convertToFwd();
    }
    return mRep.fwd;
  }
  GPUd() const TrackParFwd& asFwd() const
  {
    assert(isFwd() && "GlobalTrack::asFwd const on barrel-frame");
    return mRep.fwd;
  }

  /// Rewrite the storage in the requested frame, using the existing
  /// toBarrel/toFwd conversion routines. No-op if already there.
  void convertToBarrel()
  {
    if (isBarrel()) {
      return;
    }
    TrackPar tmp;
    mRep.fwd.toBarrelTrackPar(tmp);
    mRep.fwd.~TrackParFwd();
    ::new (&mRep.barrel) TrackPar(tmp);
    mFrame = Frame::Barrel;
  }
  void convertToFwd()
  {
    if (isFwd()) {
      return;
    }
    TrackParFwd tmp;
    mRep.barrel.toFwdTrackPar(tmp);
    mRep.barrel.~TrackPar();
    ::new (&mRep.fwd) TrackParFwd(tmp);
    mFrame = Frame::Fwd;
  }

 private:
  union Rep {
    TrackPar barrel;
    TrackParFwd fwd;
    GPUd() Rep() noexcept {}
    GPUd() ~Rep() noexcept {}
  };
  Rep mRep;
  Frame mFrame;

  GPUd() void destroyActive() noexcept
  {
    if (isBarrel()) {
      mRep.barrel.~TrackPar();
    } else {
      mRep.fwd.~TrackParFwd();
    }
  }
};

/// Param + covariance flavour.
class GlobalTrackCov
{
 public:
  using value_t = float;

  GPUd() GlobalTrackCov() noexcept : mFrame(Frame::Barrel) { ::new (&mRep.barrel) TrackParCov(); }
  GPUd() GlobalTrackCov(const TrackParCov& src) noexcept : mFrame(Frame::Barrel) { ::new (&mRep.barrel) TrackParCov(src); }
  GPUd() GlobalTrackCov(const TrackParCovFwd& src) noexcept : mFrame(Frame::Fwd) { ::new (&mRep.fwd) TrackParCovFwd(src); }

  GPUd() GlobalTrackCov(const GlobalTrackCov& o) noexcept : mFrame(o.mFrame)
  {
    if (isBarrel()) {
      ::new (&mRep.barrel) TrackParCov(o.mRep.barrel);
    } else {
      ::new (&mRep.fwd) TrackParCovFwd(o.mRep.fwd);
    }
  }
  GPUd() GlobalTrackCov(GlobalTrackCov&& o) noexcept : mFrame(o.mFrame)
  {
    if (isBarrel()) {
      ::new (&mRep.barrel) TrackParCov(std::move(o.mRep.barrel));
    } else {
      ::new (&mRep.fwd) TrackParCovFwd(std::move(o.mRep.fwd));
    }
  }
  GPUhd() GlobalTrackCov& operator=(const GlobalTrackCov& o) noexcept
  {
    if (this != &o) {
      destroyActive();
      mFrame = o.mFrame;
      if (isBarrel()) {
        ::new (&mRep.barrel) TrackParCov(o.mRep.barrel);
      } else {
        ::new (&mRep.fwd) TrackParCovFwd(o.mRep.fwd);
      }
    }
    return *this;
  }
  GPUhd() GlobalTrackCov& operator=(GlobalTrackCov&& o) noexcept
  {
    destroyActive();
    mFrame = o.mFrame;
    if (isBarrel()) {
      ::new (&mRep.barrel) TrackParCov(std::move(o.mRep.barrel));
    } else {
      ::new (&mRep.fwd) TrackParCovFwd(std::move(o.mRep.fwd));
    }
    return *this;
  }

  GPUhd() GlobalTrackCov& operator=(const TrackParCov& src) noexcept
  {
    if (isBarrel()) {
      mRep.barrel = src;
    } else {
      mRep.fwd.~TrackParCovFwd();
      ::new (&mRep.barrel) TrackParCov(src);
      mFrame = Frame::Barrel;
    }
    return *this;
  }
  GPUhd() GlobalTrackCov& operator=(const TrackParCovFwd& src) noexcept
  {
    if (isFwd()) {
      mRep.fwd = src;
    } else {
      mRep.barrel.~TrackParCov();
      ::new (&mRep.fwd) TrackParCovFwd(src);
      mFrame = Frame::Fwd;
    }
    return *this;
  }

  GPUd() ~GlobalTrackCov() noexcept { destroyActive(); }

  GPUd() Frame getFrame() const noexcept { return mFrame; }
  GPUd() bool isBarrel() const noexcept { return mFrame == Frame::Barrel; }
  GPUd() bool isFwd() const noexcept { return mFrame == Frame::Fwd; }

  TrackParCov& asBarrel(bool autoConvert = false)
  {
    if (!isBarrel()) {
      assert(autoConvert && "GlobalTrackCov::asBarrel on fwd-frame without autoConvert");
      convertToBarrel();
    }
    return mRep.barrel;
  }
  GPUd() const TrackParCov& asBarrel() const
  {
    assert(isBarrel() && "GlobalTrackCov::asBarrel const on fwd-frame");
    return mRep.barrel;
  }
  TrackParCovFwd& asFwd(bool autoConvert = false)
  {
    if (!isFwd()) {
      assert(autoConvert && "GlobalTrackCov::asFwd on barrel-frame without autoConvert");
      [[maybe_unused]] auto ok = convertToFwd();
      assert(ok && "GlobalTrackCov::asFwd: barrel->fwd conversion failed");
    }
    return mRep.fwd;
  }
  GPUd() const TrackParCovFwd& asFwd() const
  {
    assert(isFwd() && "GlobalTrackCov::asFwd const on barrel-frame");
    return mRep.fwd;
  }

  /// barrel-from-fwd is unconditional.
  void convertToBarrel()
  {
    if (isBarrel()) {
      return;
    }
    TrackParCov tmp;
    mRep.fwd.toBarrelTrackParCov(tmp);
    mRep.fwd.~TrackParCovFwd();
    ::new (&mRep.barrel) TrackParCov(tmp);
    mFrame = Frame::Barrel;
  }
  /// fwd-from-barrel can fail (e.g. cov inversion). On failure the rep is
  /// left untouched and the call returns false.
  bool convertToFwd()
  {
    if (isFwd()) {
      return true;
    }
    TrackParCovFwd tmp;
    if (!mRep.barrel.toFwdTrackParCov(tmp)) {
      return false;
    }
    mRep.barrel.~TrackParCov();
    ::new (&mRep.fwd) TrackParCovFwd(tmp);
    mFrame = Frame::Fwd;
    return true;
  }

 private:
  union Rep {
    TrackParCov barrel;
    TrackParCovFwd fwd;
    GPUd() Rep() noexcept {}
    GPUd() ~Rep() noexcept {}
  };
  Rep mRep;
  Frame mFrame;

  GPUd() void destroyActive() noexcept
  {
    if (isBarrel()) {
      mRep.barrel.~TrackParCov();
    } else {
      mRep.fwd.~TrackParCovFwd();
    }
  }
};

} // namespace o2::track

#endif
