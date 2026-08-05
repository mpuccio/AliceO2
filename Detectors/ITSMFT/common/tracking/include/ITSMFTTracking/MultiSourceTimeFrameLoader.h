// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// M2b (GenericTrackingEngineMigration.md; ADR 0007): loadEvent() is the one
// atomic loading transaction, participant-count-generic. It never names
// ITS, MFT, a layer count, or a fixed source position -- that knowledge
// lives only in loadITSAndMFT()'s own thin wrapper below and in whatever
// adapter builds an AtomicLoadBinding list (today, ITSMFTLegacyParticipantSet
// ::loadBindings(), through its two plan-driven participants' loadTarget(),
// called directly by
// the combined DPL task's own trackFrame() composition, M3).
// TrackingEngine::executeEvent() must only ever be called once loadEvent()
// (or loadITSAndMFT()) reports success -- this loading boundary calls
// neither.

#ifndef ALICEO2_ITSMFT_TRACKING_MULTISOURCETIMEFRAMELOADER_H_
#define ALICEO2_ITSMFT_TRACKING_MULTISOURCETIMEFRAMELOADER_H_

#include <gsl/gsl>

#include "ITSMFTTracking/ClusterSource.h"
#include "ITSMFTTracking/SurfaceTrackingScratch.h"
#include "ITSMFTTracking/TrackingConfigParam.h"

namespace o2::itsmft::tracking
{

class MultiSourceTimeFrameLoader
{
 public:
  // Type-erased, NLayers-free per-source legacy-backfill staging target for
  // the generic atomic loading transaction below. Nested here rather than
  // free-standing, and deliberately never included by TrackingEngine.h/
  // TrackingParticipant.h: this legacy staging mechanism is private to the
  // loading boundary and must not leak into the generic engine/participant
  // public contract (ADR 0007 decision 5).
  class LoadTarget
  {
   public:
    virtual ~LoadTarget() = default;

    // Stages `source`'s legacy backfill into state private to this
    // target's own implementation, bound to its live scratch's allocator
    // identity. Touches only that private staging state -- never the live
    // scratch -- so a failure here leaves the live scratch completely
    // unchanged and safe to retry or abandon.
    virtual LoadSourcesResult stage(const ClusterSourceInput& source, SurfaceCatalogView catalog,
                                    const o2::InteractionRecord& origin) = 0;

    // Swaps the already-staged backfill into the live scratch. Same-
    // allocator swaps only -- never throws. Must only be called after
    // every stage() call in the same whole-event transaction has already
    // returned ok().
    virtual void commit() noexcept = 0;
  };

  // Concrete LoadTarget bound to one live SurfaceTrackingScratch&, mirroring
  // the stage()/commit() contract over
  // SurfaceTrackingScratch::loadNormalizedSource()'s runtime-sized
  // orderedSurfaces loop instead of a fixed NLayers one.
  // Each plan-driven participant owns one of these, bound to its own scratch.
  class LoadTargetImplSurface final : public LoadTarget
  {
   public:
    explicit LoadTargetImplSurface(SurfaceTrackingScratch& live) noexcept : mLive(live) {}

    LoadTargetImplSurface(const LoadTargetImplSurface&) = delete;
    LoadTargetImplSurface& operator=(const LoadTargetImplSurface&) = delete;
    LoadTargetImplSurface(LoadTargetImplSurface&&) = delete;
    LoadTargetImplSurface& operator=(LoadTargetImplSurface&&) = delete;

    LoadSourcesResult stage(const ClusterSourceInput& source, SurfaceCatalogView catalog,
                            const o2::InteractionRecord& origin) override;
    void commit() noexcept override;

   private:
    SurfaceTrackingScratch& mLive;
    SurfaceTrackingScratch mStaged;
  };

  // One caller-supplied element of the ordered atomic-load transaction:
  // `source` is that participant's own per-event cluster/decoder input
  // (the "input/decoder" naming); `target` is that same participant's own
  // LoadTarget (never owned by this binding -- it must outlive the
  // loadEvent() call it is passed to).
  struct AtomicLoadBinding {
    ClusterSourceInput source;
    LoadTarget& target;
  };

  // Participant-count-generic atomic loading transaction (M2b). Stages the
  // shared normalized frame once, generically, over every binding's source
  // in one loadSources() call (already source-count-agnostic -- untouched
  // by this milestone), then stages every binding's own legacy backfill
  // via its LoadTarget, in `bindings` order. Only once every stage has
  // succeeded does this commit: `frame`'s normalized frame first, then
  // every binding's target, each a no-throw operation -- so once
  // committing starts, it cannot partially fail. On any decode,
  // allocation, validation, or staging failure -- normalized staging or
  // any one binding's stage() -- this returns immediately and `frame` plus
  // every binding's live scratch are left completely unchanged, since
  // nothing has been committed yet at that point. Never calls
  // TrackingEngine::executeEvent() itself; the caller applies its own
  // whole-event reset policy on a non-ok() result. Names none of ITS, MFT,
  // a layer count, or a fixed source position.
  static LoadSourcesResult loadEvent(TimeFrame& frame, gsl::span<const AtomicLoadBinding> bindings,
                                     SurfaceCatalogView catalog, const o2::InteractionRecord& origin);

  // Thin application-adapter wrapper over loadEvent(). The fixed ITS=0/MFT=1
  // position contract lives only here, never inside the generic transaction
  // or its target interface.
  static LoadSourcesResult loadITSAndMFT(TimeFrame& frame,
                                         SurfaceTrackingScratch& itsScratch,
                                         SurfaceTrackingScratch& mftScratch,
                                         const ClusterSourceInput& itsSource,
                                         const ClusterSourceInput& mftSource,
                                         SurfaceCatalogView catalog,
                                         const o2::InteractionRecord& origin);

  // Full shared-event reset: scratch state first (so no legacy cache can
  // outlive normalized measurements), then the common owner exactly once.
  static void resetITSAndMFTEvent(TimeFrame& frame,
                                  SurfaceTrackingScratch& itsScratch,
                                  SurfaceTrackingScratch& mftScratch) noexcept;
};

} // namespace o2::itsmft::tracking

#endif
