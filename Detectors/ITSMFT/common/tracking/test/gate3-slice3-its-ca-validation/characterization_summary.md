# Gate 3 Slice 3: ITS common-CA characterization summary

Commit under test: `4a2bf6c3644deca2fe6cd75e31f691572a4d2a5d` ("ITSCAWorkflow:
add vertex-free ITS common-CA writer, driver, and o2-its-ca-tracker-workflow
executable (Slice 2C)"), branch `codex/itsmft-gate3-its-common-ca-validation`.
Package `daily-20260717-0700-local1` (pinned explicitly, requested ==
resolved). Fixture `pp-20ev-run303000-seed20260716-daily20260717`, checksums
verified before use. Full commands/paths in `manifest.json`.

**This is characterization, not a bitwise-parity gate**, per task scope.
No physics tolerance is invented; no equivalence to legacy real-vertex/CCDB
modes is claimed. The paired legacy-diamond replay in Sec.3 below is a
diagnostic, root-causing why the frozen legacy tracker's forced-diamond leg
yields zero tracks -- it is NOT an equivalent common-versus-legacy A/B
comparison: the frozen legacy static diamond has an unstamped vertex, while
the common-CA leg uses a common-only, ROF-local timestamped diamond. Do not
read Sec.3's table as evidence of equal constraints, parity, or comparison
validity between the two legs.

## 1. Build

`o2-its-ca-tracker-workflow`, `o2-its-cluster-reader-workflow`,
`o2-its-reco-workflow` (frozen legacy), `o2-mft-cluster-reader-workflow`,
`o2-mft-ca-tracker-workflow`, and the two dlopen'd DPL plugin libraries
(`libO2FrameworkAnalysisSupport`, `libO2FrameworkCCDBSupport`) all built
clean (1363/1363 ninja steps, exit 0) from this exact worktree/commit
against the pinned package.

## 2. ITS common-CA leg: namespace, output, no vertex output

- `ITSCommonCATrackerParam` namespace accepted via `--configKeyValues`
  (confirmed both positively -- the diamond replay ran to completion -- and
  negatively -- a legacy `ITSCATrackerParam.*` override is rejected before
  `updateFromString()`, see §4).
- `o2trac_its_ca.root` non-empty (54745 bytes, canonical-1), tree `o2sim`,
  branches `ITSTrack`/`ITSTrackClusIdx`/`ITSTracksROF`/`ITSTrackMCTruth`.
- No `ITSVertices`/`ITSVerticesROF` branch or `ITS/VERTICES` OutputSpec is
  ever produced by this workflow (confirmed by source inspection of
  `CATrackerSpec.cxx`/`TrackWriterSpec.cxx`, not just by absence in the
  output file) -- matches task requirement #4 exactly.

## 3. Diagnostic replay metrics (canonical-1, single-threaded) -- NOT an equivalent comparison

Same explicit `diamondPos=(0,0,0)`, `pvRes=0.05` cm CLI values passed to
both legs (`ITSCommonCATrackerParam.*` / `ITSCATrackerParam.*`
respectively -- same field names, different registered namespace), but this
does NOT make the two legs' resulting diamond vertices equivalent: the
frozen legacy leg's vertex carries no timestamp at all (root-caused below),
while the common-CA leg's diamond is a common-only, ROF-local timestamped
construction. The table below is a diagnostic pairing, not a
common-versus-legacy A/B comparison, and no claim of equal constraints or
comparison validity is made.

| Metric | ITS common-CA | ITS legacy (forced useDiamond=true) |
|---|---|---|
| Input clusters | 4103 | 4103 |
| Input ROFs | 2304 | 2304 |
| Output tracks | 203 | **0** |
| Clusters/track (mean/median/min/max) | 6.22 / 7 / 4 / 7 | n/a |
| chi2 (mean/median/min/max) | 6.81 / 5.02 / 0.263 / 60.5 | n/a |
| MC reconstructable | 142 | 142 |
| Matched | 135 | 0 |
| Efficiency | 0.951 | 0 |
| Fake tracks / rate | 14 / 0.069 | 0 / 0 |
| Clone tracks / rate | 0 / 0 | 0 / 0 |
| MC labels available/aligned | true / true | n/a (0 tracks) |
| Track content hash | `ee7f7c794d60f2362fd2564258b7887e` | n/a |

**Legacy leg produced zero tracks under this forced configuration --
root-caused, not just observed.** The legacy tracker log
(`its_legacy_diamond_replay.log`) shows: the *vertexer* stage (a separate,
always-running pass in legacy Sync mode, independent of
`ITSCATrackerParam.useDiamond`) found 434|363 tracklets and 5 vertices as
usual; but the *tracker* stage's own tracklet finding then reports **0**
tracklets found, immediately collapsing cell/neighbour/track finding to 0.
Tracing this in `Detectors/ITSMFT/ITS/tracking/src/TrackerTraits.cxx:61`:
the diamond constraint is represented as `Vertex diamondVert(Diamond,
DiamondCov, 1, 1.f)` -- a 4-argument constructor that never sets a
timestamp. `computeLayerTracklets()` gates every candidate cluster pair on
`ROFLookupTables::isVertexCompatible()`
(`Detectors/ITSMFT/ITS/tracking/include/ITStracking/ROFLookupTables.h:497-506`),
which compares the vertex's `getTimeStamp()` against each ROF's BC window
and rejects the pair if they don't overlap. A synthetic diamond vertex built
this way has no per-ROF timestamp and fails this check for every ROF in the
fixture, so tracklet finding at the tracker stage (as opposed to the
separate vertexer-stage tracklet pass) yields nothing tracker-wide,
regardless of `Diamond`/`DiamondCov` values. This is a legacy-code behavior
specific to forcing `useDiamond=true` onto Sync mode via `--configKeyValues`
-- not something this validation's scripts caused, and not something this
task is in scope to fix (no production code was touched). It is recorded
here as an honest characterization finding, not a parity result: **the
common-CA leg's diamond constraint mechanism instead derives a common-only,
ROF-local timestamp and does not have this timestamp gate, producing
sensible non-zero, non-degenerate output (203 tracks, 95.1% efficiency
against the same MC-reconstructable denominator) under the same
diamondPos/pvRes CLI values that leave the legacy leg's unstamped vertex
failing the timestamp-compatibility check for every ROF. The two legs'
diamond vertices are not equivalent constraints, and this pairing is
diagnostic, not a comparison of equal setups.**

## 4. Determinism (2 canonical single-threaded common-CA ITS replays)

`its-common-ca-canonical-1` and `its-common-ca-canonical-2`, same
configuration, fresh directories. **Every field in the extracted metrics
JSON was identical between the two runs, including the ordered
track-content hash `ee7f7c794d60f2362fd2564258b7887e`** (`diff` of the two
JSON files is empty). Raw `o2trac_its_ca.root` file bytes differ between
the two runs (`shasum -a 256`: `b17762d2...` vs `919a6cd9...`), confirmed
to be ROOT `TFile` UUID/timestamp metadata, not a content difference -- the
same distinction already established for the MFT common-CA baseline.
**No nondeterminism observed** in any emitted metric or the content hash.

## 5. Negative checks

Both run via `negative_checks_its_common_ca.sh`, invoking
`o2-its-ca-tracker-workflow` standalone (no piped input needed -- the
preflight fires inside `defineDataProcessing()`, before any device/input is
touched):

| Case | Exit | Fatal message present | Output file created |
|---|---|---|---|
| `useDiamond=false` (default) | 134 (SIGABRT via `LOGP(fatal,...)`) | yes: `"requires ITSCommonCATrackerParam.useDiamond=true"` | no |
| `--configKeyValues 'ITSCATrackerParam.useDiamond=true'` (legacy namespace) | 134 | yes: `"rejects legacy 'ITSCATrackerParam' --configKeyValues override"` | no |

Both cases fatal before configuration update/tracking/output, confirmed by
source (`ConfigPreflight.cxx`'s `applyConfigKeyValuesOrFatal()` runs
*before* `o2::conf::ConfigurableParam::updateFromString()` is ever called on
a rejected string) and by the absence of any `o2trac_its_ca.root` in either
case's directory.

## 6. MFT common-CA non-regression fence

`gate0-baseline/replay_tracking_common_ca.sh` +
`gate0-baseline/extract_metrics_common_ca.C`, unchanged, against the same
fixture:

| | Expected (per `itsmft-stageb-mft-normalized-refit` record) | Observed |
|---|---|---|
| Output tracks | 70 | 70 |
| Track content hash | `24737e73b7146bf3bd35a90a2517c527` | `24737e73b7146bf3bd35a90a2517c527` |

**Exact match. MFT non-regression fence held** -- this Gate 3 Slice 3
validation work did not disturb the MFT common-CA path.

## 7. Fixture integrity

`shasum -a 256 -c checksums.sha256` run against the fixture directory before
any replay: all 43 listed files OK. Fixture never regenerated, overwritten,
or chmod'd; all replay/metrics output went to fresh directories under
`O2-validation-artifacts/itsmft/gate3-slice3-its-ca-validation/`, never
inside the fixture directory.

## Known limits

1. **Legacy leg produces zero tracks under the forced diamond
   configuration** (§3) -- a legacy-code timestamp-compatibility gap when
   `useDiamond=true` is forced onto Sync mode via CLI, traced to source but
   not fixed (out of scope: no production code was touched). This means the
   paired-replay table above cannot show a legacy non-zero reference point
   for this exact configuration; the common-CA leg's own output (§2, §3) is
   independently non-empty and characterized regardless.
2. **Operational flakiness**: the legacy-diamond replay hung on 2 of 4
   attempts before completing (see `README.md`); root cause not identified,
   recorded as an operational note.
3. **No thread-count knob for `ITSCommonCATrackerParam`** (unlike MFT's
   `MFTCATrackerParam.nThreads`) -- "canonical single-threaded" here is
   simply the workflow's only mode, not a configured value; no parallel
   characterization run was possible for the ITS common-CA leg for this
   reason.
4. All CCDB-network-dependency, no-local-snapshot, and installed
   sim/digi-stack-provenance caveats already documented for the shared
   fixture in `gate0-baseline/README.md`/`manifest.json` apply identically
   here, since the fixture is reused unmodified.
5. Per task instruction #5, no physics tolerance/threshold is defined for
   the common-CA leg's own metrics (efficiency, fake rate, chi2 shape,
   etc.) -- the numbers in §3 are recorded as a labelled reference point for
   future comparison, not evaluated against a pass/fail bound.

**Not marking this gate complete** -- the integration owner reviews this
evidence, per `AgentCoordination.md` §12 ("An agent may not declare a gate
passed based only on compilation").
