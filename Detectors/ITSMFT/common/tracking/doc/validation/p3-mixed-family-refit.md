# P3 mixed-family refit validation

## Contract exercised

`testDetectorAdapterCompatibility` builds a three-surface
Cylinder-to-Disk-to-Cylinder seed with source mapping `[0, 1, 0]`. The generic
refit succeeds through the shared propagator conversion path. The same fixture
with a foreign source at the Disk position fails before candidate mutation.

Existing `testCommonTrack` coverage retains the typed-output boundary: a track
with references spanning the requested detector/source is rejected as
`MixedDetector` without changing `TimeFrame` owners.

## Validation

The full affected build was rebuilt after the refit-function signature
change. The normal macro-on configuration remains blocked by the known
unregistered ROOT macro, so the configured macro-off suite ran serially:
**87/87** ITS/MFT tests passed with no `Not Run` entries. The focused mixed
refit and refit-boundary guards passed as part of that run.

The canonical fixture manifest verified **43/43** entries before and after
replay. Pinned standalone and disconnected-combined replays preserved the
parent candidate products field-for-field:

| Output | Tracks | Content hash |
| --- | ---: | --- |
| ITS standalone | 184 | `e6da9d94faed581d5bff044993698e30` |
| MFT standalone | 66 | `32555b198d9b094f3f3600ec619cd2e2` |
| ITS combined | 184 | `e6da9d94faed581d5bff044993698e30` |
| MFT combined | 94 | `96f4c632b7e0111501a63660774480ef` |

The retained parent replay comparison reported exact initialized ITS output
and exact MFT output under its established float-projection contract (2,904
standalone and 4,136 combined projected values, each with zero delta).
Mixed tracking remains synthetic-only in this slice; detector-specific DPL
publication is deliberately not attempted.

No accelerator result is claimed without a CUDA or HIP toolchain.
