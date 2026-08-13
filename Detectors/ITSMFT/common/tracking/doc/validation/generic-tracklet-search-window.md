# Generic tracklet search window validation

- Date: 2026-08-13
- Package: `daily-20260717-0700-local1`
- Source/build: `triplet-tracking-rnd-scratch` and its reusable build
- Artifacts: `/private/tmp/tracklet-window-replay`

## Change under characterization

The two family-shaped search-window structs and their variant were replaced by
one data-only `TrackletSearchWindow`. Both leaves now fill a fixed
two-coordinate prediction and packed covariance. The cylinder target is the
mean radial coordinate of its surface interval; the interval is treated as a
uniform uncertainty and contributes `tanLambda^2 DeltaR^2 / 12` to the
longitudinal variance. Its periodic azimuth residual is wrapped before the
same strict covariance-normalized ellipse used for disk `(x,y)` residuals.

This is an approved physics migration. It is not a claim of equivalence to
the former cylinder rectangular `z`/`phi` selection.

## Validation

Focused numerical tests cover target-radius mean projection, radial-spread
variance, wrapped phi residuals, disk projected coordinates/covariance, and
the data-only-window/source guard. The full serial ITS/MFT CTest gate passed.
Whitespace and pinned formatting checks passed.

The fixed fixture checksum manifest passed 43/43 before and after replay.
Canonical standalone and opt-in combined replays completed with these results:

| Product | Tracks | Content hash |
|---|---:|---|
| ITS standalone | 189 | `6343211326990c75370a76b06aad5840` |
| ITS combined | 189 | `6343211326990c75370a76b06aad5840` |
| MFT standalone | 66 | `32555b198d9b094f3f3600ec619cd2e2` |
| MFT combined | 94 | `96f4c632b7e0111501a63660774480ef` |

ITS standalone and combined are byte-identical on this fixture. MFT remains
unchanged from the preceding candidate; its standalone/combined population
difference is the established common-policy composition effect. The new
cylinder ellipse changes ITS from the preceding 184-track candidate to 189.
