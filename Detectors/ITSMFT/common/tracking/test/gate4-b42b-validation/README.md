# B4.2b CommonTrack publication A/B contract

The opt-in common-CA routes retain legacy publication as the default. ITS
CommonTrack output is required to equal legacy output field-for-field.

MFT `TrackMFT` persists forward states and covariances as `Double_t`, while
the fixed CommonTrack contract deliberately owns their canonical `float`
representation. The opt-in MFT route therefore has this explicit rule for
every serialized forward state, covariance, and fitted chi2 value:

```text
commonOutputDouble == double(static_cast<float>(legacyOutputDouble))
```

All other MFT products remain exact: track count and order, cluster IDs,
ROF payload and ranges, seed pattern/scalars, MC labels (including fake
bits), timestamps, and source associations. The temporary MFT publication
sidecar retains only publication facts not represented by CommonTrack,
including the float canonicalized `outParam` chi2 required by `TrackMFT`.

`compare_common_ca_output.C` enforces this contract on persisted ROOT
products. It compares all split track leaves, cluster-index vectors, seed
patterns, ROF fields, and optional MC labels; it reports the maximum
absolute and relative legacy-to-float deltas for the MFT projected fields.

The accepted replay baselines are:

- ITS legacy and opt-in: 203 tracks / `ee7f7c794d60f2362fd2564258b7887e`
- MFT legacy: 70 tracks / `24737e73b7146bf3bd35a90a2517c527`
- MFT opt-in float output: 70 tracks / `1a858f59a722891799c1dacfdaafdf71`
