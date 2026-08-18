# TimeFrame measurement ownership

Status: historical validation. The current loader resets and fills the
configured TimeFrame directly; failures leave it empty. TimeFrame still owns
its per-surface measurements, while adapters retain raw inputs, publication
sidecars, and the storage behind non-owning timing views.

The reusable macro-off build completed after the follow-up correction. Its
serial ITS/MFT CTest gate passed 87/87 tests (split into 44 and 43 tests only
to stay below the command-runner window). Diff and formatting checks passed.

The canonical fixture manifest passed 43/43 checksums before and after the
replay campaign. Standalone ITS reproduced 184 tracks with hash
`e6da9d94faed581d5bff044993698e30`; standalone MFT reproduced 66 with hash
`32555b198d9b094f3f3600ec619cd2e2`. The combined replay reproduced ITS
unchanged and MFT's established 94-track combined product with hash
`96f4c632b7e0111501a63660774480ef`. ITS standalone and combined products
match field-for-field. MFT standalone and combined intentionally do not: the
combined graph uses the common scalar policy and its 94-track candidate was
already the accepted P3 combined baseline.

Artifacts are in `/private/tmp/itsmft-fold-measurement-store-replay/`.
