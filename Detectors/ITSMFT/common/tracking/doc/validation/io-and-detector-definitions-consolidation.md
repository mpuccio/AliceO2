# I/O and detector-definition consolidation

`IOUtils.h` now owns the atomic multi-source loading declarations and
`IOUtils.cxx` owns their implementation. The former loader header and source
were deleted; no forwarding identity remains.

`ITSMFTDetectorDefinitions.h` owns the ITS/MFT static surface catalogs and
their nominal material defaults. The former catalog and material-default
headers were deleted without changing catalog values or material constants.

Validation used the reusable macro-off build. The complete serial ITS/MFT
CTest suite passed 88/88 tests, in two serial invocations of 47 and 41 tests.
The tracking library and ITS, MFT, and combined CA workflows built cleanly.
`git diff --check` and the pinned `git clang-format --diff` gate passed.

No replay was run: this slice changes only source/header organization.
