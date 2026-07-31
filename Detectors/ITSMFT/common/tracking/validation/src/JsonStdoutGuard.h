// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// Gate 4 B1 Slice 1 correction: with --format json,
// nominal-geometry-validator.cxx's stdout must be exactly one JSON
// document, on both the success and every typed-failure path. Without
// this guard, GeometryManager::loadGeometry()/GeometryTGeo's own
// FairLogger console sink (and GeometryTGeo's singleton destructor, which
// logs again during static teardown after main() returns) write [INFO]
// lines to stdout around the JSON payload -- observed directly in the
// Slice 1 real-geometry validation run.

#ifndef ALICEO2_ITSMFT_VALIDATION_JSONSTDOUTGUARD_H_
#define ALICEO2_ITSMFT_VALIDATION_JSONSTDOUTGUARD_H_

#include <string>

namespace o2::itsmft::validation
{

// Redirects the process's stdout file descriptor (fd 1) to stderr (fd 2)
// for the remainder of the process's life, starting at construction --
// deliberately never restored (see JsonStdoutGuard.cxx) -- so anything a
// workload writes via std::cout, C stdio, or a library's own console sink
// lands on stderr instead of stdout, visible rather than suppressed.
// Diagnostics are redirected, never discarded.
//
// emit() then writes exactly one payload line to the *real* original
// stdout, via a duplicated file descriptor captured before the redirect,
// bypassing std::cout entirely (now bound to stderr).
//
// If the underlying dup()/dup2() calls fail, ok() reports false and
// emit() falls back to writing via std::cout un-redirected, so a payload
// is still produced (degraded, not silently dropped).
class JsonStdoutGuard
{
 public:
  JsonStdoutGuard();
  ~JsonStdoutGuard();

  JsonStdoutGuard(const JsonStdoutGuard&) = delete;
  JsonStdoutGuard& operator=(const JsonStdoutGuard&) = delete;

  bool ok() const noexcept { return mOriginalStdoutFd >= 0; }

  void emit(const std::string& payload) const;

 private:
  int mOriginalStdoutFd{-1};
};

} // namespace o2::itsmft::validation

#endif /* ALICEO2_ITSMFT_VALIDATION_JSONSTDOUTGUARD_H_ */
