// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// With --format json, nominal-geometry-validator.cxx must write exactly one
// JSON document on success and every typed-failure path. This guard keeps
// GeometryTGeo/FairLogger diagnostics, including late singleton-destructor
// logs, off stdout.

#ifndef ALICEO2_ITSMFT_VALIDATION_JSONSTDOUTGUARD_H_
#define ALICEO2_ITSMFT_VALIDATION_JSONSTDOUTGUARD_H_

#include <string>

namespace o2::itsmft::validation
{

// Redirects stdout (fd 1) to stderr (fd 2) for the rest of the process.
// Output from std::cout, C stdio, and library console sinks is redirected,
// not discarded; see JsonStdoutGuard.cxx for why stdout is not restored.
//
// emit() writes one payload line to the original stdout through a duplicated
// descriptor captured before redirection, bypassing std::cout.
//
// If dup()/dup2() fails, ok() is false and emit() falls back to std::cout;
// the payload is still produced, though diagnostics may share stdout.
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
