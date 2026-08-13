// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include "JsonStdoutGuard.h"

#include <iostream>

#include <unistd.h>

namespace o2::itsmft::validation
{

JsonStdoutGuard::JsonStdoutGuard()
{
  const int duplicated = ::dup(STDOUT_FILENO);
  if (duplicated < 0) {
    return;
  }
  if (::dup2(STDERR_FILENO, STDOUT_FILENO) < 0) {
    ::close(duplicated);
    return;
  }
  mOriginalStdoutFd = duplicated;
}

JsonStdoutGuard::~JsonStdoutGuard()
{
  // Keep fd 1 bound to stderr: GeometryTGeo's singleton destructor logs
  // during static teardown after guarded call sites have emitted JSON.
  // Restoring stdout would send that late noise to the real stdout.
  // Close only the duplicated descriptor owned by this instance.
  if (mOriginalStdoutFd >= 0) {
    ::close(mOriginalStdoutFd);
  }
}

void JsonStdoutGuard::emit(const std::string& payload) const
{
  const std::string line = payload + "\n";
  if (!ok()) {
    std::cout << line;
    std::cout.flush();
    return;
  }
  const char* data = line.data();
  size_t remaining = line.size();
  while (remaining > 0) {
    const ssize_t written = ::write(mOriginalStdoutFd, data, remaining);
    if (written <= 0) {
      break;
    }
    data += written;
    remaining -= static_cast<size_t>(written);
  }
}

} // namespace o2::itsmft::validation
