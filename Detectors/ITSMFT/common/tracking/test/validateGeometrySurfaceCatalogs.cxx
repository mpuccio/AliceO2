// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

#include <exception>
#include <iostream>

#include "GeometrySurfaceCatalogValidation.h"

int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cerr << "usage: " << argv[0] << " <geometry-prefix>\n";
    return 2;
  }
  try {
    o2::itsmft::tracking::test::validateGeometrySurfaceCatalogs(argv[1]);
  } catch (const std::exception& error) {
    std::cerr << "geometry surface catalog validation failed: " << error.what() << '\n';
    return 1;
  } catch (...) {
    std::cerr << "geometry surface catalog validation failed: unknown exception\n";
    return 1;
  }
  return 0;
}
