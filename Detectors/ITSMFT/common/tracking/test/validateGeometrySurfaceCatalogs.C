// Compile-only ROOT wrapper for the shared validation implementation. The
// staged executable is the supported real-geometry validation protocol.

#include "GeometrySurfaceCatalogValidation.h"

void validateGeometrySurfaceCatalogs(const std::string& geometryPrefix)
{
  o2::itsmft::tracking::test::validateGeometrySurfaceCatalogs(geometryPrefix);
}
