#define BOOST_TEST_MODULE ITSMFT measurement ownership guard
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include <boost/test/unit_test.hpp>

namespace
{
std::string readFile(const std::filesystem::path& path)
{
  std::ifstream input{path};
  return {std::istreambuf_iterator<char>{input}, {}};
}
} // namespace

BOOST_AUTO_TEST_CASE(ProductionUsesTimeFrameMeasurementsAndDescriptorBackfill)
{
  const auto root = std::filesystem::path{__FILE__}.parent_path().parent_path();
  for (const auto& relative : {"include/ITSMFTTracking", "src"}) {
    for (const auto& entry : std::filesystem::recursive_directory_iterator{root / relative}) {
      if (!entry.is_regular_file() || (entry.path().extension() != ".h" && entry.path().extension() != ".cxx")) {
        continue;
      }
      const auto source = readFile(entry.path());
      BOOST_CHECK(source.find("MultiSourceFrame") == std::string::npos);
      BOOST_CHECK(source.find("SourceMetadata") == std::string::npos);
    }
  }

  const auto loader = readFile(root / "src/MultiSourceTimeFrameLoader.cxx");
  BOOST_CHECK(loader.find("owner->detector") == std::string::npos);
  BOOST_CHECK(loader.find("catalog.getSurface(surface).kind") != std::string::npos);
}
