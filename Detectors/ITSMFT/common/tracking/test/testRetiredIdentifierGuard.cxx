// Source guard for the retired temporary dispatch vocabulary.

#define BOOST_TEST_MODULE ITSMFT RetiredIdentifierGuard
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace
{
namespace fs = std::filesystem;

fs::path sourceRoot()
{
  auto root = fs::path{__FILE__};
  for (int level = 0; level < 4; ++level) {
    root = root.parent_path();
  }
  return root;
}

std::string read(const fs::path& path)
{
  std::ifstream input{path};
  BOOST_REQUIRE_MESSAGE(input.good(), "cannot inspect " << path.string());
  return {std::istreambuf_iterator<char>{input}, {}};
}

bool sourceFile(const fs::path& path)
{
  const auto extension = path.extension().string();
  return extension == ".h" || extension == ".hpp" || extension == ".c" || extension == ".cc" ||
         extension == ".cxx" || extension == ".cpp" || extension == ".cu" || extension == ".cmake" ||
         extension == ".txt";
}

BOOST_AUTO_TEST_CASE(NoRetiredDispatchVocabularyRemainsInITSMFTSources)
{
  // Assemble spellings from fragments so this guard cannot exempt itself.
  const std::string retiredFamily = std::string{"Trans"} + "ition" + "Policy";
  const std::string retiredTag = retiredFamily + "Tag";
  const std::string retiredTemplateName = std::string{"For"} + "Policy";
  const std::string retiredLowerName = std::string{"transition"} + "Policy";
  const std::string retiredBindingName = std::string{"bindTransition"} + "Policy" + "Params";
  const std::vector<std::string> forbidden{retiredFamily, retiredTag, retiredTemplateName,
                                           retiredLowerName, retiredBindingName};
  const std::vector<std::string> retiredFiles{retiredFamily + ".h", retiredFamily + "Binding.h",
                                              retiredFamily + "Dispatch.h", retiredFamily + "Operations.h",
                                              retiredFamily + "State.h"};

  const auto root = sourceRoot();
  BOOST_REQUIRE_MESSAGE(fs::is_directory(root), "cannot find " << root.string());
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file() || !sourceFile(entry.path())) {
      continue;
    }
    const auto name = entry.path().filename().string();
    for (const auto& filename : retiredFiles) {
      BOOST_CHECK_MESSAGE(name != filename, "retired header filename remains: " << entry.path().string());
    }
    const auto text = read(entry.path());
    for (const auto& spelling : forbidden) {
      BOOST_CHECK_MESSAGE(text.find(spelling) == std::string::npos,
                          "retired spelling remains in " << entry.path().string() << ": " << spelling);
    }
  }
}

BOOST_AUTO_TEST_CASE(NoRetiredPrimitiveHeadersRemain)
{
  const auto root = sourceRoot();
  BOOST_REQUIRE_MESSAGE(fs::is_directory(root), "cannot find " << root.string());

  const auto includeDirectory = root / "Detectors" / "ITSMFT" / "common" / "tracking" / "include" / "ITSMFTTracking";
  const std::vector<std::string> retiredHeaders{
    std::string{"SurfaceMeasurement"} + "Index.h",
    std::string{"Layer"} + "Mask.h",
    std::string{"ClockTimingPublicationView"} + ".h",
  };
  for (const auto& header : retiredHeaders) {
    BOOST_CHECK_MESSAGE(!fs::exists(includeDirectory / header),
                        "retired primitive header remains: " << (includeDirectory / header).string());
  }
}

} // namespace
