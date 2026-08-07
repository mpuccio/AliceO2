// Gate 4 L9: standalone workflows compose the same non-owning generic
// components as the combined workflow; no interface or float failure bridge
// remains in the common-CA path.

#define BOOST_TEST_MODULE ITSMFT L9 standalone composition guard
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>

namespace
{
namespace fs = std::filesystem;

std::string readFile(const fs::path& path)
{
  std::ifstream input{path};
  return {std::istreambuf_iterator<char>{input}, {}};
}

bool isSource(const fs::path& path)
{
  static constexpr std::array<std::string_view, 3> suffixes{".h", ".cxx", ".txt"};
  const auto suffix = path.extension().string();
  for (const auto allowed : suffixes) {
    if (suffix == allowed) {
      return true;
    }
  }
  return false;
}

void checkNoLegacyInterface(const fs::path& root)
{
  BOOST_REQUIRE(fs::is_directory(root));
  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file() || !isSource(entry.path())) {
      continue;
    }
    const auto source = readFile(entry.path());
    BOOST_CHECK_MESSAGE(source.find("TrackingInterface") == std::string::npos,
                        "retired TrackingInterface spelling in " << entry.path());
    BOOST_CHECK_MESSAGE(source.find("kDroppedTimeFrameResult") == std::string::npos,
                        "retired float drop sentinel in " << entry.path());
    BOOST_CHECK_MESSAGE(source.find("isDroppedTimeFrame") == std::string::npos,
                        "retired float drop predicate in " << entry.path());
  }
}
} // namespace

BOOST_AUTO_TEST_CASE(CommonAndStandaloneSourcesHaveNoInterfaceBridge)
{
  const fs::path trackingRoot{fs::path{__FILE__}.parent_path().parent_path()};
  checkNoLegacyInterface(trackingRoot / "include");
  checkNoLegacyInterface(trackingRoot / "src");
  checkNoLegacyInterface(trackingRoot.parent_path().parent_path() / "ITS/workflow-ca/include");
  checkNoLegacyInterface(trackingRoot.parent_path().parent_path() / "ITS/workflow-ca/src");
  checkNoLegacyInterface(trackingRoot.parent_path().parent_path() / "MFT/workflow/include");
  checkNoLegacyInterface(trackingRoot.parent_path().parent_path() / "MFT/workflow/src");
  BOOST_CHECK(!fs::exists(trackingRoot / "include/ITSMFTTracking/TrackingInterface.h"));
  BOOST_CHECK(!fs::exists(trackingRoot / "src/TrackingInterface.cxx"));
}

BOOST_AUTO_TEST_CASE(StandaloneWorkflowsUseDirectComposition)
{
  const fs::path trackingRoot{fs::path{__FILE__}.parent_path().parent_path()};
  const auto itsRoot = trackingRoot.parent_path().parent_path() / "ITS/workflow-ca";
  const auto mftRoot = trackingRoot.parent_path().parent_path() / "MFT/workflow";
  for (const auto& root : {itsRoot, mftRoot}) {
    const auto header = readFile(root / "include" / (root.filename() == "workflow-ca" ? "ITSCAWorkflow/CATrackerSpec.h" : "MFTWorkflow/CATrackerSpec.h"));
    const auto source = readFile(root / "src/CATrackerSpec.cxx");
    const auto text = header + source;
    BOOST_CHECK(text.find("TimeFrame mFrame") != std::string::npos);
    BOOST_CHECK(source.find("MultiSourceTimeFrameLoader::load") != std::string::npos);
    BOOST_CHECK(source.find("mTracker->run(mFrame") != std::string::npos);
    BOOST_CHECK(source.find("mFrame.resetEvent()") != std::string::npos);
    BOOST_CHECK(source.find("stage") != std::string::npos);
  }
}
