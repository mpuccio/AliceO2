/// \file CheckVertexer.C
/// \brief Simple macro to check ITSU tracks

#if !defined(__CLING__) || defined(__ROOTCLING__)
  #include <array>
  #include <fstream>
  #include <sstream>
  #include <string>
  #include <vector>
  #include "TString.h"
  #include "ITSReconstruction/CA/Event.h"
  #include "ITSReconstruction/CA/vertexer/Vertexer.h"
#endif


constexpr int PrimaryVertexLayerId = -1;
constexpr int EventLabelsSeparator = -1;
using namespace o2::ITS::CA;

std::vector<o2::ITS::CA::Event> loadEventData(const std::string& fileName)
{
  std::vector<Event> events;
  std::ifstream inputStream;
  std::string line, unusedVariable;
  int layerId, monteCarlo;
  int clusterId = -1;
  float xCoordinate, yCoordinate, zCoordinate;
  inputStream.open(fileName);
  while (std::getline(inputStream, line)) {
    std::istringstream inputStringStream(line);
    if (inputStringStream >> layerId >> xCoordinate >> yCoordinate >> zCoordinate) {
      if (layerId == PrimaryVertexLayerId) {
        if (clusterId != 0) {
          events.emplace_back(events.size());
        }
        events.back().addPrimaryVertex(xCoordinate, yCoordinate, zCoordinate);
        clusterId = 0;
      } else {
        if (inputStringStream >> unusedVariable >> unusedVariable >> unusedVariable >> unusedVariable >> monteCarlo) {
          events.back().pushClusterToLayer(layerId, clusterId, xCoordinate, yCoordinate, zCoordinate, monteCarlo);
          ++clusterId;
        }
      }
    }
  }
  Printf(" evt vector size: %lu\n", events.size());
  return events;
}


void CheckVertexer(const std::string& fname ="data.txt")
{
  std::vector<Event> events = loadEventData(fname);
  Vertexer vertexer(events.back());
  // vertexer.debugVertexerData();
  vertexer.initialize();
  Printf(" ---------------- initialised ---------------");
  vertexer.debugVertexerData();
  vertexer.printIndexTables();
  vertexer.computeTriplets();
}
