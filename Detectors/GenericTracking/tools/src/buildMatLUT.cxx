// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file buildMatLUT.cxx
/// \brief Build the tracking material LUT from an o2sim geometry for the MFT
///        or the ALICE 3 barrel tracker.

#include "GenericTracking/Tools/Material.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Detectors/TRKFT3.h"

#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/MatLayerCylSet.h"

#include <TGeoNavigator.h>

#include <boost/program_options.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace bpo = boost::program_options;

namespace
{

struct LUTLayer {
  float rMin;
  float rMax;
  float zMax;
  float dz;
  float drPhi;
};

/// The MFT's region, which O2's MatLUT leaves out: its layers below r = 5 cm
/// end at |z| = 20 cm, the beam pipe's at 50 cm,
/// and the MFT's faces lie at z -45 to -78 cm, r 2.3 to 17.3 cm. An O2 layer
/// spans z symmetrically, so these also cover the ITS and whatever a track
/// crosses from the vertex. The pipe's wall gets its own band; from 2.2 cm
/// out, bands of 1 cm with 0.25 cm along z keep each disc's material between
/// its two faces (1.42 cm apart) and out of the gaps between discs, and
/// 0.5 cm in r-phi resolves the half-discs' ladders, which cover each face
/// only in part.
inline constexpr float kMftZMax = 82.F;
inline constexpr auto kMftLayers = [] {
  std::array<LUTLayer, 20> layers{};
  layers[0] = {.rMin = 0.F, .rMax = 1.8F, .zMax = kMftZMax, .dz = 2.F, .drPhi = 2.F};
  layers[1] = {.rMin = 1.8F, .rMax = 1.95F, .zMax = kMftZMax, .dz = 0.5F, .drPhi = 0.5F};
  layers[2] = {.rMin = 1.95F, .rMax = 2.2F, .zMax = kMftZMax, .dz = 0.5F, .drPhi = 1.F};
  for (std::size_t i = 0; i < 16; ++i) {
    const float r = i == 0 ? 2.2F : 2.F + static_cast<float>(i);
    layers[3 + i] = {.rMin = r, .rMax = 3.F + static_cast<float>(i), .zMax = kMftZMax, .dz = 0.25F, .drPhi = 0.5F};
  }
  layers[19] = {.rMin = 18.F, .rMax = 20.F, .zMax = kMftZMax, .dz = 0.5F, .drPhi = 1.F};
  return layers;
}();

/// The ALICE 3 barrel tracker's region.
///
/// Its radial bands are the layout's own: one over each layer's chips -- its
/// silicon envelope and the metal stack on one face of its sensors
/// (`tools::withStacks`) -- and one over each gap between layers, as measured
/// from the loaded geometry. A LUT layer holds its material evenly across its
/// radial span, so only bands that end where the layers end keep a layer's
/// material apart from the gap's -- which is how the layout splits it, the
/// gap's going to its passive scatterer (TRKFT3.h). A band over the silicon
/// alone would spread a VD layer's stack, most of its x/X0, across the gap to
/// the next layer, and the fit would take only part of its kink at that layer.
/// A band over a layer of flat staves is as deep as their hull (0.5-1.4 cm: a
/// flat stave's edges lie further out than its middle), so it is cut into
/// equal parts no thinner than a chip (its silicon and stack,
/// \a chipThickness): at one azimuth a stave is a chip's thickness deep (and
/// the parts' cells narrow enough in r-phi that a tilted stave stays within
/// two parts), and a track at a small angle to it (|eta| 2.6, 8 deg) turns a
/// radial spread into ten times that along z. Along z, 1 cm, to |z| = 225 cm,
/// which holds the barrel (OT to 128.4 cm) and FT3's discs (to 220.6 cm): a
/// disc lies across the radial bands its radii span, its silicon and supports
/// in the 1 cm cells of its z. In r-phi, 0.05 cm for the VD's layers and the
/// gaps between them, 0.2 cm in the gap from the VD to the ML, 0.3 cm for a
/// stave layer, which resolves the staves' overlaps, and 0.5 cm in the other
/// gaps; but every band the petals' radial side walls cross (\a walls,
/// measured from the geometry) has cells half as wide as a wall is thick. A
/// wider cell spreads the wall over every track through it, and the tracks
/// through the petal gaps read a multiple of the geometry's x/X0. A wall lies
/// across two or three cells whatever their phase (bins start at phi 0); at
/// half its thickness the spread beyond its faces is at most half of it.
///
/// A gap is first cut at the faces of the cylinders around the beam line in it
/// (\a faces, measured from the geometry: the iris, the cold plates, the VD's
/// outer wall, the beam pipe, the service and support shells), so that each
/// keeps its material at its own radius: a band from the beam line would put
/// half of the iris wall at the vertex end of the step to VD L0, a band across
/// the middle of the step from L2 to the ML half of the VD's outer wall at L2.
/// Each part is then cut into equal bands of at most 1 cm in r. Beyond the
/// barrel's length the gaps hold FT3's supports, thin cylinders along z at
/// r 39.4, 71, 74.5 and 78-80 cm; a band as wide as its gap would spread each
/// over the gap in r, onto forward tracks that never reach it. A track
/// crossing the whole gap, as a barrel track does, reads the same either way.
inline constexpr float kTrkZMax = 225.F;
inline constexpr float kTrkRMax = 83.F; ///< beyond the OT's carbon support, and FT3's discs (69.2 cm)
inline constexpr float kTrkGapDr = 1.F; ///< the widest band in a gap

/// The bands of the barrel's layers, \a envelopes with their stacks; the discs' follow none.
std::vector<LUTLayer> trkLayers(std::span<const o2::gtrk::LayerEnvelope> envelopes, const o2::gtrk::tools::SideWalls& walls, std::span<const double> faces, std::span<const double> chipThickness, float rMaxCM, float zMaxCM)
{
  std::vector<LUTLayer> layers;
  const auto cell = [&](float inner, float outer, float drPhi) {
    const auto low = static_cast<double>(inner);
    const auto high = static_cast<double>(outer);
    const bool crossed = walls.thickness > 0. && low < walls.rMax && high > walls.rMin;
    return crossed ? std::min(drPhi, static_cast<float>(0.5 * walls.thickness)) : drPhi;
  };
  const auto addPart = [&](float inner, float outer, float drPhi) {
    const auto bands = static_cast<int>(std::ceil((outer - inner) / kTrkGapDr));
    for (int band = 0; band < bands; ++band) {
      const float low = band == 0 ? inner : inner + ((outer - inner) * static_cast<float>(band) / static_cast<float>(bands));
      const float high = band + 1 == bands ? outer : inner + ((outer - inner) * static_cast<float>(band + 1) / static_cast<float>(bands));
      layers.push_back({.rMin = low, .rMax = high, .zMax = zMaxCM, .dz = 1.F, .drPhi = cell(low, high, drPhi)});
    }
  };
  // a face closer than this to a band's edge makes no band of its own
  constexpr float kFaceResolution = 1e-4F;
  const auto addGap = [&](float inner, float outer, float drPhi) {
    float low = inner;
    for (const double face : faces) {
      const auto f = static_cast<float>(face);
      if (f > low + kFaceResolution && f < outer - kFaceResolution) {
        addPart(low, f, drPhi);
        low = f;
      }
    }
    addPart(low, outer, drPhi);
  };
  // the VD's layers, the layout's first barrel layers
  constexpr std::size_t kVDLayers = o2::gtrk::detectors::trkft3::kVDRadius.size();
  float previous = 0.F;
  for (std::size_t layer = 0; layer < envelopes.size(); ++layer) {
    const bool vd = layer < kVDLayers;
    const auto rMin = static_cast<float>(envelopes[layer].rMin);
    const auto rMax = static_cast<float>(envelopes[layer].rMax);
    float gapPhi = 0.5F; // the gaps between the ML's and the OT's layers
    if (vd) {
      gapPhi = 0.05F; // to and between the VD's layers
    } else if (layer == kVDLayers) {
      gapPhi = 0.2F; // from the VD to the ML
    }
    addGap(previous, rMin, gapPhi);
    // equal parts no thinner than a chip; a flat stave meets the radial
    // direction at up to alpha (cos alpha = rMin / rMax), so over a cell of
    // width w its radius changes by w tan(alpha): cells no wider than half a
    // chip's thickness over tan(alpha) keep it within a part and a half, as
    // the side walls' cells are half a wall
    const auto thickness = static_cast<float>(chipThickness[layer]);
    const int parts = std::max(1, static_cast<int>(std::floor((rMax - rMin) / thickness)));
    float drPhi = vd ? 0.05F : 0.3F;
    if (parts > 1) {
      const float tanAlpha = std::sqrt((rMax * rMax) - (rMin * rMin)) / rMin;
      drPhi = std::min(drPhi, 0.5F * thickness / tanAlpha);
    }
    for (int part = 0; part < parts; ++part) {
      const float low = part == 0 ? rMin : rMin + ((rMax - rMin) * static_cast<float>(part) / static_cast<float>(parts));
      const float high = part + 1 == parts ? rMax : rMin + ((rMax - rMin) * static_cast<float>(part + 1) / static_cast<float>(parts));
      layers.push_back({.rMin = low, .rMax = high, .zMax = zMaxCM, .dz = 1.F, .drPhi = cell(low, high, drPhi)});
    }
    previous = rMax;
  }
  addGap(previous, rMaxCM, 0.5F);
  return layers;
}

/// One cell's material from the loaded geometry: \a trials x \a trials lines
/// across it, \a trials azimuths times \a trials offsets along z. Each line
/// runs along the cell's r-z diagonal, from (rMin, zLow + s) to
/// (rMax, zHigh + s), and wraps at the cell's upper z face back to its lower
/// one, so it stays inside the cell. A plane at one z and a cylinder at one r
/// inside the cell are each crossed exactly once by every line, at the share
/// of its length their thickness takes of the cell's dz or dr: the cell gets
/// their volume average whatever the plane's z. O2's fill
/// (MatLayerCyl::populateFromTGeo) samples along nearly radial lines at
/// \a trials fixed z, which a thin plane such as a disc either lies along or
/// misses. The offsets s are the
/// centres of \a trials equal parts of dz; the cell keeps O2's mean density
/// and mean 1/X0 along the lines.
void populateCell(o2::base::MatLayerCyl& layer, int ip, int iz, int trials, TGeoNavigator* nav)
{
  const double rMin = layer.getRMin();
  const double rMax = layer.getRMax();
  const double dr = rMax - rMin;
  const double zLow = layer.getZBinMin(iz);
  const double dz = layer.getDZ();
  const double phiLow = layer.getPhiBinMin(ip);
  const double dPhi = layer.getDPhi();
  double rhoLength = 0.;
  double x2X0 = 0.;
  double length = 0.;
  const auto add = [&](double cs, double sn, double r0, double z0, double r1, double z1) {
    const auto bud = o2::base::GeometryManager::meanMaterialBudget(static_cast<float>(r0 * cs), static_cast<float>(r0 * sn), static_cast<float>(z0), static_cast<float>(r1 * cs), static_cast<float>(r1 * sn), static_cast<float>(z1), nav);
    const double segment = bud.length;
    if (segment > 0.) {
      rhoLength += segment * static_cast<double>(bud.meanRho);
      x2X0 += static_cast<double>(bud.meanX2X0); // the segment's x/X0
      length += segment;
    }
  };
  for (int iPhi = 0; iPhi < trials; ++iPhi) {
    const double phi = phiLow + ((iPhi + 0.5) * dPhi / trials);
    const double cs = std::cos(phi);
    const double sn = std::sin(phi);
    for (int iShift = 0; iShift < trials; ++iShift) {
      const double shift = (iShift + 0.5) / trials;
      // The line reaches the upper z face at this radius, and goes on from the lower one.
      const double rWrap = rMin + ((1. - shift) * dr);
      add(cs, sn, rMin, zLow + (shift * dz), rWrap, zLow + dz);
      add(cs, sn, rWrap, zLow, rMax, zLow + (shift * dz));
    }
  }
  if (length > 0.) {
    auto& cell = layer.getCellPhiBin(ip, iz);
    cell.meanRho = static_cast<float>(rhoLength / length);
    cell.meanX2X0 = static_cast<float>(x2X0 / length); // the mean 1/X0, as O2 stores it
  }
}

/// Every cell of \a lut by populateCell, on \a threads threads, each with its
/// own TGeo navigator (`tools::forEachOnNavigators`).
void populate(o2::base::MatLayerCylSet& lut, int trials, int threads)
{
  struct Cell {
    int layer;
    int ip;
    int iz;
  };
  std::vector<Cell> cells;
  for (int layer = 0; layer < lut.getNLayers(); ++layer) {
    const auto& cyl = lut.getLayer(layer);
    for (int iz = 0; iz < cyl.getNZBins(); ++iz) {
      for (int ip = 0; ip < cyl.getNPhiBins(); ++ip) {
        cells.push_back({.layer = layer, .ip = ip, .iz = iz});
      }
    }
  }
  o2::gtrk::tools::forEachOnNavigators(cells.size(), threads, [&](std::size_t i, TGeoNavigator* nav) { populateCell(lut.getLayer(cells[i].layer), cells[i].ip, cells[i].iz, trials, nav); });
  lut.finalizeStructures();
}

/// The TRK's bands from the loaded geometry, its stacks, side walls and
/// cylinder faces printed (the LUT's provenance in its log).
std::vector<LUTLayer> trkBands(int threads, float rMaxCM, float zMaxCM)
{
  using TRKFT3 = o2::gtrk::detectors::TRKFT3;
  const o2::gtrk::Geometry<TRKFT3> geometry(o2::gtrk::alice::placeChecked<TRKFT3>());
  constexpr int kBarrel = o2::gtrk::Detector<TRKFT3>::nBarrelLayers;
  const auto specs = o2::gtrk::Geometry<TRKFT3>::layerSpecs();
  const auto stack = o2::gtrk::tools::trkStackThickness(geometry.sensors(), specs, kBarrel);
  const auto chips = o2::gtrk::tools::withStacks(o2::gtrk::PlacedEnvelopes<TRKFT3>::all(), geometry.sensors(), specs, stack);
  std::cout << "metal stacks on the silicon, inner/outer face [um]:";
  for (const auto& t : stack) {
    std::cout << ' ' << t.inner * 1e4 << '/' << t.outer * 1e4;
  }
  std::cout << '\n';
  const auto walls = o2::gtrk::tools::trkSideWalls();
  std::cout << "VD side walls: " << walls.thickness * 1e4 << " um thick, r " << walls.rMin << " .. " << walls.rMax << " cm\n";
  const auto faces = o2::gtrk::tools::beamLineCylinderFaces(rMaxCM, zMaxCM, 1., threads);
  std::cout << "cylinder faces around the beam line [cm]:";
  for (const double face : faces) {
    std::cout << ' ' << face;
  }
  std::cout << '\n';
  std::vector<double> chipThickness;
  for (int layer = 0; layer < kBarrel; ++layer) {
    const auto l = static_cast<std::size_t>(layer);
    chipThickness.push_back((2. * specs[l].halfThickness) + stack[l].inner + stack[l].outer);
  }
  return trkLayers(std::span{chips}.first(static_cast<std::size_t>(kBarrel)), walls, faces, chipThickness, rMaxCM, zMaxCM);
}

int run(int argc, char** argv)
{
  std::string output;
  std::string geometryFile;
  std::string detector;
  float trkRMaxCM = kTrkRMax;
  float trkZMaxCM = kTrkZMax;
  int trials = 10;
  int threads = -1;
  float maxRelativeDifference = 0.05F;

  bpo::options_description opts("Build the GTRK material-budget LUT");
  opts.add_options()                                                                                           //
    ("help,h", "print this help")                                                                              //
    ("output,o", bpo::value(&output)->default_value("gtrk_matlut.root"), "output LUT file")                    //
    ("trials,n", bpo::value(&trials)->default_value(10), "lines per cell axis")                                //
    ("detector,d", bpo::value(&detector)->required(),                                                          //
     "MFT: the MFT's region; TRKFT3: the ALICE 3 tracker's, from --geometry (the ITS takes O2's CCDB MatLUT)") //
    ("geometry,g", bpo::value(&geometryFile)->required(),                                                      //
     "an o2sim geometry file (o2sim_geometry-aligned.root)")                                                   //
    ("threads,j", bpo::value(&threads)->default_value(-1),                                                     //
     "fill threads; negative uses NTHREADS_MATBUD")                                                            //
    ("trk-rmax-cm", bpo::value(&trkRMaxCM)->default_value(kTrkRMax),                                           //
     "TRKFT3 radial coverage in cm; extend beyond outer TOF for TOF matching")                                 //
    ("trk-zmax-cm", bpo::value(&trkZMaxCM)->default_value(kTrkZMax),                                           //
     "TRKFT3 longitudinal coverage in cm; extend to the TOF acceptance")                                       //
    ("max-relative-difference", bpo::value(&maxRelativeDifference)->default_value(0.05F),                      //
     "maximum relative difference when merging adjacent phi slices");

  bpo::variables_map vm;
  try {
    bpo::store(bpo::parse_command_line(argc, argv, opts), vm);
    if (vm.contains("help")) {
      std::cout << opts << '\n';
      return 0;
    }
    bpo::notify(vm);
  } catch (const bpo::error& e) {
    std::cerr << "error: " << e.what() << "\n\n"
              << opts << '\n';
    return 1;
  }
  if (detector != "MFT" && detector != "TRKFT3") {
    std::cerr << "error: --detector is MFT or TRKFT3" << (detector == "ITS" ? " (the ITS takes O2's CCDB MatLUT)" : "") << "\n\n"
              << opts << '\n';
    return 1;
  }
  if (!std::isfinite(trkRMaxCM) || !std::isfinite(trkZMaxCM) || trkRMaxCM < kTrkRMax || trkZMaxCM < kTrkZMax) {
    throw std::invalid_argument("TRKFT3 material coverage can only extend the default tracker region");
  }
  threads = o2::gtrk::tools::fillThreads(threads);

  o2::base::GeometryManager::loadGeometry(geometryFile);
  std::span<const LUTLayer> layers;
  std::vector<LUTLayer> trk;
  float radialCheck = 0.F;
  if (detector == "MFT") {
    layers = kMftLayers;
    radialCheck = 20.F;
  } else {
    trk = trkBands(threads, trkRMaxCM, trkZMaxCM);
    layers = trk;
    radialCheck = trkRMaxCM;
  }

  o2::base::MatLayerCylSet lut;
  std::cout << "material LUT layers:\n";
  for (std::size_t i = 0; i < layers.size(); ++i) {
    const auto& layer = layers[i];
    std::cout << "  " << std::setw(2) << i << ": r " << std::setw(4) << layer.rMin << " .. " << std::setw(4) << layer.rMax << " cm, |z| < " << layer.zMax << " cm, dz " << layer.dz << " cm, drphi " << layer.drPhi << " cm\n";
    lut.addLayer(layer.rMin, layer.rMax, layer.zMax, layer.dz, layer.drPhi);
  }

  std::cout << "\npopulating from TGeo with " << trials << " x " << trials << " diagonals per cell\n";
  populate(lut, trials, threads);
  lut.optimizePhiSlices(maxRelativeDifference);
  lut.flatten();
  lut.writeToFile(output);

  const std::unique_ptr<o2::base::MatLayerCylSet> loaded{o2::base::MatLayerCylSet::loadFromFile(output)};
  if (!loaded) {
    std::cerr << "failed to load the written LUT " << output << '\n';
    return 1;
  }
  const auto radial = loaded->getMatBudget(0.F, 0.F, 0.F, radialCheck, 0.F, 0.F);

  std::cout << "\nwrote " << output << '\n'
            << "  flattened size  " << static_cast<double>(lut.getFlatBufferSize()) / (1024. * 1024.) << " MiB\n"
            << "  radial check    x/X0 = " << radial.meanX2X0 << ", rho*L = " << radial.getXRho() << " g/cm^2 over " << radial.length << " cm\n";
  return 0;
}

} // namespace

int main(int argc, char** argv)
{
  try {
    return run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
