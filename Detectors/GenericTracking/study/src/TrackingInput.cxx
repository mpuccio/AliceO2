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

#include "GenericTracking/Study/TrackingInput.h"

#include "GenericTracking/Study/Digest.h"
#include "GenericTracking/Workflow/RecoParam.h"

#include "CommonUtils/ConfigurableParam.h"

#include <TClass.h>
#include <TDataMember.h>
#include <TDataType.h>
#include <TList.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstring>
#include <format>
#include <fstream>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

namespace
{

constexpr std::array<char, 8> kMagic{'G', 'T', 'R', 'K', 'T', 'I', 'N', '1'};

/// Bumped whenever the layout below or a column's type changes.
constexpr std::uint32_t kVersion = 9U;

/// The header: the blocks' element counts and the beam line the rows were
/// placed about. The blocks follow: the readout frames, every measurement
/// column in schema order, the grid offsets, the vertices, the measurement
/// depths.
struct Header {
  std::array<char, 8> magic = kMagic;
  std::uint32_t version = kVersion;
  std::uint32_t columns = kNMeasurementColumns;
  std::uint64_t rowBytes = kMeasurementRowBytes;
  std::uint64_t rows = 0U;
  std::uint64_t offsets = 0U;
  std::uint64_t vertices = 0U;
  std::uint64_t frameLayers = 0U;
  std::uint64_t depthLayers = 0U;
  double beamX = 0.;
  double beamY = 0.;
  double beamZ = 0.;
  double beamSigmaX = 0.;
  double beamSigmaY = 0.;
  double beamSigmaZ = 0.;
};
static_assert(std::is_trivially_copyable_v<Header> && sizeof(Header) == 112U, "the header has no padding");
static_assert(std::is_trivially_copyable_v<VertexSeed> && std::is_trivially_copyable_v<FrameTiming>);

constexpr std::array<char, 8> kTruthMagic{'G', 'T', 'R', 'K', 'T', 'R', 'U', '1'};

/// Bumped whenever the truth's layout below or a field's type changes.
constexpr std::uint32_t kTruthVersion = 2U;

/// The truth's header; then the rows' cluster indices, particles
/// (`TruthRow::kMaxParticles` per row), counts, layers and first-label flags,
/// then the particles' identities, pT, primary flags, eta, phi, origins, PDG
/// codes, mothers' PDG codes, production processes and vertices, by identity.
struct TruthHeader {
  std::array<char, 8> magic = kTruthMagic;
  std::uint32_t version = kTruthVersion;
  std::uint32_t maxParticles = TruthRow::kMaxParticles;
  std::uint64_t rows = 0U;
  std::uint64_t particles = 0U;
};
static_assert(std::has_unique_object_representations_v<TruthHeader>);

/// A stream's view of an object's bytes.
[[nodiscard]] const char* bytesOf(const void* data) { return static_cast<const char*>(data); }
[[nodiscard]] char* bytesOf(void* data) { return static_cast<char*>(data); }

template <typename T>
std::uint64_t digestOf(const T* data, std::size_t count)
{
  Digest digest;
  if constexpr (std::is_same_v<T, VertexSeed>) {
    for (std::size_t index = 0; index < count; ++index) {
      digest.add(data[index]);
    }
  } else {
    digest.column(data, count);
  }
  return digest.result();
}

template <typename T>
void writeBlock(std::ofstream& out, const T* data, std::size_t count)
{
  const std::uint64_t elements = count;
  const std::uint64_t digest = digestOf(data, count);
  out.write(bytesOf(&elements), sizeof(elements));
  out.write(bytesOf(&digest), sizeof(digest));
  if (count != 0U) {
    out.write(bytesOf(data), static_cast<std::streamsize>(count * sizeof(T)));
  }
}

template <typename T>
void readBlock(std::ifstream& in, T* data, std::size_t count, const char* what)
{
  std::uint64_t elements = 0U;
  std::uint64_t digest = 0U;
  in.read(bytesOf(&elements), sizeof(elements));
  in.read(bytesOf(&digest), sizeof(digest));
  if (!in || elements != count) {
    throw std::runtime_error(std::string{"tracking input: bad block header for "} + what);
  }
  if (count != 0U) {
    in.read(bytesOf(data), static_cast<std::streamsize>(count * sizeof(T)));
  }
  if (!in || digestOf(data, count) != digest) {
    throw std::runtime_error(std::string{"tracking input: short or corrupt "} + what);
  }
}

} // namespace

void writeTrackingInput(const std::string& path, const MeasurementView& view, std::size_t rows, std::span<const std::uint32_t> offsets,
                        std::span<const VertexSeed> vertices, const workflow::BeamLine& beamLine, std::span<const FrameTiming> frames, std::span<const double> depths)
{
  std::ofstream out{path, std::ios::binary | std::ios::trunc};
  if (!out) {
    throw std::runtime_error("tracking input: cannot create " + path);
  }
  const Header header{.rows = rows,
                      .offsets = offsets.size(),
                      .vertices = vertices.size(),
                      .frameLayers = frames.size(),
                      .depthLayers = depths.size(),
                      .beamX = beamLine.x,
                      .beamY = beamLine.y,
                      .beamZ = beamLine.z,
                      .beamSigmaX = beamLine.sigmaX,
                      .beamSigmaY = beamLine.sigmaY,
                      .beamSigmaZ = beamLine.sigmaZ};
  out.write(bytesOf(&header), sizeof(header));
  writeBlock(out, frames.data(), frames.size());
  view.apply([&](const auto*... columns) { (writeBlock(out, columns, rows), ...); });
  writeBlock(out, offsets.data(), offsets.size());
  writeBlock(out, vertices.data(), vertices.size());
  writeBlock(out, depths.data(), depths.size());
  out.close();
  if (!out) {
    throw std::runtime_error("tracking input: cannot write " + path);
  }
}

void writeParameters(const std::string& path)
{
  o2::conf::ConfigurableParam::writeINI(path);

  // Each floating-point member of every detector's RecoSettings, by its
  // parameter's section, named as ConfigurableParam names it, in the shortest
  // form that reads back to the same bits.
  std::map<std::string, std::map<std::string, std::string>> exact;
  const auto add = [&exact](const std::string& section, const RecoSettings& param) {
    auto& values = exact[section];
    auto* members = TClass::GetClass(typeid(RecoSettings))->GetListOfDataMembers();
    for (int i = 0; i < members->GetEntries(); ++i) {
      const auto* member = dynamic_cast<const TDataMember*>(members->At(i));
      const auto* type = member == nullptr ? nullptr : member->GetDataType();
      if (type == nullptr || (member->Property() & kIsStatic) != 0 || (type->GetType() != kDouble_t && type->GetType() != kFloat_t)) {
        continue;
      }
      if (member->GetArrayDim() > 1) {
        throw std::runtime_error(std::format("parameters: {}.{} has more than one dimension", section, member->GetName()));
      }
      const int size = member->GetArrayDim() == 1 ? member->GetMaxIndex(0) : 1;
      const auto bytes = std::span{bytesOf(&param), sizeof(param)}.subspan(static_cast<std::size_t>(member->GetOffset()));
      const auto value = [&]<typename T>(T, int index) {
        T number;
        std::memcpy(&number, bytes.subspan(static_cast<std::size_t>(index) * sizeof(T), sizeof(T)).data(), sizeof(T));
        return std::format("{}", number);
      };
      for (int index = 0; index < size; ++index) {
        const auto name = size > 1 ? std::format("{}[{}]", member->GetName(), index) : std::string{member->GetName()};
        values[name] = type->GetType() == kDouble_t ? value(0., index) : value(0.F, index);
      }
    }
  };
  const auto& its = RecoParam<o2::detectors::DetID::ITS>::Instance();
  const auto& mft = RecoParam<o2::detectors::DetID::MFT>::Instance();
  const auto& trk = RecoParam<o2::detectors::DetID::TRK>::Instance();
  add(its.getName(), its);
  add(mft.getName(), mft);
  add(trk.getName(), trk);

  std::vector<std::string> lines;
  {
    std::ifstream in{path};
    if (!in) {
      throw std::runtime_error("parameters: cannot read " + path);
    }
    const std::map<std::string, std::string>* section = nullptr;
    for (std::string line; std::getline(in, line);) {
      if (!line.empty() && line.front() == '[') {
        const auto found = exact.find(line.substr(1, line.size() - 2));
        section = found == exact.end() ? nullptr : &found->second;
      } else if (const auto equals = line.find('='); section != nullptr && equals != std::string::npos) {
        if (const auto value = section->find(line.substr(0, equals)); value != section->end()) {
          line = line.substr(0, equals + 1) + value->second;
        }
      }
      lines.push_back(std::move(line));
    }
  }
  std::ofstream out{path, std::ios::trunc};
  for (const auto& line : lines) {
    out << line << '\n';
  }
  if (!out) {
    throw std::runtime_error("parameters: cannot write " + path);
  }
}

TrackingInput readTrackingInput(const std::string& path)
{
  std::ifstream in{path, std::ios::binary};
  if (!in) {
    throw std::runtime_error("tracking input: cannot open " + path);
  }
  Header header;
  in.read(bytesOf(&header), sizeof(header));
  if (!in || header.magic != kMagic) {
    throw std::runtime_error("tracking input: " + path + " is not a tracking input");
  }
  if (header.version != kVersion || header.columns != kNMeasurementColumns || header.rowBytes != kMeasurementRowBytes) {
    throw std::runtime_error("tracking input: " + path + " was written with another layout; dump it again");
  }
  TrackingInput input;
  input.beamLine = {.x = header.beamX, .y = header.beamY, .z = header.beamZ, .sigmaX = header.beamSigmaX, .sigmaY = header.beamSigmaY, .sigmaZ = header.beamSigmaZ};
  input.frames.resize(header.frameLayers);
  readBlock(in, input.frames.data(), input.frames.size(), "the readout frames");
  resizeMeasurementColumns(input.columns, header.rows);
  input.columns.apply([&](auto&... columns) { (readBlock(in, columns.data(), header.rows, "a measurement column"), ...); });
  input.offsets.resize(header.offsets);
  readBlock(in, input.offsets.data(), input.offsets.size(), "the grid offsets");
  input.vertices.resize(header.vertices);
  readBlock(in, input.vertices.data(), input.vertices.size(), "the vertices");
  input.depths.resize(header.depthLayers);
  readBlock(in, input.depths.data(), input.depths.size(), "the measurement depths");
  return input;
}

void writeTrackingTruth(const std::string& path, std::span<const std::uint32_t> clusterIndex, std::span<const TruthRow> rows, const ParticleTable& particles)
{
  if (clusterIndex.size() != rows.size()) {
    throw std::runtime_error(std::format("tracking truth: {} rows for {} cluster indices", rows.size(), clusterIndex.size()));
  }
  // Field by field: a TruthRow has padding, which neither the file nor its
  // digest may hold.
  std::vector<std::uint64_t> rowParticles;
  std::vector<std::uint8_t> counts;
  std::vector<LayerID> layers;
  std::vector<std::uint8_t> firstLabelCounts;
  rowParticles.reserve(rows.size() * TruthRow::kMaxParticles);
  counts.reserve(rows.size());
  layers.reserve(rows.size());
  firstLabelCounts.reserve(rows.size());
  for (const auto& row : rows) {
    rowParticles.insert(rowParticles.end(), row.particles.begin(), row.particles.end());
    counts.push_back(row.count);
    layers.push_back(row.layer);
    firstLabelCounts.push_back(row.firstLabelCounts ? 1U : 0U);
  }
  // By identity, so that the same truth gives the same file.
  std::vector<std::uint64_t> ids;
  ids.reserve(particles.size());
  for (const auto& [id, particle] : particles) {
    ids.push_back(id);
  }
  std::ranges::sort(ids);
  std::vector<float> pt;
  std::vector<std::uint8_t> primary;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<std::uint8_t> origin;
  std::vector<std::int32_t> pdg;
  std::vector<std::int32_t> motherPdg;
  std::vector<std::uint8_t> process;
  std::vector<float> vx;
  std::vector<float> vy;
  std::vector<float> vz;
  for (const auto id : ids) {
    const auto& particle = particles.at(id);
    pt.push_back(particle.pt);
    primary.push_back(particle.primary ? 1U : 0U);
    eta.push_back(particle.eta);
    phi.push_back(particle.phi);
    origin.push_back(static_cast<std::uint8_t>(particle.origin));
    pdg.push_back(particle.pdg);
    motherPdg.push_back(particle.motherPdg);
    process.push_back(particle.process);
    vx.push_back(particle.vx);
    vy.push_back(particle.vy);
    vz.push_back(particle.vz);
  }

  std::ofstream out{path, std::ios::binary | std::ios::trunc};
  if (!out) {
    throw std::runtime_error("tracking truth: cannot create " + path);
  }
  const TruthHeader header{.rows = rows.size(), .particles = ids.size()};
  out.write(bytesOf(&header), sizeof(header));
  writeBlock(out, clusterIndex.data(), clusterIndex.size());
  writeBlock(out, rowParticles.data(), rowParticles.size());
  writeBlock(out, counts.data(), counts.size());
  writeBlock(out, layers.data(), layers.size());
  writeBlock(out, firstLabelCounts.data(), firstLabelCounts.size());
  writeBlock(out, ids.data(), ids.size());
  writeBlock(out, pt.data(), pt.size());
  writeBlock(out, primary.data(), primary.size());
  writeBlock(out, eta.data(), eta.size());
  writeBlock(out, phi.data(), phi.size());
  writeBlock(out, origin.data(), origin.size());
  writeBlock(out, pdg.data(), pdg.size());
  writeBlock(out, motherPdg.data(), motherPdg.size());
  writeBlock(out, process.data(), process.size());
  writeBlock(out, vx.data(), vx.size());
  writeBlock(out, vy.data(), vy.size());
  writeBlock(out, vz.data(), vz.size());
  out.close();
  if (!out) {
    throw std::runtime_error("tracking truth: cannot write " + path);
  }
}

TrackingTruth readTrackingTruth(const std::string& path, const TrackingInput& input)
{
  std::ifstream in{path, std::ios::binary};
  if (!in) {
    throw std::runtime_error("tracking truth: cannot open " + path);
  }
  TruthHeader header;
  in.read(bytesOf(&header), sizeof(header));
  if (!in || header.magic != kTruthMagic) {
    throw std::runtime_error("tracking truth: " + path + " is not a tracking truth");
  }
  if (header.version != kTruthVersion || header.maxParticles != TruthRow::kMaxParticles) {
    throw std::runtime_error("tracking truth: " + path + " was written with another layout; dump it again");
  }
  const auto rows = static_cast<std::size_t>(header.rows);
  if (rows != input.columns.clusterIndex.size()) {
    throw std::runtime_error(std::format("tracking truth: {} has {} rows, the tracking input {}", path, rows, input.columns.clusterIndex.size()));
  }
  std::vector<std::uint32_t> clusterIndex(rows);
  std::vector<std::uint64_t> rowParticles(rows * TruthRow::kMaxParticles);
  std::vector<std::uint8_t> counts(rows);
  std::vector<LayerID> layers(rows);
  std::vector<std::uint8_t> firstLabelCounts(rows);
  readBlock(in, clusterIndex.data(), clusterIndex.size(), "the truth's cluster indices");
  readBlock(in, rowParticles.data(), rowParticles.size(), "the truth's row particles");
  readBlock(in, counts.data(), counts.size(), "the truth's particle counts");
  readBlock(in, layers.data(), layers.size(), "the truth's layers");
  readBlock(in, firstLabelCounts.data(), firstLabelCounts.size(), "the truth's first-label flags");
  if (!std::ranges::equal(clusterIndex, input.columns.clusterIndex)) {
    throw std::runtime_error("tracking truth: " + path + " is not the tracking input's: their cluster indices differ");
  }

  const auto particles = static_cast<std::size_t>(header.particles);
  std::vector<std::uint64_t> ids(particles);
  std::vector<float> pt(particles);
  std::vector<std::uint8_t> primary(particles);
  std::vector<float> eta(particles);
  std::vector<float> phi(particles);
  std::vector<std::uint8_t> origin(particles);
  std::vector<std::int32_t> pdg(particles);
  std::vector<std::int32_t> motherPdg(particles);
  std::vector<std::uint8_t> process(particles);
  std::vector<float> vx(particles);
  std::vector<float> vy(particles);
  std::vector<float> vz(particles);
  readBlock(in, ids.data(), ids.size(), "the particles' identities");
  readBlock(in, pt.data(), pt.size(), "the particles' pT");
  readBlock(in, primary.data(), primary.size(), "the particles' primary flags");
  readBlock(in, eta.data(), eta.size(), "the particles' eta");
  readBlock(in, phi.data(), phi.size(), "the particles' phi");
  readBlock(in, origin.data(), origin.size(), "the particles' origins");
  readBlock(in, pdg.data(), pdg.size(), "the particles' PDG codes");
  readBlock(in, motherPdg.data(), motherPdg.size(), "the particles' mothers' PDG codes");
  readBlock(in, process.data(), process.size(), "the particles' production processes");
  readBlock(in, vx.data(), vx.size(), "the particles' production x");
  readBlock(in, vy.data(), vy.size(), "the particles' production y");
  readBlock(in, vz.data(), vz.size(), "the particles' production z");
  if (std::ranges::any_of(origin, [](std::uint8_t value) { return value >= kParticleOrigins; })) {
    throw std::runtime_error("tracking truth: " + path + " has a particle of no known origin");
  }

  TrackingTruth truth;
  truth.rows.resize(rows);
  for (std::size_t row = 0; row < rows; ++row) {
    auto& entry = truth.rows[row];
    std::copy_n(rowParticles.begin() + static_cast<std::ptrdiff_t>(row * TruthRow::kMaxParticles), TruthRow::kMaxParticles, entry.particles.begin());
    if (counts[row] > TruthRow::kMaxParticles) {
      throw std::runtime_error("tracking truth: " + path + " has a row of more than TruthRow::kMaxParticles particles");
    }
    entry.count = counts[row];
    entry.layer = layers[row];
    entry.r = input.columns.r[row];
    entry.z = input.columns.z[row];
    entry.firstLabelCounts = firstLabelCounts[row] != 0U;
  }
  truth.particles.reserve(particles);
  for (std::size_t index = 0; index < particles; ++index) {
    truth.particles[ids[index]] = ParticleTruth{.pt = pt[index],
                                                .primary = primary[index] != 0U,
                                                .eta = eta[index],
                                                .phi = phi[index],
                                                .origin = static_cast<ParticleOrigin>(origin[index]),
                                                .pdg = pdg[index],
                                                .motherPdg = motherPdg[index],
                                                .process = process[index],
                                                .vx = vx[index],
                                                .vy = vy[index],
                                                .vz = vz[index]};
  }
  return truth;
}

} // namespace o2::gtrk::study
