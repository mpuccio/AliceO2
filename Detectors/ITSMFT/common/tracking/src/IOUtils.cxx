// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.

#include "ITSMFTTracking/IOUtils.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <vector>

#include "ITSMFTTracking/SurfaceMask.h"
#include "ITSMFTTracking/TimeFrame.h"
#include "Framework/Logger.h"
#include "ITSBase/GeometryTGeo.h"
#include "MFTBase/GeometryTGeo.h"
#include "MathUtils/Utils.h"

namespace
{

using o2::itsmft::ioutils::detail::addSysErrors;
using o2::itsmft::ioutils::detail::shouldApplySysErrors;

struct DecodedClusterResult {
  o2::itsmft::tracking::DecodedCluster decoded{};
  o2::itsmft::tracking::ClusterDecodeError error{o2::itsmft::tracking::ClusterDecodeError::None};

  bool ok() const noexcept { return error == o2::itsmft::tracking::ClusterDecodeError::None; }
};

template <o2::detectors::DetID::ID DetId, typename GeomT>
DecodedClusterResult decodeClusterBounded(GeomT* geom,
                                          const o2::itsmft::CompClusterExt& cluster,
                                          o2::itsmft::tracking::BoundedPatternCursor& patterns,
                                          const o2::itsmft::TopologyDictionary* dict,
                                          bool applySysErrors)
{
  using o2::itsmft::tracking::ClusterDecodeError;
  DecodedClusterResult result;
  if (dict == nullptr) {
    result.error = ClusterDecodeError::MissingDictionary;
    return result;
  }
  if (geom == nullptr) {
    result.error = ClusterDecodeError::GeometryUnavailable;
    return result;
  }

  const auto sensorID = cluster.getSensorID();
  if (!o2::itsmft::ioutils::detail::isSensorInGeometry(sensorID, geom->getSize())) {
    result.error = ClusterDecodeError::InvalidSensor;
    return result;
  }
  const int layer = geom->getLayer(sensorID);
  if (!o2::itsmft::ioutils::detail::isLayerInDetector(layer, o2::itsmft::tracking::TrackerParamRef<DetId>::nLayers())) {
    result.error = ClusterDecodeError::InvalidLayer;
    return result;
  }

  const auto clusterData = o2::itsmft::ioutils::extractClusterDataBounded(cluster, patterns, dict);
  if (!clusterData.ok()) {
    result.error = clusterData.error;
    return result;
  }
  float sigma2Row = clusterData.sig2Row;
  float sigma2Col = clusterData.sig2Col;
  if (applySysErrors && shouldApplySysErrors<DetId>()) {
    addSysErrors<DetId>(layer, sigma2Row, sigma2Col);
  }

  if constexpr (DetId == o2::detectors::DetID::ITS) {
    const auto trkXYZ = geom->getMatrixT2L(sensorID) ^ clusterData.coordinates;
    const auto gloXYZ = geom->getMatrixL2G(sensorID) * clusterData.coordinates;
    result.decoded = {{gloXYZ.x(), gloXYZ.y(), gloXYZ.z()},
                      {trkXYZ.x(), trkXYZ.y(), trkXYZ.z(), geom->getSensorRefAlpha(sensorID)},
                      {sigma2Row, 0.f, sigma2Col},
                      clusterData.shape,
                      static_cast<uint32_t>(sensorID),
                      layer};
  } else {
    if (!geom->getCacheL2G().isFilled() || geom->getCacheL2G().getSize() <= sensorID) {
      result.error = ClusterDecodeError::GeometryUnavailable;
      return result;
    }
    const auto gloXYZ = geom->getMatrixL2G(sensorID) * clusterData.coordinates;
    result.decoded = {{gloXYZ.x(), gloXYZ.y(), gloXYZ.z()}, {}, {sigma2Row, 0.f, sigma2Col}, clusterData.shape, static_cast<uint32_t>(sensorID), layer};
  }
  return result;
}

} // namespace

namespace o2::itsmft::ioutils
{

void fillMatrixCache(o2::detectors::DetID::ID detId)
{
  const auto mask = o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G);
  if (detId == o2::detectors::DetID::ITS) {
    o2::its::GeometryTGeo::Instance()->fillMatrixCache(mask);
  } else if (detId == o2::detectors::DetID::MFT) {
    o2::mft::GeometryTGeo::Instance()->fillMatrixCache(mask);
  } else {
    LOGP(fatal, "Unsupported detector id {} in fillMatrixCache", static_cast<int>(detId));
  }
}

template <o2::detectors::DetID::ID DetId>
o2::itsmft::tracking::SurfaceMeasurementDecodeResult loadClusterSurfaceMeasurement(
  const CompClusterExt& cluster, o2::itsmft::tracking::BoundedPatternCursor& patterns,
  const TopologyDictionary* dict, gsl::span<const o2::itsmft::tracking::LayerId> layerToSurface,
  o2::itsmft::tracking::ClusterSourceId source, uint32_t externalClusterIndex,
  uint32_t sourceROF, bool applySysErrors)
{
  const auto decodedResult = [&] {
    if constexpr (DetId == o2::detectors::DetID::ITS) {
      return decodeClusterBounded<DetId>(o2::its::GeometryTGeo::Instance(), cluster, patterns, dict, applySysErrors);
    } else {
      return decodeClusterBounded<DetId>(o2::mft::GeometryTGeo::Instance(), cluster, patterns, dict, applySysErrors);
    }
  }();

  o2::itsmft::tracking::SurfaceMeasurementDecodeResult result;
  if (!decodedResult.ok()) {
    result.error = decodedResult.error;
    return result;
  }
  const auto& decoded = decodedResult.decoded;
  result.layer = decoded.layer;
  if (decoded.layer < 0 || static_cast<size_t>(decoded.layer) >= layerToSurface.size()) {
    result.error = o2::itsmft::tracking::ClusterDecodeError::InvalidLayerMapping;
    return result;
  }
  result.layerMapped = true;
  const auto surface = layerToSurface[decoded.layer];
  const o2::itsmft::tracking::DetectorSensorId sensor{DetId, decoded.sensor};
  const o2::itsmft::tracking::ClusterRef clusterRef{source, externalClusterIndex};
  if constexpr (DetId == o2::detectors::DetID::ITS) {
    result.kind = o2::itsmft::tracking::SurfaceKind::Cylinder;
    result.global = o2::itsmft::tracking::makeCylinderGlobalMeasurement(decoded, sensor, surface, clusterRef, sourceROF);
    result.measurement = o2::itsmft::tracking::makeCylinderSurfaceMeasurement(decoded);
  } else {
    result.kind = o2::itsmft::tracking::SurfaceKind::Disk;
    result.global = o2::itsmft::tracking::makeDiskGlobalMeasurement(decoded, sensor, surface, clusterRef, sourceROF);
    result.measurement = o2::itsmft::tracking::makeDiskSurfaceMeasurement(decoded);
  }
  return result;
}

template o2::itsmft::tracking::SurfaceMeasurementDecodeResult loadClusterSurfaceMeasurement<o2::detectors::DetID::ITS>(
  const CompClusterExt&, o2::itsmft::tracking::BoundedPatternCursor&, const TopologyDictionary*,
  gsl::span<const o2::itsmft::tracking::LayerId>, o2::itsmft::tracking::ClusterSourceId, uint32_t, uint32_t, bool);
template o2::itsmft::tracking::SurfaceMeasurementDecodeResult loadClusterSurfaceMeasurement<o2::detectors::DetID::MFT>(
  const CompClusterExt&, o2::itsmft::tracking::BoundedPatternCursor&, const TopologyDictionary*,
  gsl::span<const o2::itsmft::tracking::LayerId>, o2::itsmft::tracking::ClusterSourceId, uint32_t, uint32_t, bool);

} // namespace o2::itsmft::ioutils

namespace o2::itsmft::tracking
{

LoadSourcesResult loadTimeFrameSources(TimeFrame& frame, gsl::span<const ClusterSourceInput> sources,
                                       SurfaceCatalogView catalog, const o2::InteractionRecord& origin)
{
  if (!frame.isConfigured()) {
    return {MultiSourceLoadError::FrameNotConfigured};
  }
  // Stage all source measurements before committing the event.
  TimeFrame stagedMeasurements;
  const auto loadResult = loadSources(stagedMeasurements, catalog, sources, origin);
  if (!loadResult.ok()) {
    return loadResult;
  }

  if (frame.getNIterations() == 0) {
    return {MultiSourceLoadError::FrameNotConfigured};
  }
  const auto layout = frame.getLayout(0);
  const auto orderedSurfaces = layout.getOrderedSurfaces();
  if (orderedSurfaces.empty()) {
    return {MultiSourceLoadError::FrameNotConfigured};
  }
  SurfaceMask configuredSurfaces;
  for (const auto surface : orderedSurfaces) {
    configuredSurfaces.set(surface);
  }
  SurfaceMask mappedSurfaces;
  for (const auto& source : sources) {
    for (const auto surface : source.layerToSurface) {
      if (!surface.isValid() || mappedSurfaces.has(surface) || !configuredSurfaces.has(surface)) {
        return {MultiSourceLoadError::InvalidLayerMapping, source.id};
      }
      if (catalog.getSurface(surface).detectorId != static_cast<uint8_t>(source.detector)) {
        return {MultiSourceLoadError::DetectorSurfaceMismatch, source.id};
      }
      mappedSurfaces.set(surface);
    }
  }
  if (mappedSurfaces != configuredSurfaces) {
    // Attribute an omitted surface only when one source owns its detector.
    for (const auto surface : orderedSurfaces) {
      if (mappedSurfaces.has(surface)) {
        continue;
      }
      ClusterSourceId owner;
      for (const auto& source : sources) {
        if (static_cast<uint8_t>(source.detector) != catalog.getSurface(surface).detectorId) {
          continue;
        }
        if (owner.isValid()) {
          return {MultiSourceLoadError::InvalidLayerMapping};
        }
        owner = source.id;
      }
      return {MultiSourceLoadError::InvalidLayerMapping, owner};
    }
    return {MultiSourceLoadError::InvalidLayerMapping};
  }

  std::vector<std::vector<int>> rofBoundaries(orderedSurfaces.size());
  std::vector<RuntimeROFViews> viewsBySurface(orderedSurfaces.size());
  std::vector<uint16_t> localLayerBySurface(orderedSurfaces.size(), 0);
  std::vector<ClusterSourceId> sourceBySurface(orderedSurfaces.size(), ClusterSourceId::invalid());
  for (std::size_t position = 0; position < orderedSurfaces.size(); ++position) {
    const auto surface = orderedSurfaces[position];
    const ClusterSourceInput* owner = nullptr;
    uint16_t localLayer = 0;
    for (const auto& source : sources) {
      const auto it = std::find(source.layerToSurface.begin(), source.layerToSurface.end(), surface);
      if (it == source.layerToSurface.end()) {
        continue;
      }
      if (owner != nullptr) {
        return {MultiSourceLoadError::InvalidLayerMapping, source.id};
      }
      owner = &source;
      localLayer = static_cast<uint16_t>(std::distance(source.layerToSurface.begin(), it));
    }
    if (owner == nullptr) {
      return {MultiSourceLoadError::InvalidLayerMapping};
    }

    const auto globals = stagedMeasurements.getGlobalMeasurements(surface);
    auto& boundaries = rofBoundaries[position];
    boundaries.assign(owner->rofs.size() + 1, 0);
    std::size_t measurement = 0;
    for (std::size_t rof = 0; rof < owner->rofs.size(); ++rof) {
      while (measurement < globals.size() && globals[measurement].sourceROF == rof) {
        if (globals[measurement].surface != surface || globals[measurement].cluster.source != owner->id) {
          return {MultiSourceLoadError::InconsistentDecoderMetadata, owner->id,
                  static_cast<uint32_t>(rof), globals[measurement].cluster.index};
        }
        ++measurement;
      }
      boundaries[rof + 1] = static_cast<int>(measurement);
    }
    if (measurement != globals.size()) {
      return {MultiSourceLoadError::InconsistentDecoderMetadata, owner->id};
    }
    viewsBySurface[position] = owner->rofViews;
    localLayerBySurface[position] = localLayer;
    sourceBySurface[position] = owner->id;
  }

  stagedMeasurements.assignLoadedEventNavigation(std::move(rofBoundaries), sources.front().rofViews,
                                                 std::move(viewsBySurface), std::move(localLayerBySurface),
                                                 std::move(sourceBySurface));
  if (!frame.commitLoadedEvent(stagedMeasurements)) {
    return {MultiSourceLoadError::OtherMalformedInput};
  }
  return {};
}

LoadSourcesResult loadTimeFrameSource(
  TimeFrame& frame,
  const ClusterDecoder& decoder,
  const o2::InteractionRecord& origin,
  const ROFTimingConfig& timing,
  gsl::span<const itsmft::CompClusterExt> clusters,
  gsl::span<const unsigned char> patterns,
  gsl::span<const o2::itsmft::ROFRecord> rofs,
  const itsmft::TopologyDictionary* dictionary,
  const dataformats::MCTruthContainer<MCCompLabel>* labels,
  o2::detectors::DetID::ID detector,
  gsl::span<const LayerId> layerToSurface,
  SurfaceCatalogView catalog,
  bool applySysErrors)
{
  constexpr ClusterSourceId sourceId{0};
  if (detector != o2::detectors::DetID::ITS && detector != o2::detectors::DetID::MFT) {
    return {MultiSourceLoadError::UnsupportedDetector, sourceId};
  }
  if (catalog.surfaces == nullptr || catalog.nSurfaces == 0) {
    return {MultiSourceLoadError::SurfaceCatalogNotConfigured, sourceId};
  }
  ClusterSourceInput source;
  source.id = sourceId;
  source.detector = detector;
  source.clusters = clusters;
  source.patterns = patterns;
  source.rofs = rofs;
  source.dictionary = dictionary;
  source.labels = labels;
  source.layerToSurface = layerToSurface;
  source.timing = timing;
  source.decoder = &decoder;
  source.applySysErrors = applySysErrors;
  source.rofViews = frame.getROFViews();
  return loadTimeFrameSources(frame, gsl::span<const ClusterSourceInput>{&source, 1}, catalog, origin);
}

namespace
{
bool isSupportedDetector(o2::detectors::DetID::ID det) noexcept
{
  return det == o2::detectors::DetID::ITS || det == o2::detectors::DetID::MFT;
}

bool covariance2DIsPositiveSemidefinite(float varianceFirst, float covariance,
                                        float varianceSecond) noexcept
{
  if (!std::isfinite(varianceFirst) || !std::isfinite(covariance) ||
      !std::isfinite(varianceSecond) || varianceFirst < 0.f || varianceSecond < 0.f) {
    return false;
  }
  const double diagonalProduct = static_cast<double>(varianceFirst) * varianceSecond;
  const double covarianceSquared = static_cast<double>(covariance) * covariance;
  const double tolerance = 16. * std::numeric_limits<float>::epsilon() *
                           std::max(diagonalProduct, covarianceSquared);
  return diagonalProduct - covarianceSquared >= -tolerance;
}

bool globalCovarianceIsPositiveSemidefinite(const GlobalCovariance3F& covariance) noexcept
{
  if (!covariance2DIsPositiveSemidefinite(covariance.xx, covariance.xy, covariance.yy) ||
      !covariance2DIsPositiveSemidefinite(covariance.xx, covariance.xz, covariance.zz) ||
      !covariance2DIsPositiveSemidefinite(covariance.yy, covariance.yz, covariance.zz)) {
    return false;
  }
  const double determinant =
    static_cast<double>(covariance.xx) * covariance.yy * covariance.zz +
    2. * static_cast<double>(covariance.xy) * covariance.xz * covariance.yz -
    static_cast<double>(covariance.xx) * covariance.yz * covariance.yz -
    static_cast<double>(covariance.yy) * covariance.xz * covariance.xz -
    static_cast<double>(covariance.zz) * covariance.xy * covariance.xy;
  const double scale = std::max({std::abs(static_cast<double>(covariance.xx) * covariance.yy * covariance.zz),
                                 std::abs(2. * static_cast<double>(covariance.xy) * covariance.xz * covariance.yz),
                                 std::abs(static_cast<double>(covariance.xx) * covariance.yz * covariance.yz),
                                 std::abs(static_cast<double>(covariance.yy) * covariance.xz * covariance.xz),
                                 std::abs(static_cast<double>(covariance.zz) * covariance.xy * covariance.xy)});
  return std::isfinite(determinant) &&
         determinant >= -32. * std::numeric_limits<float>::epsilon() * scale;
}

bool decodedMeasurementIsValid(const GlobalMeasurement& global,
                               const SurfaceMeasurement& local) noexcept
{
  return std::isfinite(global.position.x) && std::isfinite(global.position.y) &&
         std::isfinite(global.position.z) &&
         globalCovarianceIsPositiveSemidefinite(global.covariance) &&
         std::isfinite(local.frame.q) && std::isfinite(local.frame.u) &&
         std::isfinite(local.frame.v) && std::isfinite(local.frame.frameAngle) &&
         covariance2DIsPositiveSemidefinite(local.covariance.uu,
                                            local.covariance.uv,
                                            local.covariance.vv);
}

MultiSourceLoadError mapDecodeError(ClusterDecodeError error) noexcept
{
  switch (error) {
    case ClusterDecodeError::None:
      return MultiSourceLoadError::None;
    case ClusterDecodeError::MissingDictionary:
      return MultiSourceLoadError::MissingDictionary;
    case ClusterDecodeError::TruncatedExplicitPattern:
      return MultiSourceLoadError::TruncatedExplicitPattern;
    case ClusterDecodeError::MalformedExplicitPattern:
      return MultiSourceLoadError::MalformedExplicitPattern;
    case ClusterDecodeError::InvalidPatternId:
      return MultiSourceLoadError::InvalidPatternId;
    case ClusterDecodeError::InvalidSensor:
      return MultiSourceLoadError::InvalidSensor;
    case ClusterDecodeError::InvalidLayer:
      return MultiSourceLoadError::InvalidDecodedLayer;
    case ClusterDecodeError::InvalidLayerMapping:
      return MultiSourceLoadError::InvalidLayerMapping;
    case ClusterDecodeError::GeometryUnavailable:
      return MultiSourceLoadError::GeometryUnavailable;
    case ClusterDecodeError::OtherMalformedInput:
      return MultiSourceLoadError::OtherMalformedInput;
  }
  return MultiSourceLoadError::OtherMalformedInput;
}
} // namespace

LoadSourcesResult loadSources(TimeFrame& frame,
                              const SurfaceCatalogView& catalog,
                              gsl::span<const ClusterSourceInput> sources,
                              const o2::InteractionRecord& origin)
{
  const auto nSources = static_cast<uint32_t>(sources.size());

  std::vector<bool> seen(nSources, false);
  for (const auto& src : sources) {
    if (!src.id.isValid() || src.id.value() >= nSources) {
      return {MultiSourceLoadError::NonDenseSourceIds, src.id};
    }
    if (seen[src.id.value()]) {
      return {MultiSourceLoadError::DuplicateSourceId, src.id};
    }
    seen[src.id.value()] = true;
    if (!isSupportedDetector(src.detector)) {
      return {MultiSourceLoadError::UnsupportedDetector, src.id};
    }
    if (src.decoder == nullptr) {
      return {MultiSourceLoadError::MissingDecoder, src.id};
    }
    if (!src.clusters.empty() && src.dictionary == nullptr) {
      return {MultiSourceLoadError::MissingDictionary, src.id, 0, 0};
    }
  }

  std::vector<std::vector<GlobalMeasurement>> perSurfaceGlobal(catalog.nSurfaces);
  std::vector<std::vector<SurfaceMeasurement>> perSurface(catalog.nSurfaces);
  std::vector<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*> labelSources(nSources, nullptr);

  for (const auto& src : sources) {
    const auto srcIdx = src.id.value();
    labelSources[srcIdx] = src.labels;

    int64_t expectedNext = 0;
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto& rof = src.rofs[r];
      const int64_t first = rof.getFirstEntry();
      const int64_t n = rof.getNEntries();
      if (n < 0 || first != expectedNext) {
        return {MultiSourceLoadError::InvalidROFRange, src.id, r};
      }
      expectedNext = first + n;
      if (expectedNext > static_cast<int64_t>(src.clusters.size())) {
        return {MultiSourceLoadError::InvalidROFRange, src.id, r};
      }
    }
    if (expectedNext != static_cast<int64_t>(src.clusters.size())) {
      return {MultiSourceLoadError::InvalidROFRange, src.id, static_cast<uint32_t>(src.rofs.size())};
    }

    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto built = computeROFIntervalBC(src.rofs[r].getBCData(), origin, src.timing, r);
      if (!built.ok()) {
        return LoadSourcesResult{.error = MultiSourceLoadError::TimingError, .source = src.id, .rof = r, .timingDetail = built.error};
      }
    }

    src.decoder->prepare();
    BoundedPatternCursor patterns{src.patterns};
    for (uint32_t r = 0; r < src.rofs.size(); ++r) {
      const auto& rof = src.rofs[r];
      const auto firstEntry = rof.getFirstEntry();
      const auto nEntries = rof.getNEntries();
      for (int32_t clusterId = firstEntry; clusterId < firstEntry + nEntries; ++clusterId) {
        const auto& cluster = src.clusters[clusterId];
        const auto externalIndex = static_cast<uint32_t>(clusterId);
        const auto decoded = src.decoder->decode(cluster, patterns, src.dictionary, src.layerToSurface,
                                                 src.id, externalIndex, r, src.applySysErrors);
        if (!decoded.ok()) {
          return {mapDecodeError(decoded.error), src.id, r, externalIndex};
        }
        if (!decoded.layerMapped || decoded.layer < 0 ||
            static_cast<size_t>(decoded.layer) >= src.layerToSurface.size()) {
          return {MultiSourceLoadError::InvalidLayerMapping, src.id, r, externalIndex};
        }
        const auto expectedSurface = src.layerToSurface[decoded.layer];
        if (!expectedSurface.isValid() || expectedSurface.value() >= catalog.nSurfaces) {
          return {MultiSourceLoadError::InvalidLayerMapping, src.id, r, externalIndex};
        }
        auto global = decoded.global;
        if (global.surface != expectedSurface ||
            global.cluster.source != src.id ||
            global.cluster.index != externalIndex ||
            global.sourceROF != r ||
            global.sensor.detector != static_cast<uint32_t>(src.detector)) {
          return {MultiSourceLoadError::InconsistentDecoderMetadata, src.id, r, externalIndex};
        }
        const auto& surfaceDescriptor = catalog.getSurface(expectedSurface);
        if (surfaceDescriptor.detectorId != static_cast<uint8_t>(src.detector)) {
          return {MultiSourceLoadError::DetectorSurfaceMismatch, src.id, r, externalIndex};
        }
        if (surfaceDescriptor.kind != decoded.kind) {
          return {MultiSourceLoadError::SurfaceKindMismatch, src.id, r, externalIndex};
        }
        if (!decodedMeasurementIsValid(global, decoded.measurement)) {
          return {MultiSourceLoadError::OtherMalformedInput, src.id, r, externalIndex};
        }
        global.radius = std::hypot(global.position.x, global.position.y);
        perSurfaceGlobal[expectedSurface.value()].push_back(global);
        perSurface[expectedSurface.value()].push_back(decoded.measurement);
      }
    }
    if (!patterns.empty()) {
      return {MultiSourceLoadError::TrailingPatternData, src.id,
              static_cast<uint32_t>(src.rofs.size()), static_cast<uint32_t>(src.clusters.size())};
    }
  }

  frame.assignLoadedMeasurements(std::move(perSurfaceGlobal), std::move(perSurface),
                                 std::move(labelSources));
  return {};
}

namespace
{
std::string formatLoadSourcesResult(const char* label, const LoadSourcesResult& result)
{
  return std::format("{}: error={} source={} rof={} clusterIndex={} timingDetail={}",
                     label, static_cast<int>(result.error), result.source.value(),
                     result.rof, result.clusterIndex, static_cast<int>(result.timingDetail));
}
} // namespace

RecoverableLoadFailure::RecoverableLoadFailure(const LoadSourcesResult& result)
  : std::runtime_error(formatLoadSourcesResult("TimeFrame loading boundary: recoverable data failure", result)),
    mResult(result)
{
}

TimeFrameLoadException::TimeFrameLoadException(TimeFrameLoadFailureReason reason, std::string message)
  : std::runtime_error(std::move(message)), mReason(reason)
{
}

TimeFrameLoadException::TimeFrameLoadException(const LoadSourcesResult& result)
  : std::runtime_error(formatLoadSourcesResult("TimeFrame loading boundary: structural failure", result)),
    mReason(TimeFrameLoadFailureReason::LoadSourcesFailure),
    mLoadResult(result)
{
}

bool isRecoverableLoadError(MultiSourceLoadError error, TimingBuildError timingDetail) noexcept
{
  switch (error) {
    case MultiSourceLoadError::InvalidROFRange:
    case MultiSourceLoadError::TruncatedExplicitPattern:
    case MultiSourceLoadError::MalformedExplicitPattern:
    case MultiSourceLoadError::InvalidPatternId:
    case MultiSourceLoadError::InvalidSensor:
    case MultiSourceLoadError::InvalidDecodedLayer:
    case MultiSourceLoadError::OtherMalformedInput:
    case MultiSourceLoadError::TrailingPatternData:
      return true;
    case MultiSourceLoadError::TimingError:
      return timingDetail == TimingBuildError::Overflow;
    case MultiSourceLoadError::None:
    case MultiSourceLoadError::NonDenseSourceIds:
    case MultiSourceLoadError::DuplicateSourceId:
    case MultiSourceLoadError::UnsupportedDetector:
    case MultiSourceLoadError::MissingDecoder:
    case MultiSourceLoadError::InvalidLayerMapping:
    case MultiSourceLoadError::DetectorSurfaceMismatch:
    case MultiSourceLoadError::InconsistentDecoderMetadata:
    case MultiSourceLoadError::SurfaceKindMismatch:
    case MultiSourceLoadError::SurfaceCatalogNotConfigured:
    case MultiSourceLoadError::SurfaceCatalogStale:
    case MultiSourceLoadError::MissingDictionary:
    case MultiSourceLoadError::GeometryUnavailable:
    case MultiSourceLoadError::FrameNotConfigured:
      return false;
  }
  return false;
}

} // namespace o2::itsmft::tracking
