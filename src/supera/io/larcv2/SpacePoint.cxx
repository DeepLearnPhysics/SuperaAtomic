#include "SpacePoint.h"

#include <map>
#include <set>

namespace supera::io::larcv2 {
namespace {
float measurement(double value) {
  require_measurement(value);
  if (std::isnan(value))
    return -1.F; // Legacy space-point field absence; never used in shared records.
  if (std::abs(value) > std::numeric_limits<float>::max())
    throw std::out_of_range("Space-point measurement exceeds LArCV2 float range");
  return static_cast<float>(value);
}

float identity(std::uint64_t value) {
  const float output = static_cast<float>(value);
  if (static_cast<long double>(output) != static_cast<long double>(value))
    throw std::out_of_range("Space-point identity cannot be represented exactly as LArCV2 float");
  return output;
}

bool per_plane(const std::string &field) {
  return field.rfind("hit_", 0) == 0;
}
} // namespace

std::vector<std::string> spacepoint_labels(const SpacePointOutput &configuration) {
  const std::set<std::string> supported{"",         "chi2",    "occupancy",  "nhits",
                                        "cryo",     "tpc",     "hit_charge", "hit_amp",
                                        "hit_time", "hit_rms", "hit_mult",   "hit_key"};
  if (configuration.name.empty() || configuration.planes == 0 || configuration.fields.empty())
    throw std::invalid_argument("Empty space-point output configuration");
  std::set<std::string> seen;
  std::vector<std::string> labels;
  for (const auto &field : configuration.fields) {
    if (!supported.count(field) || !seen.insert(field).second)
      throw std::invalid_argument("Unknown or duplicate space-point output field: " + field);
    const auto prefix = configuration.name + (field.empty() ? "" : "_" + field);
    if (per_plane(field)) {
      for (unsigned int plane = 0; plane < configuration.planes; ++plane)
        labels.push_back(prefix + std::to_string(plane));
    } else {
      labels.push_back(prefix);
    }
  }
  return labels;
}

std::vector<larcv::SparseTensor3D> project_spacepoints(const SpacePointVoxels &collection,
                                                       const SpacePointOutput &configuration) {
  const auto labels = spacepoint_labels(configuration);
  if (collection.name != configuration.name || collection.planes != configuration.planes)
    throw std::invalid_argument("Space-point output name/plane mismatch");

  // Validate metadata before handing it to the legacy unchecked constructor.
  const auto &grid = collection.grid;
  std::uint64_t size = 1;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    if (!std::isfinite(grid.minimum_cm[axis]) || !std::isfinite(grid.maximum_cm[axis]) ||
        grid.maximum_cm[axis] <= grid.minimum_cm[axis] ||
        !std::isfinite(grid.maximum_cm[axis] - grid.minimum_cm[axis]) || grid.counts[axis] == 0 ||
        grid.counts[axis] > ((std::uint64_t{1} << 53) - 1) ||
        grid.counts[axis] >
            static_cast<std::uint64_t>(std::numeric_limits<RecordIndex>::max()) / size)
      throw std::invalid_argument("Invalid space-point output grid");
    if ((grid.maximum_cm[axis] - grid.minimum_cm[axis]) / grid.counts[axis] <= 0.)
      throw std::invalid_argument("Space-point output grid width underflow");
    size *= grid.counts[axis];
  }
  larcv::Voxel3DMeta meta;
  meta.set(grid.minimum_cm[0], grid.minimum_cm[1], grid.minimum_cm[2], grid.maximum_cm[0],
           grid.maximum_cm[1], grid.maximum_cm[2], grid.counts[0], grid.counts[1], grid.counts[2]);
  std::vector<larcv::SparseTensor3D> output(labels.size());
  for (auto &tensor : output) {
    tensor.meta(meta);
    tensor.reserve(collection.voxels.size());
  }

  // Legacy concatenation assigns earlier products offsets from later product sizes.
  // Retain this convention only here. Shared references remain exact product/key pairs.
  std::map<std::string, std::uint64_t> offsets;
  std::map<std::string, std::uint64_t> source_sizes;
  std::uint64_t total = 0;
  for (auto source = collection.sources.rbegin(); source != collection.sources.rend(); ++source) {
    if (!source->associated_hits)
      continue;
    for (auto hits = source->associated_hits->collections.rbegin();
         hits != source->associated_hits->collections.rend(); ++hits) {
      if (hits->source_product.empty())
        throw std::invalid_argument("Missing hit product identity");
      if (offsets.count(hits->source_product)) {
        if (source_sizes.at(hits->source_product) != hits->source_size)
          throw std::invalid_argument("Inconsistent size for repeated hit product");
        continue;
      }
      source_sizes.emplace(hits->source_product, hits->source_size);
      offsets.emplace(hits->source_product, total);
      if (hits->source_size >
          static_cast<std::uint64_t>(std::numeric_limits<RecordIndex>::max()) - total)
        throw std::overflow_error("Hit product offsets exceed supported range");
      total += hits->source_size;
    }
  }

  RecordIndex previous = missing_index;
  std::vector<const WireHit *> plane_hits(collection.planes);
  std::vector<RecordIndex> plane_keys(collection.planes);
  for (const auto &voxel : collection.voxels) {
    require_index(voxel.id, size);
    if (voxel.id <= previous)
      throw std::invalid_argument("Space-point output voxels must be sorted and unique");
    previous = voxel.id;
    require_index(voxel.source, collection.sources.size());
    const auto &source = collection.sources[voxel.source];
    require_index(voxel.point, source.points.size());
    const auto &point = source.points[voxel.point];

    // Derive every attribute from the same selected point, resetting absent planes.
    std::fill(plane_hits.begin(), plane_hits.end(), nullptr);
    std::fill(plane_keys.begin(), plane_keys.end(), missing_index);
    std::optional<WireAddress> location;
    if (source.associated_hits) {
      const auto &data = *source.associated_hits;
      const auto begin = data.offsets.at(voxel.point), end = data.offsets.at(voxel.point + 1);
      if (end < begin || end > data.links.size() || end - begin != voxel.hit_count)
        throw std::invalid_argument("Malformed selected point associations");
      for (auto link = begin; link < end; ++link) {
        const auto &ref = data.links[link];
        require_index(ref.collection, data.collections.size());
        const auto &hits = data.collections[ref.collection];
        require_index(ref.index, hits.hits.size());
        const auto &hit = hits.hits[ref.index];
        if (!hit.wire || hit.wire->plane >= collection.planes || plane_hits[hit.wire->plane])
          throw std::invalid_argument("Invalid or repeated selected hit plane");
        if (location &&
            (location->cryostat != hit.wire->cryostat || location->tpc != hit.wire->tpc))
          throw std::invalid_argument("Mixed selected hit TPCs");
        location = hit.wire;
        const auto key = hits.source_keys.at(ref.index);
        require_index(key, hits.source_size);
        plane_keys[hit.wire->plane] = record_index(offsets.at(hits.source_product) + key);
        plane_hits[hit.wire->plane] = &hit;
      }
    }

    std::size_t product = 0;
    for (const auto &field : configuration.fields) {
      if (per_plane(field)) {
        if (!source.associated_hits)
          throw std::invalid_argument("Requested wire output without associated hits");
        for (unsigned int plane = 0; plane < collection.planes; ++plane) {
          const auto *hit = plane_hits[plane];
          float value = -1.F;
          if (hit) {
            if (field == "hit_charge")
              value = measurement(hit->integral_adc_tick);
            else if (field == "hit_amp")
              value = measurement(hit->amplitude_adc);
            else if (field == "hit_time")
              value = measurement(hit->peak_tick);
            else if (field == "hit_rms")
              value = measurement(hit->rms_tick);
            else if (field == "hit_mult")
              value = measurement(hit->multiplicity);
            else
              value = identity(plane_keys[plane]);
          }
          output[product++].emplace(voxel.id, value, false);
        }
      } else {
        float value = -1.F;
        if (field.empty())
          value = measurement(voxel.charge);
        else if (field == "chi2")
          value = measurement(point.chi_square);
        else if (field == "occupancy")
          value = identity(voxel.occupancy);
        else if (field == "nhits") {
          if (!source.associated_hits)
            throw std::invalid_argument("Requested hit count without associated hits");
          value = identity(voxel.hit_count);
        } else if (location)
          value = identity(field == "cryo" ? location->cryostat : location->tpc);
        else if (!source.associated_hits)
          throw std::invalid_argument("Requested detector identity without associated hits");
        output[product++].emplace(voxel.id, value, false);
      }
    }
  }
  return output;
}
} // namespace supera::io::larcv2
