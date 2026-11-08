#include "supera/process/SelectSpacePoints.h"

#include <algorithm>
#include <unordered_map>

namespace supera {
namespace {
void validate_source(const SpacePointCollection &source, unsigned int planes, bool require_hits) {
  if (!source.charge || source.charge->size() != source.points.size())
    throw std::invalid_argument("Space-point selection requires a complete charge column");
  if (!source.associated_hits) {
    if (require_hits)
      throw std::invalid_argument("Hit-count ranking requires associated hits");
    return;
  }

  // Check the complete graph once, before any representative can be selected.
  const auto &data = *source.associated_hits;
  if (data.offsets.size() != source.points.size() + 1 || data.offsets.front() != 0 ||
      data.offsets.back() != data.links.size())
    throw std::invalid_argument("Malformed space-point association offsets");
  for (const auto &collection : data.collections) {
    if (collection.source_keys.size() != collection.hits.size())
      throw std::invalid_argument("Hit identities and records differ in size");
    for (auto key : collection.source_keys)
      require_index(key, collection.source_size);
  }
  std::vector<bool> seen(planes);
  for (std::size_t point = 0; point < source.points.size(); ++point) {
    const auto begin = data.offsets[point], end = data.offsets[point + 1];
    if (end < begin || end > data.links.size())
      throw std::invalid_argument("Nonmonotonic space-point association offsets");
    if (end - begin > planes)
      throw std::invalid_argument("More associated hits than configured planes");
    std::fill(seen.begin(), seen.end(), false);
    std::optional<WireAddress> location;
    for (auto link = begin; link < end; ++link) {
      const auto &ref = data.links[link];
      require_index(ref.collection, data.collections.size());
      const auto &hits = data.collections[ref.collection].hits;
      require_index(ref.index, hits.size());
      const auto &wire = hits[ref.index].wire;
      if (!wire || wire->plane >= planes || seen[wire->plane])
        throw std::invalid_argument("Invalid or repeated hit plane");
      if (location && (wire->cryostat != location->cryostat || wire->tpc != location->tpc))
        throw std::invalid_argument("Mixed cryostat/TPC in space-point associations");
      seen[wire->plane] = true;
      location = wire;
    }
  }
}
} // namespace

SpacePointVoxels select_spacepoints(std::string name, std::vector<SpacePointCollection> sources,
                                    const FixedGrid &grid, const SpacePointSelection &selection) {
  if (name.empty() || !std::isfinite(selection.minimum_charge) ||
      !std::isfinite(selection.maximum_charge) ||
      selection.minimum_charge > selection.maximum_charge || selection.planes == 0)
    throw std::invalid_argument("Invalid space-point selection configuration");
  for (auto shift : selection.shift_cm)
    if (!std::isfinite(shift))
      throw std::invalid_argument("Nonfinite space-point shift");

  SpacePointVoxels output{
      std::move(name), grid.metadata(), selection.planes, std::move(sources), {}, {}};
  std::unordered_map<RecordIndex, std::size_t> occupied;
  for (std::size_t source_index = 0; source_index < output.sources.size(); ++source_index) {
    const auto &source = output.sources[source_index];
    validate_source(source, selection.planes, selection.prefer_more_hits);
    for (std::size_t point = 0; point < source.points.size(); ++point) {
      auto position = source.points[point].position_cm;
      double charge = source.charge->at(point);
      bool finite = std::isfinite(charge);
      if (selection.float_charge && finite) {
        if (std::abs(charge) > std::numeric_limits<float>::max())
          finite = false;
        else
          charge = static_cast<float>(charge);
      }
      for (std::size_t axis = 0; axis < 3; ++axis) {
        position[axis] += selection.shift_cm[axis];
        finite = finite && std::isfinite(position[axis]);
      }
      if (!finite) {
        ++output.report.nonfinite;
        continue;
      }
      if (charge < selection.minimum_charge || charge > selection.maximum_charge) {
        ++output.report.charge_cut;
        continue;
      }
      const auto id = grid.id(position);
      if (id == missing_index) {
        ++output.report.outside_grid;
        continue;
      }

      // A voxel always refers to one whole point, never a mixture of old and new fields.
      const auto count = source.associated_hits ? source.associated_hits->offsets[point + 1] -
                                                      source.associated_hits->offsets[point]
                                                : 0;
      SpacePointVoxel candidate{id,    record_index(source_index), record_index(point), 1, count,
                                charge};
      const auto inserted = occupied.emplace(id, output.voxels.size());
      if (inserted.second) {
        output.voxels.push_back(candidate);
        continue;
      }
      auto &previous = output.voxels[inserted.first->second];
      ++previous.occupancy;
      ++output.report.collisions;
      if (selection.prefer_more_hits && count != previous.hit_count) {
        if (count < previous.hit_count)
          continue;
      } else if (charge < previous.charge) {
        continue;
      }
      candidate.occupancy = previous.occupancy;
      previous = candidate;
      ++output.report.replacements;
    }
  }

  // Sort once, instead of inserting repeatedly into a sorted sparse tensor.
  std::sort(output.voxels.begin(), output.voxels.end(),
            [](const auto &left, const auto &right) { return left.id < right.id; });
  return output;
}
} // namespace supera
