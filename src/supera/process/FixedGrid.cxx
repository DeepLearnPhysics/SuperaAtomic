#include "supera/process/FixedGrid.h"

#include <limits>
#include <stdexcept>

namespace supera {
namespace {

// Axis counts must be exactly representable by the legacy double index arithmetic.
constexpr std::uint64_t max_axis_count = (std::uint64_t{1} << 53) - 1;

void require_positive(double value) {
  if (!std::isfinite(value) || value <= 0.)
    throw std::invalid_argument("Grid lengths and widths must be finite and positive");
}

RecordIndex canonical_id(VoxelID_t value) {
  return value == kINVALID_VOXELID ? missing_index : record_index(value);
}

} // namespace

FixedGrid::FixedGrid(Grid3D metadata) : metadata_(metadata) {
  // Check bounds and counts before the legacy indexer can divide or multiply them.
  std::uint64_t total = 1;
  const auto maximum = static_cast<std::uint64_t>(std::numeric_limits<RecordIndex>::max());
  for (std::size_t axis = 0; axis < 3; ++axis) {
    const auto low = metadata.minimum_cm[axis];
    const auto high = metadata.maximum_cm[axis];
    if (!std::isfinite(low) || !std::isfinite(high))
      throw std::invalid_argument("Grid bounds must be finite");
    require_positive(high - low);

    const auto count = metadata.counts[axis];
    if (count == 0 || count > max_axis_count || count > maximum / total ||
        count > std::numeric_limits<std::size_t>::max())
      throw std::invalid_argument("Grid counts exceed supported indexing range");
    total *= count;
    require_positive((high - low) / static_cast<double>(count));
  }
  if (total >= std::numeric_limits<VoxelID_t>::max() ||
      total > std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument("Grid exceeds legacy indexer capacity");

  // Keep the established coordinate arithmetic; translate its sentinel only here.
  const auto &low = metadata.minimum_cm;
  const auto &high = metadata.maximum_cm;
  const auto &count = metadata.counts;
  indexer_.set(low[0], low[1], low[2], high[0], high[1], high[2], count[0], count[1], count[2]);
  size_ = record_index(total);
}

FixedGrid FixedGrid::from_lengths(const std::array<double, 3> &minimum_cm,
                                  const std::array<double, 3> &length_cm,
                                  const std::array<double, 3> &requested_width_cm) {
  Grid3D metadata;
  metadata.minimum_cm = minimum_cm;

  // Validate the quotient before converting to an integer; truncation is intentional.
  for (std::size_t axis = 0; axis < 3; ++axis) {
    require_positive(length_cm[axis]);
    require_positive(requested_width_cm[axis]);
    const auto quotient = length_cm[axis] / requested_width_cm[axis];
    if (!std::isfinite(quotient) || quotient < 1. || quotient > static_cast<double>(max_axis_count))
      throw std::invalid_argument("Requested voxel size produces unsupported grid counts");
    metadata.counts[axis] = static_cast<std::uint64_t>(quotient);
    metadata.maximum_cm[axis] = minimum_cm[axis] + length_cm[axis];
  }

  return FixedGrid(metadata);
}

std::array<double, 3> FixedGrid::widths_cm() const {
  return {indexer_.size_voxel_x(), indexer_.size_voxel_y(), indexer_.size_voxel_z()};
}

RecordIndex FixedGrid::id(const std::array<double, 3> &position_cm) const {
  // A missing measurement is distinct from a valid point outside the configured box.
  for (const auto value : position_cm)
    if (!std::isfinite(value))
      throw std::invalid_argument("Voxel assignment requires a finite position");

  return canonical_id(indexer_.id(position_cm[0], position_cm[1], position_cm[2]));
}

RecordIndex FixedGrid::index(const std::array<RecordIndex, 3> &coordinates) const {
  for (std::size_t axis = 0; axis < 3; ++axis)
    require_index(coordinates[axis], metadata_.counts[axis]);

  return canonical_id(indexer_.index(coordinates[0], coordinates[1], coordinates[2]));
}

std::array<double, 3> FixedGrid::position(RecordIndex id) const {
  require_index(id, size_);
  const auto point = indexer_.position(static_cast<VoxelID_t>(id));
  return {point.x, point.y, point.z};
}

RecordIndex FixedGrid::shift(RecordIndex origin, const std::array<int, 3> &displacement) const {
  require_index(origin, size_);
  return canonical_id(indexer_.shift(static_cast<VoxelID_t>(origin), displacement[0],
                                     displacement[1], displacement[2]));
}

} // namespace supera
