#ifndef SUPERA_DATA_GRID3D_H
#define SUPERA_DATA_GRID3D_H

#include "supera/data/Validity.h"

#include <array>
#include <cstdint>

namespace supera {

/// Detector-independent Cartesian grid, in centimeters, with x-fastest indexing.
/// Both outer faces are included; the upper face belongs to the last voxel.
/// Counts are authoritative. Effective widths are (maximum - minimum) / count.
/// A default record is unspecified and must be validated before use.
struct Grid3D {
  std::array<double, 3> minimum_cm{missing_value, missing_value, missing_value};
  std::array<double, 3> maximum_cm{missing_value, missing_value, missing_value};
  std::array<std::uint64_t, 3> counts{};
};

} // namespace supera
#endif
