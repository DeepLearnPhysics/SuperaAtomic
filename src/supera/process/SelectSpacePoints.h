#ifndef SUPERA_PROCESS_SELECTSPACEPOINTS_H
#define SUPERA_PROCESS_SELECTSPACEPOINTS_H

#include "supera/data/SpacePointVoxels.h"
#include "supera/process/FixedGrid.h"

namespace supera {

struct SpacePointSelection {
  double minimum_charge{0.};
  double maximum_charge{std::numeric_limits<double>::max()};
  std::array<double, 3> shift_cm{};
  bool prefer_more_hits{true};
  bool float_charge{false}; // Explicit legacy precision for cuts and tie-breaking.
  unsigned int planes{3};
};

/// Select by hit count (optional), then charge, then later source/point on exact ties.
/// Invalid measurements are counted and rejected. Malformed associations throw;
/// they never partly replace an existing representative. Inputs are owned by the result.
SpacePointVoxels select_spacepoints(std::string name, std::vector<SpacePointCollection> sources,
                                    const FixedGrid &grid, const SpacePointSelection &selection);

} // namespace supera
#endif
