#ifndef SUPERA_DATA_SPACEPOINTVOXELS_H
#define SUPERA_DATA_SPACEPOINTVOXELS_H

#include "supera/data/Grid3D.h"
#include "supera/data/SpacePoint.h"

namespace supera {

/// One representative, with provenance into the owned input collections.
struct SpacePointVoxel {
  RecordIndex id{missing_index};
  RecordIndex source{missing_index};
  RecordIndex point{missing_index};
  std::uint64_t occupancy{};
  std::uint64_t hit_count{};
  double charge{missing_value};
};

struct SpacePointSelectionReport {
  std::uint64_t nonfinite{}, charge_cut{}, outside_grid{}, collisions{}, replacements{};
};

/// Sorted, unique voxel IDs. Sources remain available for traceable representative attributes.
struct SpacePointVoxels {
  std::string name;
  Grid3D grid;
  unsigned int planes{};
  std::vector<SpacePointCollection> sources;
  std::vector<SpacePointVoxel> voxels;
  SpacePointSelectionReport report;
};

} // namespace supera
#endif
