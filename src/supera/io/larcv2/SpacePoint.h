#ifndef SUPERA_IO_LARCV2_SPACEPOINT_H
#define SUPERA_IO_LARCV2_SPACEPOINT_H

#include "supera/data/SpacePointVoxels.h"
#include "supera/io/larcv2/Writer.h"

#include "larcv/core/DataFormat/Voxel3D.h"

namespace supera::io::larcv2 {
std::vector<std::string> spacepoint_labels(const SpacePointOutput &configuration);
std::vector<larcv::SparseTensor3D> project_spacepoints(const SpacePointVoxels &collection,
                                                       const SpacePointOutput &configuration);
} // namespace supera::io::larcv2
#endif
