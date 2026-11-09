#ifndef SUPERA_IO_LARCV2_DEPOSIT_H
#define SUPERA_IO_LARCV2_DEPOSIT_H

#include "larcv/core/DataFormat/Voxel3D.h"
#include "supera/data/EventRecord.h"

#include <memory>

namespace supera::io::larcv2 {
struct DepositProducts {
  larcv::ClusterVoxel3D clusters;
  larcv::SparseTensor3D tensor;
};

DepositProducts project_deposits(const DepositCollection &input);

/// Exact input contributions and attribution evidence, independent of float voxel projection.
class DepositArchive {
public:
  explicit DepositArchive(const std::string &filename);
  ~DepositArchive();
  void write(const EventRecord &event);
  void close();

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};
} // namespace supera::io::larcv2
#endif
