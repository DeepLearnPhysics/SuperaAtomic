#include "supera/algorithm/LArTPCMLReco3D.h"
#include "supera/data/ImageMeta3D.h"
#include <iostream>
#include <stdexcept>

void require(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}

int main() {
  supera::ImageMeta3D grid;
  grid.set(0, 0, 0, 10, 10, 10, 10, 10, 10);
  require(grid.valid(), "invalid grid");
  require(grid.id(1.2, 2.3, 3.4) == grid.index(1, 2, 3), "coordinate/index mismatch");
  require(grid.id(-0.1, 0, 0) == supera::kINVALID_VOXELID, "outside point accepted");
  // Characterize the upstream inclusive upper boundary; do not silently change it.
  require(grid.id(10, 10, 10) == grid.index(9, 9, 9), "upper-bound behavior changed");
  supera::EDep first;
  first.x = 1.2;
  first.y = 2.3;
  first.z = 3.4;
  first.t = 0;
  first.e = 2;
  first.dedx = 1;
  auto second = first;
  second.x = 1.8;
  second.e = 3;
  auto voxels = grid.edep2voxelset({first, second});
  require(voxels.size() == 1, "deposits not combined");
  require(voxels.as_vector().front().value() == 5, "energy not conserved");
  supera::LArTPCMLReco3D labeler;
  std::cout << "Atomic linked; grid boundaries and deposit aggregation OK\n";
}
