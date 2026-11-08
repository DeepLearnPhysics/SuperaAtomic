#include "supera/process/FixedGrid.h"

int main() {
  const auto grid = supera::FixedGrid::from_lengths({0., 0., 0.}, {2., 3., 4.}, {1., 1., 1.});
  return grid.size() == 24 && grid.id({2., 3., 4.}) == 23 &&
                 grid.id({-1., 0., 0.}) == supera::missing_index
             ? 0
             : 1;
}
