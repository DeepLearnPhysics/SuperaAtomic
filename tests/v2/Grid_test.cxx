#include "supera/data/ImageMeta3D.h"

#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

} // namespace

int main() try {
  // The ICARUS grid has more than 2^31 voxels; even a zero shift must keep its ID.
  supera::ImageMeta3D grid;
  grid.set(-412.89, -181.86, -894.951, -412.89 + 1843.2, -181.86 + 1843.2, -894.951 + 1843.2, 6144,
           6144, 6144);
  constexpr supera::VoxelID_t high_id = 3774873600ULL;
  require(grid.index(0, 0, 100) == high_id, "Worked high-ID example");
  require(grid.shift(high_id, 0, 0, 0) == high_id, "Zero shift narrowed the origin ID");
  require(grid.shift(high_id, 1, 2, 3) == grid.index(1, 2, 103), "Positive high-ID shift");
  require(grid.shift(grid.index(1, 2, 103), -1, -2, -3) == high_id, "Negative high-ID shift");

  // Axis boundaries must be checked independently, without wrapping to another row.
  const auto invalid = supera::kINVALID_VOXELID;
  const auto last = grid.index(6143, 6143, 6143);
  require(grid.shift(last, 0, 0, 0) == last, "Last voxel identity");
  require(grid.shift(last, 1, 0, 0) == invalid, "Positive x boundary");
  require(grid.shift(last, 0, 1, 0) == invalid, "Positive y boundary");
  require(grid.shift(last, 0, 0, 1) == invalid, "Positive z boundary");
  require(grid.shift(0, -1, 0, 0) == invalid, "Negative x boundary");
  require(grid.shift(0, 0, -1, 0) == invalid, "Negative y boundary");
  require(grid.shift(0, 0, 0, -1) == invalid, "Negative z boundary");
  require(grid.shift(last, -6143, -6143, -6143) == 0, "Exact displacement to minimum");
  require(grid.shift(0, 6143, 6143, 6143) == last, "Exact displacement to maximum");

  // Extreme displacements and invalid IDs must never overflow intermediate arithmetic.
  require(grid.shift(high_id, std::numeric_limits<int>::min(), 0, 0) == invalid, "INT_MIN shift");
  require(grid.shift(high_id, std::numeric_limits<int>::max(), 0, 0) == invalid, "INT_MAX shift");
  require(grid.shift(last + 1, 0, 0, 0) == invalid, "Out-of-range origin");
  require(grid.shift(invalid, 0, 0, 0) == invalid, "Invalid origin sentinel");

  // Non-cubic grids expose accidental use of another axis's dimension.
  supera::ImageMeta3D small;
  small.set(0., 0., 0., 2., 3., 4., 2, 3, 4);
  for (std::size_t z = 0; z < 4; ++z)
    for (std::size_t y = 0; y < 3; ++y)
      for (std::size_t x = 0; x < 2; ++x) {
        const auto origin = small.index(x, y, z);
        for (int dz = -1; dz <= 1; ++dz)
          for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
              const int next_x = static_cast<int>(x) + dx;
              const int next_y = static_cast<int>(y) + dy;
              const int next_z = static_cast<int>(z) + dz;
              const bool inside = next_x >= 0 && next_x < 2 && next_y >= 0 && next_y < 3 &&
                                  next_z >= 0 && next_z < 4;
              const auto expected =
                  inside ? static_cast<supera::VoxelID_t>(next_x + 2 * next_y + 6 * next_z)
                         : invalid;
              require(small.shift(origin, dx, dy, dz) == expected, "Small-grid neighbor identity");
            }
      }

  // Match the existing coordinate-query contract for an uninitialized grid.
  bool rejected = false;
  try {
    supera::ImageMeta3D{}.shift(0, 0, 0, 0);
  } catch (const std::exception &) {
    rejected = true;
  }
  require(rejected, "Uninitialized grid accepted");

  return 0;

} catch (const std::exception &error) {

  std::cerr << error.what() << '\n';
  return 1;
}
