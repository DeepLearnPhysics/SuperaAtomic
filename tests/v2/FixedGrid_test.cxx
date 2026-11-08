#include "supera/process/FixedGrid.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}

template <class Exception, class Action> void rejects(Action action) {
  try {
    action();
  } catch (const Exception &) {
    return;
  }
  throw std::runtime_error("Invalid grid operation accepted");
}
} // namespace

int main() try {
  using namespace supera;

  // Match the frozen ICARUS configuration without copying metadata from the output.
  const auto grid =
      FixedGrid::from_lengths({-412.89, -181.86, -894.951}, {1843.2, 1843.2, 1843.2}, {.3, .3, .3});
  require(grid.metadata().counts == std::array<std::uint64_t, 3>{6144, 6144, 6144},
          "ICARUS counts");
  require(grid.size() == 231928233984LL, "64-bit grid size");
  require(grid.id(grid.metadata().minimum_cm) == 0 &&
              grid.id(grid.metadata().maximum_cm) == grid.size() - 1,
          "Both outer faces are included");
  require(grid.shift(3774873600LL, {0, 0, 0}) == 3774873600LL, "Large origin remains exact");
  require(grid.shift(0, {-1, 0, 0}) == missing_index, "Outside neighbor uses canonical absence");
  require(grid.shift(0, {std::numeric_limits<int>::min(), 0, 0}) == missing_index,
          "Extreme displacement cannot overflow");

  // Probe all outer faces on each axis, including adjacent representable values.
  for (std::size_t axis = 0; axis < 3; ++axis) {
    auto point = grid.metadata().minimum_cm;
    point[axis] = std::nextafter(point[axis], -std::numeric_limits<double>::infinity());
    require(grid.id(point) == missing_index, "Just below minimum");
    point = grid.metadata().maximum_cm;
    point[axis] = std::nextafter(point[axis], std::numeric_limits<double>::infinity());
    require(grid.id(point) == missing_index, "Just above maximum");
    point[axis] =
        std::nextafter(grid.metadata().maximum_cm[axis], grid.metadata().minimum_cm[axis]);
    require(grid.id(point) >= 0, "Just inside maximum");
  }

  // Non-cubic, exactly representable geometry makes the flattening oracle independent.
  const FixedGrid small({{0., -3., 10.}, {2., 3., 18.}, {2, 3, 4}});
  require(small.widths_cm() == std::array<double, 3>{1., 2., 2.}, "Effective voxel widths");
  for (RecordIndex z = 0; z < 4; ++z)
    for (RecordIndex y = 0; y < 3; ++y)
      for (RecordIndex x = 0; x < 2; ++x) {
        const auto expected = x + 2 * y + 6 * z;
        const std::array<double, 3> center{static_cast<double>(x) + .5,
                                           -2. + 2 * static_cast<double>(y),
                                           11. + 2 * static_cast<double>(z)};
        require(small.index({x, y, z}) == expected && small.id(center) == expected &&
                    small.position(expected) == center,
                "X-fastest indexing and voxel centers");
      }
  require(small.id({1., -3., 10.}) == 1, "Internal face selects the higher bin");
  const auto truncated = FixedGrid::from_lengths({0., 0., 0.}, {1., 1., 1.}, {.3, .3, .3});
  require(truncated.metadata().counts[0] == 3 && truncated.widths_cm()[0] == 1. / 3.,
          "Requested widths determine truncated counts, not effective widths");

  // Missing/invalid coordinates are errors, not physical out-of-volume decisions.
  rejects<std::invalid_argument>([&] { grid.id({missing_value, 0., 0.}); });
  rejects<std::invalid_argument>([&] {
    grid.id({0., std::numeric_limits<double>::infinity(), 0.});
  });
  rejects<std::out_of_range>([&] { grid.position(missing_index); });
  rejects<std::out_of_range>([&] { grid.position(grid.size()); });
  rejects<std::out_of_range>([&] { grid.shift(-2, {0, 0, 0}); });
  rejects<std::out_of_range>([&] { small.index({2, 0, 0}); });
  rejects<std::out_of_range>([&] { small.index({0, -1, 0}); });

  // Invalid metadata must fail before division, casts, or legacy count multiplication.
  rejects<std::invalid_argument>([] { FixedGrid invalid(Grid3D{}); });
  const auto reject_change = [&](auto change) {
    auto metadata = small.metadata();
    change(metadata);
    rejects<std::invalid_argument>([&] { FixedGrid invalid(metadata); });
  };
  reject_change([](auto &m) { m.counts[0] = 0; });
  reject_change([](auto &m) { m.counts[0] = std::numeric_limits<std::uint64_t>::max(); });
  reject_change([](auto &m) { m.counts = {1ULL << 32, 1ULL << 32, 1}; });
  reject_change([](auto &m) { m.maximum_cm[1] = m.minimum_cm[1]; });
  reject_change([](auto &m) { m.maximum_cm[1] = m.minimum_cm[1] - 1.; });
  reject_change([](auto &m) { m.maximum_cm[1] = missing_value; });
  reject_change([](auto &m) {
    m.minimum_cm[0] = -std::numeric_limits<double>::max();
    m.maximum_cm[0] = std::numeric_limits<double>::max();
  });
  reject_change([](auto &m) {
    m.minimum_cm[0] = 0.;
    m.maximum_cm[0] = std::numeric_limits<double>::denorm_min();
    m.counts[0] = 2;
  });
  for (const auto width : {0., -1., missing_value, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::denorm_min(), 2.})
    rejects<std::invalid_argument>([&] {
      FixedGrid::from_lengths({0., 0., 0.}, {1., 1., 1.}, {width, 1., 1.});
    });

  return 0;

} catch (const std::exception &error) {

  std::cerr << error.what() << '\n';
  return 1;
}
