#include "supera/process/SelectSpacePoints.h"

#include <iostream>
#include <stdexcept>

namespace {
void require(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}

template <class F> void rejects(F action) {
  try {
    action();
  } catch (const std::exception &) {
    return;
  }
  throw std::runtime_error("Malformed point input accepted");
}

supera::SpacePointCollection sample() {
  supera::SpacePointCollection source;
  source.source = "fixture";
  source.charge = std::vector<double>{10., 100., 20., 20., 5.};
  for (int i = 0; i < 5; ++i)
    source.points.push_back({{i == 4 ? 1.5 : .5, .5, .5}, double(i), i});
  auto &data = source.associated_hits.emplace();
  data.offsets = {0, 3, 5, 8, 11, 12};
  supera::WireHitCollection hits{"product", "hits", 3, {0, 1, 2}, {}};
  for (unsigned int plane = 0; plane < 3; ++plane) {
    supera::WireHit hit;
    hit.wire = supera::WireAddress{0, 1, plane, 10 + plane};
    hit.integral_adc_tick = 30. + plane;
    hits.hits.push_back(hit);
  }
  data.collections.push_back(hits);
  for (int count : {3, 2, 3, 3, 1})
    for (int plane = 0; plane < count; ++plane)
      data.links.push_back({0, plane});
  return source;
}
} // namespace

int main() try {
  using namespace supera;
  const auto grid = FixedGrid::from_lengths({0, 0, 0}, {2, 2, 2}, {1, 1, 1});
  SpacePointSelection options;
  options.float_charge = true;
  auto source = sample();
  const auto selected = select_spacepoints("points", {source}, grid, options);
  require(selected.voxels.size() == 2 && selected.voxels[0].point == 3 &&
              selected.voxels[0].charge == 20. && selected.voxels[0].occupancy == 4 &&
              selected.report.collisions == 3 && selected.report.replacements == 2,
          "Hit-count, charge, and last exact tie ranking");
  require(selected.sources[0].points.size() == 5 && source.points.size() == 5,
          "All raw points survive selection");

  // Ranking precision is explicit; it is not dictated by an output class.
  source.charge->at(2) = 20. + 1.e-8;
  require(select_spacepoints("p", {source}, grid, options).voxels[0].point == 3,
          "Float-precision exact tie");
  options.float_charge = false;
  require(select_spacepoints("p", {source}, grid, options).voxels[0].point == 2,
          "Double-precision ranking");
  options.prefer_more_hits = false;
  require(select_spacepoints("p", {source}, grid, options).voxels[0].point == 1,
          "Charge-only policy");
  source.associated_hits.reset();
  require(select_spacepoints("p", {source}, grid, options).voxels.size() == 2,
          "Points-only rasterization");
  options.prefer_more_hits = true;
  rejects([&] { select_spacepoints("p", {source}, grid, options); });

  // Selection is ordered across source collections, not hash-map iteration order.
  options = {};
  source = sample();
  const auto multiple = select_spacepoints("p", {source, source}, grid, options);
  require(multiple.voxels[0].source == 1 && multiple.voxels[0].occupancy == 8,
          "Cross-source collisions");
  source.points[0].position_cm[0] = missing_value;
  source.charge->at(1) = -1.;
  source.points[2].position_cm[0] = 3.;
  source.points[3].position_cm[0] = 2.;
  const auto filtered = select_spacepoints("p", {source}, grid, options);
  require(filtered.report.nonfinite == 1 && filtered.report.charge_cut == 1 &&
              filtered.report.outside_grid == 1 && filtered.voxels.size() == 1,
          "Nonfinite, cut, outside, and inclusive upper boundary");
  options.shift_cm = {1., 0., 0.};
  require(select_spacepoints("p", {sample()}, grid, options).voxels[0].id == 1, "Configured shift");
  SpacePointCollection empty;
  empty.charge.emplace();
  empty.associated_hits.emplace().offsets = {0};
  require(select_spacepoints("p", {empty}, grid, options).voxels.empty(), "Empty input");

  // Reject malformed graphs before a higher-ranked candidate can partially replace a voxel.
  options = {};
  for (int corruption = 0; corruption < 12; ++corruption) {
    source = sample();
    auto &data = *source.associated_hits;
    if (corruption == 0)
      source.charge.reset();
    if (corruption == 1)
      data.offsets.pop_back();
    if (corruption == 2)
      data.offsets[2] = 2;
    if (corruption == 3)
      data.links[0].collection = -1;
    if (corruption == 4)
      data.links[0].index = 99;
    if (corruption == 5)
      data.collections[0].hits[0].wire.reset();
    if (corruption == 6)
      data.collections[0].hits[0].wire->plane = 1;
    if (corruption == 7)
      data.collections[0].hits[0].wire->tpc = 3;
    if (corruption == 8)
      data.collections[0].source_keys[0] = 9;
    if (corruption == 9)
      data.offsets[1] = 4;
    if (corruption == 10)
      data.collections[0].hits[0].wire->plane = 3;
    if (corruption == 11)
      data.collections[0].source_keys.pop_back();
    rejects([&] { select_spacepoints("p", {source}, grid, options); });
  }
  options.minimum_charge = 2.;
  options.maximum_charge = 1.;
  rejects([&] { select_spacepoints("p", {}, grid, options); });
  options = {};
  options.shift_cm[0] = INFINITY;
  rejects([&] { select_spacepoints("p", {}, grid, options); });
  return 0;

} catch (const std::exception &error) {

  std::cerr << error.what() << '\n';
  return 1;
}
