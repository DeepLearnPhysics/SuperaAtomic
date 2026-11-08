#include "supera/io/larcv2/SpacePoint.h"

#include "TFile.h"
#include "TTree.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

#include <filesystem>
#include <iostream>

namespace {
void require(bool value) {
  if (!value)
    throw std::runtime_error("Space-point writer contract failed");
}

template <class F> void rejects(F action) {
  try {
    action();
  } catch (const std::exception &) {
    return;
  }
  throw std::runtime_error("Invalid sparse projection accepted");
}

supera::SpacePointVoxels sample() {
  supera::SpacePointVoxels result;
  result.name = "points";
  result.grid = {{0, 0, 0}, {2, 2, 2}, {2, 2, 2}};
  result.planes = 3;
  supera::SpacePointCollection source;
  source.points.push_back({{.5, .5, .5}, 4., 5});
  auto &data = source.associated_hits.emplace();
  data.offsets = {0, 1};
  data.links = {{0, 0}};
  supera::WireHit hit;
  hit.wire = supera::WireAddress{1, 2, 1, 3};
  hit.integral_adc_tick = 42.;
  hit.amplitude_adc = 7.;
  hit.peak_tick = 8.;
  hit.rms_tick = 2.;
  hit.multiplicity = 1;
  data.collections.push_back({"id", "tag", 10, {9}, {hit}});
  result.sources.push_back(source);
  result.voxels.push_back({0, 0, 0, 2, 1, 12.5});
  return result;
}
} // namespace

int main() try {
  using namespace supera::io::larcv2;
  const SpacePointOutput configuration{"points", 3};
  auto source = sample();
  const auto products = project_spacepoints(source, configuration);
  require(products.size() == 10 && products[0].as_vector()[0].value() == 12.5F &&
          products[1].as_vector()[0].value() == 4.F && products[2].as_vector()[0].value() == 1.F &&
          products[3].as_vector()[0].value() == 2.F && products[4].as_vector()[0].value() == -1.F &&
          products[5].as_vector()[0].value() == 42.F && products[8].as_vector()[0].value() == 9.F);
  const SpacePointOutput extras{
      "points", 3, {"occupancy", "nhits", "hit_amp", "hit_time", "hit_rms", "hit_mult"}};
  const auto extra = project_spacepoints(source, extras);
  require(extra[0].as_vector()[0].value() == 2.F && extra[1].as_vector()[0].value() == 1.F &&
          extra[3].as_vector()[0].value() == 7.F && extra[6].as_vector()[0].value() == 8.F &&
          extra[9].as_vector()[0].value() == 2.F && extra[12].as_vector()[0].value() == 1.F);

  // Keep legacy reverse offsets at the writer boundary, without changing shared keys.
  source.sources.push_back(source.sources.front());
  source.sources.back().associated_hits->collections[0].source_product = "second";
  require(project_spacepoints(source, configuration)[8].as_vector()[0].value() == 19.F);
  source.sources.back().associated_hits->collections[0].source_product = "id";
  require(project_spacepoints(source, configuration)[8].as_vector()[0].value() == 9.F);
  source.sources.back().associated_hits->collections[0].source_size = 11;
  rejects([&] { project_spacepoints(source, configuration); });
  source = sample();
  source.grid.counts[0] = 0;
  rejects([&] { project_spacepoints(source, configuration); });
  source = sample();
  source.sources[0].associated_hits->collections[0].hits[0].integral_adc_tick =
      supera::missing_value;
  require(project_spacepoints(source, configuration)[5].as_vector()[0].value() == -1.F);
  source = sample();
  source.sources[0].associated_hits->offsets[1] = 0;
  source.sources[0].associated_hits->links.clear();
  source.voxels[0].hit_count = 0;
  require(project_spacepoints(source, configuration)[2].as_vector()[0].value() == -1.F);
  source = sample();

  // Float output must never silently round an identity or overflow a measurement.
  source.sources[0].associated_hits->collections[0].source_keys[0] = 16777217;
  source.sources[0].associated_hits->collections[0].source_size = 16777218;
  rejects([&] { project_spacepoints(source, configuration); });
  source = sample();
  source.voxels[0].charge = 1.e100;
  rejects([&] { project_spacepoints(source, configuration); });
  source = sample();
  source.voxels.push_back(source.voxels[0]);
  rejects([&] { project_spacepoints(source, configuration); });
  source = sample();
  source.sources[0].associated_hits.reset();
  rejects([&] { project_spacepoints(source, configuration); });
  require(project_spacepoints(source, {"points", 3, {"", "chi2"}}).size() == 2);
  rejects([&] { spacepoint_labels({"points", 3, {"bad"}}); });
  rejects([&] { spacepoint_labels({"points", 3, {"", ""}}); });

  // Real persistence checks nonempty -> empty -> nonempty event clearing.
  const auto path = std::filesystem::current_path() / "spacepoint-writer-test.root";
  {
    Writer writer(path.string(), {}, {}, {}, {configuration});
    for (std::uint64_t event = 1; event <= 3; ++event) {
      supera::EventRecord record;
      record.id = {1, 2, event};
      record.spacepoints = {sample()};
      if (event == 2)
        record.spacepoints[0].voxels.clear();
      writer.write(record);
    }
    writer.close();
  }
  TFile file(path.c_str(), "READ");
  for (const auto &label : spacepoint_labels(configuration)) {
    auto *tree = file.Get<TTree>(("sparse3d_" + label + "_tree").c_str());
    require(tree && tree->GetEntries() == 3);
    larcv::EventSparseTensor3D storage;
    auto *record = &storage;
    tree->SetBranchAddress(("sparse3d_" + label + "_branch").c_str(), &record);
    for (int entry = 0; entry < 3; ++entry) {
      require(tree->GetEntry(entry) > 0 && record->event() == std::size_t(entry + 1));
      require(record->size() == std::size_t(entry == 1 ? 0 : 1) && record->meta().size() == 8);
    }
    tree->ResetBranchAddresses();
  }
  file.Close();
  std::filesystem::remove(path);
  return 0;

} catch (const std::exception &error) {

  std::cerr << error.what() << '\n';
  return 1;
}
