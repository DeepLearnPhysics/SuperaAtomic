#include "supera/process/FixedGrid.h"

#include "TFile.h"
#include "TTree.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}

void check_metadata(const supera::FixedGrid &grid, const larcv::Voxel3DMeta &reference) {
  const auto &metadata = grid.metadata();
  require(metadata.minimum_cm ==
                  std::array<double, 3>{reference.min_x(), reference.min_y(), reference.min_z()} &&
              metadata.maximum_cm ==
                  std::array<double, 3>{reference.max_x(), reference.max_y(), reference.max_z()} &&
              metadata.counts == std::array<std::uint64_t, 3>{reference.num_voxel_x(),
                                                              reference.num_voxel_y(),
                                                              reference.num_voxel_z()} &&
              grid.widths_cm() == std::array<double, 3>{reference.size_voxel_x(),
                                                        reference.size_voxel_y(),
                                                        reference.size_voxel_z()} &&
              static_cast<std::uint64_t>(grid.size()) == reference.size(),
          "Fixed-grid metadata differs from LArCV2");
}
} // namespace

int main(int argc, char **argv) try {
  // Construct from the frozen configuration, independently of reference metadata.
  const auto grid = supera::FixedGrid::from_lengths({-412.89, -181.86, -894.951},
                                                    {1843.2, 1843.2, 1843.2}, {.3, .3, .3});
  larcv::Voxel3DMeta legacy;
  legacy.set(-412.89, -181.86, -894.951, -412.89 + 1843.2, -181.86 + 1843.2, -894.951 + 1843.2,
             6144, 6144, 6144);
  check_metadata(grid, legacy);

  // Compare boundary probes with the installed implementation, not a copied formula.
  std::size_t probes = 0;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    const auto low = grid.metadata().minimum_cm[axis];
    const auto high = grid.metadata().maximum_cm[axis];
    for (const auto coordinate : {std::nextafter(low, -INFINITY), low, std::nextafter(low, high),
                                  low + grid.widths_cm()[axis], std::nextafter(high, low), high,
                                  std::nextafter(high, INFINITY)}) {
      auto point = grid.metadata().minimum_cm;
      point[axis] = coordinate;
      const auto expected = legacy.id(point[0], point[1], point[2]);
      require(grid.id(point) == (expected == larcv::kINVALID_VOXELID
                                     ? supera::missing_index
                                     : supera::record_index(expected)),
              "Boundary indexing differs from installed LArCV2");
      ++probes;
    }
  }

  // Optional immutable files exercise all twenty space-point products and event identities.
  std::size_t product_events = 0;
  const std::array<unsigned int, 10> events{373, 385, 409, 433, 445, 457, 493, 529, 541, 565};
  for (int argument = 1; argument < argc; ++argument) {
    TFile file(argv[argument], "READ");
    require(!file.IsZombie(), "Cannot read reference file");
    for (const auto *cryo : {"E", "W"})
      for (const auto *suffix : {"", "_chi2", "_cryo", "_tpc", "_hit_charge0", "_hit_charge1",
                                 "_hit_charge2", "_hit_key0", "_hit_key1", "_hit_key2"}) {
        const auto name = std::string("sparse3d_reco_cryo") + cryo + suffix;
        auto *tree = file.Get<TTree>((name + "_tree").c_str());
        require(tree && tree->GetEntries() == 10, "Missing reference tree or events");
        larcv::EventSparseTensor3D storage;
        auto *record = &storage;
        require(tree->SetBranchAddress((name + "_branch").c_str(), &record) >= 0,
                "Cannot bind reference branch");
        for (std::size_t entry = 0; entry < events.size(); ++entry) {
          require(tree->GetEntry(entry) > 0 && record->run() == 13171 && record->subrun() == 1 &&
                      record->event() == events[entry],
                  "Reference event identity differs");
          check_metadata(grid, record->meta());
          ++product_events;

          // Sample stored voxel centers across the full ordered ID range.
          const auto &voxels = record->as_vector();
          const auto stride = std::max<std::size_t>(1, voxels.size() / 32);
          for (std::size_t index = 0; index < voxels.size(); index += stride) {
            const auto id = voxels[index].id();
            const auto point = record->meta().position(id);
            require(grid.id({point.x, point.y, point.z}) == supera::record_index(id),
                    "Stored voxel center maps to a different ID");
            ++probes;
          }
        }
        tree->ResetBranchAddresses();
      }
  }
  std::cout << "PASS: " << product_events << " reference product-events, " << probes
            << " coordinate probes\n";
  return 0;

} catch (const std::exception &error) {

  std::cerr << error.what() << '\n';
  return 1;
}
