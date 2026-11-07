#include "supera/io/larcv2/Writer.h"

#include "TFile.h"
#include "TTree.h"
#include "larcv/core/DataFormat/EventFlash.h"

#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>

namespace {
template <class F> void rejects(F &&action) {
  try {
    action();
  } catch (const std::invalid_argument &) {
    return;
  }
  throw std::runtime_error("Expected invalid_argument");
}
} // namespace

int main() {
  using supera::io::larcv2::Writer;
  const std::string filename = "writer-contract.root";
  rejects([&] { Writer writer("", {"east"}); });
  rejects([&] { Writer writer(filename, {}); });
  rejects([&] { Writer writer(filename, {""}); });
  rejects([&] { Writer writer(filename, {"east", "east"}); });
  {
    Writer writer(filename, {"east"});
    supera::EventRecord event;
    event.id = {1, 1, 1};
    rejects([&] { writer.write(event); });
    event.optical.push_back({"west", {}, {}, {}});
    rejects([&] { writer.write(event); });
    event.optical.front().name = "east";
    event.optical.front().flashes.emplace_back();
    event.optical.front().flashes.front().id = supera::missing_index;
    rejects([&] { writer.write(event); });
    event.optical.front().flashes.front().id = 0;
    event.id.event = std::numeric_limits<std::uint64_t>::max();
    rejects([&] { writer.write(event); });
    event.id.event = 1;
    auto &flash = event.optical.front().flashes.front();
    flash.id = 7;
    flash.time_us = -12.5;
    flash.time_width_us = 0.25;
    flash.absolute_time = 123456.5;
    flash.frame = 77;
    flash.pe = {0., 2.5, 7.5};
    flash.in_beam_frame = true;
    flash.on_beam_time = -3;
    flash.fast_to_total = 0.125;
    flash.x_cm = -10.;
    flash.x_width_cm = 1.;
    flash.y_cm = -20.;
    flash.y_width_cm = 2.;
    flash.z_cm = 30.;
    flash.z_width_cm = 3.;
    flash.wire_centers = {4., 5.};
    flash.wire_widths = {0.4, 0.5};
    auto infinite = event;
    infinite.optical.front().flashes.front().y_cm = std::numeric_limits<double>::infinity();
    rejects([&] { writer.write(infinite); });
    writer.write(event);
    event.id.event = 2;
    event.optical.front().flashes.clear();
    writer.write(event);
    event.id.event = 3;
    supera::OpticalFlash missing;
    missing.id = 0;
    missing.y_cm = 12.;
    event.optical.front().flashes = {missing};
    writer.write(event);
    writer.close();
    writer.close();
    bool rejected = false;
    try {
      writer.write(event);
    } catch (const std::logic_error &) {
      rejected = true;
    }
    if (!rejected)
      throw std::runtime_error("Write after close accepted");
  }
  {
    TFile file(filename.c_str(), "READ");
    auto *tree = file.Get<TTree>("opflash_east_tree");
    if (!tree || tree->GetEntries() != 3)
      throw std::runtime_error("Missing writer events");
    larcv::EventFlash storage;
    auto *record = &storage;
    tree->SetBranchAddress("opflash_east_branch", &record);
    if (tree->GetEntry(0) <= 0 || record->as_vector().size() != 1)
      throw std::runtime_error("Missing first flash");
    const auto &f = record->as_vector().front();
    if (record->run() != 1 || record->subrun() != 1 || record->event() != 1 || f.id() != 7 ||
        f.time() != -12.5 || f.timeWidth() != 0.25 || f.absTime() != 123456.5 || f.frame() != 77 ||
        f.PEPerOpDet() != std::vector<double>({0., 2.5, 7.5}) || !f.inBeamFrame() ||
        f.onBeamTime() != -3 || f.fastToTotal() != 0.125 || f.xCenter() != -10. ||
        f.xWidth() != 1. || f.yCenter() != -20. || f.yWidth() != 2. || f.zCenter() != 30. ||
        f.zWidth() != 3. || f.wireCenters() != std::vector<double>({4., 5.}) ||
        f.wireWidths() != std::vector<double>({0.4, 0.5}))
      throw std::runtime_error("Writer changed a field or event identity");
    if (tree->GetEntry(1) <= 0 || record->event() != 2 || !record->as_vector().empty())
      throw std::runtime_error("Writer retained stale flashes");
    if (tree->GetEntry(2) <= 0 || record->as_vector().size() != 1)
      throw std::runtime_error("Missing normalized flash");
    const auto &missing = record->as_vector().front();
    if (missing.xCenter() != std::numeric_limits<double>::max() ||
        missing.xWidth() != std::numeric_limits<double>::max() || missing.yCenter() != 12. ||
        !std::isnan(missing.time()) || !std::isnan(missing.zCenter()))
      throw std::runtime_error("Incorrect optical missing-value projection");
    tree->ResetBranchAddresses();
  }
  std::filesystem::remove(filename);
}
