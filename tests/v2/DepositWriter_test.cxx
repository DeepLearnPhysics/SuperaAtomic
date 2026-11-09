#include "supera/io/larcv2/Deposit.h"

#include "TFile.h"
#include "TTree.h"

#include <cstdio>
#include <stdexcept>

namespace {
void check(bool value) {
  if (!value)
    throw std::runtime_error("Deposit output contract failed");
}
} // namespace

int main() {
  using namespace supera;
  DepositCollection input;
  input.name = "energy";
  input.transport = "truth";
  input.particle_count = 2;
  input.grid = Grid3D{{0., 0., 0.}, {2., 2., 2.}, {2, 2, 2}};
  input.voxels = {{missing_index, 0, 1.}, {0, 0, 2.}, {1, 1, 4.}};
  auto output = io::larcv2::project_deposits(input);
  check(output.clusters.as_vector().size() == 3);
  check(output.clusters.as_vector()[2].as_vector()[0].value() == 1.F);
  check(output.tensor.as_vector()[0].value() == 3.F && output.tensor.as_vector()[1].value() == 4.F);
  auto duplicate = input;
  duplicate.voxels.push_back(duplicate.voxels.back());
  bool duplicate_rejected = false;
  try {
    io::larcv2::project_deposits(duplicate);
  } catch (const std::invalid_argument &) {
    duplicate_rejected = true;
  }
  check(duplicate_rejected);
  input.kind = DepositKind::ReadoutIDE;
  bool rejected = false;
  try {
    io::larcv2::project_deposits(input);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected);

  // An event with no records still preserves source presence and accounting.
  EventRecord event;
  event.id = {1, 2, 3};
  event.deposits = {input};
  const char *filename = "deposit-writer-test.root";
  {
    io::larcv2::DepositArchive archive(filename);
    archive.write(event);
    Deposit d;
    d.raw_track_id = -7;
    d.track_id = 7;
    d.reassigned = true;
    d.energy_mev = 0.123456789;
    d.time_ns = 123.456789;
    event.id.event = 4;
    event.deposits[0].deposits = {d};
    archive.write(event);
    archive.close();
    archive.close();
  }
  TFile file(filename, "READ");
  auto *rows = file.Get<TTree>("deposits");
  check(rows->GetEntries() == 1 && file.Get<TTree>("events")->GetEntries() == 2);
  double energy{};
  rows->SetBranchAddress("energy_mev", &energy);
  rows->GetEntry(0);
  check(energy == 0.123456789);
  rows->ResetBranchAddresses();
  file.Close();
  std::remove(filename);
}
