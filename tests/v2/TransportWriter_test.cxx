#include "supera/io/larcv2/Transport.h"

#include "TFile.h"
#include "TTree.h"

#include <cstdio>
#include <stdexcept>

namespace {
void check(bool condition) {
  if (!condition)
    throw std::runtime_error("Transport output contract failed");
}
} // namespace

int main() {
  using namespace supera;
  TransportCollection input;
  input.name = "transport";
  TransportParticle primary;
  primary.track_id = 1;
  primary.mother_id = 0;
  primary.pdg = -13;
  primary.primary = true;
  primary.primary_ancestor = 0;
  primary.trajectory.push_back({{1., 2., 3.}, 4., {0., 0., 2.}, 2.1});
  input.particles.push_back(primary);
  auto child = primary;
  child.track_id = 5;
  child.mother_id = 99;
  child.primary = false;
  child.available_parent = 0;
  child.trajectory.clear();
  input.particles.push_back(child);
  const auto result = io::larcv2::project_transport(input);
  check(result.size() == 2 && result[0].parent_id() == 0);
  check(result[0].energy_init() == 2100. && result[0].pz() == 2000.);
  check(result[1].parent_id() == larcv::kINVALID_INSTANCEID);
  check(result[1].parent_track_id() == 99 && result[1].ancestor_track_id() == 1);
  check(result[1].parent_position().t() == io::larcv2::unavailable_transport_coordinate);
  check(std::isnan(result[1].energy_init()));

  // A finite but unrepresentable time is retained in the archive, never cast to an integer.
  input.particles[0].trajectory[0].time_ns = 1.e20;
  check(io::larcv2::project_transport(input)[0].position().t() ==
        io::larcv2::unavailable_transport_coordinate);
  input.particles[0].track_id = -1;
  bool rejected = false;
  try {
    io::larcv2::project_transport(input);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected);
  input.particles[0].track_id = 1;

  // Empty events still have collection metadata. Serialized trajectories retain original units.
  EventRecord event;
  event.id = {1, 2, 3};
  event.transport = {input};
  const char *filename = "transport-writer-test.root";
  {
    io::larcv2::TransportArchive archive(filename);
    archive.write(event);
    event.id.event = 4;
    event.transport[0].particles.clear();
    archive.write(event);
    archive.close();
    archive.close();
  }
  TFile file(filename, "READ");
  check(file.Get<TTree>("transport")->GetEntries() == 2);
  check(file.Get<TTree>("events")->GetEntries() == 2);
  std::vector<double> *trajectory = nullptr;
  auto *tree = file.Get<TTree>("transport");
  tree->SetBranchAddress("trajectory", &trajectory);
  tree->GetEntry(0);
  check(trajectory->at(3) == 1.e20 && trajectory->at(7) == 2.1);
  tree->ResetBranchAddresses();
  file.Close();
  std::remove(filename);
}
