#include "TFile.h"
#include "TTree.h"
#include "larcv/core/DataFormat/EventNeutrino.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "supera/io/larcv2/Writer.h"

#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

void require(bool condition) {
  if (!condition)
    throw std::runtime_error("Generator writer contract failed");
}

template <class F> void rejects(F action) {
  try {
    action();
  } catch (const std::invalid_argument &) {
    return;
  }
  throw std::runtime_error("Unsupported generator projection accepted");
}

} // namespace

int main() try {
  // Use hand-chosen records to test serialization independently of native extraction.
  supera::GeneratorParticle incoming;
  incoming.track_id = -7;
  incoming.pdg = 14;
  incoming.process = "primary";
  incoming.initial = supera::GeneratorStep{{1., 2., 3.}, 4., {0.25, 0.5, 0.75}, 1.25};
  supera::GeneratorNeutrino nu;
  nu.incoming = incoming;
  nu.lepton_track_id = 8;
  nu.current = 1;
  nu.mode = 2;
  nu.interaction_type = 1001;
  nu.target = 1000180400;
  nu.nucleon = 2212;
  nu.quark = 2;
  nu.hadronic_mass_gev = 1.5;
  nu.bjorken_x = 0.25;
  nu.inelasticity = 0.5;
  nu.momentum_transfer_gev2 = 0.75;
  nu.theta_rad = 1.25;
  nu.trajectory = {*incoming.initial, {{5., 6., 7.}, 8., {1., 2., 3.}, 4.}};
  auto positron = incoming;
  positron.pdg = -11;
  positron.status = 1;
  positron.track_id = -9;
  auto electron = positron;
  electron.pdg = 11;
  supera::GeneratorRecord truth{4, 1, nu, {positron, electron}};
  supera::EventRecord event;
  event.id = {1, 2, 3};
  event.generator.push_back({"mpv",
                             {{0, "missing", false, {}},
                              {1, "empty", true, {}},
                              {2, "generator", true, {truth, truth}}}});

  // Reject invalid events before any manager mutation; then write nonempty/empty events.
  const std::string filename = "generator-writer.root";
  using supera::io::larcv2::Writer;
  rejects([&] { Writer writer(filename, {}, {}, {""}); });
  rejects([&] { Writer writer(filename, {}, {}, {"mpv", "mpv"}); });
  {
    Writer writer(filename, {}, {}, {"mpv"});
    auto invalid = event;
    invalid.generator.clear();
    rejects([&] { writer.write(invalid); });
    invalid = event;
    invalid.generator[0].name = "other";
    rejects([&] { writer.write(invalid); });
    invalid = event;
    invalid.generator[0].sources[2].records[0].particles[0].status = 0;
    rejects([&] { writer.write(invalid); });
    invalid = event;
    invalid.generator[0].sources[2].records[0].neutrino->mode = 100000;
    rejects([&] { writer.write(invalid); });
    invalid = event;
    invalid.generator[0].sources[2].records[0].particles[0].initial.reset();
    rejects([&] { writer.write(invalid); });
    invalid = event;
    invalid.generator[0].sources[2].records[0].neutrino->theta_rad.reset();
    rejects([&] { writer.write(invalid); });
    invalid = event;
    invalid.generator[0].sources[2].present = false;
    rejects([&] { writer.write(invalid); });

    writer.write(event);
    event.id.event = 4;
    event.generator[0].sources[2].records.clear();
    writer.write(event);
    writer.close();
  }

  // Both output products must retain legacy units, defaults, ordering, and IDs.
  {
    TFile file(filename.c_str(), "READ");
    auto *particles = file.Get<TTree>("particle_mpv_tree");
    auto *neutrinos = file.Get<TTree>("neutrino_mpv_tree");
    require(particles && neutrinos && particles->GetEntries() == 2 && neutrinos->GetEntries() == 2);
    larcv::EventParticle particle_storage;
    larcv::EventNeutrino neutrino_storage;
    auto *pr = &particle_storage;
    auto *nr = &neutrino_storage;
    particles->SetBranchAddress("particle_mpv_branch", &pr);
    neutrinos->SetBranchAddress("neutrino_mpv_branch", &nr);
    require(particles->GetEntry(0) > 0 && neutrinos->GetEntry(0) > 0);
    require(pr->as_vector().size() == 6 && nr->as_vector().size() == 2);
    require(pr->run() == 1 && pr->subrun() == 2 && pr->event() == 3 && nr->event() == 3);
    for (std::size_t i = 0; i < pr->as_vector().size(); ++i) {
      const auto &p = pr->as_vector()[i];
      require(p.id() == larcv::kINVALID_INDEX && p.mct_index() == 1 && p.px() == 250. &&
              p.py() == 500. && p.pz() == 750. && p.energy_init() == 1.25);
      require(p.x() == 1. && p.y() == 2. && p.z() == 3. && p.t() == 4. &&
              p.creation_process() == "primary");
      require(p.parent_track_id() == p.track_id() && p.ancestor_track_id() == p.track_id() &&
              p.parent_pdg_code() == p.pdg_code() && p.ancestor_pdg_code() == p.pdg_code());
      require(p.mcst_index() == larcv::kINVALID_INDEX && p.energy_deposit() == 0.);
    }
    require(pr->as_vector()[0].track_id() == static_cast<unsigned int>(-7) &&
            pr->as_vector()[0].shape() == larcv::kShapeUnknown &&
            pr->as_vector()[0].nu_current_type() == 1 &&
            pr->as_vector()[0].nu_interaction_type() == 1001);
    require(pr->as_vector()[1].shape() == larcv::kShapeTrack &&
            pr->as_vector()[2].shape() == larcv::kShapeShower);
    const auto &n = nr->as_vector()[0];
    require(n.id() == larcv::kINVALID_INDEX && nr->as_vector()[1].id() == larcv::kINVALID_INDEX &&
            n.mct_index() == 1 && n.nu_track_id() == static_cast<unsigned int>(-7) &&
            n.lepton_track_id() == 8);
    require(n.current_type() == 1 && n.interaction_mode() == 2 && n.interaction_type() == 1001 &&
            n.target() == 1000180400 && n.nucleon() == 2212 && n.quark() == 2);
    require(n.hadronic_invariant_mass() == 1.5 && n.bjorken_x() == 0.25 &&
            n.inelasticity() == 0.5 && n.momentum_transfer() == 0.75 && n.theta() == 1.25);
    require(n.pdg_code() == 14 && n.px() == 0.25 && n.py() == 0.5 && n.pz() == 0.75 &&
            n.energy_init() == 1.25 && n.creation_process() == "primary");
    require(
        n.x() == 1. && n.y() == 2. && n.z() == 3. && n.t() == 4. &&
        n.traj_x() == std::vector<double>{1., 5.} && n.traj_y() == std::vector<double>{2., 6.} &&
        n.traj_z() == std::vector<double>{3., 7.} && n.traj_t() == std::vector<double>{4., 8.} &&
        n.traj_px() == std::vector<double>{0.25, 1.} &&
        n.traj_py() == std::vector<double>{0.5, 2.} &&
        n.traj_pz() == std::vector<double>{0.75, 3.} &&
        n.traj_e() == std::vector<double>{1.25, 4.});
    require(n.lepton_pdg_code() == 0 && n.energy_transfer() == 0. &&
            n.momentum_transfer_mag() == 0.);

    // A valid empty event clears both lists without changing event alignment.
    require(particles->GetEntry(1) > 0 && neutrinos->GetEntry(1) > 0);
    require(pr->event() == 4 && nr->event() == 4 && pr->as_vector().empty() &&
            nr->as_vector().empty());
    particles->ResetBranchAddresses();
    neutrinos->ResetBranchAddresses();
  }
  std::filesystem::remove(filename);

  return 0;

} catch (const std::exception &error) {

  std::cerr << error.what() << '\n';
  return 1;
}
