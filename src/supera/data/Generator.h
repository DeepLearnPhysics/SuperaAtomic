#ifndef SUPERA_DATA_GENERATOR_H
#define SUPERA_DATA_GENERATOR_H

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace supera {

/// Source generator units, before any legacy output scaling.
struct GeneratorStep {
  std::array<double, 3> position_cm{};
  double time_ns{};
  std::array<double, 3> momentum_gev_c{};
  double energy_gev{};
};

/// Generator identity is scoped to the containing source and interaction.
/// This is not a transported track, a semantic label, or reconstructed ancestry.
struct GeneratorParticle {
  std::optional<std::uint64_t> particle_index; // Absent for the separate incoming-neutrino copy.
  std::int32_t track_id{};
  int pdg{};
  int status{};
  std::string process;

  std::optional<GeneratorStep> initial;
};

struct GeneratorNeutrino {
  GeneratorParticle incoming;
  std::int32_t lepton_track_id{};
  int current{};
  int mode{};
  int interaction_type{};
  int target{};
  int nucleon{};
  int quark{};

  double hadronic_mass_gev{};
  double bjorken_x{};
  double inelasticity{};
  double momentum_transfer_gev2{};
  std::optional<double> theta_rad;
  std::vector<GeneratorStep> trajectory;
};

/// One producer-defined MCTruth grouping, potentially containing many particles.
/// This does not assert a single physical interaction or a common vertex.
struct GeneratorRecord {
  std::uint64_t truth_index{};
  int origin{}; // Native simb origin code; no inferred particle origin.

  std::optional<GeneratorNeutrino> neutrino;
  std::vector<GeneratorParticle> particles;
};

struct GeneratorSource {
  std::uint64_t source_index{}; // Configured slot, independent of availability.
  std::string source;
  bool present{false};

  std::vector<GeneratorRecord> records;
};

struct GeneratorCollection {
  std::string name;
  std::vector<GeneratorSource> sources;
};

} // namespace supera

#endif
