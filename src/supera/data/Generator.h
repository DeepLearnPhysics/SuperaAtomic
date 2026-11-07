#ifndef SUPERA_DATA_GENERATOR_H
#define SUPERA_DATA_GENERATOR_H

#include "supera/data/Validity.h"

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace supera {

/// Source generator units, before any legacy output scaling.
struct GeneratorStep {
  std::array<double, 3> position_cm{missing_value, missing_value, missing_value};
  double time_ns{missing_value};
  std::array<double, 3> momentum_gev_c{missing_value, missing_value, missing_value};
  double energy_gev{missing_value};
};

/// Generator identity is scoped to the containing source and interaction.
/// This is not a transported track, a semantic label, or reconstructed ancestry.
struct GeneratorParticle {
  RecordIndex particle_index{missing_index}; // Absent for the separate incoming-neutrino copy.
  std::optional<std::int32_t> track_id;
  int pdg{unspecified_pdg};
  std::optional<int> status;
  std::string process;

  std::optional<GeneratorStep> initial;
};

struct GeneratorNeutrino {
  GeneratorParticle incoming;
  std::optional<std::int32_t> lepton_track_id;
  std::optional<int> current;
  std::optional<int> mode;
  std::optional<int> interaction_type;
  std::optional<int> target;
  std::optional<int> nucleon;
  std::optional<int> quark;

  double hadronic_mass_gev{missing_value};
  double bjorken_x{missing_value};
  double inelasticity{missing_value};
  double momentum_transfer_gev2{missing_value};
  double theta_rad{missing_value};
  std::vector<GeneratorStep> trajectory;
};

/// One producer-defined MCTruth grouping, potentially containing many particles.
/// This does not assert a single physical interaction or a common vertex.
struct GeneratorRecord {
  RecordIndex truth_index{missing_index};
  int origin{}; // Native simb origin code; no inferred particle origin.

  std::optional<GeneratorNeutrino> neutrino;
  std::vector<GeneratorParticle> particles;
};

struct GeneratorSource {
  RecordIndex source_index{missing_index}; // Configured slot, independent of availability.
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
