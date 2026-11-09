#ifndef SUPERA_DATA_TRANSPORT_H
#define SUPERA_DATA_TRANSPORT_H

#include "supera/data/Generator.h"

#include <map>
#include <set>

namespace supera {

/// Explicit source association; track numbers are never generator particle indices.
struct TruthLink {
  std::string product;
  RecordIndex truth_index{missing_index};
  RecordIndex generated_index{missing_index};
  int origin{};
};

/// Preserved source representation, not a claim that simulation saved every step.
enum class TransportStorage { Retained, Dropped };
enum class AncestryState { Primary, Direct, Bridged, Unresolved, Conflict, Cycle };

struct TransportTrajectoryProcess {
  RecordIndex step{missing_index};
  unsigned char source_code{};
  std::string name;
};

struct TransportParticle {
  std::string product;
  RecordIndex source_index{missing_index};
  std::optional<std::int32_t> track_id;
  std::optional<std::int32_t> mother_id;
  std::vector<std::int32_t> daughter_ids;
  int pdg{unspecified_pdg};
  std::optional<int> status;
  std::string process;
  std::string end_process;
  TransportStorage storage{TransportStorage::Retained};
  double mass_gev{missing_value};
  double weight{missing_value};
  std::array<double, 3> polarization{};
  std::array<double, 4> generator_vertex{}; // Native generator convention, not transport cm/ns.
  std::optional<int> rescatter;
  std::vector<TransportTrajectoryProcess> trajectory_processes;
  std::vector<GeneratorStep> trajectory;
  std::optional<TruthLink> truth;

  // Derived references never overwrite the original source relationship.
  bool primary{false};
  RecordIndex direct_parent{missing_index};
  RecordIndex available_parent{missing_index};
  RecordIndex primary_ancestor{missing_index};
  AncestryState ancestry{AncestryState::Unresolved};
};

/// One configured track-ID namespace. Separate overlay namespaces into collections.
/// Suppression entries relate an ancestor to descendants, not immediate children.
struct TransportCollection {
  std::string name;
  std::vector<TransportParticle> particles;
  std::map<std::int32_t, std::set<std::int32_t>> suppression;
  std::set<std::int32_t> root_mother_ids{0};
};

} // namespace supera
#endif
