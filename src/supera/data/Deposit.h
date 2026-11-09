#ifndef SUPERA_DATA_DEPOSIT_H
#define SUPERA_DATA_DEPOSIT_H

#include "supera/data/Grid3D.h"

#include <optional>
#include <string>
#include <vector>

namespace supera {
enum class DepositKind { EnergyPoint, ReadoutIDE, EnergyStep };

struct DepositSegment {
  std::array<double, 3> start_cm{missing_value, missing_value, missing_value};
  std::array<double, 3> end_cm{missing_value, missing_value, missing_value};
  double start_ns{missing_value}, end_ns{missing_value};
  std::optional<std::int32_t> pdg;
  double photons{missing_value};
  double scintillation_yield_ratio{missing_value};
};
enum class DepositAttribution {
  Exact,
  Reassigned,
  Ancestor,
  Unresolved,
  Conflict,
  EvidenceAncestor
};

struct Deposit {
  RecordIndex source_index{missing_index};
  RecordIndex channel_index{missing_index};
  RecordIndex ide_index{missing_index};
  std::optional<std::uint32_t> channel;
  std::optional<std::uint32_t> tdc; // Electronics clock count, never a simulation time.
  std::int32_t raw_track_id{};
  std::optional<std::int32_t> raw_original_id;
  std::optional<std::int32_t> track_id;
  std::optional<std::int32_t> original_id;
  bool reassigned{false};
  std::array<double, 3> position_cm{missing_value, missing_value, missing_value};
  double time_ns{missing_value};
  double energy_mev{missing_value};
  double electrons{missing_value};
  std::optional<DepositSegment> segment; // Full simulation steps only.

  // Exact creator, upstream assigned identity, and attributed ancestor stay distinct.
  RecordIndex original_particle{missing_index};
  RecordIndex assigned_particle{missing_index};
  RecordIndex attributed_particle{missing_index};
  RecordIndex primary_ancestor{missing_index};
  DepositAttribution attribution{DepositAttribution::Unresolved};
};

struct DepositVoxel {
  RecordIndex particle{missing_index}; // Missing denotes the separate unassigned cluster.
  RecordIndex voxel{missing_index};
  double energy_mev{};
};

struct DepositReport {
  std::array<std::uint64_t, 6> counts{};
  std::array<double, 6> energy_mev{};
  std::uint64_t invalid_energy{}, invalid_position{}, outside_grid{};
  double outside_energy_mev{}, invalid_position_energy_mev{};
};

struct DepositCollection {
  std::string name;
  std::string product;
  std::string transport;
  DepositKind kind{DepositKind::EnergyPoint};
  std::vector<Deposit> deposits;
  std::vector<std::uint32_t> channel_inventory; // Source order, including empty channels.
  std::optional<Grid3D> grid;
  std::vector<DepositVoxel> voxels;
  RecordIndex particle_count{};
  DepositReport report;
};
} // namespace supera
#endif
