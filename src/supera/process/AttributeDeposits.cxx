#include "supera/process/AttributeDeposits.h"

#include "supera/process/FixedGrid.h"

#include <map>
#include <set>
#include <unordered_map>

namespace supera {
ReassignmentEvidence reassignment_evidence(const std::vector<DepositCollection> &sources,
                                           const std::string &transport) {
  ReassignmentEvidence result;
  for (const auto &source : sources) {
    if (source.kind == DepositKind::EnergyPoint || source.transport != transport)
      continue;
    for (const auto &deposit : source.deposits)
      if (deposit.reassigned && deposit.original_id && deposit.track_id &&
          *deposit.original_id != *deposit.track_id)
        result[*deposit.original_id].insert(*deposit.track_id);
  }
  return result;
}

void attribute_deposits(DepositCollection &collection, const TransportCollection &transport,
                        const ReassignmentEvidence &readout) {
  if (collection.transport != transport.name)
    throw std::invalid_argument("Deposit transport namespace mismatch");
  if (collection.grid && collection.kind == DepositKind::ReadoutIDE)
    throw std::invalid_argument("Readout IDEs cannot be projected as physical deposited energy");
  std::optional<FixedGrid> grid;
  if (collection.grid)
    grid.emplace(*collection.grid);
  collection.particle_count = record_index(transport.particles.size());
  collection.report = {};
  collection.voxels.clear();

  // Sparse lookup preserves exact source identities without offset-sized arrays.
  std::unordered_map<std::int32_t, RecordIndex> indices;
  for (std::size_t i = 0; i < transport.particles.size(); ++i) {
    const auto &p = transport.particles[i];
    if (!p.track_id || !indices.emplace(*p.track_id, record_index(i)).second)
      throw std::invalid_argument("Missing or duplicate transport identity during attribution");
    if (p.primary_ancestor != missing_index)
      require_index(p.primary_ancestor, transport.particles.size());
  }
  const auto lookup = [&](std::optional<std::int32_t> id) {
    if (!id)
      return missing_index;
    const auto found = indices.find(*id);
    return found == indices.end() ? missing_index : found->second;
  };

  // Cache available primary evidence for suppressed IDs once per event.
  std::unordered_map<std::int32_t, std::set<RecordIndex>> roots;
  for (const auto &[ancestor, descendants] : transport.suppression) {
    const auto index = lookup(ancestor);
    if (index == missing_index)
      continue;
    const auto root = transport.particles[index].primary_ancestor;
    if (root == missing_index)
      continue;
    for (auto descendant : descendants)
      roots[descendant].insert(root);
  }
  std::map<std::pair<RecordIndex, RecordIndex>, double> voxels;
  for (auto &deposit : collection.deposits) {
    deposit.original_particle = lookup(deposit.original_id);
    deposit.assigned_particle = lookup(deposit.track_id);
    deposit.attributed_particle = deposit.primary_ancestor = missing_index;
    deposit.attribution = DepositAttribution::Unresolved;
    std::set<RecordIndex> evidence;
    bool conflict = false;
    for (auto index : {deposit.original_particle, deposit.assigned_particle}) {
      if (index == missing_index)
        continue;
      const auto &p = transport.particles[index];
      conflict =
          conflict || p.ancestry == AncestryState::Conflict || p.ancestry == AncestryState::Cycle;
      if (p.primary_ancestor != missing_index)
        evidence.insert(p.primary_ancestor);
    }
    if (deposit.original_id) {
      const auto found = roots.find(*deposit.original_id);
      if (found != roots.end())
        evidence.insert(found->second.begin(), found->second.end());
    }
    if (deposit.track_id) {
      const auto found = roots.find(*deposit.track_id);
      if (found != roots.end())
        evidence.insert(found->second.begin(), found->second.end());
    }

    // Cross-product evidence may recover affiliation even when the suppression map
    // stores NoParticleId. Keep this path distinct from both direct identity and that map.
    std::set<RecordIndex> readout_candidates;
    if (deposit.original_id) {
      const auto found = readout.find(*deposit.original_id);
      if (found != readout.end()) {
        for (auto ancestor : found->second) {
          const auto index = lookup(ancestor);
          if (index == missing_index)
            continue;
          const auto &p = transport.particles[index];
          conflict = conflict || p.ancestry == AncestryState::Conflict ||
                     p.ancestry == AncestryState::Cycle;
          if (p.primary_ancestor != missing_index) {
            evidence.insert(p.primary_ancestor);
            readout_candidates.insert(index);
          }
        }
      }
    }

    // Exact creator wins when consistent. Fallback identity is explicitly weaker evidence.
    if (conflict || evidence.size() > 1) {
      deposit.attribution = DepositAttribution::Conflict;
    } else {
      if (evidence.size() == 1)
        deposit.primary_ancestor = *evidence.begin();
      if (deposit.original_particle != missing_index) {
        deposit.attributed_particle = deposit.original_particle;
        deposit.attribution = DepositAttribution::Exact;
      } else if (deposit.assigned_particle != missing_index) {
        deposit.attributed_particle = deposit.assigned_particle;
        deposit.attribution = DepositAttribution::Reassigned;
      } else if (!readout_candidates.empty()) {
        deposit.attributed_particle =
            readout_candidates.size() == 1 ? *readout_candidates.begin() : deposit.primary_ancestor;
        deposit.attribution = DepositAttribution::EvidenceAncestor;
      } else if (deposit.primary_ancestor != missing_index) {
        deposit.attributed_particle = deposit.primary_ancestor;
        deposit.attribution = DepositAttribution::Ancestor;
      }
    }

    // Invalid energy is counted and retained; it never contaminates a sum with NaN.
    const auto state = static_cast<std::size_t>(deposit.attribution);
    ++collection.report.counts[state];
    if (!std::isfinite(deposit.energy_mev) || deposit.energy_mev < 0.) {
      ++collection.report.invalid_energy;
      continue;
    }
    collection.report.energy_mev[state] += deposit.energy_mev;
    require_measurement(collection.report.energy_mev[state]);
    const auto &x = deposit.position_cm;
    if (!std::isfinite(x[0]) || !std::isfinite(x[1]) || !std::isfinite(x[2])) {
      ++collection.report.invalid_position;
      collection.report.invalid_position_energy_mev += deposit.energy_mev;
      require_measurement(collection.report.invalid_position_energy_mev);
      continue;
    }
    if (!grid)
      continue;
    const auto voxel = grid->id(x);
    if (voxel == missing_index) {
      ++collection.report.outside_grid;
      collection.report.outside_energy_mev += deposit.energy_mev;
      require_measurement(collection.report.outside_energy_mev);
      continue;
    }
    auto &energy = voxels[{deposit.attributed_particle, voxel}];
    energy += deposit.energy_mev;
    require_measurement(energy);
  }
  for (const auto &[key, energy] : voxels)
    collection.voxels.push_back({key.first, key.second, energy});
}
} // namespace supera
