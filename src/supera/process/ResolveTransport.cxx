#include "supera/process/ResolveTransport.h"

#include <algorithm>
#include <unordered_map>

namespace supera {
void resolve_transport(TransportCollection &collection) {
  auto &particles = collection.particles;
  std::unordered_map<std::int32_t, RecordIndex> indices;
  std::unordered_map<std::int32_t, std::set<RecordIndex>> candidates;

  // Build sparse identities, independent of collection order and overlay offset size.
  for (std::size_t i = 0; i < particles.size(); ++i) {
    auto &particle = particles[i];
    if (!particle.track_id || !indices.emplace(*particle.track_id, record_index(i)).second)
      throw std::invalid_argument("Transport particle identity is missing or duplicated");
    particle.direct_parent = particle.available_parent = particle.primary_ancestor = missing_index;
    particle.ancestry = AncestryState::Unresolved;
    particle.primary = particle.process == "primary" && particle.mother_id &&
                       collection.root_mother_ids.count(*particle.mother_id);
  }

  // Reverse only map links whose ancestor has a persisted particle record.
  for (const auto &[ancestor, descendants] : collection.suppression) {
    const auto found = indices.find(ancestor);
    if (found == indices.end())
      continue;
    for (auto descendant : descendants)
      candidates[descendant].insert(found->second);
  }

  // Direct parents take precedence. A suppression link must never replace one.
  for (std::size_t i = 0; i < particles.size(); ++i) {
    auto &particle = particles[i];
    if (particle.primary) {
      particle.primary_ancestor = record_index(i);
      particle.ancestry = AncestryState::Primary;
      continue;
    }
    if (particle.mother_id) {
      const auto found = indices.find(*particle.mother_id);
      if (found != indices.end()) {
        particle.direct_parent = particle.available_parent = found->second;
        particle.ancestry = AncestryState::Direct;
        continue;
      }
    }

    // Map entries may name several generations in one chain. Retain only candidates
    // not proven to be ancestors of another candidate; unrelated alternatives conflict.
    auto evidence = candidates[*particle.track_id];
    if (particle.mother_id) {
      const auto &parent_evidence = candidates[*particle.mother_id];
      evidence.insert(parent_evidence.begin(), parent_evidence.end());
    }
    evidence.erase(record_index(i));
    const auto reaches = [&](RecordIndex descendant, RecordIndex ancestor) {
      std::set<RecordIndex> seen;
      std::vector<RecordIndex> pending{descendant};
      while (!pending.empty()) {
        const auto next = pending.back();
        pending.pop_back();
        if (next == ancestor)
          return true;
        if (!seen.insert(next).second)
          continue;
        const auto &record = particles[next];
        if (record.mother_id && !record.primary) {
          const auto mother = indices.find(*record.mother_id);
          if (mother != indices.end())
            pending.push_back(mother->second);
        }
        for (auto candidate : candidates[*record.track_id])
          pending.push_back(candidate);
      }
      return false;
    };
    auto nearest = evidence;
    for (auto ancestor : evidence)
      for (auto descendant : evidence)
        if (ancestor != descendant && reaches(descendant, ancestor) &&
            !reaches(ancestor, descendant))
          nearest.erase(ancestor);
    evidence = std::move(nearest);
    if (evidence.size() == 1) {
      particle.available_parent = *evidence.begin();
      particle.ancestry = AncestryState::Bridged;
    } else if (evidence.size() > 1) {
      particle.ancestry = AncestryState::Conflict;
    }
  }

  // Memoized iterative walks avoid recursion depth and quadratic chain traversal.
  std::vector<unsigned char> visited(particles.size());
  for (std::size_t start = 0; start < particles.size(); ++start) {
    if (visited[start] == 2)
      continue;
    std::vector<RecordIndex> path;
    auto current = record_index(start);
    RecordIndex root = missing_index;
    AncestryState failure = AncestryState::Unresolved;
    while (current != missing_index) {
      auto &particle = particles[current];
      if (visited[current] == 1) {
        failure = AncestryState::Cycle;
        break;
      }
      if (visited[current] == 2 || particle.primary) {
        root = particle.primary_ancestor;
        failure = particle.ancestry;
        break;
      }
      visited[current] = 1;
      path.push_back(current);
      failure = particle.ancestry;
      current = particle.available_parent;
    }
    for (auto index : path) {
      auto &particle = particles[index];
      particle.primary_ancestor = root;
      if (root == missing_index &&
          (failure == AncestryState::Cycle || failure == AncestryState::Conflict))
        particle.ancestry = failure;
      visited[index] = 2;
    }
    if (particles[start].primary)
      visited[start] = 2;
  }

  // Check independent evidence against the chosen chain; preserve contradictory raw links.
  std::vector<std::vector<RecordIndex>> children(particles.size());
  std::vector<RecordIndex> conflicts;
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const auto &particle = particles[i];
    if (particle.available_parent != missing_index)
      children[particle.available_parent].push_back(record_index(i));
    const auto root = particle.primary_ancestor;
    if (root == missing_index)
      continue;
    bool conflict = false;
    for (auto candidate : candidates[*particle.track_id]) {
      const auto other_root = particles[candidate].primary_ancestor;
      conflict = conflict || (other_root != missing_index && other_root != root);
    }
    const auto &origin = particles[root].truth;
    if (particle.truth && origin)
      conflict = conflict || particle.truth->product != origin->product ||
                 particle.truth->truth_index != origin->truth_index;
    if (conflict)
      conflicts.push_back(record_index(i));
  }
  while (!conflicts.empty()) {
    const auto index = conflicts.back();
    conflicts.pop_back();
    auto &particle = particles[index];
    if (particle.ancestry == AncestryState::Conflict && particle.primary_ancestor == missing_index)
      continue;
    particle.ancestry = AncestryState::Conflict;
    particle.primary_ancestor = missing_index;
    conflicts.insert(conflicts.end(), children[index].begin(), children[index].end());
  }
}
} // namespace supera
