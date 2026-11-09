#include "supera/process/ResolveTransport.h"

#include <iostream>
#include <stdexcept>

void check(bool condition) {
  if (!condition)
    throw std::runtime_error("Transport ancestry contract failed");
}

supera::TransportParticle particle(int id, std::optional<int> parent, bool primary = false) {
  supera::TransportParticle result;
  result.track_id = id;
  result.mother_id = parent;
  result.process = primary ? "primary" : "secondary";
  return result;
}

int main() try {
  using namespace supera;
  TransportCollection input;
  input.particles = {particle(3, 99), particle(1, 0, true), particle(4, 3)};
  input.suppression[1] = {3, 99};
  resolve_transport(input);
  check(input.particles.size() == 3 && input.particles[1].primary);
  check(input.particles[0].mother_id == 99 && input.particles[0].direct_parent == missing_index);
  check(input.particles[0].available_parent == 1 && input.particles[0].primary_ancestor == 1);
  check(input.particles[0].ancestry == AncestryState::Bridged);
  check(input.particles[2].direct_parent == 0 && input.particles[2].primary_ancestor == 1);

  // Repeated resolution must discard stale derived state.
  input.suppression.clear();
  resolve_transport(input);
  check(input.particles[0].primary_ancestor == missing_index);
  check(input.particles[2].primary_ancestor == missing_index);

  // Conflicting evidence must not silently choose a primary.
  input.particles.push_back(particle(2, 0, true));
  input.suppression[1] = {3};
  input.suppression[2] = {99};
  resolve_transport(input);
  check(input.particles[0].ancestry == AncestryState::Conflict);
  check(input.particles[2].ancestry == AncestryState::Conflict);

  // Nested map candidates describe one chain; choose its nearest available member.
  input.particles = {particle(1, 0, true), particle(2, 1), particle(3, 99)};
  input.suppression = {{1, {2, 3, 99}}, {2, {99}}};
  resolve_transport(input);
  check(input.particles[2].available_parent == 1 && input.particles[2].primary_ancestor == 0);

  // Conflicting map evidence must also be detected when a direct parent exists.
  input.particles = {particle(1, 0, true), particle(2, 0, true), particle(3, 1), particle(4, 3)};
  input.suppression = {{2, {3}}};
  resolve_transport(input);
  check(input.particles[2].direct_parent == 0);
  check(input.particles[2].ancestry == AncestryState::Conflict);
  check(input.particles[3].primary_ancestor == missing_index);
  input.suppression.clear();
  input.particles[0].truth = TruthLink{"a", 0, 0, 1};
  input.particles[2].truth = TruthLink{"b", 0, -1, 1};
  resolve_transport(input);
  check(input.particles[2].ancestry == AncestryState::Conflict);

  // Cycles include self-parenting secondaries; neither becomes a primary.
  input.particles = {particle(1, 2), particle(2, 1), particle(3, 3)};
  input.suppression.clear();
  resolve_transport(input);
  for (const auto &p : input.particles)
    check(p.ancestry == AncestryState::Cycle && !p.primary);

  // Large and negative identities are literal keys, not vector offsets or abs(ID).
  input.root_mother_ids = {10000000};
  input.particles = {particle(10000001, 10000000, true), particle(-3, 10000001)};
  resolve_transport(input);
  check(input.particles[1].primary_ancestor == 0);
  input.particles.push_back(input.particles[0]);
  bool rejected = false;
  try {
    resolve_transport(input);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected);
  input.particles = {TransportParticle{}};
  rejected = false;
  try {
    resolve_transport(input);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected);
  input.particles.clear();
  resolve_transport(input);

  // Long chains use an iterative memoized walk rather than stack recursion.
  input.root_mother_ids = {0};
  input.particles.push_back(particle(1, 0, true));
  for (int i = 2; i <= 20000; ++i)
    input.particles.push_back(particle(i, i - 1));
  resolve_transport(input);
  check(input.particles.back().primary_ancestor == 0);

} catch (const std::exception &error) {

  std::cerr << error.what() << '\n';
  return 1;
}
