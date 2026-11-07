#include "supera/process/SelectGenerator.h"

#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition) {
  if (!condition)
    throw std::runtime_error("Generator selection contract failed");
}

template <class F> void rejects(F action) {
  try {
    action();
  } catch (const std::invalid_argument &) {
    return;
  }
  throw std::runtime_error("Missing generator kinematics accepted");
}

} // namespace

int main() try {
  // Construct independent sources with stable indices and multiple origins/statuses.
  supera::GeneratorParticle particle;
  particle.particle_index = 8;
  particle.track_id = -7;
  particle.status = 1;
  particle.initial.emplace();
  supera::GeneratorRecord truth;
  truth.truth_index = 3;
  truth.origin = 1;
  truth.particles = {particle};
  particle.particle_index = 9;
  particle.status = 0;
  particle.initial.reset();
  truth.particles.push_back(particle);
  supera::GeneratorCollection input{
      "mpv", {{0, "missing", false, {}}, {1, "empty", true, {}}, {2, "generator", true, {truth}}}};

  // Selection keeps original identities and does not mutate its input.
  auto selected = supera::select_generator(input, 1);
  const auto &saved = selected.sources[2].records[0];
  require(saved.truth_index == 3 && saved.particles.size() == 1 &&
          saved.particles[0].particle_index == 8);
  require(saved.particles[0].track_id == -7 && input.sources[2].records[0].particles.size() == 2);
  require(!selected.sources[0].present && selected.sources[1].present &&
          selected.sources[1].records.empty());
  require(supera::select_generator(input, 0).sources[2].records.empty());
  require(supera::select_generator({}).sources.empty());

  // A missing source cannot contain records, and selected kinematics must exist.
  auto invalid = input;
  invalid.sources[0].records = {truth};
  rejects([&] { supera::select_generator(invalid); });
  invalid = input;
  invalid.sources[2].records[0].particles[0].initial.reset();
  rejects([&] { supera::select_generator(invalid); });
  require(supera::select_generator(invalid, 2).sources[2].records.empty());

  // Incoming-neutrino requirements apply only to retained interactions.
  invalid = input;
  invalid.sources[2].records[0].neutrino.emplace();
  rejects([&] { supera::select_generator(invalid); });
  auto &nu = invalid.sources[2].records[0].neutrino;
  nu->incoming.initial.emplace();
  rejects([&] { supera::select_generator(invalid); });
  nu->theta_rad = 0.;
  require(supera::select_generator(invalid).sources[2].records[0].neutrino.has_value());

  return 0;

} catch (const std::exception &error) {

  std::cerr << error.what() << '\n';
  return 1;
}
