#include "Generator.h"

#include <limits>
#include <stdexcept>

namespace supera::io::larcv2 {
namespace {

unsigned int legacy_track(std::optional<std::int32_t> value) {
  return value ? static_cast<unsigned int>(*value) : std::numeric_limits<unsigned int>::max();
}

int legacy_code(std::optional<int> value) {
  return value.value_or(-1);
}

double legacy_kinematic(double value) {
  if (std::isinf(value))
    throw std::invalid_argument("Noncanonical infinite generator measurement");
  return value;
}

short legacy_short(std::optional<int> input) {
  if (!input)
    return -1;
  const auto value = *input;
  if (value < std::numeric_limits<short>::min() || value > std::numeric_limits<short>::max())
    throw std::invalid_argument("Generator interaction code exceeds LArCV2 short range");

  return static_cast<short>(value);
}

void validate_step(const GeneratorStep &step) {
  for (const auto value : step.position_cm)
    require_measurement(value);
  for (const auto value : step.momentum_gev_c)
    require_measurement(value);
  require_measurement(step.time_ns);
  require_measurement(step.energy_gev);
}

double legacy_momentum(double value) {
  const auto scaled = value * 1.e3;
  require_measurement(scaled);
  return scaled;
}

larcv::Particle legacy_particle(const GeneratorParticle &input, larcv::MCTIndex_t source) {
  if (!input.initial)
    throw std::invalid_argument("Generator particle has no initial kinematics");

  // LArCV2 uses unsigned IDs; retain the legacy conversion of signed generator IDs.
  // This conversion is confined to compatibility output, never the shared identity.
  larcv::Particle result;
  result.mct_index(source);
  result.track_id(legacy_track(input.track_id));
  result.pdg_code(input.pdg);
  result.creation_process(input.process);

  // Preserve the historical mixed units: momentum scales, energy does not.
  const auto &step = *input.initial;
  validate_step(step);
  const auto &p = step.momentum_gev_c;
  const auto &x = step.position_cm;
  result.momentum(legacy_momentum(p[0]), legacy_momentum(p[1]), legacy_momentum(p[2]));
  result.position(x[0], x[1], x[2], step.time_ns);
  result.energy_init(step.energy_gev);

  // These are legacy placeholders, not extracted generator ancestry.
  result.parent_track_id(legacy_track(input.track_id));
  result.parent_pdg_code(input.pdg);
  result.ancestor_track_id(legacy_track(input.track_id));
  result.ancestor_pdg_code(input.pdg);

  return result;
}

} // namespace

LegacyGeneratorProducts project_legacy_generator(const GeneratorCollection &input) {
  LegacyGeneratorProducts result;
  std::size_t present_source_index = 0;

  for (const auto &source : input.sources) {
    require_identity(source.source_index);
    // Reproduce the old counter, which advances for present-empty sources only.
    if (!source.present) {
      if (!source.records.empty())
        throw std::invalid_argument("Missing generator source contains interactions");
      continue;
    }
    if (present_source_index >= std::numeric_limits<larcv::MCTIndex_t>::max())
      throw std::invalid_argument("Generator source index exceeds LArCV2 range");
    const auto legacy_index = static_cast<larcv::MCTIndex_t>(present_source_index++);

    for (const auto &truth : source.records) {
      require_identity(truth.truth_index);
      // Append the incoming neutrino before the selected final-state particles.
      if (truth.neutrino) {
        const auto &input_neutrino = *truth.neutrino;
        if (!std::isfinite(input_neutrino.theta_rad))
          throw std::invalid_argument("Generator neutrino has no defined theta");
        auto particle = legacy_particle(input_neutrino.incoming, legacy_index);
        particle.nu_current_type(legacy_short(input_neutrino.current));
        particle.nu_interaction_type(legacy_short(input_neutrino.interaction_type));
        result.particles.push_back(std::move(particle));

        // Keep interaction metadata in the separate Neutrino product.
        larcv::Neutrino neutrino;
        neutrino.mct_index(legacy_index);
        neutrino.nu_track_id(legacy_track(input_neutrino.incoming.track_id));
        neutrino.lepton_track_id(legacy_track(input_neutrino.lepton_track_id));
        neutrino.current_type(legacy_short(input_neutrino.current));
        neutrino.interaction_mode(legacy_short(input_neutrino.mode));
        neutrino.interaction_type(legacy_short(input_neutrino.interaction_type));
        neutrino.target(legacy_code(input_neutrino.target));
        neutrino.nucleon(legacy_code(input_neutrino.nucleon));
        neutrino.quark(legacy_code(input_neutrino.quark));
        neutrino.hadronic_invariant_mass(legacy_kinematic(input_neutrino.hadronic_mass_gev));
        neutrino.bjorken_x(legacy_kinematic(input_neutrino.bjorken_x));
        neutrino.inelasticity(legacy_kinematic(input_neutrino.inelasticity));
        neutrino.momentum_transfer(legacy_kinematic(input_neutrino.momentum_transfer_gev2));
        neutrino.theta(input_neutrino.theta_rad);

        // Unlike Particle, the legacy Neutrino momentum remains in GeV/c.
        const auto &incoming = input_neutrino.incoming;
        const auto &step = *incoming.initial;
        const auto &p = step.momentum_gev_c;
        const auto &x = step.position_cm;
        neutrino.pdg_code(incoming.pdg);
        neutrino.creation_process(incoming.process);
        neutrino.momentum(p[0], p[1], p[2]);
        neutrino.position(x[0], x[1], x[2], step.time_ns);
        neutrino.energy_init(step.energy_gev);

        for (const auto &point : input_neutrino.trajectory) {
          validate_step(point);
          const auto &position = point.position_cm;
          const auto &momentum = point.momentum_gev_c;
          neutrino.add_trajectory_point(position[0], position[1], position[2], point.time_ns,
                                        momentum[0], momentum[1], momentum[2], point.energy_gev);
        }
        result.neutrinos.push_back(std::move(neutrino));
      }

      // Selection belongs to shared processing; reject accidentally unselected input.
      for (const auto &input_particle : truth.particles) {
        require_identity(input_particle.particle_index);
        if (input_particle.status != 1)
          throw std::invalid_argument("Generator truth must be selected before serialization");
        auto particle = legacy_particle(input_particle, legacy_index);
        particle.shape(input_particle.pdg == 11 || input_particle.pdg == 22 ? larcv::kShapeShower
                                                                            : larcv::kShapeTrack);
        result.particles.push_back(std::move(particle));
      }
    }
  }

  return result;
}

} // namespace supera::io::larcv2
