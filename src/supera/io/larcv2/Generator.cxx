#include "Generator.h"

#include <limits>
#include <stdexcept>

namespace supera::io::larcv2 {
namespace {

short legacy_short(int value) {
  if (value < std::numeric_limits<short>::min() || value > std::numeric_limits<short>::max())
    throw std::invalid_argument("Generator interaction code exceeds LArCV2 short range");

  return static_cast<short>(value);
}

larcv::Particle legacy_particle(const GeneratorParticle &input, larcv::MCTIndex_t source) {
  if (!input.initial)
    throw std::invalid_argument("Generator particle has no initial kinematics");

  // LArCV2 uses unsigned IDs; retain the legacy conversion of signed generator IDs.
  // This conversion is confined to compatibility output, never the shared identity.
  larcv::Particle result;
  result.mct_index(source);
  result.track_id(static_cast<unsigned int>(input.track_id));
  result.pdg_code(input.pdg);
  result.creation_process(input.process);

  // Preserve the historical mixed units: momentum scales, energy does not.
  const auto &step = *input.initial;
  const auto &p = step.momentum_gev_c;
  const auto &x = step.position_cm;
  result.momentum(p[0] * 1.e3, p[1] * 1.e3, p[2] * 1.e3);
  result.position(x[0], x[1], x[2], step.time_ns);
  result.energy_init(step.energy_gev);

  // These are legacy placeholders, not extracted generator ancestry.
  result.parent_track_id(static_cast<unsigned int>(input.track_id));
  result.parent_pdg_code(input.pdg);
  result.ancestor_track_id(static_cast<unsigned int>(input.track_id));
  result.ancestor_pdg_code(input.pdg);

  return result;
}

} // namespace

LegacyGeneratorProducts project_legacy_generator(const GeneratorCollection &input) {
  LegacyGeneratorProducts result;
  std::size_t present_source_index = 0;

  for (const auto &source : input.sources) {
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
      // Append the incoming neutrino before the selected final-state particles.
      if (truth.neutrino) {
        const auto &input_neutrino = *truth.neutrino;
        if (!input_neutrino.theta_rad)
          throw std::invalid_argument("Generator neutrino has no defined theta");
        auto particle = legacy_particle(input_neutrino.incoming, legacy_index);
        particle.nu_current_type(legacy_short(input_neutrino.current));
        particle.nu_interaction_type(legacy_short(input_neutrino.interaction_type));
        result.particles.push_back(std::move(particle));

        // Keep interaction metadata in the separate Neutrino product.
        larcv::Neutrino neutrino;
        neutrino.mct_index(legacy_index);
        neutrino.nu_track_id(static_cast<unsigned int>(input_neutrino.incoming.track_id));
        neutrino.lepton_track_id(static_cast<unsigned int>(input_neutrino.lepton_track_id));
        neutrino.current_type(legacy_short(input_neutrino.current));
        neutrino.interaction_mode(legacy_short(input_neutrino.mode));
        neutrino.interaction_type(legacy_short(input_neutrino.interaction_type));
        neutrino.target(input_neutrino.target);
        neutrino.nucleon(input_neutrino.nucleon);
        neutrino.quark(input_neutrino.quark);
        neutrino.hadronic_invariant_mass(input_neutrino.hadronic_mass_gev);
        neutrino.bjorken_x(input_neutrino.bjorken_x);
        neutrino.inelasticity(input_neutrino.inelasticity);
        neutrino.momentum_transfer(input_neutrino.momentum_transfer_gev2);
        neutrino.theta(*input_neutrino.theta_rad);

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
          const auto &position = point.position_cm;
          const auto &momentum = point.momentum_gev_c;
          neutrino.add_trajectory_point(position[0], position[1], position[2], point.time_ns,
                                        momentum[0], momentum[1], momentum[2], point.energy_gev);
        }
        result.neutrinos.push_back(std::move(neutrino));
      }

      // Selection belongs to shared processing; reject accidentally unselected input.
      for (const auto &input_particle : truth.particles) {
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
