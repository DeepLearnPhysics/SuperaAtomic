#ifndef SUPERA_PROCESS_SELECTGENERATOR_H
#define SUPERA_PROCESS_SELECTGENERATOR_H

#include "supera/data/Generator.h"

#include <algorithm>
#include <stdexcept>

namespace supera {

/// Select final-state generator particles while preserving source identities.
/// An absent origin accepts all origins; zero is a real origin when supplied.
/// Missing sources remain explicit, including after all records are removed.
inline GeneratorCollection select_generator(GeneratorCollection input,
                                            std::optional<int> origin = {}) {
  for (auto &source : input.sources) {
    // Missing products cannot simultaneously contain interaction records.
    if (!source.present && !source.records.empty())
      throw std::invalid_argument("Missing generator source contains records");

    // Filter records before requiring any particle kinematics.
    auto &records = source.records;
    records.erase(
        std::remove_if(records.begin(), records.end(),
                       [&](const auto &truth) { return origin && truth.origin != *origin; }),
        records.end());

    for (auto &truth : records) {
      // Preserve the original source particle indices after status selection.
      auto &particles = truth.particles;
      particles.erase(std::remove_if(particles.begin(), particles.end(),
                                     [](const auto &particle) { return particle.status != 1; }),
                      particles.end());

      // Fail explicitly instead of reading trajectory point zero out of bounds.
      for (const auto &particle : particles)
        if (!particle.initial)
          throw std::invalid_argument(
              "Selected generator particle has no trajectory: " + source.source +
              " truth=" + std::to_string(truth.truth_index) +
              " particle=" + std::to_string(particle.particle_index.value_or(0)));
      if (truth.neutrino && (!truth.neutrino->incoming.initial || !truth.neutrino->theta_rad))
        throw std::invalid_argument("Generator neutrino/lepton has no trajectory: " +
                                    source.source);
    }
  }

  return input;
}

} // namespace supera

#endif
