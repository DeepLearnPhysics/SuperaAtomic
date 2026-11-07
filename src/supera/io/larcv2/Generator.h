#ifndef SUPERA_IO_LARCV2_GENERATOR_H
#define SUPERA_IO_LARCV2_GENERATOR_H

#include "larcv/core/DataFormat/Neutrino.h"
#include "larcv/core/DataFormat/Particle.h"
#include "supera/data/Generator.h"

namespace supera::io::larcv2 {

struct LegacyGeneratorProducts {
  std::vector<larcv::Particle> particles;
  std::vector<larcv::Neutrino> neutrinos;
};

/// Explicitly reproduce pinned SuperaMCTruth output, including documented quirks.
LegacyGeneratorProducts project_legacy_generator(const GeneratorCollection &input);

} // namespace supera::io::larcv2

#endif
