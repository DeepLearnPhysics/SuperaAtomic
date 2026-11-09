#ifndef SUPERA_IO_LARCV2_TRANSPORT_H
#define SUPERA_IO_LARCV2_TRANSPORT_H

#include "supera/data/EventRecord.h"

#include "larcv/core/DataFormat/Particle.h"

#include <memory>

namespace supera::io::larcv2 {

/// Safe, finite writer-only sentinel: the legacy Vertex API cannot accept NaN.
inline constexpr double unavailable_transport_coordinate = -9.e12;
std::vector<larcv::Particle> project_transport(const TransportCollection &input);

/// Companion ROOT tables preserve fields that LArCV2 Particle cannot represent.
/// This is a versioned migration archive, not the proposed future output format.
class TransportArchive {
public:
  explicit TransportArchive(const std::string &filename);
  ~TransportArchive();
  void write(const EventRecord &event);
  void close();

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace supera::io::larcv2
#endif
