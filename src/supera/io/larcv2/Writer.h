#ifndef SUPERA_IO_LARCV2_WRITER_H
#define SUPERA_IO_LARCV2_WRITER_H

#include "supera/io/EventWriter.h"

#include <memory>
#include <string>
#include <vector>

namespace supera::io::larcv2 {
/// Compatibility projection of optical and supported auxiliary records into the existing LArCV2
/// schema. Collection provenance/conventions are not representable in that schema. This class is
/// serial and must not be called concurrently.
class Writer final : public EventWriter {
public:
  Writer(const std::string &filename, const std::vector<std::string> &collections,
         const std::vector<std::string> &auxiliary_collections = {},
         const std::vector<std::string> &generator_collections = {});
  ~Writer() override;
  Writer(const Writer &) = delete;
  Writer &operator=(const Writer &) = delete;
  void write(const EventRecord &event) override;
  void close() override;

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};
} // namespace supera::io::larcv2
#endif
