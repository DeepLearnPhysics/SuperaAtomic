#ifndef SUPERA_DATA_EVENTRECORD_H
#define SUPERA_DATA_EVENTRECORD_H

#include "supera/data/Optical.h"

#include <cstdint>
#include <vector>

namespace supera {
struct EventID {
  std::uint64_t run{};
  std::uint64_t subrun{};
  std::uint64_t event{};
};

/// Initial v2 event envelope. This stage carries optical data only.
/// Existing particle EventInput/EventOutput remain intact until their port.
/// The envelope owns its records and contains no framework pointers or services.
struct EventRecord {
  EventID id;
  std::vector<OpticalFlashCollection> optical;
};
} // namespace supera
#endif
