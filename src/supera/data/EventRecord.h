#ifndef SUPERA_DATA_EVENTRECORD_H
#define SUPERA_DATA_EVENTRECORD_H

#include "supera/data/AuxDet.h"
#include "supera/data/Deposit.h"
#include "supera/data/Generator.h"
#include "supera/data/Optical.h"
#include "supera/data/SpacePointVoxels.h"
#include "supera/data/Transport.h"

#include <cstdint>
#include <vector>

namespace supera {
struct EventID {
  std::uint64_t run{};
  std::uint64_t subrun{};
  std::uint64_t event{};
};

/// Initial v2 event envelope. Carries the currently migrated input products.
/// Existing particle EventInput/EventOutput remain intact until their port.
/// The envelope owns its records and contains no framework pointers or services.
struct EventRecord {
  EventID id;
  std::vector<OpticalFlashCollection> optical;
  std::vector<AuxDetSpacePointCollection> auxiliary;
  std::vector<GeneratorCollection> generator;
  std::vector<TransportCollection> transport;
  std::vector<DepositCollection> deposits;
  std::vector<SpacePointVoxels> spacepoints;
};
} // namespace supera
#endif
