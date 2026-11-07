#ifndef SUPERA_DATA_AUXDET_H
#define SUPERA_DATA_AUXDET_H

#include "supera/data/Measurement.h"

#include "supera/data/Validity.h"

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace supera {
enum class AuxDetTechnology { Unspecified, Scintillator };
enum class AuxDetRole { CosmicTagging, BeamTiming, MuonTagging, Veto };
/// A source contract, not a detector role or an output format.
enum class AuxDetSourceSchema { Generic, SbnCrtHit };

struct AuxDetConventions {
  std::string detector;
  std::string coordinate_frame;

  AuxDetTechnology technology{AuxDetTechnology::Unspecified};
  std::vector<AuxDetRole> roles;
  AuxDetSourceSchema source_schema{AuxDetSourceSchema::Generic};
};

struct AuxDetSignal {
  double value{missing_value};
  SignalUnit unit{SignalUnit::Unspecified};
};

/// Preserve integer clocks exactly. Source-defined readings are not normalized
/// timestamps: a name such as ts0_ns does not establish its origin or epoch.
struct AuxDetTimeReading {
  std::string name;
  std::variant<std::monostate, std::uint64_t, double> value;
  TimeUnit unit{TimeUnit::Unspecified};
  TimeReference reference{TimeReference::Unspecified};
  std::string clock;
  std::optional<double> uncertainty;
};

struct AuxDetPosition {
  std::array<double, 3> cm{missing_value, missing_value, missing_value};
  std::optional<std::array<double, 3>> uncertainty_cm;
};

struct AuxDetChannelSignal {
  std::int64_t channel{};
  AuxDetSignal signal;
};

/// Hardware identifiers are source-scoped, not geometry indices. Empty readout
/// groups are representable; do not lose them by flattening into channel hits.
struct AuxDetReadout {
  std::uint64_t id{};
  std::vector<AuxDetChannelSignal> channels;
};

/// A reconstructed localized measurement, not a scintillator slab, raw trigger,
/// track, or cosmic-particle label. Missing estimates stay explicitly absent.
struct AuxDetSpacePoint {
  RecordIndex id{missing_index};

  std::optional<AuxDetPosition> position;
  std::optional<AuxDetSignal> signal;

  std::vector<AuxDetTimeReading> times;

  std::vector<std::uint64_t> readout_ids;
  std::vector<AuxDetReadout> readouts;

  std::optional<std::int64_t> plane;
  std::string region;
  std::optional<bool> complete;
};

struct AuxDetSpacePointCollection {
  std::string name;
  std::string source;
  AuxDetConventions conventions;
  std::vector<AuxDetSpacePoint> points;
};
} // namespace supera
#endif
