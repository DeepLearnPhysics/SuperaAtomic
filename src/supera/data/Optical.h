#ifndef SUPERA_DATA_OPTICAL_H
#define SUPERA_DATA_OPTICAL_H

#include "supera/data/Measurement.h"

#include <cstdint>
#include <string>
#include <vector>

namespace supera {

/// How positions in the PE vector are interpreted; unknown must stay explicit.
enum class OpticalIndexBasis { Unspecified, OpticalChannel, OpticalDetector };

/// Collection-level conventions, independent of a particular detector framework.
/// Empty strings mean unspecified, not a default detector or coordinate system.
/// Extend this typed contract when a real adapter requires another convention.
struct OpticalConventions {
  std::string detector;
  std::string region;
  std::string coordinate_frame;

  OpticalIndexBasis pe_index_basis{OpticalIndexBasis::Unspecified};
  TimeUnit absolute_time_unit{TimeUnit::Unspecified};
  std::string absolute_time_reference;
};

/// Existing reconstructed flash attributes, with no validity cuts or recalibration.
/// Values may contain upstream sentinels/defaults. A stored value is not by itself
/// evidence that a quantity was measured. Vectors retain source ordering.
struct OpticalFlash {
  std::uint64_t id{}; // Local to its named collection, not a global event identity.

  double time_us{}; // Trigger-relative, as supplied by the input adapter.
  double time_width_us{};
  double absolute_time{}; // Interpretation belongs to OpticalConventions.
  unsigned int frame{};

  std::vector<double> pe;
  bool in_beam_frame{};
  int on_beam_time{}; // Keep the source integer; do not narrow to bool.
  double fast_to_total{};

  double x_cm{};
  double x_width_cm{};
  double y_cm{};
  double y_width_cm{};
  double z_cm{};
  double z_width_cm{};

  std::vector<double> wire_centers;
  std::vector<double> wire_widths;
};

struct OpticalFlashCollection {
  std::string name;   // Stable logical name, distinct from native source provenance.
  std::string source; // Opaque provenance; algorithms must not parse framework tags.
  OpticalConventions conventions;
  std::vector<OpticalFlash> flashes;
};

} // namespace supera
#endif
