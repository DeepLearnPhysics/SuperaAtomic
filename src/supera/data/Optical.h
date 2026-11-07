#ifndef SUPERA_DATA_OPTICAL_H
#define SUPERA_DATA_OPTICAL_H

#include "supera/data/Measurement.h"

#include "supera/data/Validity.h"

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

/// Reconstructed flash attributes with canonical missing measurements.
/// Vectors retain source ordering; no calibration or physics cuts are applied.
struct OpticalFlash {
  RecordIndex id{missing_index}; // Local to its named collection, not a global event identity.

  double time_us{missing_value}; // Trigger-relative, as supplied by the input adapter.
  double time_width_us{missing_value};
  double absolute_time{missing_value}; // Interpretation belongs to OpticalConventions.
  unsigned int frame{};

  std::vector<double> pe;
  bool in_beam_frame{};
  int on_beam_time{}; // Keep the source integer; do not narrow to bool.
  double fast_to_total{missing_value};

  double x_cm{missing_value};
  double x_width_cm{missing_value};
  double y_cm{missing_value};
  double y_width_cm{missing_value};
  double z_cm{missing_value};
  double z_width_cm{missing_value};

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
