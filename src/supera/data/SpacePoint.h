#ifndef SUPERA_DATA_SPACEPOINT_H
#define SUPERA_DATA_SPACEPOINT_H

#include "supera/data/Validity.h"

#include <array>
#include <optional>
#include <string>
#include <vector>

namespace supera {

struct WireAddress {
  unsigned int cryostat{}, tpc{}, plane{}, wire{};
};

/// Reconstructed wire signal. Times are ticks; amplitudes are ADC; integrals ADC*tick.
/// View/signal codes retain the source convention, not inferred semantic labels.
struct WireHit {
  std::optional<unsigned int> channel;
  std::optional<WireAddress> wire;
  int view_code{}, signal_type_code{};
  double start_tick{missing_value}, end_tick{missing_value};
  double peak_tick{missing_value}, peak_error_tick{missing_value}, rms_tick{missing_value};
  double amplitude_adc{missing_value}, amplitude_error_adc{missing_value};
  double roi_sum_adc{missing_value}, hit_sum_adc{missing_value};
  double integral_adc_tick{missing_value}, integral_error_adc_tick{missing_value};
  double goodness_of_fit{missing_value};
  int multiplicity{};
  std::optional<int> local_index, degrees_of_freedom;
};

/// Product identity is scoped to the containing event. Original keys survive compaction.
struct WireHitCollection {
  std::string source_product;
  std::string source;
  std::uint64_t source_size{}; // Full native size, not just the referenced subset.
  std::vector<RecordIndex> source_keys;
  std::vector<WireHit> hits;
};

struct HitReference {
  RecordIndex collection{missing_index};
  RecordIndex index{missing_index};
};

struct SpacePoint {
  std::array<double, 3> position_cm{missing_value, missing_value, missing_value};
  double chi_square{missing_value};
  std::optional<int> source_id;
};

/// CSR adjacency: links for point i occupy [offsets[i], offsets[i+1]).
/// A present empty table differs from associations that were not requested.
struct SpacePointHits {
  std::vector<std::uint64_t> offsets;
  std::vector<HitReference> links;
  std::vector<WireHitCollection> collections;
};

/// Point index is its original collection offset; extraction does not select/reorder.
/// Optional columns, when present, have exactly points.size() entries.
struct SpacePointCollection {
  std::string source;
  std::vector<SpacePoint> points;
  std::optional<std::vector<std::array<double, 6>>> covariance_cm2; // xx, xy, yy, xz, yz, zz
  std::optional<std::vector<double>> charge;
  std::optional<std::vector<double>> charge_asymmetry;
  std::string charge_convention; // Source-defined units/estimator, not inferred from a slot.
  std::optional<SpacePointHits> associated_hits;
};

} // namespace supera
#endif
