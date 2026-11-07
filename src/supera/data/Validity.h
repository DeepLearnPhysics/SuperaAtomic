#ifndef SUPERA_DATA_VALIDITY_H
#define SUPERA_DATA_VALIDITY_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace supera {

using RecordIndex = std::int64_t;
inline constexpr RecordIndex missing_index = -1;
inline constexpr double missing_value = std::numeric_limits<double>::quiet_NaN();
inline constexpr int unspecified_pdg = 0;

/// Only collection indices use -1; native identifiers have their own domains.
inline RecordIndex record_index(std::uint64_t value) {
  if (value > static_cast<std::uint64_t>(std::numeric_limits<RecordIndex>::max()))
    throw std::out_of_range("Collection index exceeds signed 64-bit range");
  return static_cast<RecordIndex>(value);
}

/// Source-specific finite sentinels must be recognized by the input adapter.
inline double finite_or_missing(double value) {
  return std::isfinite(value) ? value : missing_value;
}

/// Shared records permit NaN, but never infinity, as a missing measurement.
inline void require_measurement(double value) {
  if (std::isinf(value))
    throw std::invalid_argument("Infinite measurement violates the shared validity contract");
}

inline void require_identity(RecordIndex value) {
  if (value < 0)
    throw std::invalid_argument("Required record identity is missing or malformed");
}

inline void require_index(RecordIndex value, std::uint64_t size) {
  if (value < 0 || static_cast<std::uint64_t>(value) >= size)
    throw std::out_of_range("Invalid collection reference");
}

} // namespace supera
#endif
