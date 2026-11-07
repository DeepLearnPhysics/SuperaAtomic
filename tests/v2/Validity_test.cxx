#include "supera/data/Validity.h"
#include <stdexcept>

namespace {
template <class Exception, class Action> void rejects(Action action) {
  try {
    action();
  } catch (const Exception &) {
    return;
  }
  throw std::runtime_error("Invalid value accepted");
}
} // namespace

int main() {
  using namespace supera;
  if (!std::isnan(missing_value) || missing_index != -1 || record_index(65535) != 65535 ||
      finite_or_missing(-1.) != -1. || !std::signbit(finite_or_missing(-0.)) ||
      !std::isnan(finite_or_missing(std::numeric_limits<double>::infinity())))
    throw std::runtime_error("Canonical validity contract");

  // Reject overflow and malformed references, including missing required indices.
  for (const auto value : {missing_index, RecordIndex{3}}) {
    bool rejected = false;
    try {
      require_index(value, 3);
    } catch (const std::out_of_range &) {
      rejected = true;
    }
    if (!rejected)
      throw std::runtime_error("Invalid reference accepted");
  }
  require_index(2, 3);
  // Signed index limits and valid zero must not be confused with absence.
  const auto largest = std::numeric_limits<RecordIndex>::max();
  if (record_index(static_cast<std::uint64_t>(largest)) != largest || record_index(0) != 0)
    throw std::runtime_error("Index boundary conversion");
  require_identity(0);
  require_identity(largest);
  require_measurement(missing_value);
  require_measurement(-1.);
  rejects<std::out_of_range>([] { record_index(std::numeric_limits<std::uint64_t>::max()); });
  rejects<std::invalid_argument>([] { require_identity(missing_index); });
  rejects<std::invalid_argument>([] { require_identity(-2); });
  rejects<std::invalid_argument>(
      [] { require_measurement(std::numeric_limits<double>::infinity()); });
  rejects<std::invalid_argument>(
      [] { require_measurement(-std::numeric_limits<double>::infinity()); });
}
