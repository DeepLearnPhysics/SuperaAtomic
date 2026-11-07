#include "AuxDet.h"

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace supera::io::larcv2 {
namespace {
float exact_float(double value) {
  if (std::isfinite(value) && std::abs(value) > std::numeric_limits<float>::max())
    throw std::invalid_argument("Auxiliary value exceeds LArCV2 float range");

  const auto narrowed = static_cast<float>(value);
  if (!std::isnan(value) && static_cast<double>(narrowed) != value)
    throw std::invalid_argument("Auxiliary value loses precision in LArCV2");

  return narrowed;
}

std::uint8_t board_id(std::uint64_t id) {
  if (id > std::numeric_limits<std::uint8_t>::max())
    throw std::invalid_argument("Auxiliary board ID exceeds LArCV2 range");

  return static_cast<std::uint8_t>(id);
}

int integer(std::int64_t value) {
  if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
    throw std::invalid_argument("Auxiliary integer exceeds LArCV2 range");

  return static_cast<int>(value);
}

float pe(const AuxDetSignal &signal) {
  if (signal.unit != SignalUnit::Photoelectron)
    throw std::invalid_argument("LArCV2 CRT projection requires PE, not ADC or energy");

  return exact_float(signal.value);
}

template <class T>
T time_value(const AuxDetSpacePoint &point, const std::string &name, TimeUnit unit) {
  const AuxDetTimeReading *found = nullptr;

  for (const auto &reading : point.times) {
    if (reading.name != name)
      continue;
    if (found)
      throw std::invalid_argument("Duplicate auxiliary time reading");
    found = &reading;
  }

  if (!found || found->unit != unit || found->reference != TimeReference::SourceDefined ||
      found->clock != "sbn_crt_source" || found->uncertainty ||
      !std::holds_alternative<T>(found->value))
    throw std::invalid_argument("Unsupported auxiliary clock contract for LArCV2");

  return std::get<T>(found->value);
}
} // namespace

std::vector<larcv::CRTHit> project_auxiliary(const AuxDetSpacePointCollection &collection) {
  // Require an explicitly supported source contract before projecting any points.
  if (collection.conventions.source_schema != AuxDetSourceSchema::SbnCrtHit)
    throw std::invalid_argument("No LArCV2 projection for this auxiliary source schema");

  std::vector<larcv::CRTHit> result;
  result.reserve(collection.points.size());

  for (const auto &point : collection.points) {
    // Reject missing estimates and information the legacy schema cannot retain.
    if (point.id >= larcv::kINVALID_INSTANCEID || !point.position ||
        !point.position->uncertainty_cm || !point.signal || !point.plane || point.complete ||
        point.times.size() != 5)
      throw std::invalid_argument("Incomplete or unsupported auxiliary record for LArCV2");

    // Set every LArCV scalar explicitly: its default constructor does not do so.
    larcv::CRTHit hit;
    hit.id(point.id);
    hit.peshit(pe(*point.signal));

    // Preserve raw clock fields; do not reinterpret their origins.
    hit.ts0_s(time_value<std::uint64_t>(point, "ts0_s", TimeUnit::Second));
    hit.ts0_s_corr(time_value<double>(point, "ts0_s_corr", TimeUnit::Second));
    hit.ts0_ns(time_value<double>(point, "ts0_ns", TimeUnit::Nanosecond));
    hit.ts0_ns_corr(time_value<double>(point, "ts0_ns_corr", TimeUnit::Nanosecond));
    hit.ts1_ns(time_value<double>(point, "ts1_ns", TimeUnit::Nanosecond));

    // Transfer source location annotations and reconstructed spatial estimates.
    hit.plane(integer(*point.plane));
    hit.tagger(point.region);
    const auto &xyz = point.position->cm;
    const auto &error = *point.position->uncertainty_cm;
    hit.x_pos(exact_float(xyz[0]));
    hit.y_pos(exact_float(xyz[1]));
    hit.z_pos(exact_float(xyz[2]));
    hit.x_err(exact_float(error[0]));
    hit.y_err(exact_float(error[1]));
    hit.z_err(exact_float(error[2]));

    // Validate hardware ranges while preserving board order and channel groups.
    std::vector<std::uint8_t> boards;
    for (const auto id : point.readout_ids)
      boards.push_back(board_id(id));
    hit.feb_id(boards);

    std::map<std::uint8_t, std::vector<std::pair<int, float>>> signals;
    for (const auto &readout : point.readouts) {
      auto [entry, inserted] =
          signals.emplace(board_id(readout.id), std::vector<std::pair<int, float>>{});
      if (!inserted)
        throw std::invalid_argument("Duplicate auxiliary readout group");
      for (const auto &channel : readout.channels)
        entry->second.emplace_back(integer(channel.channel), pe(channel.signal));
    }
    hit.pesmap(signals);

    result.push_back(std::move(hit));
  }

  return result;
}
} // namespace supera::io::larcv2
