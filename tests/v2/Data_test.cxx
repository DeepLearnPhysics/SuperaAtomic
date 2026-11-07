#include "supera/data/EventRecord.h"

#include <iostream>
#include <stdexcept>

int main() try {
  supera::EventRecord first;
  first.id = {1, 2, 3};
  supera::OpticalFlashCollection collection;
  collection.name = "east";
  collection.source = "native-product";
  collection.conventions.detector = "example";
  collection.conventions.pe_index_basis = supera::OpticalIndexBasis::OpticalChannel;
  supera::OpticalFlash flash;
  flash.pe = {1., 2.};
  collection.flashes.push_back(flash);
  first.optical.push_back(collection);
  auto copy = first;
  first.optical.front().flashes.front().pe.clear();
  if (copy.optical.front().flashes.front().pe.size() != 2 || copy.id.event != 3)
    throw std::runtime_error("Event record must own its values");
  supera::AuxDetSpacePointCollection auxiliary;
  auxiliary.name = "beam";
  auxiliary.conventions.technology = supera::AuxDetTechnology::Scintillator;
  auxiliary.conventions.roles = {supera::AuxDetRole::BeamTiming, supera::AuxDetRole::Veto};
  supera::AuxDetSpacePoint point;
  if (point.position || point.signal || point.plane || point.complete)
    throw std::runtime_error("Unavailable auxiliary estimates must stay absent");
  point.readout_ids = {1234};
  point.times.push_back({"counter",
                         (std::uint64_t{1} << 53) + 1,
                         supera::TimeUnit::Tick,
                         supera::TimeReference::SourceDefined,
                         "test-clock",
                         {}});
  auxiliary.points.push_back(point);
  first.auxiliary.push_back(auxiliary);
  auto auxiliary_copy = first;
  first.auxiliary.front().points.front().readout_ids.clear();
  const auto &saved = auxiliary_copy.auxiliary.front();
  if (saved.points.front().readout_ids != std::vector<std::uint64_t>{1234} ||
      std::get<std::uint64_t>(saved.points.front().times.front().value) !=
          (std::uint64_t{1} << 53) + 1 ||
      saved.conventions.roles.size() != 2)
    throw std::runtime_error("Auxiliary records must own values and preserve integer clocks");
  const supera::AuxDetConventions unknown_auxiliary;
  if (unknown_auxiliary.technology != supera::AuxDetTechnology::Unspecified ||
      unknown_auxiliary.source_schema != supera::AuxDetSourceSchema::Generic ||
      !unknown_auxiliary.roles.empty())
    throw std::runtime_error("No implicit auxiliary role or source schema");
  const supera::EventRecord empty;
  if (!empty.optical.empty() || !empty.auxiliary.empty())
    throw std::runtime_error("New event contains stale data");
  const supera::OpticalConventions unknown;
  if (unknown.pe_index_basis != supera::OpticalIndexBasis::Unspecified ||
      unknown.absolute_time_unit != supera::TimeUnit::Unspecified || !unknown.detector.empty())
    throw std::runtime_error("Unspecified detector conventions must remain unknown");

} catch (const std::exception &error) {

  std::cerr << error.what() << "\n";
  return 1;
}
