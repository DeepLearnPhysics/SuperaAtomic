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
  const supera::EventRecord empty;
  if (!empty.optical.empty())
    throw std::runtime_error("New event contains stale data");
  const supera::OpticalConventions unknown;
  if (unknown.pe_index_basis != supera::OpticalIndexBasis::Unspecified ||
      unknown.absolute_time_unit != supera::TimeUnit::Unspecified || !unknown.detector.empty())
    throw std::runtime_error("Unspecified detector conventions must remain unknown");
} catch (const std::exception &error) {
  std::cerr << error.what() << "\n";
  return 1;
}
